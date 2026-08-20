#!/usr/bin/env Rscript
# pc_tissue_map.R -- open the PCA black box: link each PC to tissue openness and
#                    decompose which PCs build the positive vs negative SHAP bins.
#
# Runs on the AWS instance (needs the model .rds + the bundled openness atlas).
# Everything it reads is either model internals or the PUBLIC atlas; the two CSVs
# it writes are aggregates over PCs/tissues -- no per-sample / per-patient data.
#
# Two exact computations (both stages of the model are linear):
#
#  1. PC <-> tissue alignment.  For PC k, the bin-loading vector L[,k] is a
#     genome-wide spatial pattern. Its Pearson correlation with tissue t's
#     openness vector C_t over the shared bins says how much that PC's pattern
#     resembles tissue t's open chromatin. We orient by sign(beta_k) so a
#     POSITIVE value means: in the direction this PC pushes predicted FF UP, it
#     aligns with tissue t's openness.
#
#  2. Per-PC contribution to signed-SHAP bin sets.  The effective per-bin weight
#     decomposes exactly:  w_bin = sum_k L[bin,k]*beta_k.  So contrib[bin,k] =
#     L[bin,k]*beta_k is PC k's exact contribution to that bin's weight. Averaged
#     over the top-N most-positive (FF-up) and top-N most-negative (FF-down) bins,
#     it shows which PCs assemble each SHAP direction.
#
# Usage:
#   Rscript pc_tissue_map.R --model ff_model.rds \
#       --atlas analysis/fetal_fraction/reference/ff_openness_atlas_hg19_50kb.csv.gz \
#       --topn 2000 --outdir pc_tissue_out
#
# Outputs (under --outdir):
#   pc_tissue_corr.csv      PC x tissue matrix: pc, beta, var_frac, best_tissue,
#                           then corr_<tissue> columns (beta-oriented)
#   pc_shap_contrib.csv     per-PC: pc, beta, pos_contrib, neg_contrib,
#                           pos_frac, neg_frac (share of the set's total |w|)
#   pc_tissue_meta.json     shapes, topn, n bins matched

suppressMessages({library(data.table); library(jsonlite)})

getarg <- function(flag, default = NA) {
  a <- commandArgs(TRUE); i <- match(flag, a)
  if (is.na(i) || i == length(a)) return(default); a[i + 1]
}
model_path <- getarg("--model")
atlas_path <- getarg("--atlas",
  "analysis/fetal_fraction/reference/ff_openness_atlas_hg19_50kb.csv.gz")
topn       <- as.integer(getarg("--topn", "2000"))
outdir     <- getarg("--outdir", "pc_tissue_out")
prefix     <- getarg("--prefix", "C_")
stopifnot(!is.na(model_path))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- model ------------------------------------------------------------------
obj <- readRDS(model_path)
if (is.null(obj$bin_loadings) || is.null(obj$beta))
  stop("pc_tissue_map.R requires a PCA+lm model (fields $bin_loadings and $beta)")
L <- as.data.frame(obj$bin_loadings, stringsAsFactors = FALSE)
bin_name <- as.character(L[["bin_name"]])
pc_cols  <- grep("^PC[0-9]+$", names(L), value = TRUE)
Lmat <- as.matrix(L[, pc_cols, drop = FALSE]); rownames(Lmat) <- bin_name
beta <- obj$beta[pc_cols]
if (any(is.na(beta))) stop("beta missing PC coefficients present in bin_loadings")
# PC variance fraction (if pc_sdev available), else NA
var_frac <- rep(NA_real_, length(pc_cols))
if (!is.null(obj$pc_sdev)) {
  sd <- obj$pc_sdev
  if (length(sd) >= length(pc_cols)) {
    v <- sd[seq_along(pc_cols)]^2
    var_frac <- v / sum(sd^2)
  }
}
cat(sprintf("[model] %d bins x %d PCs\n", length(bin_name), length(pc_cols)))

# ---- reformat bin names to atlas keys: chr1_550000_600000 -> chr1:550000-600000
to_key <- function(f) {
  m <- regmatches(f, regexec("^(chr[0-9XYM]+)[:_.]([0-9]+)[:_.-]([0-9]+)$", f))[[1]]
  if (length(m) == 4) return(paste0(m[2], ":", m[3], "-", m[4]))
  NA_character_
}
keys <- vapply(bin_name, to_key, character(1))
ok <- !is.na(keys)
Lmat <- Lmat[ok, , drop = FALSE]; keys <- keys[ok]
rownames(Lmat) <- keys

# ---- atlas ------------------------------------------------------------------
# gz-robust read: fread needs R.utils for .gz, so decompress via a connection
# when needed (keeps the dependency surface at data.table + jsonlite).
read_atlas <- function(path) {
  if (grepl("\\.gz$", path)) {
    if (requireNamespace("R.utils", quietly = TRUE))
      return(fread(path, showProgress = FALSE))
    con <- gzfile(path, "rt"); on.exit(close(con))
    return(fread(text = readLines(con), showProgress = FALSE))
  }
  fread(path, showProgress = FALSE)
}
atlas <- read_atlas(atlas_path)
key_col <- if ("key" %in% names(atlas)) "key" else names(atlas)[4]
Ccols <- grep(paste0("^", prefix), names(atlas), value = TRUE)
# drop reference-contrast pseudo-tissues if present
Ccols <- Ccols[!grepl("(_specific|_tumor|_immune|background|signal)$", Ccols)]
tissues <- sub(paste0("^", prefix), "", Ccols)
atlas_keys <- as.character(atlas[[key_col]])
Cmat <- as.matrix(atlas[, ..Ccols]); rownames(Cmat) <- atlas_keys
colnames(Cmat) <- tissues

# ---- align on shared bins ---------------------------------------------------
common <- intersect(rownames(Lmat), rownames(Cmat))
if (length(common) < 10) stop(sprintf("only %d shared bins -- key format mismatch?", length(common)))
Lc <- Lmat[common, , drop = FALSE]
Cc <- Cmat[common, , drop = FALSE]
cat(sprintf("[align] %d shared bins; %d tissues\n", length(common), length(tissues)))

# ============================================================================
# 1. PC x tissue correlation (beta-oriented)
# ============================================================================
# raw Pearson corr of each PC loading vector vs each tissue openness vector
# pairwise.complete.obs: single-cell layers (e.g. Wang MFI W_) leave ~7% of grid
# bins NA where windows failed the hg38->hg19 lift; without this every corr collapses
# to NA->0. Bulk C_ atlases are fully covered, so this is a no-op for them.
corr <- suppressWarnings(cor(Lc, Cc, use = "pairwise.complete.obs"))   # npc x ntissue
corr[is.na(corr)] <- 0
# orient so + means "FF-up direction of this PC aligns with tissue openness"
corr_or <- corr * sign(beta)
best_idx <- max.col(abs(corr_or), ties.method = "first")
best_tissue <- tissues[best_idx]

corr_df <- data.frame(pc = pc_cols, beta = as.numeric(beta),
                      var_frac = var_frac, best_tissue = best_tissue,
                      best_corr = corr_or[cbind(seq_along(pc_cols), best_idx)],
                      stringsAsFactors = FALSE)
cor_block <- as.data.frame(corr_or); names(cor_block) <- paste0("corr_", tissues)
corr_df <- cbind(corr_df, cor_block)
fwrite(corr_df, file.path(outdir, "pc_tissue_corr.csv"))

# ============================================================================
# 2. per-PC contribution to top pos / neg SHAP bins
# ============================================================================
w <- as.numeric(Lmat %*% beta); names(w) <- rownames(Lmat)   # effective weight
ord_pos <- names(sort(w, decreasing = TRUE))
ord_neg <- names(sort(w, decreasing = FALSE))
pos_set <- head(ord_pos[w[ord_pos] > 0], topn)
neg_set <- head(ord_neg[w[ord_neg] < 0], topn)
# contribution matrix rows = bins, cols = PC:  L[bin,k]*beta_k
contrib_pos <- sweep(Lmat[pos_set, , drop = FALSE], 2, beta, "*")
contrib_neg <- sweep(Lmat[neg_set, , drop = FALSE], 2, beta, "*")
pos_c <- colMeans(contrib_pos)     # mean PC contribution across the pos set
neg_c <- colMeans(contrib_neg)
# share of each set's total signed weight carried by each PC
pos_frac <- colSums(contrib_pos) / sum(w[pos_set])
neg_frac <- colSums(contrib_neg) / sum(w[neg_set])
contrib_df <- data.frame(pc = pc_cols, beta = as.numeric(beta),
                         pos_contrib = as.numeric(pos_c),
                         neg_contrib = as.numeric(neg_c),
                         pos_frac = as.numeric(pos_frac),
                         neg_frac = as.numeric(neg_frac),
                         best_tissue = best_tissue,
                         stringsAsFactors = FALSE)
fwrite(contrib_df, file.path(outdir, "pc_shap_contrib.csv"))

writeLines(toJSON(list(n_pcs = length(pc_cols), n_tissues = length(tissues),
                       n_shared_bins = length(common), topn = topn,
                       n_pos = length(pos_set), n_neg = length(neg_set),
                       tissues = tissues), auto_unbox = TRUE, pretty = TRUE),
           file.path(outdir, "pc_tissue_meta.json"))

cat(sprintf("[write] %s/{pc_tissue_corr.csv, pc_shap_contrib.csv, pc_tissue_meta.json}\n", outdir))
cat(sprintf("[pos] top-%d FF-up bins: leading PC = %s (contrib %.3g)\n",
            length(pos_set), pc_cols[which.max(pos_c)], max(pos_c)))
cat(sprintf("[neg] top-%d FF-down bins: leading PC = %s (contrib %.3g)\n",
            length(neg_set), pc_cols[which.min(neg_c)], min(neg_c)))
