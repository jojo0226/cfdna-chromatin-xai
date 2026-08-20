#!/usr/bin/env Rscript
# pc_origin_map.R -- maternal vs fetal ORIGIN marking for the FF explainer, using the
# Wang et al. 2026 maternal-fetal interface (MFI) snATAC atlas (W_ layer).
#
# The bulk C_ atlas cannot tell fetal from maternal chromatin; the Wang MFI atlas can,
# because every cell type is origin-labelled. This script answers two questions the
# PC<->tissue map cannot:
#
#  1. PER-PC ORIGIN SCORE. For PC k with loading L[,k]:
#        origin_score_k = ( corr(L[,k], openness_fetal)
#                         - corr(L[,k], openness_maternal) ) * sign(beta_k)
#     where openness_fetal / openness_maternal are the clean-origin specificity
#     composites in the W_ layer. > 0 => this PC's FF-UP direction opens FETAL
#     chromatin (expected for a fetal-fraction signal); < 0 => maternal. A bin-shuffle
#     null (n_perm) gives z and empirical p per PC.
#
#  2. PER-BIN ORIGIN CALL + 2x2 test. For each top-N FF-up and FF-down bin (by the
#     effective weight w_bin = sum_k L[bin,k]*beta_k), assign the dominant Wang cell
#     type (argmax of row-centered specificity over the W_ columns) and read its origin
#     label. The prediction the whole tool rests on:
#        FF-up bins   -> predominantly FETAL cell types
#        FF-down bins -> predominantly MATERNAL cell types
#     Reported as a FF-direction x origin contingency with a Fisher exact test.
#
# Reuses the same model contract as pc_tissue_map.R (obj$bin_loadings + obj$beta) and
# the same key reformatting. Runs on the synthetic demo model now; on the real AWS
# PC-100 model when available (same deferral as pc_tissue_map).
#
# Usage:
#   Rscript pc_origin_map.R --model synth_pca_model.rds \
#       --atlas reference/wang_mfi_openness_hg19_50kb.csv.gz \
#       --manifest reference/manifest_wang_mfi.json \
#       --topn 2000 --n-perm 2000 --outdir _wang_demo
#
# Outputs (in --outdir):
#   pc_origin_score.csv   pc, beta, corr_fetal, corr_maternal, origin_score, z, p, call
#   bin_origin_calls.csv  key, direction (up/down), w, dominant_celltype, origin, purity
#   origin_contingency.csv  2x2 counts + Fisher p + odds ratio
#   pc_origin_meta.json   parameters + summary

suppressWarnings(suppressMessages({library(data.table); library(jsonlite)}))

getarg <- function(flag, default = NA) {
  a <- commandArgs(TRUE); i <- match(flag, a)
  if (is.na(i) || i == length(a)) return(default)
  a[i + 1]
}
model_path <- getarg("--model")
atlas_path <- getarg("--atlas", "reference/wang_mfi_openness_hg19_50kb.csv.gz")
manifest_p <- getarg("--manifest", "reference/manifest_wang_mfi.json")
topn       <- as.integer(getarg("--topn", "2000"))
n_perm     <- as.integer(getarg("--n-perm", "2000"))
seed       <- as.integer(getarg("--seed", "0"))
outdir     <- getarg("--outdir", "pc_origin_out")
stopifnot(!is.na(model_path))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
set.seed(seed)

# ---- model ------------------------------------------------------------------
obj <- readRDS(model_path)
if (is.null(obj$bin_loadings) || is.null(obj$beta))
  stop("pc_origin_map.R requires a PCA+lm model (fields $bin_loadings and $beta)")
L <- obj$bin_loadings
pc_cols <- grep("^PC[0-9]+$", colnames(L), value = TRUE)
bin_name <- if ("bin_name" %in% colnames(L)) as.character(L[, "bin_name"]) else rownames(L)
Lmat <- as.matrix(L[, pc_cols, drop = FALSE]); rownames(Lmat) <- bin_name
storage.mode(Lmat) <- "double"
beta <- obj$beta[pc_cols]
if (any(is.na(beta))) stop("beta missing PC coefficients present in bin_loadings")
cat(sprintf("[model] %d bins x %d PCs\n", nrow(Lmat), length(pc_cols)))

to_key <- function(f) {
  m <- regmatches(f, regexec("^(chr[0-9XYM]+)[:_.]([0-9]+)[:_.-]([0-9]+)$", f))[[1]]
  if (length(m) == 4) return(paste0(m[2], ":", m[3], "-", m[4]))
  NA_character_
}
keys <- vapply(bin_name, to_key, character(1))
ok <- !is.na(keys); Lmat <- Lmat[ok, , drop = FALSE]; rownames(Lmat) <- keys[ok]

# ---- atlas (W_ columns + origin composites) ---------------------------------
read_atlas <- function(path) {
  if (grepl("\\.gz$", path)) { con <- gzfile(path, "rt"); on.exit(close(con))
    return(fread(text = readLines(con), showProgress = FALSE)) }
  fread(path, showProgress = FALSE)
}
atlas <- read_atlas(atlas_path)
key_col <- if ("key" %in% names(atlas)) "key" else names(atlas)[4]
akeys <- as.character(atlas[[key_col]])
Wcols <- grep("^W_", names(atlas), value = TRUE)
Wmat <- as.matrix(atlas[, ..Wcols]); rownames(Wmat) <- akeys
colnames(Wmat) <- sub("^W_", "", Wcols)
fetal_v <- atlas[["openness_fetal"]]; mat_v <- atlas[["openness_maternal"]]
names(fetal_v) <- akeys; names(mat_v) <- akeys

manifest <- fromJSON(manifest_p)
origin_of <- vapply(colnames(Wmat), function(ct) manifest$celltypes[[ct]]$origin, character(1))
purity_of <- vapply(colnames(Wmat), function(ct) manifest$celltypes[[ct]]$origin_purity, numeric(1))

common <- intersect(rownames(Lmat), rownames(Wmat))
Lc <- Lmat[common, , drop = FALSE]
Wc <- Wmat[common, , drop = FALSE]
fc <- fetal_v[common]; mc <- mat_v[common]
cat(sprintf("[align] %d shared bins; %d cell types\n", length(common), ncol(Wc)))

# ---- 1. per-PC origin score + bin-shuffle null (vectorized) -----------------
# Restrict to bins where BOTH composites are present (unmapped grid bins are NA in
# fetal AND maternal), so every correlation runs on one common complete-bin subset.
# Then corr(L[,k], v) = (1/(n-1)) * <zscore(L[,k]), zscore(v)>, so a whole perm batch
# is two matrix multiplies (npc x n) %*% (n x n_perm) -- seconds, not minutes.
zscore <- function(x) { x <- x - mean(x); s <- sqrt(sum(x^2)); if (s == 0) x else x / s }
cb <- which(!is.na(fc) & !is.na(mc))
n <- length(cb)
ZL <- apply(Lc[cb, , drop = FALSE], 2, zscore)     # n x npc, unit-norm columns
zf <- zscore(fc[cb]); zm <- zscore(mc[cb])         # unit-norm composite vectors
cf <- as.numeric(crossprod(ZL, zf))                # npc: corr(L[,k], fetal)
cm <- as.numeric(crossprod(ZL, zm))                # npc: corr(L[,k], maternal)
score <- (cf - cm) * sign(beta)

# null: permute composite bin order; unit-norm is preserved under permutation
Pf <- matrix(0, n, n_perm); Pm <- matrix(0, n, n_perm)
for (b in seq_len(n_perm)) { Pf[, b] <- zf[sample.int(n)]; Pm[, b] <- zm[sample.int(n)] }
NF <- crossprod(ZL, Pf); NM <- crossprod(ZL, Pm)   # npc x n_perm
null <- t((NF - NM) * sign(beta))                  # n_perm x npc
mu <- colMeans(null); sdv <- apply(null, 2, sd) + 1e-12
z <- (score - mu) / sdv
# per-column centering: sweep mu/score across columns (matrix - vector would recycle
# the length-npc vector down each n_perm-length column and misalign the comparison)
null_dev <- abs(sweep(null, 2, mu, "-"))               # n_perm x npc
obs_dev  <- abs(score - mu)                            # npc
ge <- sweep(null_dev, 2, obs_dev, FUN = ">=")          # n_perm x npc logical
p_emp <- (colSums(ge) + 1) / (n_perm + 1)
call <- ifelse(score > 0, "fetal", "maternal")
pc_origin <- data.table(pc = pc_cols, beta = as.numeric(beta),
                        corr_fetal = cf, corr_maternal = cm,
                        origin_score = score, z = z, p = p_emp, call = call)
fwrite(pc_origin, file.path(outdir, "pc_origin_score.csv"))

# ---- 2. per-bin origin call + Fisher ----------------------------------------
w <- as.numeric(Lmat %*% beta); names(w) <- rownames(Lmat)
up  <- names(sort(w[w > 0], decreasing = TRUE))[seq_len(min(topn, sum(w > 0)))]
dn  <- names(sort(w[w < 0]))[seq_len(min(topn, sum(w < 0)))]
# dominant cell type = argmax row-centered specificity over W_ columns
spec <- Wc - rowMeans(Wc, na.rm = TRUE)
dom_idx <- max.col(replace(spec, is.na(spec), -Inf), ties.method = "first")
dom_ct <- colnames(Wc)[dom_idx]; names(dom_ct) <- rownames(Wc)

mk <- function(bins, dir) {
  bins <- intersect(bins, rownames(Wc))
  dt <- data.table(bin_key = bins, direction = dir, w = w[bins],
                   dominant_celltype = dom_ct[bins],
                   origin = origin_of[dom_ct[bins]],
                   purity = purity_of[dom_ct[bins]])
  setnames(dt, "bin_key", "key")   # rename after construction (key= is reserved in data.table())
  dt[]
}
bin_calls <- rbind(mk(up, "up"), mk(dn, "down"))
fwrite(bin_calls, file.path(outdir, "bin_origin_calls.csv"))

# 2x2 (direction x origin), restrict to Fetal/Maternal (drop Unknown)
bc <- bin_calls[origin %in% c("Fetal", "Maternal")]
tab <- table(factor(bc$direction, c("up", "down")),
             factor(bc$origin, c("Fetal", "Maternal")))
ft <- fisher.test(tab)
cont <- as.data.table(as.data.frame.matrix(tab), keep.rownames = "direction")
cont[, fisher_p := ft$p.value][, odds_ratio := unname(ft$estimate)]
fwrite(cont, file.path(outdir, "origin_contingency.csv"))

meta <- list(model = model_path, atlas = atlas_path, n_shared_bins = length(common),
             topn = topn, n_perm = n_perm, seed = seed,
             n_pc_fetal = sum(call == "fetal"), n_pc_maternal = sum(call == "maternal"),
             fisher_p = ft$p.value, odds_ratio = unname(ft$estimate),
             contingency = as.data.frame.matrix(tab))
write_json(meta, file.path(outdir, "pc_origin_meta.json"), auto_unbox = TRUE, pretty = TRUE)

cat(sprintf("[origin] up-bins fetal=%d maternal=%d | down-bins fetal=%d maternal=%d\n",
            tab["up","Fetal"], tab["up","Maternal"], tab["down","Fetal"], tab["down","Maternal"]))
cat(sprintf("[fisher] OR=%.3g  p=%.3g\n", unname(ft$estimate), ft$p.value))
top <- pc_origin[order(-abs(origin_score))][1:5]
cat("[top PCs by |origin_score|]\n"); print(top[, .(pc, beta, origin_score = round(origin_score,3), z = round(z,1), p, call)])
cat(sprintf("[write] %s/{pc_origin_score.csv, bin_origin_calls.csv, origin_contingency.csv, pc_origin_meta.json}\n", outdir))
