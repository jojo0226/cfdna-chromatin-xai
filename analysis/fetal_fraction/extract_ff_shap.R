#!/usr/bin/env Rscript
# extract_ff_shap.R -- per-bin importance from a fetal-fraction model, emitted as
# the bin-keyed importance table that run_ff_tissue_track.py consumes.
#
# RUNS IN PLACE on the AWS instance. Reads the (encrypted-at-rest) model .rds and,
# optionally, the per-sample bin matrix that was fed into the model's PCA. The
# OUTPUT is a per-BIN aggregate (an effective weight, optionally x a column
# statistic) -- NOT per-patient data -- so it is safe to copy off the instance.
#
# ---------------------------------------------------------------------------
# TWO MODEL SHAPES ARE SUPPORTED (auto-detected):
#
# (A) PCA + linear regression  [the NIPT/seqFF model: a list with
#     $bin_loadings, $beta, $pc_sdev, $bin_params, $train_dt, ...]
#
#     The pipeline is two stacked LINEAR maps, so it collapses EXACTLY into a
#     single linear model in bin space:
#         score_k(sample) = sum_bin  loadings[bin,k] * (x_bin - xbar_bin)
#         FF(sample)       = b0 + sum_k beta_k * score_k
#                          = b0 + sum_bin  w_bin * (x_bin - xbar_bin)
#     where the EFFECTIVE PER-BIN WEIGHT is
#         w_bin = sum_k  loadings[bin,k] * beta_k       (matrix-vector product).
#
#     SHAP is then closed form in bin space (NOT on the PCs -- a PC is a
#     whole-genome linear combination and has no location):
#         phi_bin(sample) = w_bin * (x_bin - xbar_bin)
#     Global importance   = mean_i |phi_bin| = |w_bin| * mean_i |x_bin - xbar_bin|.
#     With no bin matrix, |w_bin| alone IS the model-level importance track and
#     needs nothing but the model object.
#
# (B) glmnet / cv.glmnet directly in bin space  [fallback]
#         phi_ij = beta_j * (x_ij - xbar_j);  importance = |beta_j| * mean|x-xbar|.
#
# Usage:
#   Rscript extract_ff_shap.R \
#       --model   ff_model.rds \
#       --x       bin_matrix.csv   # optional; rows=samples, cols=bins (the PCA input)
#       --lambda  min              # glmnet only: "min" | "1se" | numeric
#       --out     ff_shap_importance.csv
#
# --x is the per-sample bin coverage matrix you fed into the PCA (same bins as
# $bin_loadings$bin_name). It stays on the instance; only the aggregate leaves.
# If omitted, importance = |w_bin| (model-only, exact up to the per-bin scale).

suppressWarnings(suppressMessages({
  ok_glmnet <- requireNamespace("glmnet", quietly = TRUE)
  if (ok_glmnet) library(glmnet)
}))

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NA) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[i + 1]
}
model_path <- getarg("--model")
x_path     <- getarg("--x", NA)
lambda_sel <- getarg("--lambda", "min")
out_path   <- getarg("--out", "ff_shap_importance.csv")
stopifnot(!is.na(model_path))

obj <- readRDS(model_path)

# ===========================================================================
# Detect the model shape.
# ===========================================================================
is_pca_lm <- is.list(obj) && !is.null(obj$bin_loadings) && !is.null(obj$beta)

# result: a data.frame with columns feature (bin name) + weight (effective coef)
if (is_pca_lm) {
  cat("[detect] PCA + linear-regression model\n")
  L <- obj$bin_loadings                         # data.table: bin_name + PC1..PCk
  bin_name <- as.character(L[["bin_name"]])
  pc_cols  <- grep("^PC[0-9]+$", names(L), value = TRUE)
  Lmat <- as.matrix(as.data.frame(L)[, pc_cols, drop = FALSE])  # nbin x npc

  beta <- obj$beta                              # named: (Intercept), PC1..PCk
  bpc  <- beta[pc_cols]                          # align to the loadings' PCs
  if (any(is.na(bpc))) stop("beta is missing some PC coefficients present in bin_loadings")

  # effective per-bin weight: w = L %*% beta_pc
  w <- as.numeric(Lmat %*% bpc)
  names(w) <- bin_name
  cat(sprintf("[compose] %d bins x %d PCs -> per-bin weights; |w| range [%.3g, %.3g]\n",
              length(bin_name), length(pc_cols), min(abs(w)), max(abs(w))))

  imp_mode <- "pca_weight_only"
  imp_vec  <- abs(w)

  if (!is.na(x_path) && file.exists(x_path)) {
    X <- as.matrix(read.csv(x_path, row.names = 1, check.names = FALSE))
    common <- intersect(bin_name, colnames(X))
    cat(sprintf("[X] %d samples x %d bins; %d/%d model bins present\n",
                nrow(X), ncol(X), length(common), length(bin_name)))
    xbar <- colMeans(X[, common, drop = FALSE])
    center <- sweep(X[, common, drop = FALSE], 2, xbar, "-")
    mad_j  <- colMeans(abs(center))             # mean_i |x_ij - xbar_j|
    imp_vec <- setNames(rep(0, length(w)), bin_name)
    imp_vec[common] <- abs(w[common]) * mad_j   # exact global linear-SHAP
    imp_mode <- "pca_linear_shap_exact"
  }
  imp <- data.frame(feature = names(imp_vec),
                    mean_abs_shap = as.numeric(imp_vec),
                    stringsAsFactors = FALSE)

} else {
  cat("[detect] glmnet / linear model in bin space (fallback path)\n")
  if (inherits(obj, c("glmnet", "cv.glmnet"))) {
    fit <- obj; feats <- NULL
  } else {
    fit <- obj$model
    feats <- obj$features_used; if (is.null(feats)) feats <- obj$predictors
  }
  lam <- NULL
  if (lambda_sel == "min" && !is.null(fit$lambda.min)) lam <- fit$lambda.min
  else if (lambda_sel == "1se" && !is.null(fit$lambda.1se)) lam <- fit$lambda.1se
  else if (!lambda_sel %in% c("min", "1se")) lam <- as.numeric(lambda_sel)
  co <- as.matrix(coef(fit, s = lam))
  beta <- co[-1, 1, drop = TRUE]
  feat_names <- if (!is.null(feats)) feats else rownames(co)[-1]
  if (length(feat_names) != length(beta)) feat_names <- names(beta)
  names(beta) <- feat_names
  nz <- beta[beta != 0]
  cat(sprintf("[glmnet] %d features, %d non-zero (lambda=%s)\n",
              length(beta), length(nz), format(lam)))
  imp_mode <- "coef_only"; imp_vec <- abs(nz); nm <- names(nz)
  if (!is.na(x_path) && file.exists(x_path)) {
    X <- as.matrix(read.csv(x_path, row.names = 1, check.names = FALSE))
    common <- intersect(names(nz), colnames(X))
    xbar <- colMeans(X[, common, drop = FALSE])
    mad_j <- colMeans(abs(sweep(X[, common, drop = FALSE], 2, xbar, "-")))
    imp_vec <- abs(nz[common]) * mad_j; nm <- common
    imp_mode <- "linear_shap_exact"
  }
  imp <- data.frame(feature = nm, mean_abs_shap = as.numeric(imp_vec),
                    stringsAsFactors = FALSE)
}
cat(sprintf("[importance] mode=%s, %d features scored\n", imp_mode, nrow(imp)))

# ---- map feature names -> hg19 50kb bin keys (chrN:start-end) ----------------
# Coverage/PCA bins carry a genomic location; any non-genomic feature is dropped.
to_key <- function(f) {
  # chr1:50000-100000 | chr1:50000:100000 | chr1_50000_100000 | chr1.50000.100000
  m <- regmatches(f, regexec("^(chr[0-9XYM]+)[:_.]([0-9]+)[:_.-]([0-9]+)$", f))[[1]]
  if (length(m) == 4) return(paste0(m[2], ":", m[3], "-", m[4]))
  NA_character_
}
imp$key <- vapply(imp$feature, to_key, character(1))
mapped <- imp[!is.na(imp$key), c("key", "mean_abs_shap")]
mapped <- aggregate(mean_abs_shap ~ key, data = mapped, FUN = max)
cat(sprintf("[map] %d/%d features -> genomic bins (%d dropped: non-genomic)\n",
            nrow(mapped), nrow(imp), nrow(imp) - nrow(mapped)))

write.csv(mapped, out_path, row.names = FALSE)
cat(sprintf("[write] %s  (feed to run_ff_tissue_track.py --importance)\n", out_path))
