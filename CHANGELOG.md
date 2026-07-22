# Changelog

All notable changes to **cfdna-chromatin-xai** are documented here.
The format loosely follows [Keep a Changelog](https://keepachangelog.com/).

## [0.5.2] - 2026-07-22
### Fixed
- **FF SHAP extractor now handles the real NIPT model shape (PCA + linear
  regression), not just glmnet.** The production fetal-fraction model is a
  two-stage linear pipeline: a PCA over ~52k genomic 50 kb bins followed by an
  `lm` on the top 100 PCs. Because both stages are linear, they collapse
  *exactly* into a single linear model in bin space, with per-bin effective
  weight `w_bin = Σ_k loadings[bin,k] · beta_k`. `extract_ff_shap.R` now
  auto-detects this shape (`$bin_loadings` + `$beta`) and composes the bin-space
  weights directly; the previous glmnet path is kept as a fallback.
- SHAP is computed in **bin space**, never on PCs — a principal component has no
  genomic coordinate, so per-PC attributions cannot be mapped to the chromatin
  atlas. Verified: reconstructed prediction and Shapley sum agree to ~1e-14, and
  `|w_bin|` vs per-sample `mean|phi_bin|` correlate at r = 0.997.

### Changed
- `FF_IMPLEMENTATION.md` Step 1 rewritten for the PCA+lm model, with the
  weight-composition derivation, the data-egress table (model `.rds` and bin
  matrix stay on the instance; only the per-bin aggregate CSV leaves), and a
  validation section.

## [0.5.1] - 2026-07-22
### Added
- 11-tissue chromatin openness atlas (hg19, 50 kb bins): placenta, liver,
  endothelial, monocyte, B cell, CD4, CD8, NK, neutrophil, CMP (CD34+ common
  myeloid progenitor), K562 control.

## [0.5.0] - 2026-07-21
### Added
- Fetal-fraction (seqFF++) analysis track: per-bin importance → tissue-of-origin
  enrichment via openness-matched null model.
