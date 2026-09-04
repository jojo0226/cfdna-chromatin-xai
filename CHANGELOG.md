# Changelog

All notable changes to **cfdna-chromatin-xai** are documented here.
The format loosely follows [Keep a Changelog](https://keepachangelog.com/).

## [0.7.0] - 2026-09-04
### Added
- **New packaged `cfdna_chromatin.chromhmm` module** — a build-aware, user-facing
  façade over the internal `references`/`engine`/`liftover` layers that annotates a
  region's importance table (SHAP or differential) against the bundled ChromHMM
  18-state reference epigenomes and reports the state/group composition of the
  top-|importance| bins. Handles hg19 inputs three ways: (1) per-analysis auto-lift
  hg19→hg38, (2) pre-lift the references hg38→hg19 for a whole cohort
  (`scripts/prelift_chromhmm_hg19.py`), or (3) drop in native Roadmap/full-stack
  hg19 tracks. Key functions: `annotate_importance_bins()` (one-call driver),
  `state_composition()`, `composition_matrix()`, `multi_cutoff_composition()`,
  `plot_composition_heatmaps()`.
- **New `scripts/run_chromhmm.py`** — CLI over the module. Auto-detects key/importance/
  direction columns, ranks by |importance|, and writes per-tissue composition CSVs +
  annotated bins + `chromhmm_meta.json`. Flags: `--from-build {hg38,hg19}`, `--chain`
  (forward UCSC chain), `--back-chain` (reverse chain — enables offline round-trip QC),
  `--top-ns 1000,2000,5000` + `--plot` for a multi-cutoff composition heatmap, and
  `--fasta-dir` for matched-null state enrichment.
- **New `scripts/run_fullstack_chromhmm.py`** — cell-type-agnostic alternative that
  annotates against the Ernst-lab **universal (full-stack) 100-state** ChromHMM model
  (native hg19 and hg38, no per-cell-type choice), collapsing the 100 states into ~16
  functional groups. Answers "what functional chromatin are the top regions in?"
  without committing to any single reference epigenome.
- **Liftover round-trip fix** — `liftover.lift_regions()` gains a `back_chain_path`
  parameter so reverse-direction round-trip QC uses a local chain instead of forcing a
  network auto-download (which fails offline). Threaded through `prepare_query()` /
  `annotate_importance_bins()` and the `--back-chain` CLI flag.

## [0.6.2] - 2026-08-18
### Added
- **New `run_ff_tissue_track.py` lollipop layout** (`--style lollipop`,
  `--top-k`) — horizontal dot-and-stem cell-of-origin plot over the top-K
  tissues (default 20, placenta always kept) sorted by |spec_z|. Tissue labels
  sit on the y-axis so they never rotate, one colored dot per top-N cutoff makes
  cutoff stability legible at a glance, and the near-zero tail is dropped with a
  "+N more near zero" caption instead of crowding the axis. Purely additive —
  the grouped-`bars` layout remains the default, so existing invocations are
  unchanged.
- **New `analysis/fetal_fraction/plot_ff_tissue_dumbbell.py`** — POS-vs-NEG
  cell-of-origin dumbbell. Reads only the `ff_tissue_results.csv` the track
  already writes for the positive-SHAP (FF-raising) and negative-SHAP
  (FF-lowering) bin sets — no recomputation, no atlas needed. One row per tissue,
  two dots (POS `spec_z` vs NEG `spec_z`) joined by a line whose length *is* the
  fetal↔maternal contrast; rows sorted by POS−NEG so fetal tissues
  (placenta/trophoblast) collect at the top and maternal hematopoietic tissues
  (monocyte/neutrophil/CMP) at the bottom. `--topn`, `--top-k`, `--metric`.
- **Raw 73-tissue atlas source BEDs now shipped** — `data/access_ff/` (72
  public-ENCODE hg19 DNase-seq narrowPeaks in the 4-column
  `chrom,start,end,signalValue` form + `encode_accessions.json` tissue→accession
  manifest, 164 replicate files) is now un-gitignored and included, and the 73rd
  track ships under `data/neutrophil_dnase/`. The curated73 atlas is now fully
  reconstructable offline from the shipped inputs, not just re-downloadable. The
  bulky `_raw` staging dirs and the unused `histone_ff/` tracks remain excluded
  (still re-fetchable via `scripts/fetch_ff_atlas_beds.py`).

### Changed
- **Tissue-track label legibility for the 73-tissue atlas** — bar-layout tissue
  labels are now tilted (55° for >16 tissues, 45° otherwise) instead of vertical,
  with an enlarged size-adaptive font (12–14 pt), wider per-tissue spacing, and
  deeper bottom margins so names no longer collapse to ~7 pt at 73 tissues. New
  `--tick-fs` override and `--layout {auto,row,col,separate}` control.

### Fixed
- **Duplicate top-N cutoffs collapsed** — when fewer bins are available than a
  requested cutoff, that cutoff now collapses onto the previous identical one and
  is skipped (with a printed reason) instead of drawing two identical bar/dot
  groups.

### Added (earlier this release)
- **Placental snATAC openness layer** — `analysis/fetal_fraction/reference/
  ff_placenta_atac_hg19_50kb.csv.gz` (+ `manifest_placenta_atac.json`,
  `PLACENTA_ATAC_LAYER.md`), built by
  `analysis/fetal_fraction/build_placenta_atac_atlas.py` from the Wang et al. 2024
  pregnancy snATAC atlas (Nat Genet, GSE247036; hg38→hg19 lifted, 283,406 peaks).
  57,633-bin hg19 grid, six cell-type columns (`A_STB/A_vCTB/A_EVT/A_Endothelial/
  A_Erythroid/A_Fibroblast`, z across bins) + `placenta_atac_spec`. Kept as a
  **separate layer** — never merged into the DNase `C_` columns — for cross-modal
  validation of the placental signal (Spearman 0.91–0.93 vs `C_placenta`) and for
  PC alignment via the existing `pc_tissue_map.R --prefix A_` (no code change to the
  PC tool; it is atlas-agnostic).
- **`ATLAS_BUILD_LOG.md` + `ff_atlas_accessions.csv`** in
  `analysis/fetal_fraction/reference/` — full provenance trail for the 73-tissue
  atlas (per-tissue ENCODE accessions, biosample classifications, rep-selection
  ranking, 50 kb binning) and the 164-row machine-readable accession table.

## [0.6.1] - 2026-08-17
### Added
- **`scripts/fetch_ff_atlas_beds.py`** — one-command, public-ENCODE-only driver
  that reproduces the 73-tissue atlas's source BEDs. It re-runs the same hg19
  DNase-seq narrowPeak query the atlas was built from, applies the same
  ≤3-replicate-per-biosample ranking and the same curation (classification==
  `tissue` + the fixed ADDITIONS list of immune / myeloid-progenitor / cfDNA-
  relevant biosamples), downloads the narrowPeaks, converts them to the shipped
  4-column form (`chrom,start,end,signalValue`, hg19 autosomes) into
  `data/access_ff/`, and writes `encode_accessions.json` (raw biosample name →
  file accessions, joinable to the atlas `C_<tissue>` columns). `--dry-run`,
  `--resume`, `--out` supported. **No S3, no AWS, no credentials.** The
  72-biosample ENCODE panel is recovered by re-querying the public portal; the
  73rd track (neutrophil) has no hg19 DNase-seq on ENCODE and is supplied from
  local BEDs (reported as a known gap by the script).
- **`scripts/restore_ff_atlas_beds.py` — DEPRECATED tombstone.** The prior
  release shipped this as an "S3-restore" driver on the mistaken belief that the
  atlas accessions survived only in a private archived-S3 index. That was wrong:
  the atlas is built entirely from public ENCODE. The file now exits with a
  pointer to `fetch_ff_atlas_beds.py` and keeps the corrected provenance on
  record.

## [0.6.0] - 2026-08-17
### Added
- **`--marks` / `--composite` on `build_openness_atlas.py`** — the per-tissue
  combined-openness track `C_<tissue>` is now configurable instead of always
  averaging all four active histone marks + DNase. `--marks` selects which
  histone marks feed `H_<tissue>` (e.g. `H3K4me3` alone, or
  `H3K4me3,H3K27ac`); `--composite {dnase,histone,histone+dnase}` picks whether
  `C_` is built from accessibility only, selected histone marks only, or both
  (default `histone+dnase`). A tissue missing the requested modality falls back
  to whatever it has, so every tissue still gets a `C_` column on the same
  z-scale. Same two params are on the `build_openness_atlas()` API.

### Changed
- **Atlas renamed `curated72` → `curated73`** everywhere (files
  `ff_openness_atlas_hg19_50kb_curated73.csv.gz`, `manifest_curated73.json`, and
  all code/doc references). The atlas genuinely has 73 tissue columns (72 curated
  + real neutrophil accessibility added in 0.5.7); the `72` name was stale.
- **`plot_pc_tissue.py` panel A x-axis** — width now scales with tissue count and
  tick labels rotate 90° with per-count font sizing, so all 73 tissue names are
  legible and non-overlapping.

## [0.5.7] - 2026-08-10
### Added
- **Real neutrophil accessibility column in both openness atlases.** ENCODE has
  no hg19 neutrophil DNase-seq, so the 11-tissue atlas previously carried
  neutrophil as histone-only and the curated-72 atlas substituted the CD34+
  common-myeloid-progenitor (CMP) as a proxy. Three user-supplied hg19
  neutrophil accessibility peak sets (fetal-lung neutrophil + two GMP/bone-marrow
  neutrophil, 4-col BED chrom/start/end/signalValue, ~146k peaks each) replace
  those proxies with real data:
  - **curated-72 -> curated-73:** new `C_neutrophil` (signalValue-weighted zmean
    coverage), 73rd tissue.
  - **11-tissue:** new `D_neutrophil` / `n_D_neutrophil` / `ACC_neutrophil`, and
    `C_neutrophil` upgraded from histone-only (`z(H)`) to the full
    `mean(z(H), z(ACC))` composite used by every other tissue.
- **`add_neutrophil_dnase.py`** reproducible integration script (bins the three
  BEDs onto the 50 kb grid with the atlas's exact signalValue-weighted covered-
  fraction + per-replicate z-score + zmean aggregation, recomputes the composite
  columns, and updates both manifests). Source BEDs shipped under
  `data/neutrophil_dnase/`.
- **Manifests updated:** `manifest_curated73.json` `neutrophil_note` marked
  RESOLVED (real column now present; CMP retained as the distinct myeloid-
  progenitor axis); `ff_reference_manifest.json` neutrophil tissue-note records
  the real accessibility provenance.
### Validation
- Lineage sanity check on the shipped column (pan-subtracted specificity): most
  similar to **monocyte (rho=0.71)** and **CMP (rho=0.66)** (both myeloid),
  then lymphoid (B/CD8/CD4/NK rho~0.49-0.53), and **most distinct from placenta
  (0.25) / trophoblast (0.14)**; brain ~0. Correct granulocyte lineage signature
  (`analysis/fetal_fraction/neutrophil_atlas_validation.png`).

## [0.5.6] - 2026-08-04
### Fixed
- **`run_ff_tissue_track.py` x-axis labels no longer collide on large atlases.**
  The `_plot()` panel width now scales with the number of tissues
  (~0.34 in/tissue/panel), rotates labels to vertical past 16 tissues, and
  shrinks tick font past 24/40 tissues -- so the curated-72 atlas renders
  legibly instead of overprinting 72 names in a fixed 12.5-in canvas.
### Added
- **`--max-tissues N` CLI flag** on `run_ff_tissue_track.py`: cap the x-axis to
  the N most informative tissues (always including placenta, then largest
  `|spec_z|` across cutoffs). Use with the 73-tissue atlas to keep a compact,
  presentation-ready figure, e.g. `--max-tissues 20`. Default: show all.
- **`--layout {auto,row,col,separate}` CLI flag** on `run_ff_tissue_track.py`
  controlling how the cell-of-origin (A) and proportion (B) panels are arranged:
  `row` = A|B side by side; `col` = A stacked over B, **each panel labeled with
  its own tissue-name row** so A is readable without tracing down to B (half the
  width of `row`; recommended for the 73-tissue atlas); `separate` = two
  standalone PNGs (`<stem>_A.png` / `<stem>_B.png`); `auto` = row for <=16
  tissues else col. Default: auto.
- **Default top-N cutoffs now `500,1000,2000,5000`** (was `500,1000,2000`) for
  both `run_ff_tissue_track.py --topns` and `ff_tissue_proportion.DEF_TOPNS`,
  matching the four signed-SHAP selection tiers in `select_signed_shap.py`, so
  the top-5000 series is drawn by default instead of silently missing.
- **Collapsed-cutoff guard:** a cutoff larger than the number of available
  importance bins (e.g. top-5000 with only 1,500 bins) previously produced a
  bar group identical to the largest real cutoff. Such duplicate cutoffs are
  now skipped with an explanatory `[skip]` message instead of double-drawing.

## [0.5.5] - 2026-07-28
### Added
- **Repressive (H3K27me3) axis as a separate cell-of-origin test (Step 5).**
  `build_repressive_atlas.py` builds a 57,633-bin × 39-tissue hg19 atlas of
  z-mean H3K27me3 coverage (`R_<tissue>`) plus pan-tissue-subtracted
  specificity (`repr_spec_<tissue>`) from ENCODE narrowPeak, covering every
  cfDNA-critical tissue (placenta ×3, trophoblast, chorion, monocyte,
  B/CD4/CD8/NK, liver, CMP, HUVEC). It is a **fully separate axis** — never
  merged into the `C_` openness columns — so a genuine origin can be required
  to be active-enriched **and** repressive-depleted.
- **`plot_repressive.py`** emits two separate figures: (1) tissue clustering on
  the H3K27me3 specificity profile (hematopoietic and placenta blocks recover
  correctly), and (2) the two-sided test — per-tissue Spearman ρ between
  active-specificity and repressive-specificity, where ρ<0 validates the test.
  cfDNA-critical tissues are the most anti-correlated (placenta ρ≈−0.16,
  CD8/CD4/CMP ρ≈−0.17…−0.21). Ships `repr_twosided_rho.csv`.
- FF guide gains **Step 5** documenting the axis, the commands, and why the
  comparison is done on *specificity* (raw z is confounded by gene density at
  50 kb, ρ≈+0.2; pan-subtraction reveals the expected anti-correlation).

## [0.5.4] - 2026-07-28
### Added
- **PC ↔ tissue interpretation step (Step 4).** `pc_tissue_map.R` opens the PCA
  black box: it correlates each principal component's genome-wide loading
  pattern `L[:,k]` against every tissue's open chromatin (β-oriented, so the
  sign tracks the FF-raising direction) and decomposes the exact per-bin weight
  `w_bin = Σ_k L[bin,k]·β_k` into per-PC contributions over the top positive
  (FF-up) and negative (FF-down) SHAP bin sets. Writes `pc_tissue_corr.csv`,
  `pc_shap_contrib.csv`, `pc_tissue_meta.json`.
- **`plot_pc_tissue.py`.** Renders the two-panel figure — (A) PC × tissue
  alignment heatmap, (B) per-PC contribution to each SHAP direction — ranking
  PCs by total signed contribution. Documented as Step 4 in FF_IMPLEMENTATION.md.
  Reads only the model object and the public atlas; no patient data.

## [0.5.3] - 2026-07-28
### Added
- **Signed SHAP direction.** `extract_ff_shap.R` now emits a `signed_shap`
  column alongside `mean_abs_shap`, carrying the sign of the effective per-bin
  weight `w_bin` (`>0` raises predicted FF = fetal/placental-leaning; `<0`
  lowers it = maternal-leaning). Threaded through both the PCA+lm and the
  glmnet fallback paths; duplicate keys collapse to the largest-magnitude
  signed value. The magnitude column is unchanged, so the tissue track stays
  backward-compatible.
- **`select_signed_shap.py`.** Splits the importance table into most-positive
  and most-negative bin sets per cutoff (default 5000/2000/1000/500), writing
  `pos_top{N}.csv` / `neg_top{N}.csv` and a `signed_selection_manifest.json`.
  Each directed set feeds Step 2 separately so the cell-of-origin readout is
  resolved per SHAP direction. Documented as Step 1b in FF_IMPLEMENTATION.md.

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
