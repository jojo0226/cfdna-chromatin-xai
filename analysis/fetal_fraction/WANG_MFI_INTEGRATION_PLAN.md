# Integrating the Wang MFI scATAC atlas into the FF explainer
### PC↔cell-type correlation + maternal/fetal origin marking

**Atlas:** Wang et al. 2026, *Single-cell spatiotemporal dissection of the human maternal-fetal interface*, Nature 653:167-179 (doi:10.1038/s41586-026-10316-x). Source file: `wang_mfi_atlas_hg19.tar.gz` (~21 MB), hg19.

---

## 0. Why this atlas is the right one

The existing PC↔tissue layer (`pc_tissue_map.R`) correlates each PCA loading vector against a **bulk multi-tissue** openness atlas (73 `C_` columns). That atlas has no placenta cell-type resolution and no maternal/fetal axis — so it can say "PC3 looks placental" only because one bulk placenta track happens to be in it. The Wang atlas replaces that coarse signal with **33 maternal-fetal-interface cell types, each pre-labelled Fetal or Maternal**. This turns the mechanism step from "which tissue" into two sharper questions:

1. **Cell-of-origin at MFI resolution** — does a PC's FF-up pattern align with trophoblast (SCT/VCT/EVT), decidual stroma, maternal immune, or fetal endothelium?
2. **Origin marking** — is a PC (and each top-SHAP bin) driven by **fetal** or **maternal** chromatin? This is the direct biological test of the FF sign convention: FF-up bins *should* resolve to fetal/placental cell types, FF-down to maternal.

It slots into the existing `w_b = Σ_k L[b,k]·β_k` → PC↔tissue machinery with no change to the model or the SHAP core — it's a new **reference layer**, exactly like the placenta-ATAC layer added in v0.6.2.

---

## 1. What's in the box (verified)

| item | detail |
|---|---|
| `wang_mfi_win50kb_hg19.h5ad` | 33 cell types × 56,734 windows |
| `wang_mfi_win5kb_hg19.h5ad`  | 33 cell types × 567,239 windows |
| `celltype_manifest.csv` / `celltype_metadata.csv` | origin, origin_purity, nuclei counts, QC flags |
| `X` | raw fragment-end counts, uint32 (dense) |
| `var` | `chrom`, `start` (hg19); `end = start + uns['bin_size_bp']`; `hg38_start` = source window |

**Origin composition:** 17 Fetal, 15 Maternal, 1 Unknown (cDC).
**Clean fetal (purity ≥ 0.9):** SCT, SCT_A, SCT_B, proSCT, VCT, cycling_VCT, avVCT, EVT, iEVT, proEVT, eEVT, pEVT, GC.
**Clean maternal (purity ≥ 0.9):** Epi, uDSC, DSC_A, DSC_B, eS, CD16_M, LEC, cilated.
**Ambiguous (purity < 0.75, exclude from origin composites):** pvSMC 0.46, cDC 0.51, FB 0.55, Ery 0.61, B 0.66, dNK 0.67, T 0.69.
**Low-n (< 500 nuclei, downweight/flag):** aEC, eEVT, cDC, B, cilated, pEVT, Ery (Ery covers only 41% of 5 kb windows — drop).

---

## 2. The grid problem, and how we solve it

The FF model bins are a clean 50 kb hg19 grid (starts ≡ 0 mod 50000; 57,633 bins). The Wang 50 kb windows are **phase-offset** — only the start edge was lifted from hg38, so their starts land at 0, 1, 2, 8, 11 … mod 50000. Consequences (measured):

- **Exact (chrom,start) match: 420 / 57,633 = 0.7%.** Unusable as a key join.
- **Midpoint re-binning** (floor(window_midpoint / 50000)·50000) onto the FF grid: **53,506 / 57,633 = 92.8% of FF bins covered.** This is the join.

**Recommended build path — aggregate from the 5 kb matrix, not the 50 kb one.** The 50 kb matrix is an exact 10:1 per-chromosome sum of the 5 kb matrix (README), so nothing is lost by starting at 5 kb; and re-summing 5 kb windows onto the *FF* 50 kb grid by midpoint gives cleaner, phase-correct boundaries than re-binning the already-offset 50 kb windows. Each 5 kb window's midpoint → FF bin; sum counts per (cell type, FF bin). This yields a Wang matrix natively on the FF grid.

---

## 3. Normalization (counts → comparable openness)

Raw counts span 1000× in depth, so correlation against loadings would be dominated by depth, not accessibility. Pipeline per cell type, matching how the existing `C_` openness columns behave:

1. **Depth-normalize:** CPM = count / total_fragments(celltype) × 1e6 (per cell-type total, from metadata).
2. **Stabilize:** `log1p(CPM)`.
3. **Make cross-cell-type comparable:** per-bin, keep the value; for the specificity view (matching `ff_tissue_proportion.py`), `spec_ct = openness_ct − row-mean over cell types`.
4. **For the PC correlation:** z-score each cell-type column across bins so all 33 enter the Pearson correlation on equal footing (the same scale-free treatment the loadings get).

Output table `wang_mfi_openness_hg19_50kb.csv.gz` with columns `chrom,start,end,key, W_<celltype>×33` — deliberately prefixed `W_` (not `C_`) so it's a distinct, switchable atlas the existing scripts can consume via `--atlas`.

---

## 4. PC ↔ cell-type correlation (reuse existing math)

Feed the new atlas straight into `pc_tissue_map.R` (it already takes `--atlas` and `--prefix`):

```
Rscript pc_tissue_map.R --model <model.rds> \
  --atlas reference/wang_mfi_openness_hg19_50kb.csv.gz --prefix W_ --outdir <out>
```

Per PC *k*, on the shared bins:
- `align_{k,ct} = corr(L[·,k], openness_ct) × sign(β_k)` — β-oriented so **red = FF-up** across the whole matrix (identical convention to the current heatmap).
- `best_tissue` per PC → its dominant MFI cell type.
- Per-PC contribution `Contrib_{b,k} = L[b,k]·β_k` averaged over the FF-up / FF-down bin sets — which PCs build each direction, now mapped to MFI cell types.

No new math — this is the existing PC↔tissue decomposition pointed at a better reference.

---

## 5. Origin marking (the new axis) — maternal vs fetal

This is what the bulk atlas cannot do. Two levels:

**(a) Per-PC origin score.** Build two composite openness tracks from the *clean-origin* cell types only (purity ≥ 0.9; ambiguous and Unknown excluded):
- `openness_fetal`  = mean over 13 clean fetal cell types
- `openness_maternal` = mean over 8 clean maternal cell types

Then per PC:
- `origin_score_k = [corr(L[·,k], openness_fetal) − corr(L[·,k], openness_maternal)] × sign(β_k)`

Positive ⇒ this PC's **FF-up direction opens fetal chromatin** (expected for the fetal-fraction signal); negative ⇒ maternal. A permutation null (shuffle bins, same as the enrichment step's `n_perm`) gives a z and p per PC.

**(b) Per-bin origin call.** For each top-N FF-up and FF-down bin, assign the dominant Wang cell type (argmax specificity) and read its origin label. The prediction the whole tool rests on:
- **FF-up bins → predominantly Fetal cell types** (trophoblast/placental).
- **FF-down bins → predominantly Maternal** (decidua, maternal immune).

Report as a 2×2 (FF-direction × origin) contingency with a Fisher test, plus a stacked-bar "origin composition of top-N bins." This is a direct, external validation of the FF sign convention using an independent single-cell atlas.

---

## 6. Deliverables

**Reference layer**
- `reference/wang_mfi_openness_hg19_50kb.csv.gz` — 57,633 bins × 33 `W_` cell-type columns (grid-aligned, normalized).
- `reference/wang_mfi_origin_composites.csv.gz` — `openness_fetal`, `openness_maternal` per bin.
- `reference/manifest_wang_mfi.json` — provenance, cell-type→origin map, purity, QC flags, grid-mapping stats (92.8% coverage), normalization recipe.

**Scripts (new / extended)**
- `build_wang_mfi_atlas.py` — extract h5ad → 5 kb → FF-grid re-bin → normalize → write layer + composites + manifest.
- `pc_origin_map.R` (or a `--origin` flag on `pc_tissue_map.R`) — per-PC origin score + null; per-bin origin call + Fisher.
- Reuse `pc_tissue_map.R`, `run_ff_tissue_track.py`, `plot_pc_tissue.py`, `plot_ff_tissue_dumbbell.py` unchanged via `--atlas`/`--prefix W_`.

**Figures**
- PC × Wang-cell-type heatmap (β-oriented, rows = PC·best-MFI-cell-type).
- Per-PC origin-score bar (fetal red / maternal blue), with null band.
- Origin composition of FF-up vs FF-down bins (stacked bar) + 2×2 Fisher panel.
- Dumbbell: FF-up vs FF-down cell-type contrast at MFI resolution.

**Docs**
- `WANG_MFI_LAYER.md` — what/why/how, caveats, interpretation guide (mirrors `PLACENTA_ATAC_LAYER.md`).

---

## 7. Validation & caveats (build into the manifest and doc)

- **Grid coverage:** 92.8% of FF bins mapped; the 7.2% unmapped (failed lift / phase gaps) are flagged NA, not zero — they must be dropped from correlations, not counted as closed.
- **Depth confounding:** verify post-normalization that per-PC alignment does **not** correlate with cell-type sequencing depth (guards against a depth artifact masquerading as origin signal).
- **Origin composites use clean types only** — 7 ambiguous types (incl. dNK, FB, pvSMC) are excluded from fetal/maternal means and reported separately; including them would blur the very axis we're testing.
- **Low-n types flagged** (Ery dropped; aEC/eEVT/cDC/B/cilated/pEVT downweighted or annotated).
- **Resolution honesty:** correlations are at 50 kb because the loadings are 50 kb; note the README's dilution caveat and treat cell-type *rankings* as more reliable than absolute correlation magnitudes.
- **Positive control:** the trophoblast lineage (SCT/VCT/EVT) should top the FF-up alignment and the fetal origin score; if it doesn't, the grid mapping or normalization is wrong.
- **Cross-check** against the existing placenta-ATAC layer (v0.6.2): the Wang fetal composite should track the placenta-ATAC track (concordance sanity check).

---

## 8. Sequencing of the work

1. `build_wang_mfi_atlas.py`: h5ad → FF-grid layer + origin composites + manifest. *(quantitative — I run it)*
2. Run `pc_tissue_map.R --atlas W_` on the synthetic demo model → PC × MFI-cell-type heatmap. *(the real AWS-encrypted model is deferred as before; synthetic first)*
3. `pc_origin_map.R`: per-PC origin score + per-bin origin Fisher → origin figures.
4. Depth-confound + placenta-ATAC concordance checks.
5. `WANG_MFI_LAYER.md` + fold layer, scripts, figures into a **v0.6.3** bundle.

Steps 1–4 are runnable now on the synthetic model; re-running on the real PC-100 model happens whenever the AWS model is available (same deferral as the current pc_tissue_map real-model item).
