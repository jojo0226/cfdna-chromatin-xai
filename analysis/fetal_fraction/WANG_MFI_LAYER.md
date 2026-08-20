# Wang MFI snATAC layer (`W_`) — maternal-vs-fetal origin marking

The `W_` layer adds a **maternal-fetal interface (MFI) single-nucleus ATAC**
reference to the FF explainer, on top of the bulk DNase openness atlas (`C_`) and
the Wang-2024 placental snATAC layer (`A_`). Its purpose is the one thing the bulk
atlas cannot do: **tell fetal chromatin from maternal chromatin**, because every
cell type in it carries an origin label.

## Source

Wang et al. **Nature 653:167-179 (2026)**, doi:10.1038/s41586-026-10316-x
(cell.ucsf.edu/snPlacenta). Pseudobulk snATAC over **33 MFI cell types**, 23
donors, 3 gestational stages, 4.39 B fragments. Delivered as
`wang_mfi_atlas_hg19.tar.gz` (win5kb + win50kb `.h5ad`, metadata, README).

The 33 cell types span the fetal placenta (syncytio-/villous-/extravillous
trophoblast: SCT, VCT, EVT and subtypes), fetal stroma/endothelium/immune (FB,
fEC, HB, GC), and the maternal decidua/immune compartment (DSC, uDSC, eS, Epi,
LEC, CD14/CD16 monocytes, T, dNK, B).

## Build (`build_wang_mfi_atlas.py`)

```
python build_wang_mfi_atlas.py \
    --tar  ~/Downloads/wang_mfi_atlas_hg19.tar.gz \
    --grid reference/ff_openness_atlas_hg19_50kb_curated73.csv.gz \
    --out  reference/wang_mfi_openness_hg19_50kb.csv.gz
```

**Grid join (the key subtlety).** The Wang windows are hg19 only on their *start*
edge (lifted from hg38), so their starts are phase-offset from the clean 50 kb FF
grid. Exact `(chrom,start)` matching recovers **0.7 %** of bins; re-binning by
window **midpoint** — `floor((start+2500)/50000)*50000` — recovers **93.2 %**
(53,727 / 57,633). We aggregate from the **5 kb** matrix (the 50 kb matrix is an
exact 10:1 sum of it, so nothing is lost) onto the FF grid by midpoint. Unmapped
bins stay `NaN`, never 0.

**Normalization** (raw fragment-end counts → comparable openness), per cell type:
CPM (÷ that cell type's `total_fragments`) → `log1p` → **z-score across bins**.
The `W_<celltype>` columns are left z-scored, matching the `C_`/`A_` layers so
Pearson correlations and PC alignment are on the same footing.

**Origin composites.** Two extra columns, `openness_fetal` and
`openness_maternal`, are built from **clean-origin** cell types only
(`origin_purity ≥ 0.90`, Unknown excluded → 13 fetal, 8 maternal). Critically
they are built from **row-centered specificity** (`z − per-bin mean across cell
types`), not raw z-score. At 50 kb every cell type shares a large
"open-everywhere" component (the README's dilution caveat); raw-z composites
correlate **+0.91** and cannot separate origin, whereas specificity composites
correlate **−0.83** and form a genuine fetal-vs-maternal axis (trophoblast
+0.6..0.76 fetal / −0.5..−0.72 maternal).

Output: `reference/wang_mfi_openness_hg19_50kb.csv.gz` (57,633 bins × 39 cols:
chrom, start, end, key, 33 `W_` columns, `openness_fetal`, `openness_maternal`)
plus `reference/manifest_wang_mfi.json` (provenance, per-cell-type origin/purity/
QC, grid stats, recipe).

## PC correlation (`pc_tissue_map.R --prefix W_`)

Runs unchanged against the `W_` layer (one edit: `cor(..., use =
"pairwise.complete.obs")`, needed because the single-cell layer has ~7 % NaN
bins; a no-op for the fully-covered bulk atlas). On the synthetic demo model
(PC3 seeded placental β=+0.6, PC7 seeded monocyte β=−0.5):

- **PC3** (FF-up) → syncytiotrophoblast: SCT_B/SCT/SCT_A/proSCT, β-oriented corr **+0.81..0.83**
- **PC7** (FF-down) → maternal monocytes: CD16_M/CD14_M **−0.71/−0.70**

## Origin marking (`pc_origin_map.R`) — new

```
Rscript pc_origin_map.R --model MODEL.rds \
    --atlas reference/wang_mfi_openness_hg19_50kb.csv.gz \
    --manifest reference/manifest_wang_mfi.json \
    --topn 2000 --n-perm 2000 --outdir OUT
```

Two outputs the tissue map cannot give:

1. **Per-PC origin score** `[corr(L[,k],fetal) − corr(L[,k],maternal)] · sign(β_k)`
   with a vectorized bin-shuffle null (n_perm). `>0` = this PC's FF-up direction
   opens fetal chromatin (expected for a fetal-fraction signal); `<0` = maternal.
   → `pc_origin_score.csv` (pc, beta, corr_fetal, corr_maternal, origin_score, z,
   p, call).

2. **Per-bin origin call + 2×2 Fisher.** Each top-N FF-up / FF-down bin gets its
   dominant Wang cell type (argmax row-centered specificity) and that type's
   origin label, then a FF-direction × origin contingency + Fisher exact test.
   → `bin_origin_calls.csv`, `origin_contingency.csv`.

**Synthetic demo result** (positive control):

| | Fetal | Maternal |
|---|---|---|
| FF-up bins (w>0) | **1,727** | 272 |
| FF-down bins (w<0) | 521 | **1,434** |

Fisher **OR = 17.5, p < 1e-3**. FF-up bins are dominated by syncytiotrophoblast
(SCT_B, SCT, VCT); FF-down by maternal immune cells (T, CD14/CD16 monocytes,
dNK, B). PC3 origin_score **+0.49 (z=81, p<1e-3, fetal)**, PC7 **+0.22 (z=38,
p<1e-3)** — both seeded PCs recovered on the correct origin side. This confirms
the FF sign convention (FF-up = fetal/placental, FF-down = maternal) against an
independent origin-labelled atlas.

## Status

Validated end-to-end on the **synthetic demo model** (a seeded positive control,
not a benchmark). The identical commands run on the real AWS-encrypted PC-100 FF
model when available — that turns these into genuine origin claims, with the Wang
labels as external ground truth. Figure: `_wang_demo/pc_origin_demo.png`.
