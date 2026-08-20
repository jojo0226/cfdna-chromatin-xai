# seqFF++ tissue-of-origin track — running it on your AWS FF model

This is the fetal-fraction fork of the chromatin-XAI pipeline. It answers one
question about your trained FF model: **when the model leans on a genomic
region, whose open chromatin is that region?** For a genuine fetal-fraction
model the leading tissue should be **placenta**; a model that only tracks total
cfDNA openness will spread across everything instead.

Everything here runs **in place on the AWS instance**. The only thing that
leaves the instance is a per-*feature* importance table (one number per
genomic bin) — never per-patient data.

---

## What ships in the repo

```
analysis/fetal_fraction/
├── reference/
│   └── ff_openness_atlas_hg19_50kb.csv.gz   # 57,633 hg19 50kb bins × 11 tissues
├── ff_tissue_proportion.py                   # analysis library
├── run_ff_tissue_track.py                    # CLI (Python)
└── extract_ff_shap.R                         # importance extractor (R: PCA+lm or glmnet)
```

**Atlas tissues (11):** `placenta` (fetal signal); `liver`, `endothelial`
(maternal solid); `monocyte`, `Bcell`, `CD4`, `CD8`, `NK`, `neutrophil`
(maternal hematopoietic, mature); `CMP` (CD34+ common myeloid progenitor);
`K562` (control). Four openness tracks per tissue: `D_` DNase, `H_` histone
(mean of H3K4me3 / H3K4me1 / H3K27ac / H3K36me3), `ACC_` accessibility
composite, `C_` combined (the default the track uses). Two coverage caveats,
both handled by `C_`: **neutrophil is histone-only** (no ENCODE DNase in
hg19/GRCh38 → `C_neutrophil = z(H)`); **CMP has no H3K27ac** narrowPeak in hg19.
hg19, autosomes only. See `reference/README.md` provenance notes and
`build_openness_atlas.py` to regenerate.

---

## The three-step run

### Step 1 — extract per-bin importance from the model (R, on AWS)

**Your NIPT model is a two-stage linear pipeline: PCA → linear regression.**
The `.rds` is a list with `$bin_loadings` (52,459 genomic 50 kb bins ×
PC1…PC100), `$beta` (Intercept + PC1…PC100 regression coefficients on
`FF_Yplus`), `$bin_params`, `$pc_sdev`, `$train_dt`. The extractor auto-detects
this shape.

**Why SHAP on the PCs is the wrong move — and what to do instead.** A PC is a
whole-genome linear combination of bins; "importance of PC17" has no genomic
location, so it can't be mapped to a tissue. But because *both* stages are
linear, the pipeline collapses **exactly** into one linear model in bin space:

```
score_k(sample) = Σ_bin  loadings[bin,k] · (x_bin − x̄_bin)
FF(sample)      = b0 + Σ_k beta_k · score_k
                = b0 + Σ_bin  w_bin · (x_bin − x̄_bin)
```

where the **effective per-bin weight** is a single matrix–vector product:

```
w_bin = Σ_k  loadings[bin, k] · beta_k        # bin_loadings %*% beta[PC1..PC100]
```

`w` *is* the model's genome-wide importance track. SHAP is then closed form in
bin space (verified exact to ~1e-14 — prediction and Shapley sum both
reconstruct):

```
phi_bin(sample) = w_bin · (x_bin − x̄_bin)              # exact Shapley value
global_bin      = mean_i |phi_bin| = |w_bin| · mean_i |x_bin − x̄_bin|
```

`extract_ff_shap.R` runs `bin_loadings %*% beta`, writes the two-column table
(`key,mean_abs_shap`), and needs **nothing but the model object** for the
model-level track (`--x` optional):

```bash
Rscript analysis/fetal_fraction/extract_ff_shap.R \
    --model  /secure/ff_model.rds \
    --x      /secure/bin_matrix.csv \          # optional: the per-sample bin coverage
    --out    ff_shap_importance.csv            #           matrix you fed into the PCA
```

- **Without `--x`** (model-only, exact up to per-bin scale): importance =
  `|w_bin|`. Comes entirely from `$bin_loadings` and `$beta` — no patient data
  touched at all. This is usually all you need.
- **With `--x`** (exact global linear-SHAP): `|w_bin| · mean_i|x_bin − x̄_bin|`,
  weighting each bin by how much it actually varies in your cohort. The matrix
  is read only for per-column means and mean-abs-deviation; the output stays a
  per-bin aggregate, never per-sample. `--x` is the same bin matrix that went
  into the model's PCA (columns = `$bin_loadings$bin_name`).

**Bin-key mapping.** Your bins are pure genomic 50 kb (`chr1_550000_600000`),
so the only transform is a name reformat to the atlas key form:

```
chr1_550000_600000   →   chr1:550000-600000
```

The extractor does this automatically (it also accepts `chr:start-end`,
`chr:start:start`, `chr.start.end`) and drops anything non-genomic. Console:

```
[compose] 52459 bins x 100 PCs -> per-bin weights
[map]     52459/52459 features -> genomic bins (0 dropped: non-genomic)
```

> **Grid match — no liftover.** Your 52,459 bins are hg19 50 kb, a subset of
> the atlas's 57,633 hg19 50 kb bins, so they map 1:1. Bins the model has but
> the atlas doesn't (or vice-versa) are simply reported in the coverage line;
> `run_ff_tissue_track.py` uses the intersection. If your build/width ever
> differs, rebuild the atlas with `build_openness_atlas.py` (`--bin-size`,
> `--assembly`) — never lift your cohort features.

> **glmnet fallback.** If you ever point the extractor at a bare
> `glmnet`/`cv.glmnet` fit in bin space instead, it detects that and uses
> `phi = beta·(x−x̄)` with `--lambda min|1se|<numeric>`. The PCA+lm path above
> is what your current `.rds` triggers.

### Step 1b — split by SHAP direction (optional but recommended)

`extract_ff_shap.R` now also emits a **`signed_shap`** column carrying the sign
of the effective weight `w_bin` (the magnitude column `mean_abs_shap` is
unchanged, so the downstream tissue track is fully backward-compatible):

- **`signed_shap > 0`** — more coverage/openness in the bin raises predicted FF
  → the **fetal / placental-leaning** direction.
- **`signed_shap < 0`** — lowers predicted FF → the **maternal-leaning**
  direction.

Because per-sample SHAP is `phi_bin(i) = w_bin·(x_i,bin − x̄_bin)` and `x` is
centered, the mean signed φ is ~0 by construction; the meaningful *global*
signed quantity is `w_bin` itself (or `w_bin · mean_i|x−x̄|` with `--x`). The
two directions are biologically distinct — pooling them by `|SHAP|` mixes
"open chromatin drives FF up" with "…drives FF down" — so we split them before
the cell-of-origin test:

```bash
python analysis/fetal_fraction/select_signed_shap.py \
    --importance ff_shap_importance.csv \
    --topns 5000,2000,1000,500 \
    --outdir ff_signed_sets
```

This writes, for each cutoff N, `pos_top{N}.csv` (largest positive `w_bin`) and
`neg_top{N}.csv` (most negative), plus `signed_selection_manifest.json` with
per-cutoff counts and value ranges. Each set is then fed to Step 2 separately
(`--importance ff_signed_sets/pos_top2000.csv`, etc.) so the tissue-of-origin
readout is resolved *per direction* — the expectation for a genuine FF model is
that the **positive** set concentrates in placenta/trophoblast openness while
the **negative** set leans hematopoietic (monocyte/neutrophil).

### Step 2 — run the tissue track (Python, on AWS or anywhere)

The importance table has no patient data, so this step can run on the instance
or on your laptop. The atlas path auto-resolves to the shipped `reference/`.

```bash
python analysis/fetal_fraction/run_ff_tissue_track.py \
    --importance ff_shap_importance.csv \
    --outdir     out_ff_tissue
# columns autodetect; name them if needed:
#   --key-col key --shap-col mean_abs_shap
# other tracks: --prefix D_  (DNase)  | H_ (histone) | ACC_ (accessibility)
```

Outputs in `out_ff_tissue/`:

| file | contents |
|------|----------|
| `ff_tissue_results.csv` | per (top-N × tissue): `abs_z/abs_p` (absolute openness), `spec_z/spec_p` (tissue-specific), bins used |
| `ff_tissue_proportions.csv` | normalized tissue-of-origin share per top-N |
| `ff_tissue_meta.json` | grid coverage, matched fraction, params |
| `ff_tissue_track.png` | 2-panel figure (specificity + proportion) |

### Step 3 — read the two panels

- **Panel A / `spec_z` (cell-of-origin).** `spec_<tissue> = C_<tissue> − mean
  over tissues`, so it removes the pan-openness component. **Placenta leading
  here, significant against the matched-random-region null, is the result you
  want** — it says the model's attribution falls on fetal-specific open
  chromatin, i.e. it learned tissue-of-origin, not just "cfDNA is open here."
- **Panel B / proportion.** Normalized share of positive specificity — a
  quick "what fraction of the tissue-of-origin signal is placenta vs maternal."
- **`abs_z` (absolute openness).** If *every* tissue is high on `abs_z` but
  `spec_z` is flat, the model is riding general active chromatin — a
  pan-openness confound rather than cell-of-origin. That is the negative reading.

The console prints the verdict automatically (`placenta leads` / `placenta does
NOT lead`).

### Step 4 — open the PCA black box: which PC ↔ which tissue (optional)

Steps 1–3 collapse the PCA away and explain the model in bin space. This step
looks *inside* the latent basis: it asks which principal components carry the
fetal vs maternal signal, and what tissue biology each one encodes. Two exact
facts make this well-defined:

- Each PC *k* has a genome-wide loading pattern `L[:,k]`. Correlating it against
  a tissue's openness vector (β-oriented, so the sign tracks the FF-raising
  direction) says **which tissue that latent component looks like**.
- The per-bin weight decomposes exactly: `w_bin = Σ_k L[bin,k]·β_k`. So
  `Contrib[bin,k] = L[bin,k]·β_k` is PC *k*'s additive share of that bin's SHAP,
  and averaging it over the top FF-up vs FF-down bin sets shows **which PCs
  build each SHAP direction**.

```bash
Rscript analysis/fetal_fraction/pc_tissue_map.R \
    --model  /secure/ff_model.rds \
    --atlas  analysis/fetal_fraction/reference/ff_openness_atlas_hg19_50kb.csv.gz \
    --pos    ff_signed_sets/pos_top2000.csv \
    --neg    ff_signed_sets/neg_top2000.csv \
    --outdir pc_tissue_out
python analysis/fetal_fraction/plot_pc_tissue.py \
    --indir pc_tissue_out --topk 15 --out pc_tissue_figure.png
```

Outputs in `pc_tissue_out/`: `pc_tissue_corr.csv` (PC × tissue β-oriented
correlations + each PC's single best-aligned tissue), `pc_shap_contrib.csv`
(per-PC mean contribution to the pos/neg bin sets), `pc_tissue_meta.json`. The
figure has two panels — **A** the PC↔tissue alignment heatmap (rows = top PCs by
total signed contribution, red = this PC's FF-up direction aligns with the
tissue's openness), **B** the per-PC contribution to the FF-up (red) vs FF-down
(blue) bin sets. The expectation for a genuine FF model: the PCs that dominate
the **positive** bars align with placenta/trophoblast openness, while the PCs
dominating the **negative** bars align with hematopoietic tissues. Only the
model object and the public atlas are read — no patient data enters this step.

---

### Step 5 — repressive (H3K27me3) axis: the two-sided cell-of-origin test (optional)

Steps 1–4 explain the model through **open** chromatin (DNase). A genuine
cell-of-origin signal should also show the mirror image: where a tissue's DNA
contributes to cfDNA, that tissue's regions are open **and** free of its own
repressive mark. This step adds a fully **separate** H3K27me3 axis — it is
never merged into the `C_` openness columns — so you can require both
directions before calling a tissue the source.

Build the repressive atlas once (ENCODE hg19 H3K27me3 narrowPeak, ≤3
replicates/biosample, z-mean coverage, pan-tissue-subtracted specificity),
then plot it against the openness atlas:

```bash
python analysis/fetal_fraction/build_repressive_atlas.py \
    --probe  h3k27me3_probe.json \
    --outdir analysis/fetal_fraction/reference   # -> ff_repressive_atlas_*.csv.gz
python analysis/fetal_fraction/plot_repressive.py \
    --repressive analysis/fetal_fraction/reference/ff_repressive_atlas_hg19_50kb_H3K27me3.csv.gz \
    --openness   analysis/fetal_fraction/reference/ff_openness_atlas_hg19_50kb_curated73.csv.gz \
    --outdir     repr_figs
```

The shipped atlas covers **39 tissues** — every cfDNA-critical one (placenta
×3, trophoblast, chorion, monocyte, B/CD4/CD8/NK, liver, CMP, HUVEC). Two
figures come out:

- **`repr_fig1_clustering.png`** — tissue×tissue correlation of the H3K27me3
  *specificity* profile, hierarchically ordered. A real repressive axis
  clusters by lineage: the hematopoietic block (CMP/NK/Bcell/CD4/CD8) is tight
  and the placenta/chorion pair sits together — confirmed in the shipped run.
- **`repr_fig2_twosided.png`** (+ `repr_twosided_rho.csv`) — per-tissue
  Spearman ρ between active-specificity and repressive-specificity. **ρ < 0**
  means the two-sided test is valid for that tissue. The cfDNA-critical tissues
  are the most anti-correlated (placenta ρ≈−0.16, CD8/CD4/CMP ρ≈−0.17 to −0.21).

**Why specificity, not raw z.** At 50 kb, raw DNase and raw H3K27me3 both
concentrate in gene-rich euchromatin, so raw-vs-raw is *positively* correlated
(a gene-density baseline, ρ≈+0.2 — this is expected, not a failure).
Pan-tissue subtraction removes that shared baseline; the residual specificity
is the interpretable quantity, and only then does the expected anti-correlation
appear. As with every other step, only public ENCODE tracks and the model's
per-bin aggregate are used — no patient data.

---

## Data-egress summary

| artifact | leaves instance? | why safe |
|----------|------------------|----------|
| `ff_model.rds`, bin matrix | **no** | stay on AWS; read in place by R |
| `ff_shap_importance.csv` | yes (optional) | per-bin aggregate — `|w|` or `|w|·mean|x−x̄|`, no per-sample rows |
| track outputs / figure | yes | derived from importance + public ENCODE atlas |

No patient-level cfDNA ever moves. The importance table is a genome-indexed
vector of the same kind you would put in a methods supplement.

---

## Validation done before shipping

The library and CLI were checked against the real 57,633-bin, 11-tissue atlas:

- **PCA→bin composition is exact** — on a synthetic PCA+lm model, both the
  prediction (`b0 + Σ w·(x−x̄)`) and the Shapley sum (`Σ phi + b0`) reconstruct
  the PC-space prediction to ~1e-14; global `|w|` correlates with per-sample
  `mean|phi|` at r = 0.997.
- **Extractor on the real object shape** — a synthetic `.rds` with your exact
  fields (`$bin_loadings` 3,000 bins × PC1..PC100, `$beta`, `$bin_params`,
  `$train_dt`) is auto-detected as PCA+lm, composed to per-bin weights, and its
  `chr_start_end` names normalize to atlas keys 1:1 (0 dropped).
- **Positive control** — importance planted on placenta-specific bins recovers
  placenta as the leading tissue (`spec_z ≈ 127`, 100% of the top-500 share);
  maternal tissues correctly pushed negative.
- **Negative control** — random importance yields a non-significant leader
  (`p ≈ 0.07–0.28`), no false tissue signal.
- **Full R→Python chain** — synthetic PCA+lm `.rds` → `extract_ff_shap.R`
  (both `--x` and model-only modes) → `run_ff_tissue_track.py` runs clean,
  100% grid coverage, all four outputs written.
