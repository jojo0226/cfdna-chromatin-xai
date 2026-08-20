# cfdna-chromatin-xai — Fetal-Fraction Module: Complete Method Walkthrough

*A step-by-step companion to the math slide deck. Every formula here is taken verbatim from the tool's source (`extract_ff_shap.R`, `ff_tissue_proportion.py`, `pc_tissue_map.R`, and the plot scripts). For each stage: **concept · formula · input · output · how to interpret · script / how to run.***

> **Running the scripts.** All commands below are run from `analysis/fetal_fraction/`. Two interpreters are used: **R** scripts (`extract_ff_shap.R`, `pc_tissue_map.R`, `pc_origin_map.R`) via `Rscript`, and **Python** scripts (everything else) via `python`. The one model object every R stage reads is a PCA→lm `.rds` with a `bin_loadings` matrix (bin × PC) and a named `beta` vector — either your trained FF model, or `synth_pca_model.rds` (the seeded positive-control used for all demo numbers here). The default atlas is the DNase openness table `reference/ff_openness_atlas_hg19_50kb.csv.gz` (a curated-73 variant, `..._curated73.csv.gz`, and the `A_` placenta-ATAC, `W_` Wang-MFI, and H3K27me3 repressive tables all sit alongside it in `reference/`); swap `--atlas`/`--prefix` to use any other layer.
>
> **Two cwd conventions.** `run_ff_tissue_track.py` resolves its default atlas relative to the *script's own location*, so it works from anywhere. The R scripts (`pc_tissue_map.R`, `pc_origin_map.R`) default to a **repo-root-relative** atlas path — so either run them from the repo root, or (running from `analysis/fetal_fraction/`) pass `--atlas reference/<file>.csv.gz` explicitly, as shown in each command below.

**End-to-end, in order:**
```bash
# 1) importance: model → per-bin weight + SHAP  (Stages 3–4)
Rscript extract_ff_shap.R --model synth_pca_model.rds --out ff_shap_importance.csv
# 2) partition into FF-up / FF-down sets         (Stage 5)
python  select_signed_shap.py --importance ff_shap_importance.csv --outdir ff_signed_sets
# 3) cell-of-origin enrichment + bars/lollipop   (Stages 6–7)  — ranks by |SHAP|
python  run_ff_tissue_track.py --importance ff_shap_importance.csv --prefix C_ --outdir out_ff_tissue
# 3b) directional runs on the signed subsets (for the dumbbell)
python  run_ff_tissue_track.py --importance ff_signed_sets/pos_top2000.csv --prefix C_ --outdir out_ff_up   # FF-up
python  run_ff_tissue_track.py --importance ff_signed_sets/neg_top2000.csv --prefix C_ --outdir out_ff_dn   # FF-down
# 4) POS-vs-NEG contrast dumbbell                (Stage 8)
python  plot_ff_tissue_dumbbell.py --pos out_ff_up/ff_tissue_results.csv --neg out_ff_dn/ff_tissue_results.csv
# 5) PC ↔ tissue mechanism + figure              (Stages 9–10)
Rscript pc_tissue_map.R --model synth_pca_model.rds --atlas reference/ff_openness_atlas_hg19_50kb.csv.gz --outdir pc_tissue_out
python  plot_pc_tissue.py --indir pc_tissue_out
# 6) maternal-vs-fetal origin (Wang MFI layer)   (optional, cross-check)
Rscript pc_origin_map.R --model synth_pca_model.rds --atlas reference/wang_mfi_openness_hg19_50kb.csv.gz --manifest reference/manifest_wang_mfi.json
```

---

## 0. What this tool is

The fetal-fraction (FF) module is a **feature-agnostic explainer** wrapped around a fetal-fraction predictor built as **PCA → linear regression** (the seqFF++ NIPT model). It answers three questions in sequence, each feeding the next:

1. **Importance** — which genomic bins drive the FF prediction, and in which direction? *(SHAP, closed-form)*
2. **Explanation** — what cell types are those bins open in? *(tissue enrichment vs a null model)*
3. **Mechanism** — which latent principal components build each FF direction, and what tissue does each resemble? *(PC ↔ tissue decomposition)*

A defining constraint: the explainer reads **only** the trained model's `bin_loadings` and `beta`. No per-patient data is required for the importance and mechanism tracks, and none ever leaves the instance — only per-feature aggregates do.

The whole pipeline in one line:

$$\text{coverage } x_{i,b}\ \xrightarrow{\text{PCA}\to\text{lm}}\ \hat{FF}_i\ \xrightarrow{\text{SHAP}}\ w_b,\ \phi_{i,b}\ \xrightarrow{\text{atlas}}\ \text{cell of origin}$$

---

## 1. The feature — normalized coverage residual

**Concept.** The model does not see raw read counts. Raw coverage is confounded by GC bias and sequencing depth. The input feature is the *residual* left after regressing those technical factors out, so a non-zero value is biological signal, not an artifact.

**Formula.**
$$x_{i,b} = \text{coverage}_{i,b} - \mathbb{E}\!\left[\text{coverage}\mid \text{GC}_b,\ \text{mappability}_b\right]$$

**Input.** Whole-genome cfDNA sequencing for sample $i$, binned into 50 kb windows on hg19 autosomes (~52,000 bins).

**Output.** One corrected coverage residual $x_{i,b}$ per (sample, bin).

**How to interpret.** A positive residual means the bin is *over-represented* relative to what GC predicts (relatively protected / less fragmented); a negative residual means *under-represented*. The biologically crucial inversion: **open chromatin fragments more readily → fewer intact fragments span the bin → lower coverage.** So "open chromatin ⇒ low coverage." This is why, later, a *negative*-weight bin corresponds to *high* placental openness.

**Script / how to run.** *Produced upstream — no script in this module.* The GC/mappability-corrected residual matrix is an **input** to the explainer, generated by your NIPT coverage-preprocessing pipeline (per-bin GC + negative-binomial depth normalization). The explainer never re-derives it; it consumes the resulting `.rds` model (Stage 2). If you have the training coverage matrix as a CSV (bins × samples), you only need it for the *exact* linear-SHAP magnitude in Stage 4 (`extract_ff_shap.R --x`); the entire signed-direction pipeline runs from the model object alone.

---

## 2. The model — PCA, then linear regression

**Concept.** With ~52 k bins and far fewer training samples, regressing FF directly on bins is rank-deficient and overfits. PCA first compresses the bins into ~100 orthogonal components; the regression then runs on those components.

**Formula.**
$$\text{score}_k(i) = \sum_b L[b,k]\,(x_{i,b}-\bar{x}_b), \qquad \hat{FF}_i = b_0 + \sum_{k=1}^{100} \beta_k\,\text{score}_k(i)$$

**Input.** The centered bin matrix $(x_{i,b}-\bar x_b)$; the trained loadings $L$ (bin × PC) and regression coefficients $\beta$.

**Output.** A scalar prediction $\hat{FF}_i$ per sample, plus the two reusable model objects $L$ and $\beta$.

**How to interpret.**
- $L[b,k]$ — **loadings.** Column $k$ is a genome-wide spatial pattern (a "meta-bin"): which bins move together.
- $\beta_k$ — **coefficient.** How strongly, and in which direction, PC $k$ moves FF. The *sign* of $\beta_k$ is load-bearing throughout the rest of the pipeline.

Because both maps are linear, they collapse **exactly** into a single linear model in bin space — no approximation. That collapse is what makes SHAP closed-form.

**Script / how to run.** *Trained upstream — supplied as an `.rds` to the explainer.* The model is fit by your NIPT training code (e.g. `prcomp` on the centered bin matrix, then `lm`/`glmnet` of FF on the top PCs) and saved as a list with a `bin_loadings` matrix (rows = bins named `chr1_0_50000`, cols = `PC1..PC100`) and a named `beta` vector. To generate a drop-in **positive control** with a seeded fetal PC and maternal PC (used for every demo number below):
```r
# reproduces synth_pca_model.rds — see the header of extract_ff_shap.R
# PC3 seeded to placenta (β=+0.6), PC7 to CD14 monocyte (β=−0.5)
```
Confirm the object before running anything else:
```bash
Rscript -e 'm<-readRDS("synth_pca_model.rds"); print(dim(m$bin_loadings)); print(head(m$beta))'
```

---

## 3. Effective per-bin weight $w_b$

**Concept.** Substitute the score into $\hat{FF}$ and swap the order of summation. The two stacked linear maps become one weight per bin.

**Formula.**
$$\hat{FF}_i = b_0 + \sum_b \Big(\underbrace{\textstyle\sum_k L[b,k]\,\beta_k}_{w_b}\Big)(x_{i,b}-\bar x_b), \qquad \boxed{\,w_b = \sum_k L[b,k]\,\beta_k\,}$$

Equivalently the matrix–vector product $w = L\beta$.

**Input.** $L$ and $\beta$ from the model object — nothing else.

**Output.** One signed number $w_b$ per bin: the bin's net effect on predicted FF.

**How to interpret** *(verbatim from `extract_ff_shap.R`)*:
- $w_b > 0$ — more coverage/openness in that bin ⇒ **higher** predicted FF → **fetal / placental-leaning**.
- $w_b < 0$ — more coverage/openness there ⇒ **lower** predicted FF → **maternal-leaning**.

Because $w_b$ needs only the trained model, the entire importance track is computable with **no patient data at all**.

**Script / how to run.** `extract_ff_shap.R` (Stages 3 **and** 4 in one pass — it computes $w_b=L\beta$ and both SHAP columns).
```bash
Rscript extract_ff_shap.R --model synth_pca_model.rds --out ff_shap_importance.csv
```
Arguments: `--model` (required, the `.rds`) · `--out` (output CSV, default `ff_shap_importance.csv`) · `--lambda` (`min`/`1se`/numeric, only if the model is a `cv.glmnet`, default `min`) · `--x` (optional coverage matrix CSV — see Stage 4). Output CSV columns: `feature`, `key` (`chr1:0-50000`), `mean_abs_shap`, `signed_shap`.

---

## 4. SHAP — closed-form importance

**Concept.** For a linear model the Shapley value has an exact closed form: it is the feature's weight times its centered value. No sampling, no KernelExplainer.

**Formula (per sample).**
$$\phi_{i,b} = w_b\,(x_{i,b}-\bar x_b), \qquad \sum_b \phi_{i,b} = \hat{FF}_i - b_0 \quad(\text{local accuracy})$$

**Formula (aggregated across samples — the two emitted tracks).**
$$\text{mean\_abs\_shap}_b = |w_b|\cdot\overline{|x_{i,b}-\bar x_b|}, \qquad \text{signed\_shap}_b = w_b\cdot\underbrace{\overline{|x_{i,b}-\bar x_b|}}_{\text{mad}_b\,\ge\,0}$$

**Input.** $w_b$ (always); optionally the training coverage matrix $X$ to compute $\text{mad}_b=\overline{|x_{i,b}-\bar x_b|}$.

**Output.** Two per-bin columns: `mean_abs_shap` (magnitude) and `signed_shap` (direction).

**How to interpret.**
- **The sign identity.** Since $\text{mad}_b \ge 0$, it can never flip the sign: $\operatorname{sign}(\text{signed\_shap}_b) \equiv \operatorname{sign}(w_b)$. The aggregate signed track is anchored to the model weight, patient-independent.
- **Two operating modes.**
  - *Model-only track* (no matrix supplied): importance $= |w_b|$. Ranking bins by signed SHAP is identical to ranking by signed weight.
  - *Exact linear-SHAP track* (matrix available, `imp_mode = pca_linear_shap_exact`): magnitude $= |w_b|\cdot\text{mad}_b$ — the weight scaled by how variable that bin's coverage is across patients. The **direction is unchanged**, but *which* bin is most extreme gets reweighted: a huge-weight bin with flat coverage lands near zero; a moderate-weight bin with variable coverage can outrank it.
- **Per-sample caveat.** At the individual level, $\phi_{i,b}=w_b(x_{i,b}-\bar x_b)$ flips sign for any patient whose coverage in that bin is *below* the bin mean. Only the aggregate track has a stable sign.

**Figure — coverage / SHAP track.** Two stacked panels sharing one x-axis (genomic region, chr1→chr22). Top panel y = normalized coverage $x_{i,b}$ (the input); bottom panel y = SHAP value $\phi_{i,b}$ (the signed contribution). This is the *local* explanation for one prediction: a tall positive bar pushed that sample's FF up, a deep negative bar pushed it down. Aligning the panels shows *why* — a bin earns a large $|\phi|$ only when its coverage deviates from the mean **and** it carries a non-trivial weight.

**Script / how to run.** Same script as Stage 3, `extract_ff_shap.R` — the two SHAP columns are emitted in the same run. To get the **exact linear-SHAP magnitude** ($|w_b|\cdot\text{mad}_b$) instead of the model-only $|w_b|$, add the training coverage matrix:
```bash
# model-only track (magnitude = |w_b|), no patient data:
Rscript extract_ff_shap.R --model synth_pca_model.rds --out ff_shap_importance.csv
# exact linear-SHAP track (magnitude = |w_b|·mad_b), needs coverage matrix:
Rscript extract_ff_shap.R --model synth_pca_model.rds --x coverage_matrix.csv --out ff_shap_importance.csv
```
`--x` is a bins × samples CSV of centered/uncentered coverage residuals (the script centers per bin). When supplied, the run reports `imp_mode = pca_linear_shap_exact`; the **direction (`signed_shap` sign) is identical either way** — only the `mean_abs_shap` ranking changes. *Figure:* the coverage/SHAP panel is drawn by `plot_ff_coverage_shap.py` (edit its data block to point at your `coverage_df` + `shap_df`, then `python plot_ff_coverage_shap.py`).

---

## 5. Partition — FF-up vs FF-down bin sets

**Concept.** Split every bin by the sign of its weight, then take the top-$N$ most extreme in each direction.

**Formula.**
$$\mathcal{U}_N = \text{top-}N\{b : w_b > 0\}\ \text{(FF-up)}, \qquad \mathcal{D}_N = \text{top-}N\{b : w_b < 0\}\ \text{(FF-down)}$$
with defaults $N \in \{500, 1000, 2000, 5000\}$.

**Input.** The `signed_shap` column from Step 4.

**Output.** Two ranked bin sets per cutoff — the FF-up ("POS") and FF-down ("NEG") sets that every downstream figure consumes.

**How to interpret — why "FF-up / FF-down" and not "positive / negative SHAP".** They are the *same* partition of bins; the naming is deliberate:

| Name | Sign anchored to | Ambiguous? |
|---|---|---|
| positive / negative SHAP | arithmetic sign of a number — but *aggregate* or *per-sample* $\phi_{i,b}$? | Yes — per-sample $\phi$ flips for patients below the bin mean |
| **FF-up / FF-down** | $\operatorname{sign}(w_b)$ — the sample-independent model weight | No — names the biological direction the bin pushes FF |

"FF-up" states *what the sign means for the model output* and sidesteps the per-sample sign-flip ambiguity that "positive SHAP" invites. That is why the tool and its docs use it.

**Script / how to run.** `select_signed_shap.py`.
```bash
python select_signed_shap.py --importance ff_shap_importance.csv --outdir ff_signed_sets
```
Arguments: `--importance` (required, the CSV from Stage 3–4) · `--topns` (comma-separated cutoffs, default `5000,2000,1000,500`) · `--outdir` (default `ff_signed_sets`) · `--key-col`/`--signed-col`/`--abs-col` (override auto-detected column names). Writes one FF-up and one FF-down bin list per cutoff (e.g. `ff_up_top2000.txt`, `ff_dn_top2000.txt`). Note: `run_ff_tissue_track.py` (Stage 6) does this partition internally, so this standalone step is only needed when you want the raw bin lists as files.

---

## 6. Explanation — cell-of-origin enrichment vs a matched null

**Concept.** Ask what cell types the top SHAP bins are preferentially open in — but calibrate against random genome so "these bins are just open everywhere" cannot masquerade as a real signal.

**Formula.** For each top-$N$ set and each tissue $t$:
$$\text{spec\_obs}_t = \frac1N\!\!\sum_{b\in\text{top-}N}\!\! SP[b,t], \qquad \text{spec\_z}_t = \frac{\text{spec\_obs}_t - \mu^{\text{null}}_t}{\sigma^{\text{null}}_t}$$
where the null $(\mu^{\text{null}}_t, \sigma^{\text{null}}_t)$ comes from drawing $N$ random bins and recomputing the mean, repeated $n_{\text{perm}}$ times (default 2000). Two-sided empirical $p$ = fraction of null draws at least as extreme as observed. The composition track is
$$\text{share}_t = \frac{[\text{spec\_obs}_t]_+}{\sum_{t'}[\text{spec\_obs}_{t'}]_+}.$$

**Input.** A top-$N$ bin set (Step 5) and the tissue **specificity atlas** $SP[b,t]$ (73 tissues, per-bin z-scored openness).

**Output.** Per (cutoff, tissue): `spec_obs`, `spec_z`, `spec_p`, `abs_*` counterparts, and the normalized `share`. Written to `ff_tissue_results.csv` and `ff_tissue_proportions.csv`.

**How to interpret.** `spec_z` = how many standard deviations *more tissue-specific-open* the SHAP-selected bins are than random genome → **calibrated significance** ("is tissue $t$ real?"). `share` = tissue $t$'s slice of the total positive origin signal → **composition** ("how big is $t$?"). A tissue leading in *both* is a solid cell-of-origin call.

**Figure — tissue track (grouped bars, default).** Panel A: x = tissue, y = spec_z, grouped bars = the four cutoffs. Panel B: x = tissue, y = share. On the synthetic FF-up set, placenta dominates both (spec_z ≈ 105) — the fetal-signal confirmation.

**Script / how to run.** `run_ff_tissue_track.py` (CLI wrapper around `ff_tissue_proportion.py` — computes the enrichment **and** draws the figure).
```bash
python run_ff_tissue_track.py --importance ff_shap_importance.csv --prefix C_ --outdir out_ff_tissue
```
Key arguments: `--importance` (required) · `--atlas` (default the curated-73 DNase table; point at `reference/ff_placenta_atac_hg19_50kb.csv.gz`, the Wang MFI or repressive tables to switch layers) · `--prefix` (atlas column prefix — `C_` bulk DNase, `A_` placenta ATAC, `W_` Wang MFI; default `C_`) · `--topns` (default `500,1000,2000,5000`) · `--n-perm` (null draws, default 2000) · `--seed` (default 0) · `--outdir` (default `out_ff_tissue`). **This ranks bins by `|SHAP|` (both directions together)** — it is the main importance run. Outputs: `ff_tissue_results.csv` (per topn × tissue: `spec_obs/z/p`, `abs_*`), `ff_tissue_proportions.csv` (per-tissue `share`), `ff_tissue_meta.json`, `ff_tissue_track.png`. For **directional** enrichment (needed by the Stage 8 dumbbell), first run `select_signed_shap.py` (Stage 5), then point `--importance` at `ff_signed_sets/pos_top2000.csv` (FF-up) and `neg_top2000.csv` (FF-down) into separate `--outdir`s.

---

## 7. Explanation, new layout — lollipop *(new in v0.6.2)*

**Concept.** At 73 tissues the vertical bar labels shrink to ~7 pt. The lollipop transposes and sorts the same data so tissue names stay horizontal and legible, and shows enrichment stability across cutoffs.

**Formula.** No new statistic — it plots the same `spec_z` and `share` as Step 6.

**Input.** The `ff_tissue_results.csv` / `ff_tissue_proportions.csv` from Step 6. CLI: `--style lollipop` (plus `--top-k`, default 20; placenta always kept).

**Output.** `ff_tissue_track.png` in lollipop form.

**How to interpret.** Axes flip: y = tissue (top-K by $|\text{spec\_z}|$, labels horizontal), x = spec_z (panel A) / share (panel B). One colored dot **per top-$N$ cutoff**, stem from zero; the near-zero tail is collapsed into a "+N more near zero" note. Because the dots for all cutoffs sit on one row, you can read at a glance whether a tissue's enrichment is **stable** as the bin set widens — a robustness signal the grouped bars don't surface.

**Script / how to run.** Same script as Stage 6, `run_ff_tissue_track.py`, with `--style lollipop` (no recomputation — it re-plots the same enrichment):
```bash
python run_ff_tissue_track.py --importance ff_shap_importance.csv --prefix C_ --style lollipop --top-k 20 --outdir out_ff_tissue
```
Added arguments for this layout: `--style {bars,lollipop}` (default `bars`) · `--top-k` (tissues to show, default 20; placenta always kept) · `--layout {auto,row,col,separate}` · `--tick-fs` (tick font size override). Writes `ff_tissue_track.png` in lollipop form.

---

## 8. Contrast, new figure — POS-vs-NEG dumbbell *(new in v0.6.2)*

**Concept.** Put both FF directions on a single axis so the fetal↔maternal contrast per tissue is directly visible.

**Formula.**
$$\text{contrast}_t = \text{spec\_z}^{(\text{FF-up})}_t - \text{spec\_z}^{(\text{FF-down})}_t$$

**Input.** The *two* `ff_tissue_results.csv` files — one from the FF-up (POS) run, one from the FF-down (NEG) run. No recomputation, no atlas reload. CLI: `--pos`, `--neg`, `--top-k` (default 25).

**Output.** `ff_tissue_dumbbell.png`.

**How to interpret.** y = tissue (top-K by contrast), x = spec_z, with **two dots per tissue** — FF-up (red) and FF-down (blue) — joined by a line. The **line length is the contrast**: how strongly that tissue distinguishes fetal from maternal bins. Placenta sits far right on its POS dot (fetal); CD14-positive monocyte does the reverse (maternal hematopoietic). On the demo, placenta contrast ≈ +119, CD14 monocyte ≈ −143.

**Script / how to run.** `plot_ff_tissue_dumbbell.py`. It consumes the *two* directional `ff_tissue_results.csv` from the FF-up and FF-down runs of Stage 6 (see the "directional" note there — build them from `select_signed_shap.py` subsets):
```bash
# prerequisites: the two directional runs
python select_signed_shap.py --importance ff_shap_importance.csv --outdir ff_signed_sets
python run_ff_tissue_track.py --importance ff_signed_sets/pos_top2000.csv --prefix C_ --outdir out_ff_up
python run_ff_tissue_track.py --importance ff_signed_sets/neg_top2000.csv --prefix C_ --outdir out_ff_dn
# the dumbbell itself:
python plot_ff_tissue_dumbbell.py --pos out_ff_up/ff_tissue_results.csv --neg out_ff_dn/ff_tissue_results.csv --outdir . --out ff_tissue_dumbbell.png
```
Arguments: `--pos` / `--neg` (required, the two results CSVs) · `--metric` (default `spec_z`) · `--top-k` (default 25) · `--topn` (which cutoff row to read, default = all) · `--outdir` / `--out`. No atlas reload, no permutation — pure plotting of pre-computed enrichment.

---

## 9. Mechanism — PC ↔ tissue alignment (β-oriented)

**Concept.** The tissue track treats the model as a weight per bin, but those weights are *built* from PCs ($w_b=\sum_k L[b,k]\beta_k$). This stage decomposes the model one layer down: how much does each PC's loading pattern resemble each tissue's openness — reoriented so the answer speaks in FF direction, not the PC's arbitrary polarity.

**Formula.**
$$\text{align}_{k,t} = \operatorname{corr}\!\big(L[\cdot,k],\ \text{openness}_t\big)\times \operatorname{sign}(\beta_k)$$

**Input.** $L$ (bin × PC), $\beta$, and the tissue openness atlas.

**Output.** `pc_tissue_corr.csv` — a PC × tissue matrix of β-oriented alignments, plus each PC's `best_tissue`.

**How to interpret.** The raw correlation answers only "does PC $k$ look like tissue $t$?" — direction-agnostic. Multiplying by $\operatorname{sign}(\beta_k)$ rotates it onto the FF axis:
- **red** = the *FF-raising* direction of PC $k$ resembles tissue $t$'s open chromatin;
- **blue** = the FF-raising direction anti-resembles $t$ (that tissue's openness *lowers* FF).

Without the flip, a PC with $\beta_k<0$ that correlates +0.5 with placenta would appear "placenta-positive" while actually lowering FF. After the flip, **red consistently means FF-up across the entire matrix** — which is the only reason the colorbar can carry one honest "FF-up →" label.

*Cell color vs colorbar:* they are one encoding. The **cell color is the value itself** (a β-oriented alignment, on `RdBu_r`, symmetric about 0). The **colorbar is its legend** — the ticks give the magnitude (white = 0), and the "FF-up →" label gives the meaning of the sign. To decode a cell, match its shade against the bar.

**Script / how to run.** `pc_tissue_map.R` (computes both Stage 9 alignment **and** Stage 10 contribution in one run).
```bash
Rscript pc_tissue_map.R --model synth_pca_model.rds --atlas reference/ff_openness_atlas_hg19_50kb.csv.gz --outdir pc_tissue_out
```
Arguments: `--model` (required, the `.rds`) · `--atlas` (default is repo-root-relative — pass `reference/...` when running from this dir; supply a `W_`/`A_` table to map PCs against another layer) · `--prefix` (atlas column prefix, default `C_`) · `--topn` (bins per FF direction for the Stage-10 contribution average, default 2000) · `--outdir` (default `pc_tissue_out`). Outputs: `pc_tissue_corr.csv` (PC × tissue β-oriented alignment + each PC's `best_tissue`) and `pc_shap_contrib.csv` (Stage 10). *Note v2:* the correlation uses `use="pairwise.complete.obs"` so single-cell layers with partial-coverage bins (e.g. the Wang `W_` atlas, ~7% NaN bins) don't collapse to NA.

---

## 10. Mechanism — exact per-PC contribution to each direction

**Concept.** Because $w_b = \sum_k L[b,k]\beta_k$ is an identity, each term $L[b,k]\beta_k$ is PC $k$'s **exact** additive share of that bin's weight. Averaging those shares over the FF-up and FF-down sets shows which PCs construct each direction.

**Formula.**
$$\text{Contrib}_{b,k} = L[b,k]\,\beta_k, \qquad \text{pos\_contrib}_k = \frac{1}{|\mathcal U_N|}\!\!\sum_{b\in\mathcal U_N}\!\!\text{Contrib}_{b,k}, \qquad \text{neg\_contrib}_k = \frac{1}{|\mathcal D_N|}\!\!\sum_{b\in\mathcal D_N}\!\!\text{Contrib}_{b,k}$$

**Input.** $L$, $\beta$, and the FF-up/FF-down bin sets from Step 5.

**Output.** `pc_shap_contrib.csv` — per PC: `pos_contrib`, `neg_contrib`, and their shares of the set's total $|w|$.

**How to interpret.** The per-PC bars for a bin set **sum to that set's mean weight** — an exact decomposition, nothing left over. A PC with a long red (FF-up) bar and near-zero blue bar is a pure fetal-direction builder; one with bars in both directions contributes to both. This answers "of the total FF-up signal, how much is PC3 vs PC7 vs …?"

**Figure — PC ↔ tissue plot.** Panel A (heatmap): y = PCs (labeled `PC · best_tissue`), x = tissue, color = β-oriented alignment (red = FF-up). Panel B (bars): y = same PCs, x = mean contribution $L[b,k]\beta_k$, red = FF-up set, blue = FF-down set. **Combined reading** gives one mechanistic sentence per PC: "PC3 aligns red with placenta *and* has a long FF-up bar → **PC3 is the fetal axis**"; "PC7 aligns blue with CD14 monocyte *and* has a long FF-down bar → **maternal axis**." On the synthetic model (seeded PC3→placenta β=+0.6, PC7→CD14 β=−0.5) the figure reproduces exactly this.

**Script / how to run.** The statistic is computed by `pc_tissue_map.R` (Stage 9, `pc_shap_contrib.csv`); the two-panel figure is drawn by `plot_pc_tissue.py`.
```bash
python plot_pc_tissue.py --indir pc_tissue_out --topk 15 --out pc_tissue_figure.png
```
Arguments: `--indir` (the `pc_tissue_map.R` output dir, default `pc_tissue_out`) · `--topk` (PCs to show, default 15) · `--out` (default `pc_tissue_figure.png`). Reads `pc_tissue_corr.csv` + `pc_shap_contrib.csv` and renders both panels.

**Optional cross-check — maternal-vs-fetal origin (Wang MFI layer).** `pc_origin_map.R` uses the origin-labelled Wang 2026 MFI placenta atlas (`W_` layer, 33 cell types tagged fetal/maternal) to (1) score each PC's fetal−maternal openness oriented by $\operatorname{sign}(\beta_k)$, and (2) call each top bin's dominant cell type and run an FF-direction × origin 2×2 Fisher test — an independent confirmation of the FF-up=fetal / FF-down=maternal sign convention.
```bash
Rscript pc_origin_map.R --model synth_pca_model.rds \
  --atlas reference/wang_mfi_openness_hg19_50kb.csv.gz \
  --manifest reference/manifest_wang_mfi.json --outdir pc_origin_out
```
Arguments: `--model` (required) · `--atlas` (default the Wang MFI table) · `--manifest` (fetal/maternal origin labels, default `reference/manifest_wang_mfi.json`) · `--topn` (default 2000) · `--n-perm` (default 2000) · `--seed` (default 0) · `--outdir` (default `pc_origin_out`). Outputs: `pc_origin_score.csv`, `bin_origin_calls.csv`, `origin_contingency.csv`, `pc_origin_meta.json`. On the synthetic model: PC3 origin +0.49 (fetal, p<1e-3), PC7 maternal; FF-up vs FF-down bins Fisher **OR = 17.5, p < 1e-3**. See `WANG_MFI_LAYER.md` for the full recipe.

---

## Summary — the through-line

| Figure | x-axis | y-axis | question answered |
|---|---|---|---|
| coverage / SHAP | genomic region | coverage / SHAP value | *which bins* drive this prediction, and why |
| tissue bars A | tissue | spec_z (vs null) | are the top bins significantly open in tissue $t$? |
| tissue bars B | tissue | share (normalized) | how big is $t$'s slice of the origin signal? |
| lollipop | spec_z / share | tissue | same, legible at 73 tissues + cutoff stability |
| dumbbell | spec_z | tissue | POS vs NEG → fetal vs maternal contrast |
| PC ↔ tissue A | tissue | PC | which PC resembles which tissue, oriented to FF |
| PC ↔ tissue B | contribution | PC | which PCs build each FF direction (exact) |

The SHAP figure explains **the model's numbers**; the tissue figures translate them into **biology**; the PC figure opens the box to **the latent components**. Two devices keep every step honest: the **null model** (spec_z is enrichment over random genome, not raw openness) and the **β-orientation** (alignment is reported in FF direction, not the PC's arbitrary sign). The single identity that threads them together:

$$w_b = \sum_k L[b,k]\beta_k \ \Rightarrow\ \operatorname{sign}(\text{SHAP}_b) = \operatorname{sign}(w_b) = \text{FF direction} \ \Rightarrow\ \text{cell of origin.}$$
