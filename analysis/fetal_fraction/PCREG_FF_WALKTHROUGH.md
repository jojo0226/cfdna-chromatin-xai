# From coverage file to tissue contribution — the `pcregFF` pipeline, step by step

How a fetal-fraction (FF) prediction is produced by the PCA→linear-regression
("pcreg") model, and how the XAI layer opens it up: per-bin SHAP importance, and
the decomposition of each principal component into tissue open-chromatin
contributions.

The whole point that makes the math exact: **both stages of the model are linear**,
so the two stacked maps collapse into a *single* linear model in genomic-bin space.
Nothing below is an approximation — the SHAP values are closed form, not sampled.

Symbols used throughout:

| symbol | shape | meaning |
|--------|-------|---------|
| $x$ | $B$ | one sample's per-bin coverage feature vector (after normalization) |
| $\bar{x}$ | $B$ | training-set mean of each bin (the PCA centering vector) |
| $L$ | $B \times K$ | `bin_loadings` — PCA loadings (each column = one PC's genome-wide pattern) |
| $\beta$ | $K$ | `beta` — linear-regression coefficients on the PC scores |
| $b_0$ | scalar | regression intercept |
| $C$ | $B \times T$ | tissue openness atlas (`C_<tissue>` columns) |
| $B$ | 52,459 | model bins (production `.rds`) |
| $K$ | 100 | principal components |
| $T$ | 73 / 6 | tissues (DNase atlas) / cell types (placental ATAC layer) |

Production model object (`ff_model.rds`): `bin_params` (52,459×4), `bin_loadings`
(52,459 × [bin_name + PC1..PC100]), `pc_sdev` (100), `beta` (101 = intercept +
PC1..PC100), `train_dt` (2,477 samples × 108; contains private `FF_Yplus`). Training
fit: Pearson $r=0.854$, $R^2=0.73$, RMSE $=0.0373$.

---

## Step 0 — Input: the coverage file

The pipeline starts from a **per-bin cfDNA coverage matrix**: rows = samples, columns
= genomic bins (50 kb hg19). Entry = number of cfDNA fragments (or fragment ends)
whose midpoint falls in that bin, from the aligned BAM.

$$\text{raw}_{i,b} = \#\{\text{fragments of sample } i \text{ in bin } b\}$$

Raw counts are not comparable across bins or samples: they scale with sequencing
depth, and per-bin capture is biased by GC content and mappability. Step 1 removes
those nuisances.

---

## Step 1 — Normalize coverage → the model feature $x$

Each raw bin count is turned into a **normalized coverage residual**, stored per-bin
as `bin_params` (52,459 × 4: the per-bin normalization constants — baseline level,
GC/mappability terms, dispersion). Conceptually:

$$
x_{i,b} \;=\; \frac{\text{raw}_{i,b}/D_i}{\widehat{\mu}_b(\text{GC}_b,\,\text{map}_b)}
\;-\; 1
$$

- $D_i$ — sample depth (total usable fragments), so bins are depth-normalized.
- $\widehat{\mu}_b(\cdot)$ — the expected normalized coverage for bin $b$ given its GC
  and mappability, learned across the training cohort (the GC/negative-binomial
  correction). Dividing by it flattens the GC "smile."
- Subtracting 1 makes $x_{i,b}$ a **residual**: 0 = the bin behaves as expected;
  $>0$ = over-covered; $<0$ = under-covered relative to the GC/mappability model.

**This $x$ is the SHAP feature.** It is the quantity every importance number below is
expressed in. (In `extract_ff_shap.R` this matrix is the optional `--x` input; the
per-bin normalization lives in `bin_params`.)

> **Why openness shows up as coverage *depletion*.** In cfDNA, open/accessible
> chromatin is nucleosome-poor → fragmented more → **fewer** intact fragments map
> there. So a tissue that is highly *open* at bin $b$ tends to push $x_{i,b}$ **down**
> in samples where that tissue contributes cfDNA. Keep this sign in mind for Step 5.

---

## Step 2 — Center on the training mean

PCA operates on centered features. The centering vector $\bar{x}$ is the per-bin
training mean and is frozen into the model:

$$\tilde{x}_{i,b} = x_{i,b} - \bar{x}_b$$

$\bar{x}$ is the *only* piece of the training feature distribution the SHAP step needs
(it is the SHAP baseline — see Step 4). No per-patient data is required downstream.

---

## Step 3 — PCA scores, then linear regression on the scores

**3a. PC scores.** Project the centered sample onto each loading vector:

$$
s_{i,k} \;=\; \sum_{b=1}^{B} L_{b,k}\,\tilde{x}_{i,b}
\;=\; \sum_{b=1}^{B} L_{b,k}\,(x_{i,b}-\bar{x}_b),
\qquad k = 1,\dots,K
$$

$s_{i,k}$ is sample $i$'s coordinate on PC $k$. $L_{\cdot,k}$ (a column of
`bin_loadings`) is a **genome-wide spatial pattern** — a weighting of all 52,459 bins.
`pc_sdev` gives each PC's standard deviation across training samples; its square over
the total is the **variance fraction** each PC explains.

**3b. Regression.** Predicted FF is linear in the scores:

$$
\widehat{FF}_i \;=\; b_0 + \sum_{k=1}^{K} \beta_k\, s_{i,k}
$$

with $\beta$ = `beta`. That is the model's output. Everything so far is two linear
maps in a row.

---

## Step 4 — Collapse to one linear model, and read off SHAP

Substitute 3a into 3b and swap the order of summation:

$$
\widehat{FF}_i
= b_0 + \sum_{k} \beta_k \sum_{b} L_{b,k}(x_{i,b}-\bar{x}_b)
= b_0 + \sum_{b} \underbrace{\Big(\sum_{k} L_{b,k}\,\beta_k\Big)}_{\displaystyle w_b}\,(x_{i,b}-\bar{x}_b)
$$

So the entire pipeline is a **single linear model in bin space** with an
**effective per-bin weight**

$$
\boxed{\,w_b \;=\; \sum_{k=1}^{K} L_{b,k}\,\beta_k\,}
\qquad\Longleftrightarrow\qquad
\mathbf{w} = L\,\beta \quad (B\text{-vector} = (B\times K)\cdot(K)).
$$

This is the one line of code at the heart of the XAI layer
(`w <- as.numeric(Lmat %*% bpc)` in `extract_ff_shap.R`).

### SHAP is closed form

For a linear model $f(x)=b_0+\sum_b w_b(x_b-\bar{x}_b)$ with the interventional
(marginal) baseline $\mathbb{E}[x]=\bar{x}$, the exact Shapley value of feature $b$ for
sample $i$ is simply that feature's own term:

$$
\boxed{\,\phi_{i,b} \;=\; w_b\,(x_{i,b}-\bar{x}_b)\,}
$$

No sampling, no KernelSHAP — for a linear model the Shapley value *is* the centered
contribution. It satisfies local accuracy exactly:

$$
\sum_{b}\phi_{i,b} \;=\; \widehat{FF}_i - b_0 \;=\; \widehat{FF}_i - \mathbb{E}[\widehat{FF}].
$$

**Two importance tracks** (columns of the extractor's output CSV):

$$
\text{signed\_shap}_b = w_b \cdot \operatorname{mean}_i|x_{i,b}-\bar{x}_b|,
\qquad
\text{mean\_abs\_shap}_b = |w_b|\cdot \operatorname{mean}_i|x_{i,b}-\bar{x}_b|
$$

- The **sign** is carried entirely by $w_b$ (the mean-absolute-deviation factor is
  $\ge 0$). $w_b>0$: more coverage at bin $b$ ⇒ higher predicted FF; $w_b<0$: ⇒ lower.
- With **no** sample matrix supplied, $|w_b|$ alone is the model-level importance
  track — it needs nothing but the model object (`imp_mode = pca_weight_only`). With
  the matrix, you get the exact SHAP magnitude (`pca_linear_shap_exact`).

> **A PC has no genomic location.** SHAP is defined in *bin* space, not on the PCs — a
> PC is a whole-genome linear combination, so "the SHAP of PC$k$" has no coordinate to
> point at. The right question for a PC is Step 5.

The signed per-bin weights $w_b$ are what the tissue-of-origin step (`run_ff_tissue_
track.py`) reads: it asks whether the FF-up bins ($w_b>0$) and FF-down bins ($w_b<0$)
are enriched for each tissue's open chromatin, via the null-model permutation test.

---

## Step 5 — PC → tissue contribution (`pc_tissue_map.R`)

This is the "which biology does each PC encode" step. It has two exact computations.

### 5a. PC ↔ tissue alignment

Each loading column $L_{\cdot,k}$ is a genome-wide pattern; each tissue column
$C_{\cdot,t}$ is that tissue's open-chromatin pattern over the same bins. Their Pearson
correlation over the shared bins measures how much PC $k$'s spatial pattern *resembles*
tissue $t$'s open chromatin:

$$
\rho_{k,t} \;=\; \operatorname{corr}\!\big(L_{\cdot,k},\, C_{\cdot,t}\big)
\;=\; \frac{\sum_b (L_{b,k}-\bar L_k)(C_{b,t}-\bar C_t)}
{\sqrt{\sum_b (L_{b,k}-\bar L_k)^2}\,\sqrt{\sum_b (C_{b,t}-\bar C_t)^2}}
$$

**Orient by the sign of $\beta_k$** so the number reads in the direction the PC pushes
FF, not the arbitrary sign PCA assigned to the eigenvector:

$$
\boxed{\,\rho^{\text{orient}}_{k,t} \;=\; \operatorname{sign}(\beta_k)\,\rho_{k,t}\,}
$$

Interpretation: $\rho^{\text{orient}}_{k,t} > 0$ ⇒ *in the direction PC $k$ pushes
predicted FF up, its genome-wide pattern aligns with tissue $t$'s open chromatin.*
`pc_tissue_corr.csv` stores one row per PC: `pc, beta, var_frac, best_tissue`, then a
`corr_<tissue>` column for every tissue.

Run it once against the **DNase** atlas (`--prefix C_`, 73 tissues) and once against the
**placental ATAC** layer (`--prefix A_`, 6 cell types). A PC that aligns with placenta
in **both** assays is a cross-modally supported fetal-origin component — that is the
payoff of keeping the ATAC layer as an independent track.

### 5b. Per-PC contribution to each SHAP direction

From Step 4, each PC contributes additively to every bin's effective weight:

$$
w_b = \sum_k \underbrace{L_{b,k}\,\beta_k}_{\displaystyle \text{contrib}_{b,k}}
$$

So $\text{contrib}_{b,k}=L_{b,k}\beta_k$ is **PC $k$'s exact share of bin $b$'s weight.**
Average it over the top-$N$ most FF-up bins ($\mathcal{P}$: largest $w_b>0$) and the
top-$N$ most FF-down bins ($\mathcal{N}$: most negative $w_b$):

$$
\text{pos\_contrib}_k = \frac{1}{|\mathcal{P}|}\sum_{b\in\mathcal{P}} L_{b,k}\beta_k,
\qquad
\text{neg\_contrib}_k = \frac{1}{|\mathcal{N}|}\sum_{b\in\mathcal{N}} L_{b,k}\beta_k
$$

and, as a share of each set's total weight,
$\text{pos\_frac}_k = \big(\sum_{b\in\mathcal P}L_{b,k}\beta_k\big)/\sum_{b\in\mathcal P}w_b$.
These say **which PCs assemble the fetal-leaning vs maternal-leaning bin sets**
(`pc_shap_contrib.csv`). The leading PC of $\mathcal{P}$ is the model's main
"FF-up" axis; cross-referencing it with 5a names the tissue that axis resembles.

---

## The chain, end to end

$$
\text{BAM} \xrightarrow{\text{count}} \text{raw}_{i,b}
\xrightarrow[\text{bin\_params}]{\text{GC/depth norm}} x_{i,b}
\xrightarrow{-\,\bar x} \tilde x_{i,b}
\xrightarrow{L^\top} s_{i,k}
\xrightarrow{\beta} \widehat{FF}_i
$$

$$
\underbrace{\phi_{i,b}=w_b(x_{i,b}-\bar x_b)}_{\text{SHAP, bin space}}
\quad\Big|\quad
\underbrace{w_b=\textstyle\sum_k L_{b,k}\beta_k}_{\text{effective weight}}
\quad\Big|\quad
\underbrace{\rho^{\text{orient}}_{k,t}=\operatorname{sign}(\beta_k)\operatorname{corr}(L_{\cdot,k},C_{\cdot,t})}_{\text{PC}\to\text{tissue}}
$$

**Privacy note.** Everything that leaves the AWS instance — $w_b$, the SHAP tracks,
the PC↔tissue tables — is an aggregate over bins/PCs/tissues. The per-sample matrix
$x$ and the private label `FF_Yplus` never leave; the XAI layer reads only
`bin_loadings`, `beta`, `pc_sdev` (model internals) and the public openness atlas.

---

### Scripts implementing each step

| step | script | key output |
|------|--------|-----------|
| 1–4 SHAP | `analysis/fetal_fraction/extract_ff_shap.R` | `key, mean_abs_shap, signed_shap` |
| tissue enrichment of SHAP bins | `analysis/fetal_fraction/run_ff_tissue_track.py` | per-tissue enrichment + permutation p |
| 5 PC→tissue | `analysis/fetal_fraction/pc_tissue_map.R` | `pc_tissue_corr.csv`, `pc_shap_contrib.csv` |
| atlas (DNase) | `analysis/fetal_fraction/build_openness_atlas.py` | `ff_openness_atlas_hg19_50kb*.csv.gz` |
| atlas (placental ATAC) | `analysis/fetal_fraction/build_placenta_atac_atlas.py` | `ff_placenta_atac_hg19_50kb.csv.gz` |
