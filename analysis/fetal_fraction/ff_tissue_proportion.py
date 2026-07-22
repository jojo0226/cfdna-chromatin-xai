"""ff_tissue_proportion.py -- tissue-of-origin proportion of top-SHAP regions (seqFF++).

The fetal-fraction fork's core question: for the genomic regions a fetal-fraction
model relies on most (top |SHAP| bins), which tissues' open chromatin do they fall in,
and is that concentration real vs a random set of regions?

No blood background is subtracted -- every tissue in the panel is compared on equal
footing (placenta, liver, endothelial, monocyte, Bcell, CD4, CD8, NK, K562). Two
readouts are produced for each top-N cut (default 500 / 1000 / 2000):

  1. ABSOLUTE openness enrichment per tissue
        mean C_<tissue> over the top-N bins, where C_ is the genome-wide z-scored
        combined-openness composite (mean of z(histone), z(accessibility)). Because C_
        is centred genome-wide, a random region set averages ~0; a positive value means
        the top-SHAP regions sit in open chromatin for that tissue.
     Confound: top-SHAP bins tend to be in generally-active chromatin, so most tissues
     rise together. This says "the regions are open", not "which tissue".

  2. TISSUE-SPECIFIC enrichment (cell-of-origin readout)
        spec_<tissue> = C_<tissue> - mean_over_all_tissues(C)   [per bin]
        mean over top-N bins. This removes the pan-openness confound: a positive value
        means the top-SHAP regions are open MORE in that tissue than tissues on average
        -> a cell-of-origin signal. For a genuine FF model, placenta should lead here.

Significance for both: a matched random-region null. For each top-N cut we draw the
same number of bins uniformly at random `n_perm` times, recompute the statistic, and
report z = (obs - null_mean)/null_sd and an empirical two-sided p. The absolute
readout can optionally use an openness-matched null (draw random bins with similar
pan-tissue openness) to test specificity beyond the active-chromatin confound.

Tissue "proportion/contribution" is the normalized share of positive tissue-specific
enrichment across the panel -- a simplex that answers "what mix of tissues explains the
top-SHAP regions", the direct seqFF++ analog of a cell-type deconvolution.

Inputs
------
  - a per-bin importance table with a genomic key (chrom/start/end or 'chr:start-end'
    index) and a |SHAP| column, and
  - the FF openness atlas (build_openness_atlas.py output; hg19 50kb, C_/ACC_/D_/H_).
Both are hg19-native; no liftover of your cohort ever happens (build-native rule).
"""
import numpy as np
import pandas as pd

DEF_TOPNS = (500, 1000, 2000)


def _norm_key(s: pd.Series) -> pd.Series:
    """Normalise a bin key to 'chrom:start-end' (accepts chr:start:end too)."""
    return s.astype(str).str.replace(r"^(chr[\w]+):(\d+):(\d+)$", r"\1:\2-\3", regex=True)


def atlas_key_from_coords(chrom, start, end) -> str:
    return f"{chrom}:{int(start)}-{int(end)}"


def load_importance(imp: pd.DataFrame, shap_col=None, key_col=None) -> pd.DataFrame:
    """Return a frame with columns ['key','importance'] from a flexible importance table.

    - key: from `key_col`, or a 'key'/'bin'/'interval' column, or the index.
    - importance: from `shap_col`, or the first of
      ['mean_abs_shap','importance','abs_shap','shap'].
    """
    df = imp.copy()
    if key_col and key_col in df.columns:
        key = df[key_col]
    elif "key" in df.columns:
        key = df["key"]
    elif {"chrom", "start", "end"}.issubset(df.columns):
        key = df.apply(lambda r: atlas_key_from_coords(r["chrom"], r["start"], r["end"]),
                       axis=1)
    elif "interval" in df.columns:
        key = df["interval"]
    elif "bin" in df.columns:
        key = df["bin"]
    else:
        key = df.index.to_series()
    if shap_col is None:
        for c in ("mean_abs_shap", "importance", "abs_shap", "shap", "mean_shap"):
            if c in df.columns:
                shap_col = c
                break
    if shap_col is None:
        raise KeyError("no importance column found; pass shap_col=")
    imp_val = df[shap_col].abs().values
    return pd.DataFrame({"key": _norm_key(pd.Series(key)).values, "importance": imp_val})


def _tissue_cols(atlas: pd.DataFrame, prefix="C_"):
    """Panel tissue columns for a given track prefix, excluding reference contrasts."""
    skip_suffix = ("_specific", "infl_specific")
    skip_exact = ("background", "immune", "tumor", "signal")
    out = []
    for c in atlas.columns:
        if not c.startswith(prefix):
            continue
        name = c[len(prefix):]
        if any(name.endswith(s) or name == s for s in skip_suffix):
            continue
        if name in skip_exact:
            continue
        out.append(name)
    return out


def build_specificity(atlas: pd.DataFrame, prefix="C_", tissues=None) -> pd.DataFrame:
    """Add spec_<tissue> = C_<tissue> - row-mean over the panel tissues.

    Returns a frame indexed by atlas key with C_<t>, spec_<t> for each panel tissue.
    """
    if tissues is None:
        tissues = _tissue_cols(atlas, prefix)
    C = atlas[[f"{prefix}{t}" for t in tissues]].copy()
    C.columns = tissues
    rowmean = C.mean(axis=1)
    spec = C.sub(rowmean, axis=0)
    spec.columns = [f"spec_{t}" for t in tissues]
    out = pd.concat([C.add_prefix("C_"), spec], axis=1)
    out.index = atlas["key"].values if "key" in atlas.columns else atlas.index
    return out, tissues


def tissue_proportion(imp: pd.DataFrame, atlas: pd.DataFrame,
                      topns=DEF_TOPNS, prefix="C_", tissues=None,
                      n_perm=2000, seed=0, shap_col=None, key_col=None):
    """Core analysis: per-tissue absolute + specific enrichment of top-N SHAP bins
    vs a matched random-region null, with a normalized cell-of-origin proportion.

    Returns (results_df, proportions_df, meta):
      results_df rows = (topn, tissue) with columns
        abs_obs, abs_z, abs_p, spec_obs, spec_z, spec_p, n_bins_used
      proportions_df: index tissue, columns per topn -> share of positive spec_obs.
    """
    rng = np.random.default_rng(seed)
    S, tissues = build_specificity(atlas, prefix, tissues)
    C = S[[f"C_{t}" for t in tissues]].values      # genome-wide, ~z-scored
    SP = S[[f"spec_{t}" for t in tissues]].values   # tissue-specific
    keys = np.array(S.index)
    keypos = {k: i for i, k in enumerate(keys)}

    ii = load_importance(imp, shap_col=shap_col, key_col=key_col)
    ii = ii[ii["key"].isin(keypos)].copy()
    ii["pos"] = ii["key"].map(keypos)
    ii = ii.sort_values("importance", ascending=False)
    n_avail = len(ii)

    rows = []
    props = {}
    Ngrid = C.shape[0]
    for topn in topns:
        k = min(topn, n_avail)
        sel = ii["pos"].values[:k]
        abs_obs = np.nanmean(C[sel], axis=0)     # per tissue
        spec_obs = np.nanmean(SP[sel], axis=0)
        # matched random-region null: draw k random bins n_perm times
        abs_null = np.empty((n_perm, len(tissues)))
        spec_null = np.empty((n_perm, len(tissues)))
        for b in range(n_perm):
            r = rng.choice(Ngrid, k, replace=False)
            abs_null[b] = np.nanmean(C[r], axis=0)
            spec_null[b] = np.nanmean(SP[r], axis=0)
        abs_mu, abs_sd = abs_null.mean(0), abs_null.std(0) + 1e-12
        spec_mu, spec_sd = spec_null.mean(0), spec_null.std(0) + 1e-12
        abs_z = (abs_obs - abs_mu) / abs_sd
        spec_z = (spec_obs - spec_mu) / spec_sd
        # two-sided empirical p
        abs_p = (np.abs(abs_null - abs_mu) >= np.abs(abs_obs - abs_mu)[None]).mean(0)
        spec_p = (np.abs(spec_null - spec_mu) >= np.abs(spec_obs - spec_mu)[None]).mean(0)
        for j, t in enumerate(tissues):
            rows.append(dict(topn=topn, tissue=t,
                             abs_obs=abs_obs[j], abs_z=abs_z[j], abs_p=abs_p[j],
                             spec_obs=spec_obs[j], spec_z=spec_z[j], spec_p=spec_p[j],
                             n_bins_used=k))
        pos = np.clip(spec_obs, 0, None)
        props[topn] = pos / (pos.sum() or 1.0)
    results = pd.DataFrame(rows)
    proportions = pd.DataFrame(props, index=tissues)
    meta = dict(n_bins_in_importance=len(ii), n_bins_matched=n_avail,
                n_atlas_bins=Ngrid, tissues=tissues, prefix=prefix, n_perm=n_perm)
    return results, proportions, meta
