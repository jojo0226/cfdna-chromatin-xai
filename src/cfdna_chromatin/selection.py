"""
selection.py -- importance / selection front-end for the cfDNA explainer.

Stage 1 of the three-stage architecture (ANY per-region feature -> IMPORTANCE
-> chromatin/GO EXPLANATION). It is feature-agnostic and model-agnostic: it
turns a per-region *importance vector* into a ranked region list ready for
attribute.attribute_signal(), and answers the key interpretation question:
does the importance concentrate in signal-specific (cell-of-origin) or
background-specific chromatin, once the genome-wide openness confound is removed?

Two ways to produce importance:

    rank_by_shap(shap_df)          -- from a SHAP value matrix (samples x bins),
                                      e.g. TabPFN + KernelSHAP back-projected to
                                      genomic bins. Importance = mean|SHAP|;
                                      direction = mean SHAP (sign convention is
                                      the model's positive class).
    rank_by_differential(mat, y)   -- from a feature matrix + binary labels, via
                                      a per-bin Mann-Whitney U test. Importance =
                                      -log10(p); direction = group mean diff.

Both return a tidy DataFrame indexed by bin name with columns
(chrom, start, end, importance, direction, rank). Feed the top rows to
to_regions() -> attribute_signal(), or pass the whole frame to
compartment_importance_test() for the genomic-null-corrected readout.
"""
from __future__ import annotations
import re

import numpy as np
import pandas as pd
from scipy import stats

from . import histone as H
from . import attribute as AT

_BIN_RE = re.compile(r"(chr[\w]+)[:\-_](\d+)[\-_](\d+)")


def parse_bin(name):
    """'chr19:50000000-50500000' -> ('chr19', 50000000, 50500000)."""
    m = _BIN_RE.search(str(name))
    if not m:
        raise ValueError(f"cannot parse genomic bin from {name!r}")
    return m.group(1), int(m.group(2)), int(m.group(3))


def _coord_frame(names):
    rows = [parse_bin(n) for n in names]
    return pd.DataFrame(rows, columns=["chrom", "start", "end"], index=list(names))


def rank_by_shap(shap_df, samples_axis=0):
    """Rank bins by mean|SHAP| from a SHAP matrix.

    shap_df : DataFrame with bins on one axis and samples on the other.
              samples_axis=0 -> rows are samples, columns are bins (default,
              matches the shap_values_*_original_space.csv layout).
    Returns a DataFrame indexed by bin, sorted by importance (descending).
    """
    df = shap_df if samples_axis == 0 else shap_df.T
    importance = df.abs().mean(axis=0)
    direction = df.mean(axis=0)
    out = _coord_frame(importance.index)
    out["importance"] = importance.values
    out["direction"] = direction.values
    out = out.sort_values("importance", ascending=False)
    out["rank"] = np.arange(1, len(out) + 1)
    return out


def shap_from_model(model, X, background=None, explainer="auto",
                    nsamples=200, random_state=42):
    """Compute a per-region SHAP matrix from *your own* fitted FF model.

    This is the ingestion point for the common case where the model already
    exists (e.g. an in-house seqFF++ fetal-fraction regressor) and you only
    want the explanation stages. It returns a SHAP DataFrame (samples x
    features) laid out exactly like rank_by_shap() expects, so the whole
    Stage 1 -> Stage 3 chain runs on your model's attributions.

    Parameters
    ----------
    model : a fitted estimator. Anything with ``predict`` (regressor, e.g.
            FF in [0,1]) or ``predict_proba`` (classifier). Tree models use
            the exact TreeExplainer; linear models use LinearExplainer; anything
            else falls back to the model-agnostic KernelExplainer.
    X     : DataFrame (n_samples x n_features). Columns MUST be genomic-bin
            names parseable by parse_bin() ('chrX:start-end') so the resulting
            importance is region-anchored. Non-region feature columns are kept
            in the SHAP frame but will fail downstream coordinate parsing --
            drop them first if the model mixes region and non-region features.
    background : reference dataset for the explainer's expected value. Defaults
            to shap.kmeans(X, 10) for Kernel/Linear explainers; ignored by
            TreeExplainer. Keep it small -- KernelExplainer cost scales with it.
    explainer : 'auto' | 'tree' | 'linear' | 'kernel'. 'auto' picks tree for
            tree ensembles, linear for linear_model estimators, else kernel.

    Returns
    -------
    DataFrame (n_samples x n_features), signed SHAP values, same index/columns
    as X. Feed straight to rank_by_shap().

    Notes
    -----
    SHAP sign is in the model's output units (higher FF for a regressor). Keep
    the sign: a fetal-fraction signal is directional, and collapsing to
    mean|SHAP| too early hides which regions push FF up vs down (a sign-trap
    that collapsing to magnitude too early would introduce).
    """
    import shap as _shap

    Xv = X.values if hasattr(X, "values") else np.asarray(X)
    cols = list(X.columns) if hasattr(X, "columns") else list(range(Xv.shape[1]))
    idx = list(X.index) if hasattr(X, "index") else list(range(Xv.shape[0]))

    def _predict(data):
        if hasattr(model, "predict_proba"):
            proba = model.predict_proba(data)
            return proba[:, 1] if getattr(proba, "ndim", 1) == 2 else proba
        return model.predict(data)

    kind = explainer
    if kind == "auto":
        mod = type(model).__module__.lower()
        name = type(model).__name__.lower()
        if any(k in mod or k in name for k in
               ("xgboost", "lightgbm", "catboost", "forest", "tree", "gradientboost")):
            kind = "tree"
        elif "linear_model" in mod or name in ("lasso", "ridge", "elasticnet",
                                               "linearregression", "logisticregression"):
            kind = "linear"
        else:
            kind = "kernel"

    if kind == "tree":
        expl = _shap.TreeExplainer(model)
        sv = expl.shap_values(Xv)
    elif kind == "linear":
        # LinearExplainer needs a raw background dataset (builds an Independent
        # masker); shap.kmeans DenseData is not an accepted masker here.
        bg = Xv if background is None else (
            background.data if hasattr(background, "data") else np.asarray(background))
        expl = _shap.LinearExplainer(model, bg)
        sv = expl.shap_values(Xv)
    else:
        bg = _shap.kmeans(Xv, min(10, len(Xv))) if background is None else background
        expl = _shap.KernelExplainer(_predict, bg)
        sv = expl.shap_values(Xv, nsamples=nsamples, l1_reg="num_features(10)")

    if isinstance(sv, list):          # some explainers return per-class lists
        sv = sv[-1]
    sv = np.asarray(sv)
    if sv.ndim == 3:                  # (samples, features, classes)
        sv = sv[:, :, -1]
    return pd.DataFrame(sv, index=idx, columns=cols)


def rank_by_differential(matrix, labels, positive=None):
    """Rank bins by a per-bin two-sided Mann-Whitney U test.

    matrix : DataFrame (samples x bins).
    labels : array-like of length n_samples, exactly two distinct values.
    positive : which label is the 'positive' group for the direction sign
               (direction = mean_positive - mean_other). Defaults to the
               second sorted label.
    """
    labels = np.asarray(labels)
    groups = sorted(pd.unique(labels[~pd.isnull(labels)]))
    if len(groups) != 2:
        raise ValueError(f"need exactly 2 label groups, got {groups}")
    if positive is None:
        positive = groups[-1]
    other = [g for g in groups if g != positive][0]
    A = matrix[labels == positive]
    B = matrix[labels == other]
    imp, direction = [], []
    for b in matrix.columns:
        a, c = A[b].values, B[b].values
        try:
            _, p = stats.mannwhitneyu(a, c, alternative="two-sided")
        except ValueError:
            p = 1.0
        imp.append(-np.log10(max(p, 1e-300)))
        direction.append(np.mean(a) - np.mean(c))
    out = _coord_frame(matrix.columns)
    out["importance"] = imp
    out["direction"] = direction
    out["p"] = np.power(10.0, -np.asarray(imp))
    out = out.sort_values("importance", ascending=False)
    out["rank"] = np.arange(1, len(out) + 1)
    return out


def to_regions(ranked, top_n=None, top_frac=None):
    """Return [(chrom,start,end), ...] for the top rows of a ranked frame."""
    df = ranked
    if top_frac is not None:
        top_n = int(np.ceil(len(df) * top_frac))
    if top_n is not None:
        df = df.head(top_n)
    return list(df[["chrom", "start", "end"]].itertuples(index=False, name=None))


def _compartment_activity(ranked, hpanel, panel_name):
    """Per-bin active-mark coverage for signal / background / inflammation tissues."""
    from . import references as R2
    roles = {t: R2.tissue_role(panel_name, t) for t in R2.panel_tissues(panel_name)}
    sig = [t for t, r in roles.items() if r == "signal"]
    non_sig = [t for t, r in roles.items() if r in ("background", "inflammation")]
    rows = {}
    for name, r in ranked.iterrows():
        act = {}
        for t in sig + non_sig:
            fp = H.region_mark_fingerprint(r["chrom"], int(r["start"]), int(r["end"]),
                                           hpanel[t], marks=AT.ACTIVE_MARKS)
            act[t] = float(np.mean([fp["coverage"][m] for m in AT.ACTIVE_MARKS]))
        signal = np.mean([act[t] for t in sig]) if sig else 0.0
        background = np.mean([act[t] for t in non_sig]) if non_sig else 0.0
        rows[name] = {"signal_activity": signal, "background_activity": background,
                      "signal_specific": signal - background}
    return pd.DataFrame(rows).T


def compartment_importance_test(ranked, hpanel, panel_name, use_abs_direction=False):
    """Genomic-null-corrected interpretation of an importance vector.

    Correlates per-bin importance with signal-specific vs background-specific
    chromatin openness (Spearman). A positive correlation with signal_specific
    means the model/feature weight concentrates in the cell-of-origin (signal)
    chromatin; a positive correlation with background_specific means it
    concentrates in the hematopoietic background compartment.

    Reporting BOTH the raw-openness and the specific-openness correlations
    exposes a common confound: importance often correlates positively with
    *total* openness (any active chromatin), and only the compartment-*specific*
    contrast reveals cell-of-origin.

    Returns dict with per-target Spearman rho/p and the joined per-bin frame.
    """
    act = _compartment_activity(ranked, hpanel, panel_name)
    J = ranked.join(act)
    score = J["direction"].abs() if use_abs_direction else J["importance"]
    J = J.assign(_score=score.values)
    res = {}
    for comp in ("signal_specific", "signal_activity", "background_activity"):
        rho, p = stats.spearmanr(J["_score"], J[comp])
        res[comp] = {"rho": float(rho), "p": float(p)}
    # background_specific is the mirror of signal_specific
    res["background_specific"] = {"rho": -res["signal_specific"]["rho"],
                                  "p": res["signal_specific"]["p"]}
    res["n"] = int(len(J))
    res["panel"] = panel_name
    res["score"] = "abs_direction" if use_abs_direction else "importance"
    return res, J
