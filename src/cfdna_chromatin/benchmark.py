"""
benchmark.py -- B1 statistical-calibration harness for the enrichment engine.

B1 is the gate every downstream claim inherits: feed the engine query sets that
carry NO real signal (regions drawn from the same matched-null process as the
background) and verify the engine is calibrated --

  * empirical p-values are UNIFORM on [0,1] under the null (KS test vs Uniform);
  * the realized false-positive rate at nominal alpha matches alpha;
  * a QQ-plot of observed vs expected -log10(p) hugs the diagonal.

If B1 fails, the engine is mis-calibrated and B2/B3/B4 biology is uninterpretable.
"""
from __future__ import annotations
import numpy as np
from scipy import stats

from .engine import enrichment_test


def sample_null_query_sets(null_model, n_sets, set_size, region_len, seed=0):
    """Draw n_sets 'fake' query sets (each set_size regions of region_len bp) from
    the mappable genome with no signal -- the null-of-the-null for B1."""
    rng = np.random.default_rng(seed)
    sets = []
    for k in range(n_sets):
        # sample uniformly mappable regions of fixed length (no GC constraint here:
        # these ARE the null hypothesis -- their enrichment vs matched null should be ~0)
        regs = []
        tries = 0
        while len(regs) < set_size and tries < set_size * 50:
            tries += 1
            ci = rng.choice(len(null_model.chroms), p=null_model._w)
            c = null_model.chroms[ci]
            n = null_model.genome[c]["n"]
            if n <= region_len:
                continue
            s = int(rng.integers(0, n - region_len))
            from .genome import region_gc_nfrac
            gcf, nf, L = region_gc_nfrac(c, s, s + region_len, null_model.genome)
            if nf <= null_model.max_nfrac and not np.isnan(gcf):
                regs.append((c, s, s + region_len))
        sets.append(regs)
    return sets


def run_b1_calibration(null_model, seg_tissue, n_sets=100, set_size=50,
                       region_len=50_000, by="group", n_per=50, n_boot=500, seed=0):
    """Run B1: collect one empirical p per (null query set x label) and test uniformity.

    Returns dict with per-label KS statistic/p vs Uniform, realized FPR at alpha=0.05,
    and the raw p-value matrix for QQ plotting.
    """
    qsets = sample_null_query_sets(null_model, n_sets, set_size, region_len, seed=seed)
    p_by_label = {}
    for i, qs in enumerate(qsets):
        res, _ = enrichment_test(qs, seg_tissue, null_model, by=by,
                                 n_per=n_per, n_boot=n_boot, seed=seed + i + 1)
        for r in res:
            p_by_label.setdefault(r["label"], []).append(r["p_emp"])
    summary = {}
    for lab, ps in p_by_label.items():
        ps = np.array([p for p in ps if np.isfinite(p)])
        if len(ps) < 5:
            continue
        ks = stats.kstest(ps, "uniform")
        summary[lab] = {
            "n": int(len(ps)),
            "ks_stat": float(ks.statistic),
            "ks_p": float(ks.pvalue),
            "fpr_0.05": float(np.mean(ps < 0.05)),
            "mean_p": float(np.mean(ps)),
        }
    return {"summary": summary, "p_by_label": {k: list(map(float, v)) for k, v in p_by_label.items()}}
