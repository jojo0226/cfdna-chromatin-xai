"""
engine.py -- chromatin-state annotation and matched-null enrichment.

Given query regions (important / differential cfDNA bins), the engine:
  1. annotates each region with its ChromHMM 18-state composition + dominant state
     in one or more reference epigenomes;
  2. builds a GC / length / N-fraction-matched null background;
  3. tests state / group enrichment of the query set vs. the matched null, with
     fold-enrichment, bootstrap CIs and an empirical p-value.

All coordinates are 0-based half-open, hg38.
"""
from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np

from .references import STATES, I2S, GROUP, ACTIVE
from .genome import region_gc_nfrac


# --------------------------------------------------------------------------- #
#  Annotation                                                                  #
# --------------------------------------------------------------------------- #
def region_state_bp(chrom, start, end, seg_tc):
    """bp overlap of [start,end) with each ChromHMM state in one tissue/chrom.

    seg_tc = segmentation[chrom] for one tissue (sorted non-overlapping segments).
    Returns {state_code: bp}. Vectorised interval slice via searchsorted.
    """
    starts, ends, code = seg_tc["starts"], seg_tc["ends"], seg_tc["code"]
    lo = int(np.searchsorted(ends, start, side="right"))
    hi = int(np.searchsorted(starts, end, side="left"))
    out = {}
    for i in range(lo, hi):
        ov = min(end, int(ends[i])) - max(start, int(starts[i]))
        if ov > 0:
            k = int(code[i])
            out[k] = out.get(k, 0) + ov
    return out


def annotate_region(chrom, start, end, seg_tissue):
    """Full annotation of one region in one tissue.

    Returns dict with dominant_state, composition (state->fraction of covered bp),
    covered_bp, and active_fraction.
    """
    bp = region_state_bp(chrom, start, end, seg_tissue[chrom]) if chrom in seg_tissue else {}
    tot = sum(bp.values())
    if tot == 0:
        return {"dominant_state": None, "composition": {}, "covered_bp": 0, "active_fraction": np.nan}
    comp = {I2S[k]: v / tot for k, v in bp.items()}
    dom = max(comp, key=comp.get)
    active = sum(f for s, f in comp.items() if s in ACTIVE)
    return {"dominant_state": dom, "composition": comp, "covered_bp": tot, "active_fraction": active}


def annotate_regions(regions, seg_tissue):
    """regions: iterable of (chrom, start, end). Returns list of annotation dicts."""
    return [annotate_region(c, s, e, seg_tissue) for (c, s, e) in regions]


# --------------------------------------------------------------------------- #
#  Matched-null background                                                     #
# --------------------------------------------------------------------------- #
@dataclass
class NullModel:
    """GC / length / N-fraction-matched background sampler over the mappable genome."""
    genome: dict
    chroms: list
    max_nfrac: float = 0.10          # reject intervals with >10% N (unmappable)
    gc_tol: float = 0.02             # matched-null GC window (+/-)
    seed: int = 0
    _rng: np.random.Generator = field(default=None, repr=False)

    def __post_init__(self):
        self._rng = np.random.default_rng(self.seed)
        # sampling weights proportional to mappable length per chrom
        self._lens = np.array([self.genome[c]["n"] for c in self.chroms], dtype=float)
        self._w = self._lens / self._lens.sum()

    def _draw_one(self, length, gc_lo, gc_hi, tries=200):
        for _ in range(tries):
            ci = self._rng.choice(len(self.chroms), p=self._w)
            c = self.chroms[ci]
            n = self.genome[c]["n"]
            if n <= length:
                continue
            s = int(self._rng.integers(0, n - length))
            gcf, nf, L = region_gc_nfrac(c, s, s + length, self.genome)
            if nf <= self.max_nfrac and not np.isnan(gcf) and gc_lo <= gcf <= gc_hi:
                return (c, s, s + length)
        return None  # GC-matched draw failed within budget

    def sample_matched(self, query_regions, n_per=100):
        """For each query region draw up to n_per GC/length/N-matched background regions.

        Returns (null_regions, meta) where meta flags regions that could not be matched.
        """
        null, meta = [], []
        for (c, s, e) in query_regions:
            gcf, nf, L = region_gc_nfrac(c, s, e, self.genome)
            if np.isnan(gcf) or L <= 0:
                meta.append({"query": (c, s, e), "matched": 0, "reason": "query undefined GC"})
                continue
            lo, hi = gcf - self.gc_tol, gcf + self.gc_tol
            got = 0
            for _ in range(n_per):
                d = self._draw_one(L, lo, hi)
                if d is not None:
                    null.append(d)
                    got += 1
            meta.append({"query": (c, s, e), "matched": got, "gc": gcf, "len": L})
        return null, meta


# --------------------------------------------------------------------------- #
#  Enrichment test                                                             #
# --------------------------------------------------------------------------- #
def _dominant_counts(regions, seg_tissue, by="state"):
    """Fraction of regions whose dominant label == each state/group."""
    keys = STATES if by == "state" else sorted(set(GROUP.values()))
    counts = {k: 0 for k in keys}
    n = 0
    for (c, s, e) in regions:
        a = annotate_region(c, s, e, seg_tissue)
        dom = a["dominant_state"]
        if dom is None:
            continue
        lab = dom if by == "state" else GROUP[dom]
        counts[lab] += 1
        n += 1
    frac = {k: (counts[k] / n if n else np.nan) for k in keys}
    return counts, frac, n


def enrichment_test(query_regions, seg_tissue, null_model, by="state",
                    n_per=100, n_boot=1000, seed=0):
    """Fold-enrichment of each state/group in query vs. GC/length-matched null.

    Returns a list of per-label result dicts:
        label, query_frac, null_frac, log2_fold, ci_low, ci_high, p_emp, n_query, n_null
    p_emp is a two-sided empirical p from bootstrap resampling of the null (how often
    a null-sized draw reaches the observed |log2 fold|).
    """
    rng = np.random.default_rng(seed)
    null_regions, meta = null_model.sample_matched(query_regions, n_per=n_per)

    q_counts, q_frac, nq = _dominant_counts(query_regions, seg_tissue, by=by)
    n_counts, n_frac, nn = _dominant_counts(null_regions, seg_tissue, by=by)

    # bootstrap the null fraction to get CI + empirical p on the fold-change
    labels = list(q_frac.keys())
    null_labels = []
    for (c, s, e) in null_regions:
        dom = annotate_region(c, s, e, seg_tissue)["dominant_state"]
        null_labels.append(None if dom is None else (dom if by == "state" else GROUP[dom]))
    null_labels = np.array([x for x in null_labels if x is not None], dtype=object)

    results = []
    eps = 0.5 / max(nq, 1)  # pseudocount on the fraction scale
    for lab in labels:
        qf = q_frac[lab]
        nf = n_frac[lab]
        log2fold = np.log2((qf + eps) / (nf + eps))
        # bootstrap null fractions at query sample size
        if len(null_labels) and nq:
            boot = rng.choice(null_labels, size=(n_boot, nq), replace=True)
            bfrac = (boot == lab).mean(axis=1)
            blog2 = np.log2((bfrac + eps) / (nf + eps))
            ci_low, ci_high = np.percentile(blog2, [2.5, 97.5])
            p_emp = (np.sum(np.abs(blog2) >= abs(log2fold)) + 1) / (n_boot + 1)
        else:
            ci_low = ci_high = p_emp = np.nan
        results.append({
            "label": lab, "query_frac": qf, "null_frac": nf,
            "log2_fold": log2fold, "ci_low": ci_low, "ci_high": ci_high,
            "p_emp": p_emp, "n_query": nq, "n_null": nn,
        })
    return results, meta

