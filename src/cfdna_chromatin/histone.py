"""
histone.py -- histone-mark fingerprint layer.

An orthogonal, state-independent cross-check on the ChromHMM call. For each
region we compute, per histone mark, the fraction of the region covered by
ENCODE ChIP-seq peaks (and the peak-weighted mean signal). Six marks form a
per-region fingerprint whose expected pattern is well known:

    H3K4me3   sharp    active promoters        (TssA / TssFlnk)
    H3K4me1   broad    enhancers (poised+active)
    H3K27ac   broad    ACTIVE enhancers/promoters
    H3K36me3  broad    transcribed gene bodies (Tx)
    H3K27me3  broad    Polycomb repression     (ReprPC / Bivalent)
    H3K9me3   broad    constitutive heterochromatin (Het)

This gives a chromatin-mark readout that does not depend on the ChromHMM
segmentation, so agreement between the two is genuine corroboration rather than
circularity.

Peak files: ENCODE Histone ChIP-seq, GRCh38, one experiment per (tissue, mark),
subset to the working chromosomes. Columns: chrom, start, end, signalValue
(0-based half-open). See data/histone/manifest.json for accessions.
"""
from __future__ import annotations
import gzip
import json
import os
from collections import defaultdict

import numpy as np

# canonical mark order (promoter -> enhancer -> transcription -> repressive)
MARKS = ["H3K4me3", "H3K4me1", "H3K27ac", "H3K36me3", "H3K27me3", "H3K9me3"]

# which marks are "active" (used for an active-vs-repressive fingerprint score)
ACTIVE_MARKS = {"H3K4me3", "H3K4me1", "H3K27ac", "H3K36me3"}
REPRESSIVE_MARKS = {"H3K27me3", "H3K9me3"}


def load_peaks(bed_gz_path, chroms=None):
    """Load one narrowPeak-derived .bed.gz into per-chromosome sorted arrays.

    Expects 4 columns: chrom, start, end, signalValue. Returns
    {chrom: {"starts": int64[], "ends": int64[], "signal": float64[]}} sorted by
    start. Peaks may overlap slightly; overlap accounting handles that.
    """
    tmp = defaultdict(list)
    keep = set(chroms) if chroms is not None else None
    with gzip.open(bed_gz_path, "rt") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            c = p[0]
            if keep is not None and c not in keep:
                continue
            s, e = int(p[1]), int(p[2])
            sig = float(p[3]) if len(p) > 3 and p[3] not in ("", ".") else 1.0
            tmp[c].append((s, e, sig))
    out = {}
    for c, rows in tmp.items():
        rows.sort()
        starts = np.array([r[0] for r in rows], dtype=np.int64)
        ends = np.array([r[1] for r in rows], dtype=np.int64)
        sig = np.array([r[2] for r in rows], dtype=np.float64)
        out[c] = {"starts": starts, "ends": ends, "signal": sig}
    return out


def load_mark_panel(data_dir, chroms=None, manifest="manifest.json"):
    """Load every peak file listed in data/histone/manifest.json.

    Returns nested dict panel[tissue][mark] = peak arrays (as load_peaks).
    """
    man = json.load(open(os.path.join(data_dir, manifest)))
    panel = defaultdict(dict)
    for rec in man:
        path = os.path.join(data_dir, rec["file"])
        panel[rec["tissue"]][rec["mark"]] = load_peaks(path, chroms=chroms)
    return dict(panel)


def _peak_overlap_bp(chrom, start, end, peaks_c):
    """(covered_bp, signal-weighted covered_bp) of [start,end) with peaks on one chrom.

    Peaks may overlap; covered_bp is the union of covered positions (no double
    count), while the signal-weighted term sums signal*overlap per peak (used for
    a mean-signal readout, tolerant of small peak overlaps).
    """
    starts, ends, signal = peaks_c["starts"], peaks_c["ends"], peaks_c["signal"]
    lo = int(np.searchsorted(ends, start, side="right"))
    hi = int(np.searchsorted(starts, end, side="left"))
    if hi <= lo:
        return 0, 0.0
    # clip intervals to [start,end) then union via sort-merge for covered_bp
    segs = []
    sig_bp = 0.0
    for i in range(lo, hi):
        a = max(start, int(starts[i]))
        b = min(end, int(ends[i]))
        if b > a:
            segs.append((a, b))
            sig_bp += signal[i] * (b - a)
    if not segs:
        return 0, 0.0
    segs.sort()
    covered = 0
    cs, ce = segs[0]
    for a, b in segs[1:]:
        if a > ce:
            covered += ce - cs
            cs, ce = a, b
        else:
            ce = max(ce, b)
    covered += ce - cs
    return covered, sig_bp


def region_mark_fingerprint(chrom, start, end, tissue_panel, marks=MARKS):
    """Per-mark peak coverage fraction + mean signal for one region in one tissue.

    tissue_panel = panel[tissue] (dict mark -> peak arrays).
    Returns dict with, per mark: coverage (fraction of region under a peak) and
    mean_signal (signal-weighted mean over covered bp); plus a scalar
    active_minus_repressive fingerprint score.
    """
    L = end - start
    frac, mean_sig = {}, {}
    for m in marks:
        peaks = tissue_panel.get(m, {})
        pc = peaks.get(chrom)
        if pc is None or L <= 0:
            frac[m] = 0.0
            mean_sig[m] = 0.0
            continue
        cov, sig_bp = _peak_overlap_bp(chrom, start, end, pc)
        frac[m] = cov / L
        mean_sig[m] = (sig_bp / cov) if cov > 0 else 0.0
    act_vals = [frac[m] for m in marks if m in ACTIVE_MARKS]
    rep_vals = [frac[m] for m in marks if m in REPRESSIVE_MARKS]
    act = float(np.mean(act_vals)) if act_vals else 0.0
    rep = float(np.mean(rep_vals)) if rep_vals else 0.0
    return {"coverage": frac, "mean_signal": mean_sig,
            "active_minus_repressive": act - rep}


def annotate_regions_marks(regions, tissue_panel, marks=MARKS):
    """regions: iterable of (chrom, start, end) -> list of fingerprint dicts."""
    return [region_mark_fingerprint(c, s, e, tissue_panel, marks) for (c, s, e) in regions]


def _mean_coverage(regions, tissue_panel, marks):
    """Mean per-mark peak coverage fraction across a region set."""
    if not regions:
        return {m: np.nan for m in marks}
    acc = {m: 0.0 for m in marks}
    for (c, s, e) in regions:
        fp = region_mark_fingerprint(c, s, e, tissue_panel, marks)
        for m in marks:
            acc[m] += fp["coverage"][m]
    n = len(regions)
    return {m: acc[m] / n for m in marks}


def mark_enrichment_test(query_regions, tissue_panel, null_model, marks=MARKS,
                         n_per=100, n_boot=1000, seed=0):
    """Per-mark enrichment of query peak coverage vs. the GC/length/N-matched null.

    Mirrors engine.enrichment_test but on the continuous mark-coverage fingerprint
    rather than dominant ChromHMM state. For each mark: mean query coverage, mean
    null coverage, log2 fold-enrichment (with a coverage pseudocount), bootstrap
    CI, and a two-sided empirical p from resampling the null at the query size.

    Returns (results, meta) where results is a list of per-mark dicts:
        mark, query_cov, null_cov, log2_fold, ci_low, ci_high, p_emp, n_query, n_null
    """
    rng = np.random.default_rng(seed)
    null_regions, meta = null_model.sample_matched(query_regions, n_per=n_per)

    q_cov = _mean_coverage(query_regions, tissue_panel, marks)
    # per-null-region coverage vectors, for bootstrap
    null_fp = np.array([[region_mark_fingerprint(c, s, e, tissue_panel, marks)["coverage"][m]
                         for m in marks] for (c, s, e) in null_regions], dtype=float)
    n_null = len(null_regions)
    nq = len(query_regions)
    n_cov = {m: (null_fp[:, j].mean() if n_null else np.nan) for j, m in enumerate(marks)}

    eps = 1e-4  # coverage-scale pseudocount
    results = []
    for j, m in enumerate(marks):
        qf, nf = q_cov[m], n_cov[m]
        log2fold = np.log2((qf + eps) / (nf + eps))
        if n_null and nq:
            idx = rng.integers(0, n_null, size=(n_boot, nq))
            bmean = null_fp[idx, j].mean(axis=1)
            blog2 = np.log2((bmean + eps) / (nf + eps))
            ci_low, ci_high = np.percentile(blog2, [2.5, 97.5])
            p_emp = (np.sum(np.abs(blog2) >= abs(log2fold)) + 1) / (n_boot + 1)
        else:
            ci_low = ci_high = p_emp = np.nan
        results.append({"mark": m, "query_cov": qf, "null_cov": nf,
                        "log2_fold": log2fold, "ci_low": ci_low, "ci_high": ci_high,
                        "p_emp": p_emp, "n_query": nq, "n_null": n_null})
    return results, meta
