"""
accessibility.py -- open-chromatin (DNase-seq) layer.

A third, state-independent readout on top of ChromHMM and the histone-mark
fingerprint. For each region we compute the fraction covered by ENCODE
open-chromatin peaks (DNase-seq hypersensitive sites) and a peak-weighted mean
signal, then test enrichment of query coverage against the same
GC/length/mappability-matched null the other layers use.

Open chromatin is the most direct footprint of regulatory activity: active
promoters and enhancers are nucleosome-depleted and hypersensitive to DNase,
while heterochromatin and quiescent regions are not. In the cfDNA setting,
tissue-specific accessibility is what makes cell-of-origin inference possible
(nucleosome spacing in plasma tracks the accessibility landscape of the cells
that shed the DNA), so an accessibility axis is the natural complement to the
mark- and state-based layers.

One track per tissue (unlike the 6-mark histone panel), so the API is a thin
specialisation: load_track / load_access_panel and access_enrichment_test.

Peak files: ENCODE DNase-seq narrowPeak, GRCh38, one experiment per tissue,
subset to the working chromosomes. Columns: chrom, start, end, signalValue
(0-based half-open). See data/access/manifest.json for accessions.

NOTE: neutrophils have highly condensed chromatin and are essentially absent
from ENCODE accessibility assays; the "neutrophil" track is a CD14+ monocyte
DNase-seq experiment used as a labeled myeloid proxy (manifest: proxy_for).
"""
from __future__ import annotations
import json
import os

import numpy as np

# reuse the peak I/O and union-overlap accounting from the histone layer so the
# two open-region layers share one implementation (and one set of tests).
from .histone import load_peaks, _peak_overlap_bp

ASSAY = "DNase-seq"


def load_track(bed_gz_path, chroms=None):
    """Load one DNase-seq peak .bed.gz into per-chromosome sorted arrays.

    Thin alias of histone.load_peaks (same 4-column schema).
    """
    return load_peaks(bed_gz_path, chroms=chroms)


def load_access_panel(data_dir, chroms=None, manifest="manifest.json"):
    """Load every accessibility track listed in data/access/manifest.json.

    Returns panel[tissue] = peak arrays (as load_track). One track per tissue.
    """
    man = json.load(open(os.path.join(data_dir, manifest)))
    panel = {}
    for rec in man:
        path = os.path.join(data_dir, rec["file"])
        panel[rec["tissue"]] = load_track(path, chroms=chroms)
    return panel


def region_accessibility(chrom, start, end, track):
    """Open-chromatin coverage fraction + mean signal for one region.

    track = panel[tissue] (peak arrays). Returns {coverage, mean_signal}.
    """
    L = end - start
    pc = track.get(chrom)
    if pc is None or L <= 0:
        return {"coverage": 0.0, "mean_signal": 0.0}
    cov, sig_bp = _peak_overlap_bp(chrom, start, end, pc)
    return {"coverage": cov / L,
            "mean_signal": (sig_bp / cov) if cov > 0 else 0.0}


def annotate_regions_access(regions, track):
    """regions: iterable of (chrom, start, end) -> list of accessibility dicts."""
    return [region_accessibility(c, s, e, track) for (c, s, e) in regions]


def _mean_coverage(regions, track):
    if not regions:
        return np.nan
    return float(np.mean([region_accessibility(c, s, e, track)["coverage"]
                          for (c, s, e) in regions]))


def access_enrichment_test(query_regions, track, null_model,
                           n_per=100, n_boot=1000, seed=0):
    """Enrichment of query open-chromatin coverage vs. the matched null.

    Mirrors histone.mark_enrichment_test for the single accessibility track:
    mean query coverage, mean null coverage, log2 fold-enrichment (coverage
    pseudocount), bootstrap CI, two-sided empirical p.

    Returns (result, meta) where result is a dict:
        query_cov, null_cov, log2_fold, ci_low, ci_high, p_emp, n_query, n_null
    """
    rng = np.random.default_rng(seed)
    null_regions, meta = null_model.sample_matched(query_regions, n_per=n_per)

    q_cov = _mean_coverage(query_regions, track)
    null_cov_vec = np.array([region_accessibility(c, s, e, track)["coverage"]
                             for (c, s, e) in null_regions], dtype=float)
    n_null = len(null_regions)
    nq = len(query_regions)
    n_cov = float(null_cov_vec.mean()) if n_null else np.nan

    eps = 1e-4  # coverage-scale pseudocount
    log2fold = np.log2((q_cov + eps) / (n_cov + eps))
    if n_null and nq:
        idx = rng.integers(0, n_null, size=(n_boot, nq))
        bmean = null_cov_vec[idx].mean(axis=1)
        blog2 = np.log2((bmean + eps) / (n_cov + eps))
        ci_low, ci_high = np.percentile(blog2, [2.5, 97.5])
        p_emp = (np.sum(np.abs(blog2) >= abs(log2fold)) + 1) / (n_boot + 1)
    else:
        ci_low = ci_high = p_emp = np.nan

    result = {"query_cov": q_cov, "null_cov": n_cov, "log2_fold": log2fold,
              "ci_low": ci_low, "ci_high": ci_high, "p_emp": p_emp,
              "n_query": nq, "n_null": n_null}
    return result, meta
