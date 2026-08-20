"""Unit tests for the histone-mark fingerprint layer (synthetic peaks, no I/O)."""
import numpy as np

from cfdna_chromatin import histone as H


def _peaks(intervals):
    """Build a load_peaks-style dict for one chrom from (start,end,signal) tuples."""
    intervals = sorted(intervals)
    return {"chr1": {
        "starts": np.array([i[0] for i in intervals], dtype=np.int64),
        "ends": np.array([i[1] for i in intervals], dtype=np.int64),
        "signal": np.array([i[2] for i in intervals], dtype=np.float64),
    }}


def test_overlap_fraction_exact():
    # region [100,200); one peak covering [120,170) -> 50/100 = 0.5 coverage
    panel = {"H3K4me3": _peaks([(120, 170, 5.0)])}
    fp = H.region_mark_fingerprint("chr1", 100, 200, panel, marks=["H3K4me3"])
    assert abs(fp["coverage"]["H3K4me3"] - 0.5) < 1e-9
    assert abs(fp["mean_signal"]["H3K4me3"] - 5.0) < 1e-9


def test_overlapping_peaks_union_not_double_counted():
    # two overlapping peaks [100,160) and [140,200) over region [100,200)
    # union covers the whole region -> coverage 1.0, never >1
    panel = {"H3K4me1": _peaks([(100, 160, 1.0), (140, 200, 1.0)])}
    fp = H.region_mark_fingerprint("chr1", 100, 200, panel, marks=["H3K4me1"])
    assert abs(fp["coverage"]["H3K4me1"] - 1.0) < 1e-9


def test_no_peak_zero_coverage():
    panel = {"H3K9me3": _peaks([(1000, 2000, 3.0)])}
    fp = H.region_mark_fingerprint("chr1", 0, 500, panel, marks=["H3K9me3"])
    assert fp["coverage"]["H3K9me3"] == 0.0
    assert fp["mean_signal"]["H3K9me3"] == 0.0


def test_active_minus_repressive_sign():
    # active marks fully cover, repressive absent -> score > 0
    full = _peaks([(0, 100, 1.0)])
    empty = _peaks([(1000, 1001, 1.0)])
    panel = {m: (full if m in H.ACTIVE_MARKS else empty) for m in H.MARKS}
    fp = H.region_mark_fingerprint("chr1", 0, 100, panel)
    assert fp["active_minus_repressive"] > 0.9
