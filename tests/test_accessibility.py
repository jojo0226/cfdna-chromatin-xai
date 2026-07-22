"""Unit tests for the open-chromatin (DNase-seq) layer (synthetic peaks, no I/O)."""
import numpy as np

from cfdna_chromatin import accessibility as A


def _track(intervals):
    """Build a load_track-style dict for one chrom from (start,end,signal) tuples."""
    intervals = sorted(intervals)
    return {"chr1": {
        "starts": np.array([i[0] for i in intervals], dtype=np.int64),
        "ends": np.array([i[1] for i in intervals], dtype=np.int64),
        "signal": np.array([i[2] for i in intervals], dtype=np.float64),
    }}


def test_coverage_fraction_exact():
    # region [100,200); one DHS covering [120,170) -> 0.5 coverage
    track = _track([(120, 170, 8.0)])
    r = A.region_accessibility("chr1", 100, 200, track)
    assert abs(r["coverage"] - 0.5) < 1e-9
    assert abs(r["mean_signal"] - 8.0) < 1e-9


def test_overlapping_peaks_union_not_double_counted():
    # two overlapping DHS over region [100,200) -> union covers all -> 1.0, never >1
    track = _track([(100, 160, 1.0), (140, 200, 1.0)])
    r = A.region_accessibility("chr1", 100, 200, track)
    assert abs(r["coverage"] - 1.0) < 1e-9


def test_no_peak_zero_coverage():
    track = _track([(1000, 2000, 3.0)])
    r = A.region_accessibility("chr1", 0, 500, track)
    assert r["coverage"] == 0.0
    assert r["mean_signal"] == 0.0


def test_enrichment_open_region_positive():
    # A query region sitting fully inside a broad DHS should enrich vs a null
    # drawn from mostly closed chromatin. Use a tiny hand-built null model.
    track = _track([(0, 10_000, 5.0)])

    class _NM:
        def sample_matched(self, query_regions, n_per=100):
            # null: regions far from any peak (closed) -> ~0 coverage
            null = [("chr1", 2_000_000 + i * 1000, 2_000_500 + i * 1000)
                    for i in range(50)]
            return null, {"n_null": len(null)}

    query = [("chr1", 1000, 1500), ("chr1", 2000, 2500)]  # inside the DHS
    res, _ = A.access_enrichment_test(query, track, _NM(), n_boot=200, seed=0)
    assert res["query_cov"] > 0.9      # query fully open
    assert res["null_cov"] < 0.01      # null closed
    assert res["log2_fold"] > 2.0      # strong positive enrichment
