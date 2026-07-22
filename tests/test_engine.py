"""Fast unit tests for the chromatin engine (no network, no FASTA required).

Uses a tiny synthetic segmentation + genome so tests run in <1s and don't depend
on the large sequence files. A separate integration test (test_integration.py)
exercises the real references when they are present.
"""
import numpy as np
import pytest

from cfdna_chromatin import references as R, engine as E


def _toy_seg():
    # one chrom, two segments: [0,1000)=TssA, [1000,2000)=Quies
    return {"chrT": {
        "starts": np.array([0, 1000], dtype=np.int64),
        "ends": np.array([1000, 2000], dtype=np.int64),
        "code": np.array([R.S2I["TssA"], R.S2I["Quies"]], dtype=np.int8),
    }}


def _toy_genome():
    n = 2000
    # first 1000 bp are GC-rich, last 1000 AT-rich; no Ns
    cumGC = np.concatenate([[0], np.cumsum(np.r_[np.ones(1000), np.zeros(1000)].astype(np.int64))])
    cumN = np.concatenate([[0], np.zeros(n, dtype=np.int64)])
    return {"chrT": {"n": n, "cumGC": cumGC, "cumN": cumN}}


def test_state_bp_split():
    seg = _toy_seg()["chrT"]
    bp = E.region_state_bp("chrT", 500, 1500, seg)
    assert bp[R.S2I["TssA"]] == 500
    assert bp[R.S2I["Quies"]] == 500


def test_annotate_dominant():
    seg = _toy_seg()
    a = E.annotate_region("chrT", 0, 800, seg)
    assert a["dominant_state"] == "TssA"
    assert a["active_fraction"] == pytest.approx(1.0)


def test_annotate_empty_region():
    seg = _toy_seg()
    a = E.annotate_region("chrT", 5000, 6000, seg)  # beyond segments
    assert a["dominant_state"] is None


def test_null_model_matches_length_and_gc():
    genome = _toy_genome()
    nm = E.NullModel(genome=genome, chroms=["chrT"], seed=0, gc_tol=0.05)
    # a GC-rich query in [0,200): matched null must be length 200 and GC-rich
    null, meta = nm.sample_matched([("chrT", 0, 200)], n_per=10)
    assert all((e - s) == 200 for _, s, e in null)


def test_enrichment_promoter_positive():
    seg = _toy_seg()
    genome = _toy_genome()
    nm = E.NullModel(genome=genome, chroms=["chrT"], seed=0, gc_tol=1.0)
    # query = promoter region; expect Promoter group enriched vs null
    q = [("chrT", 0, 300)] * 10
    res, _ = E.enrichment_test(q, seg, nm, by="group", n_per=20, n_boot=200, seed=0)
    prom = [r for r in res if r["label"] == "Promoter"][0]
    assert prom["log2_fold"] > 0
