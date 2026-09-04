"""Tests for the chromhmm façade: composition, annotation, and hg19 build crossing."""
import os

import pandas as pd
import pytest

from cfdna_chromatin import chromhmm as CH
from cfdna_chromatin import references as R

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, "data")
CHROMS = ["chr19", "chr20", "chr21", "chr22"]


def _has_ref(t):
    try:
        CH.reference_bed_path(DATA, t)
        return True
    except FileNotFoundError:
        return False


pytestmark = pytest.mark.skipif(
    not _has_ref("placenta"),
    reason="ChromHMM reference BEDs not present in data/",
)


def _load_placenta():
    return R.load_segmentation(CH.reference_bed_path(DATA, "placenta"), chroms=CHROMS)


def test_reference_bed_path_accession_and_bare():
    p = CH.reference_bed_path(DATA, "placenta")
    assert p.endswith(".bed.gz") and os.path.exists(p)


def test_regions_frame_from_names():
    df = CH._regions_frame(["chr19:1000000-1050000", "chr20:2000000-2050000"])
    assert list(df.columns[:3]) == ["chrom", "start", "end"]
    assert df.iloc[0]["chrom"] == "chr19" and int(df.iloc[0]["start"]) == 1000000


def test_regions_frame_preserves_importance_columns():
    src = pd.DataFrame({"key": ["chr19:1000000-1050000"],
                        "mean_abs_shap": [3.2], "signed_shap": [-3.2]})
    df = CH._regions_frame(src)
    assert "mean_abs_shap" in df.columns and "signed_shap" in df.columns
    assert {"chrom", "start", "end"}.issubset(df.columns)


def test_state_composition_sums_to_one():
    seg = _load_placenta()
    bins = [f"chr19:{s}-{s+200000}" for s in range(1_000_000, 5_000_000, 200_000)]
    df = CH._regions_frame(bins)
    comp = CH.state_composition(df, seg, by="state")
    # fractions over bins with a defined dominant state should sum ~1
    assert abs(comp["query_frac"].sum() - 1.0) < 1e-6
    assert set(comp["label"]).issubset(set(R.STATES))


def test_annotate_adds_dominant_state():
    seg = _load_placenta()
    df = CH._regions_frame([f"chr19:{s}-{s+200000}" for s in range(1_000_000, 2_000_000, 200_000)])
    ann = CH.annotate(df, seg)
    assert "dominant_state" in ann.columns
    assert len(ann) == len(df)


def test_prepare_query_no_lift_when_same_build():
    df, qc = CH.prepare_query(["chr19:1000000-1050000"], query_build="hg38", ref_build="hg38")
    assert qc is None and len(df) == 1


@pytest.mark.skipif(
    __import__("importlib").util.find_spec("pyliftover") is None,
    reason="pyliftover not installed",
)
def test_prepare_query_hg19_lift_smoke():
    # a well-behaved region; just assert the lift path runs and QC is reported.
    # Prefer a local chain to avoid a network download; skip cleanly if the
    # chain is unavailable and pyliftover's auto-download can't reach UCSC.
    local_chain = "/tmp/hg19ToHg38.over.chain.gz"
    chain = local_chain if os.path.exists(local_chain) else None
    try:
        df, qc = CH.prepare_query(
            pd.DataFrame({"chrom": ["chr19"], "start": [1_000_000], "end": [1_050_000],
                          "mean_abs_shap": [1.0]}),
            query_build="hg19", ref_build="hg38",
            chain_path=chain, round_trip=False,
        )
    except (AttributeError, OSError, ConnectionError) as e:
        pytest.skip(f"hg19->hg38 chain unavailable offline: {e}")
    assert qc is not None and "mapping_rate" in qc
    # importance column preserved through the lift (if it mapped)
    if len(df):
        assert "mean_abs_shap" in df.columns


def _toy_comp():
    # two tissues, hand-made state_composition frames
    a = pd.DataFrame({"label": ["Quies", "TxWk"], "query_frac": [0.8, 0.2],
                      "group": ["Quiescent", "Transcription"], "n_bins": [10, 10]})
    b = pd.DataFrame({"label": ["Quies", "Tx"], "query_frac": [0.6, 0.4],
                      "group": ["Quiescent", "Transcription"], "n_bins": [10, 10]})
    return {"placenta": a, "Bcell": b}


def test_composition_matrix_shape_and_order():
    mat = CH.composition_matrix(_toy_comp(), by="state")
    # columns are tissues, rows are canonical STATES subset (nonzero), values fractions
    assert list(mat.columns) == ["placenta", "Bcell"]
    assert mat.at["Quies", "placenta"] == 0.8
    assert mat.at["Tx", "Bcell"] == 0.4
    assert mat.at["Tx", "placenta"] == 0.0          # zero-filled
    # canonical order preserved (Tx before TxWk before Quies)
    assert list(mat.index) == [s for s in R.STATES if s in set(mat.index)]


def test_multi_cutoff_composition_clamps_and_keys():
    import numpy as np
    rng = np.random.default_rng(0)
    n = 40
    df = pd.DataFrame({
        "chrom": ["chr19"] * n,
        "start": np.arange(n) * 100_000 + 1_000_000,
        "end": np.arange(n) * 100_000 + 1_050_000,
    })
    # segmentation is {chrom: {starts, ends, code}} with integer state codes;
    # one Quies segment spanning all query bins.
    seg = {"chr19": {"starts": np.array([0], dtype=np.int64),
                     "ends": np.array([60_000_000], dtype=np.int64),
                     "code": np.array([R.S2I["Quies"]], dtype=np.int8)}}
    out = CH.multi_cutoff_composition(df, {"t": seg}, cutoffs=[10, 25, 999], by="state")
    assert set(out.keys()) == {10, 25, 999}
    assert out[10]["t"]["n_bins"].iloc[0] == 10       # subset honored
    assert out[999]["t"]["n_bins"].iloc[0] == n       # clamped to len(regions)


def test_plot_composition_heatmaps_writes_png(tmp_path):
    import matplotlib
    matplotlib.use("Agg")
    mats = {1000: CH.composition_matrix(_toy_comp()),
            2000: CH.composition_matrix(_toy_comp())}
    out = tmp_path / "hm.png"
    fig, axes = CH.plot_composition_heatmaps(mats, outpath=str(out))
    assert out.exists() and out.stat().st_size > 0
    assert len(axes) == 2
