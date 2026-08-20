"""Unit tests for the importance/selection front-end (selection.py)."""
import numpy as np
import pandas as pd

from cfdna_chromatin import selection as SEL


def test_parse_bin():
    assert SEL.parse_bin("chr19:50000000-50500000") == ("chr19", 50000000, 50500000)
    assert SEL.parse_bin("chr7_157000000_157500000")[0] == "chr7"


def test_rank_by_shap_orders_by_mean_abs():
    bins = ["chr1:1000000-1500000", "chr1:2000000-2500000", "chr1:3000000-3500000"]
    # bin 2 has largest magnitude, bin 0 smallest
    mat = pd.DataFrame(
        {bins[0]: [0.01, -0.01], bins[1]: [0.5, -0.5], bins[2]: [0.1, -0.1]},
        index=["s1", "s2"],
    )
    r = SEL.rank_by_shap(mat)
    assert list(r.index) == [bins[1], bins[2], bins[0]]
    assert r.iloc[0]["rank"] == 1
    # direction ~ 0 for symmetric shap
    assert abs(r.loc[bins[1], "direction"]) < 1e-9


def test_rank_by_differential_direction():
    bins = ["chr1:1000000-1500000", "chr1:2000000-2500000"]
    # bin 0 clearly higher in group B, bin 1 flat
    rng = np.random.default_rng(0)
    A = pd.DataFrame({bins[0]: rng.normal(0, 0.01, 20), bins[1]: rng.normal(0, 0.01, 20)})
    B = pd.DataFrame({bins[0]: rng.normal(1, 0.01, 20), bins[1]: rng.normal(0, 0.01, 20)})
    mat = pd.concat([A, B], ignore_index=True)
    y = ["A"] * 20 + ["B"] * 20
    r = SEL.rank_by_differential(mat, y, positive="A")
    assert r.index[0] == bins[0]           # separable bin ranks first
    assert r.loc[bins[0], "direction"] < 0  # A < B
    assert r.loc[bins[0], "p"] < 1e-3


def test_shap_from_model_linear_region_anchored():
    # a linear FF-style regressor: only bin 1 drives the target
    bins = ["chr19:50000000-50500000", "chr19:51000000-51500000",
            "chr19:52000000-52500000"]
    rng = np.random.default_rng(3)
    X = pd.DataFrame(rng.normal(size=(40, 3)), columns=bins,
                     index=[f"s{i}" for i in range(40)])
    y = 5.0 * X[bins[1]] + 0.01 * rng.normal(size=40)   # bin 1 is the signal
    from sklearn.linear_model import LinearRegression
    model = LinearRegression().fit(X, y)

    sv = SEL.shap_from_model(model, X, explainer="auto")
    assert sv.shape == X.shape
    assert list(sv.columns) == bins                      # region-anchored, aligned
    # ranked importance must place the true signal bin first
    ranked = SEL.rank_by_shap(sv)
    assert ranked.index[0] == bins[1]
    # and it must be coordinate-parseable downstream
    reg = SEL.to_regions(ranked, top_n=1)[0]
    assert reg[0] == "chr19" and reg[2] > reg[1]


def test_to_regions_top_n_and_frac():
    bins = [f"chr1:{i}000000-{i}500000" for i in range(1, 11)]
    mat = pd.DataFrame(np.random.default_rng(1).normal(size=(4, 10)), columns=bins)
    r = SEL.rank_by_shap(mat)
    assert len(SEL.to_regions(r, top_n=3)) == 3
    assert len(SEL.to_regions(r, top_frac=0.5)) == 5
    reg = SEL.to_regions(r, top_n=1)[0]
    assert reg[0].startswith("chr") and reg[2] > reg[1]
