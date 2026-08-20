"""Unit tests for compartment attribution (attribute.py).

These exercise the generic signal / inflammation / background / control engine
against a synthetic in-memory panel, so they do not depend on any particular
application panel being registered in references.PANELS.
"""
import numpy as np

from cfdna_chromatin import attribute as AT
from cfdna_chromatin import references as R


# A neutral synthetic panel wired to the same tissue slots the fixture builds.
# Registered at import time so attribute_region(..., "test_panel") resolves.
R.PANELS.setdefault("test_panel", {
    "signal":       ["keratinocyte"],
    "background":   ["neutrophil", "Bcell"],
    "inflammation": ["monocyte"],
    "control":      ["K562"],
    "note": "Synthetic panel for unit tests only.",
})


def _panel_from_marks(coverage_by_tissue, region=("chr1", 1000, 3000)):
    """Build a minimal hpanel where each tissue's active marks are a single peak
    spanning the whole region at the requested coverage fraction (approx)."""
    c, s, e = region
    L = e - s
    hp = {}
    for tissue, cov in coverage_by_tissue.items():
        marks = {}
        # place a peak covering `cov` fraction of the region for each ACTIVE mark
        peak_len = int(cov * L)
        for m in AT.ACTIVE_MARKS:
            if peak_len > 0:
                starts = np.array([s], dtype=np.int64)
                ends = np.array([s + peak_len], dtype=np.int64)
                sig = np.array([10.0], dtype=float)
            else:
                starts = np.array([], dtype=np.int64)
                ends = np.array([], dtype=np.int64)
                sig = np.array([], dtype=float)
            marks[m] = {c: {"starts": starts, "ends": ends, "signal": sig}}
        # repressive marks empty
        for m in ("H3K27me3", "H3K9me3"):
            marks[m] = {c: {"starts": np.array([], dtype=np.int64),
                            "ends": np.array([], dtype=np.int64),
                            "signal": np.array([], dtype=float)}}
        hp[tissue] = marks
    return hp


def test_signal_specific_region_called_signal():
    # signal tissue open, everything else closed -> signal
    hp = _panel_from_marks({"keratinocyte": 0.8, "monocyte": 0.0,
                            "neutrophil": 0.0, "Bcell": 0.0, "K562": 0.0})
    r = AT.attribute_region("chr1", 1000, 3000, hp, "test_panel")
    assert r["compartment"] == "signal"
    assert r["top"] == "signal"
    assert r["log2_margin"] > 0.5


def test_shared_active_region_is_ambiguous_not_signal():
    # open everywhere (housekeeping-like) -> ambiguous, never signal
    hp = _panel_from_marks({"keratinocyte": 0.7, "monocyte": 0.7,
                            "neutrophil": 0.7, "Bcell": 0.7, "K562": 0.7})
    r = AT.attribute_region("chr1", 1000, 3000, hp, "test_panel")
    assert r["compartment"] == "ambiguous"


def test_inflammation_region_labeled_inflammation():
    hp = _panel_from_marks({"keratinocyte": 0.0, "monocyte": 0.8,
                            "neutrophil": 0.0, "Bcell": 0.0, "K562": 0.0})
    r = AT.attribute_region("chr1", 1000, 3000, hp, "test_panel")
    assert r["compartment"] == "inflammation"


def test_quiescent_region_unexplained():
    hp = _panel_from_marks({"keratinocyte": 0.0, "monocyte": 0.0,
                            "neutrophil": 0.0, "Bcell": 0.0, "K562": 0.0})
    r = AT.attribute_region("chr1", 1000, 3000, hp, "test_panel")
    assert r["compartment"] == "unexplained"
    assert r["called_open"] is False


def test_background_region_labeled_background():
    hp = _panel_from_marks({"keratinocyte": 0.0, "monocyte": 0.0,
                            "neutrophil": 0.8, "Bcell": 0.7, "K562": 0.0})
    r = AT.attribute_region("chr1", 1000, 3000, hp, "test_panel")
    assert r["compartment"] == "background"


def test_summary_fractions_sum_to_one():
    hp = _panel_from_marks({"keratinocyte": 0.8, "monocyte": 0.0,
                            "neutrophil": 0.0, "Bcell": 0.0, "K562": 0.0})
    regions = [("chr1", 1000, 3000)] * 5
    per, summ = AT.attribute_signal(regions, hp, "test_panel")
    assert summ["n"] == 5
    assert abs(sum(summ["fractions"].values()) - 1.0) < 1e-9
    assert summ["signal_fraction"] == 1.0


def test_panels_and_roles_defined():
    # the shipped application panel plus the synthetic test panel
    assert "fetal" in R.PANELS
    assert R.tissue_role("test_panel", "keratinocyte") == "signal"
    assert R.tissue_role("test_panel", "monocyte") == "inflammation"
    assert "monocyte" in R.panel_tissues("test_panel")
    # the shipped fetal panel: placenta is the cell-of-origin signal
    assert R.tissue_role("fetal", "placenta") == "signal"
