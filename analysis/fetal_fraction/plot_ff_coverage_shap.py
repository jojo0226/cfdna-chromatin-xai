#!/usr/bin/env python
"""plot_ff_coverage_shap.py -- slide-ready 2-panel coverage + SHAP figure for FF.

A presentation figure that tells the seqFF++ story in one picture:

  TOP panel   normalized sequencing coverage across a set of genomic regions,
              split into "High FF" (>30%) vs "Low FF" (<5%) samples. One zone of
              regions is DIFFERENTIAL (coverage separates the two FF groups); a
              second zone is a BASELINE CONTROL (the two groups overlap).

  BOTTOM panel the model's SHAP value for the SAME regions, on the SAME x-axis,
              vertically aligned bin-for-bin. SHAP is large/non-zero over the
              differential zone and ~0 over the baseline zone -- i.e. the model
              keys on exactly the regions where coverage carries FF information.

The two panels share the x-axis (2 rows x 1 col, sharex) so region i sits at the
same horizontal position in both, making the coverage->SHAP correspondence read
directly off the slide.

--------------------------------------------------------------------------------
HOW TO USE WITH YOUR OWN DATA
--------------------------------------------------------------------------------
Replace the MOCK DATA block below with your two DataFrames (same columns/dtypes):

  coverage_df : one row per (sample x region) -- LONG/tidy format
     sample_id : str        e.g. "S001"                unique per sample
     ff_class  : str        "High FF" or "Low FF"      exactly these two labels
     region    : str        e.g. "chr7:5,300,000"      region/bin label (x tick)
     coverage  : float      normalized coverage        e.g. GC/NB-corrected residual,
                                                        or depth / mean-depth (~1.0)

  shap_df : one row per region -- the model-level signed SHAP for that region
     region     : str       must match coverage_df.region labels
     shap_value : float     signed SHAP (w_bin, or w_bin * mean|x-xbar|)

  region_order  : list[str] left-to-right order of regions on the x-axis
  region_zones  : dict      {region_label: "differential" | "baseline"} -- used
                            only to shade/annotate the two zones; optional.

Then just call:  make_figure(coverage_df, shap_df, region_order, region_zones)

Everything below the MOCK DATA block is generic and needs no editing.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# ============================================================================
# PRESENTATION STYLE -- large fonts, high-contrast palette (edit freely)
# ============================================================================
HIGH_FF_COLOR = "#D62728"   # crimson  -- High FF (>30%)
LOW_FF_COLOR  = "#1F4E79"   # deep navy -- Low FF (<5%)
SHAP_POS_COLOR = "#C0392B"  # positive/large SHAP
SHAP_ZERO_COLOR = "#95A5A6" # near-zero SHAP (baseline)
DIFF_ZONE_SHADE = "#FDECEA" # differential zone background
BASE_ZONE_SHADE = "#EEF1F4" # baseline zone background

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 15,          # base -- readable on a projected slide
    "axes.titlesize": 19,
    "axes.labelsize": 16,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "legend.fontsize": 14,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 120,
})

# ============================================================================
# ----------------------------  MOCK DATA  -----------------------------------
# Delete/replace this whole block with your real coverage_df and shap_df.
# It exists so you can (a) see the exact expected structure and (b) run the
# script standalone to preview the figure before plugging in your data.
# ============================================================================
def _make_mock_data(seed: int = 0):
    rng = np.random.default_rng(seed)

    # Ordered regions: 4 DIFFERENTIAL regions then 4 BASELINE-control regions.
    diff_regions = ["chr7:5.30M", "chr7:5.35M", "chr6:32.1M", "chr6:32.2M"]
    base_regions = ["chr2:88.0M", "chr2:88.1M", "chr11:65.4M", "chr11:65.5M"]
    region_order = diff_regions + base_regions
    region_zones = {r: "differential" for r in diff_regions}
    region_zones.update({r: "baseline" for r in base_regions})

    # Per-region "true" mean coverage for each FF group.
    #   differential zone: High FF sits ABOVE Low FF (placental-open -> the
    #     FF-informative separation the model uses)
    #   baseline zone: the two groups coincide (no FF information)
    means = {
        #                       High FF   Low FF
        "chr7:5.30M":          (1.45,     1.00),
        "chr7:5.35M":          (1.55,     1.02),
        "chr6:32.1M":          (1.38,     1.01),
        "chr6:32.2M":          (1.30,     0.99),
        "chr2:88.0M":          (1.00,     1.00),
        "chr2:88.1M":          (1.02,     1.01),
        "chr11:65.4M":         (0.98,     0.99),
        "chr11:65.5M":         (1.01,     1.00),
    }

    n_high, n_low = 40, 40   # samples per FF group
    rows = []
    for r in region_order:
        mu_hi, mu_lo = means[r]
        for i in range(n_high):
            rows.append(("H%03d" % i, "High FF", r,
                         rng.normal(mu_hi, 0.10)))
        for i in range(n_low):
            rows.append(("L%03d" % i, "Low FF", r,
                         rng.normal(mu_lo, 0.10)))
    coverage_df = pd.DataFrame(rows,
                               columns=["sample_id", "ff_class", "region", "coverage"])

    # Model-level signed SHAP per region: large over the differential zone,
    # ~0 over the baseline zone. (Sign here is positive for illustration; in the
    # real tool the sign is w_bin and may be negative for placenta-open bins.)
    shap_vals = {
        "chr7:5.30M": 0.041, "chr7:5.35M": 0.052, "chr6:32.1M": 0.036,
        "chr6:32.2M": 0.030, "chr2:88.0M": 0.002, "chr2:88.1M": -0.001,
        "chr11:65.4M": 0.001, "chr11:65.5M": 0.000,
    }
    shap_df = pd.DataFrame({"region": region_order,
                            "shap_value": [shap_vals[r] for r in region_order]})
    return coverage_df, shap_df, region_order, region_zones


# ============================================================================
# ---------------------------  PLOTTING  -------------------------------------
# Generic below this line -- works on any coverage_df / shap_df in the format
# documented in the header.
# ============================================================================
def _summarize_coverage(coverage_df: pd.DataFrame, region_order: list[str]):
    """mean +/- 95% CI of coverage per (ff_class, region), on the given order."""
    g = (coverage_df.groupby(["ff_class", "region"])["coverage"]
         .agg(["mean", "std", "count"]).reset_index())
    g["sem"] = g["std"] / np.sqrt(g["count"].clip(lower=1))
    g["ci95"] = 1.96 * g["sem"]
    g["x"] = g["region"].map({r: i for i, r in enumerate(region_order)})
    return g.sort_values("x")


def make_figure(coverage_df: pd.DataFrame,
                shap_df: pd.DataFrame,
                region_order: list[str],
                region_zones: dict | None = None,
                out_path: str = "ff_coverage_shap.png",
                title: str = "Coverage and SHAP at FF-informative chromatin regions"):
    """Build the stacked coverage (top) + SHAP (bottom) slide figure."""
    x = np.arange(len(region_order))
    cov = _summarize_coverage(coverage_df, region_order)

    # align SHAP to the region order
    shap_map = dict(zip(shap_df["region"], shap_df["shap_value"]))
    shap_y = np.array([shap_map.get(r, np.nan) for r in region_order])

    fig, (ax_cov, ax_shap) = plt.subplots(
        2, 1, figsize=(12, 8.5), sharex=True,
        gridspec_kw={"height_ratios": [1.6, 1.0], "hspace": 0.12})

    # ---- zone background shading on BOTH panels (aligned regions) ----------
    if region_zones:
        for i, r in enumerate(region_order):
            z = region_zones.get(r)
            shade = DIFF_ZONE_SHADE if z == "differential" else (
                BASE_ZONE_SHADE if z == "baseline" else None)
            if shade:
                for ax in (ax_cov, ax_shap):
                    ax.axvspan(i - 0.5, i + 0.5, color=shade, zorder=0)

    # ---- TOP: coverage mean +/- CI per FF group ----------------------------
    for cls, color in [("High FF", HIGH_FF_COLOR), ("Low FF", LOW_FF_COLOR)]:
        sub = cov[cov["ff_class"] == cls]
        ax_cov.plot(sub["x"], sub["mean"], "-o", color=color, lw=2.6,
                    ms=9, label=cls, zorder=3)
        ax_cov.fill_between(sub["x"], sub["mean"] - sub["ci95"],
                            sub["mean"] + sub["ci95"], color=color,
                            alpha=0.20, zorder=2)
    ax_cov.set_ylabel("Normalized coverage")
    ax_cov.set_title(title, pad=12, fontweight="bold")
    ax_cov.legend(loc="center right", frameon=True, ncol=1)
    ax_cov.grid(axis="y", ls=":", alpha=0.5)

    # ---- BOTTOM: signed SHAP per region (same x) ---------------------------
    bar_colors = [SHAP_POS_COLOR if abs(v) >= 0.01 else SHAP_ZERO_COLOR
                  for v in np.nan_to_num(shap_y)]
    ax_shap.bar(x, shap_y, width=0.7, color=bar_colors, zorder=3,
                edgecolor="white", linewidth=0.6)
    ax_shap.axhline(0, color="0.3", lw=1.0)
    ax_shap.set_ylabel("SHAP value")
    ax_shap.set_xlabel("Genomic region")
    ax_shap.grid(axis="y", ls=":", alpha=0.5)

    # ---- shared x tick labels (bottom only) --------------------------------
    ax_shap.set_xticks(x)
    ax_shap.set_xticklabels(region_order, rotation=35, ha="right")
    ax_shap.set_xlim(-0.6, len(region_order) - 0.4)

    # ---- legend for zone shading ------------------------------------------
    if region_zones:
        handles = [Patch(color=DIFF_ZONE_SHADE, label="Differential zone"),
                   Patch(color=BASE_ZONE_SHADE, label="Baseline control")]
        ax_shap.legend(handles=handles, loc="upper right", frameon=True, ncol=2)

    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    print(f"[plot] wrote {out_path}")
    return fig


def _annotate_zones(ax, region_order, region_zones):
    """Draw a labeled bracket spanning each contiguous zone at the top."""
    labels = {"differential": "Differential (FF-informative)",
              "baseline": "Baseline control"}
    ymax = ax.get_ylim()[1]
    i = 0
    n = len(region_order)
    while i < n:
        z = region_zones.get(region_order[i])
        j = i
        while j + 1 < n and region_zones.get(region_order[j + 1]) == z:
            j += 1
        if z in labels:
            xmid = (i + j) / 2.0
            ax.text(xmid, ymax * 0.99, labels[z], ha="center", va="top",
                    fontsize=13, fontweight="bold", color="0.25")
        i = j + 1


if __name__ == "__main__":
    coverage_df, shap_df, region_order, region_zones = _make_mock_data()
    make_figure(coverage_df, shap_df, region_order, region_zones,
                out_path="ff_coverage_shap.png")
