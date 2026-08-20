#!/usr/bin/env python
"""plot_pc_tissue.py -- visualize the PC<->tissue map and per-PC SHAP contributions.

Reads the two aggregate tables written by pc_tissue_map.R and renders a
two-panel figure:

  A) PC x tissue alignment heatmap. Each cell is the beta-oriented Pearson
     correlation between a PC's genome-wide loading pattern and a tissue's open
     chromatin. Rows are the top-K PCs by total signed-SHAP contribution; the
     row label carries the PC's single best-aligned tissue. Diverging map
     centered at 0: red = this PC's FF-raising direction aligns with the
     tissue's openness, blue = anti-aligns.

  B) Per-PC contribution to the two SHAP directions. For the same K PCs, the
     mean exact contribution L[bin,k]*beta_k over the top-positive (FF-up) bins
     and over the top-negative (FF-down) bins. Shows which latent components
     assemble the fetal/placental-leaning vs maternal-leaning weight.

Usage:
  python plot_pc_tissue.py --indir pc_tissue_out --topk 15 --out pc_tissue_figure.png
"""
import argparse
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    from cfdna_chromatin.figstyle import apply_figure_style  # if packaged later
except Exception:
    def apply_figure_style(*a, **k):
        mpl.rcParams.update({"font.size": 8, "axes.spines.top": False,
                             "axes.spines.right": False})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", default="pc_tissue_out")
    ap.add_argument("--topk", type=int, default=15,
                    help="number of PCs to show, ranked by |pos|+|neg| SHAP contribution")
    ap.add_argument("--out", default="pc_tissue_figure.png")
    args = ap.parse_args()

    corr = pd.read_csv(os.path.join(args.indir, "pc_tissue_corr.csv"))
    contrib = pd.read_csv(os.path.join(args.indir, "pc_shap_contrib.csv"))
    tissue_cols = [c for c in corr.columns if c.startswith("corr_")]
    tissues = [c[len("corr_"):] for c in tissue_cols]

    # rank PCs by total signed-SHAP contribution magnitude
    m = contrib.copy()
    m["rank_score"] = m["pos_contrib"].abs() + m["neg_contrib"].abs()
    order = m.sort_values("rank_score", ascending=False)["pc"].head(args.topk).tolist()
    corr_i = corr.set_index("pc").loc[order]
    contrib_i = contrib.set_index("pc").loc[order]

    H = corr_i[tissue_cols].to_numpy()               # K x tissues
    row_lab = [f"{pc} · {corr_i.loc[pc, 'best_tissue']}" for pc in order]

    apply_figure_style(sizes=(8, 7, 6))
    # Panel A width scales with the tissue count so many-tissue atlases (e.g. the
    # 73-tissue panel) don't crush the x labels. ~0.22 in/tissue for panel A,
    # plus a fixed ~4 in for panel B. Cap the width so it stays slide-friendly.
    n_tissue = len(tissues)
    panelA_w = min(max(4.5, 0.22 * n_tissue), 16.0)
    fig_w = panelA_w + 4.5
    # x-label font shrinks a touch as the panel fills up
    xlab_fs = 7 if n_tissue <= 30 else (6 if n_tissue <= 55 else 5)
    fig, (axA, axB) = plt.subplots(
        1, 2, figsize=(fig_w, max(4.0, 0.32 * len(order) + 1.6)),
        gridspec_kw={"width_ratios": [panelA_w, 4.0], "wspace": 0.45})

    # ---- Panel A: PC x tissue heatmap -------------------------------------
    vmax = np.nanmax(np.abs(H)) or 1.0
    im = axA.imshow(H, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    axA.set_xticks(range(len(tissues)))
    # rotate fully vertical + anchor at the tick centre so long tissue names
    # (e.g. "esophagus squamous epithelium") stack without overlapping neighbours
    axA.set_xticklabels(tissues, rotation=90, ha="center", va="top",
                        fontsize=xlab_fs)
    axA.tick_params(axis="x", pad=2)
    axA.set_yticks(range(len(order)))
    axA.set_yticklabels(row_lab, fontsize=6)
    axA.set_title("PC \u2194 tissue openness alignment\n(\u03b2-oriented correlation)",
                  fontsize=8, loc="left")
    if H.shape[0] * H.shape[1] <= 200:
        for i in range(H.shape[0]):
            for j in range(H.shape[1]):
                v = H[i, j]
                axA.text(j, i, f"{v:.2f}".lstrip("0").replace("-0.", "-."),
                         ha="center", va="center", fontsize=5,
                         color="white" if abs(v) > 0.6 * vmax else "black")
    cb = fig.colorbar(im, ax=axA, fraction=0.046, pad=0.02)
    cb.set_label("aligns with FF-up \u2192", fontsize=6)
    cb.ax.tick_params(labelsize=6)

    # ---- Panel B: per-PC pos/neg SHAP contribution ------------------------
    y = np.arange(len(order))[::-1]                  # top PC at top
    pos = contrib_i["pos_contrib"].to_numpy()
    neg = contrib_i["neg_contrib"].to_numpy()
    axB.barh(y + 0.18, pos, height=0.36, color="#c0392b",
             label="top FF-up bins")
    axB.barh(y - 0.18, neg, height=0.36, color="#2471a3",
             label="top FF-down bins")
    axB.axvline(0, color=".3", lw=0.7)
    axB.set_yticks(y)
    axB.set_yticklabels(order, fontsize=6)
    axB.set_xlabel("mean per-PC contribution  L[bin,k]\u00b7\u03b2\u2096", fontsize=7)
    axB.set_title("Which PCs build each SHAP direction", fontsize=8, loc="left")
    axB.legend(fontsize=6, frameon=False, loc="lower right")

    fig.savefig(args.out, dpi=200, bbox_inches="tight")
    print(f"[plot] {args.out}  ({len(order)} PCs x {len(tissues)} tissues)")


if __name__ == "__main__":
    main()
