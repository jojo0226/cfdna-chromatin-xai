#!/usr/bin/env python3
"""POS-vs-NEG tissue-enrichment dumbbell plot.

The FF SHAP attribution splits into two signed sets: bins whose coverage
*raises* the fetal-fraction estimate (POS set) and bins that *lower* it (NEG
set). Run `run_ff_tissue_track.py` once on each set (writing to a pos/ and a
neg/ outdir); this script overlays the two resulting cell-of-origin profiles as
a dumbbell so the fetal <-> maternal contrast reads directly:

    one row per tissue, two dots -- POS-set spec_z and NEG-set spec_z --
    joined by a line whose length IS the contrast. Tissues are sorted by
    (POS - NEG), so genuine fetal tissues (placenta/trophoblast, high in POS,
    low in NEG) collect at the top and maternal hematopoietic tissues
    (monocyte/neutrophil/CMP, high in NEG) at the bottom.

Reads only the `ff_tissue_results.csv` that the track already writes -- no
recomputation, no atlas needed here.

Usage:
    python plot_ff_tissue_dumbbell.py --pos _fig_demo/pos --neg _fig_demo/neg \
        --outdir _fig_demo --topn 500 --top-k 25
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def _load_specz(results_dir, topn, metric):
    path = os.path.join(results_dir, "ff_tissue_results.csv")
    df = pd.read_csv(path)
    topns = sorted(df["topn"].unique())
    if topn is None:
        topn = topns[0]
    elif topn not in topns:
        raise SystemExit(
            f"[error] topn {topn} not in {path} (have {topns}); pass --topn")
    sub = df[df.topn == topn].set_index("tissue")
    if metric not in sub.columns:
        raise SystemExit(f"[error] metric '{metric}' not a column of {path}")
    return sub[metric], topn


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pos", required=True,
                    help="outdir of the POS-set track run (has ff_tissue_results.csv)")
    ap.add_argument("--neg", required=True,
                    help="outdir of the NEG-set track run")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--out", default="ff_tissue_dumbbell.png",
                    help="output PNG filename (written under --outdir)")
    ap.add_argument("--topn", type=int, default=None,
                    help="top-N cutoff to read from each results table "
                         "(default: smallest available, usually the sharpest)")
    ap.add_argument("--metric", default="spec_z",
                    help="results column to contrast (default spec_z)")
    ap.add_argument("--top-k", type=int, default=25,
                    help="show the K tissues with the largest |POS-NEG| "
                         "contrast (default 25). Placenta always included.")
    args = ap.parse_args()

    pos, tn_pos = _load_specz(args.pos, args.topn, args.metric)
    neg, tn_neg = _load_specz(args.neg, args.topn, args.metric)
    if tn_pos != tn_neg:
        print(f"[warn] pos topn={tn_pos} != neg topn={tn_neg}", file=sys.stderr)

    tissues = pos.index.intersection(neg.index)
    df = pd.DataFrame({"pos": pos.reindex(tissues), "neg": neg.reindex(tissues)})
    df["contrast"] = df["pos"] - df["neg"]

    # keep the K largest |contrast|, always retaining placenta
    keep = list(df["contrast"].abs().sort_values(ascending=False).index[:args.top_k])
    if "placenta" in df.index and "placenta" not in keep:
        keep = keep[:args.top_k - 1] + ["placenta"]
    d = df.loc[keep].sort_values("contrast", ascending=True)  # bottom -> top
    n_hidden = len(df) - len(d)

    y = np.arange(len(d))
    fig_h = max(4.5, 0.36 * len(d) + 1.6)
    fig, ax = plt.subplots(figsize=(9.5, fig_h))

    # connecting stems, colored by sign of the contrast (fetal-leaning warm,
    # maternal-leaning cool)
    for yi, (_, row) in zip(y, d.iterrows()):
        c = "#c0392b" if row["contrast"] >= 0 else "#2c6fbb"
        ax.plot([row["neg"], row["pos"]], [yi, yi], color=c, lw=1.6,
                alpha=0.6, zorder=1)
    ax.scatter(d["neg"], y, s=55, color="#2c6fbb", label="NEG set (FF-lowering)",
               zorder=3, edgecolor="white", linewidth=0.6)
    ax.scatter(d["pos"], y, s=55, color="#c0392b", label="POS set (FF-raising)",
               zorder=3, edgecolor="white", linewidth=0.6)
    ax.axvline(0, color="0.4", lw=0.8, zorder=0)

    ax.set_yticks(y)
    ax.set_yticklabels(d.index, fontsize=11)
    ax.set_ylim(-0.7, len(d) - 0.3)
    ax.set_xlabel(f"tissue-specific enrichment ({args.metric}, top-{tn_pos} bins)",
                  fontsize=12)
    ax.set_title("POS vs NEG cell-of-origin contrast (fetal \u2194 maternal)",
                 fontsize=13)
    ax.tick_params(axis="x", labelsize=10)
    ax.legend(fontsize=10, frameon=False, loc="lower right")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    cap = f"Top {len(d)} of {len(df)} shared tissues by |POS-NEG| contrast"
    if n_hidden > 0:
        cap += f"  (+{n_hidden} more, smaller contrast, omitted)"
    fig.suptitle(cap, y=0.995, fontsize=9, color="0.35")
    fig.tight_layout(rect=(0, 0, 1, 0.98))

    os.makedirs(args.outdir, exist_ok=True)
    outpng = os.path.join(args.outdir, args.out)
    fig.savefig(outpng, dpi=150)
    plt.close(fig)
    print(f"[write] {outpng}")

    # console summary
    top_fetal = d.iloc[-1]
    top_mat = d.iloc[0]
    print(f"[result] most fetal-leaning (POS>>NEG): {d.index[-1]}  "
          f"pos={top_fetal['pos']:.2f} neg={top_fetal['neg']:.2f} "
          f"contrast={top_fetal['contrast']:.2f}")
    print(f"[result] most maternal-leaning (NEG>>POS): {d.index[0]}  "
          f"pos={top_mat['pos']:.2f} neg={top_mat['neg']:.2f} "
          f"contrast={top_mat['contrast']:.2f}")


if __name__ == "__main__":
    main()
