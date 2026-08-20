#!/usr/bin/env python
"""
run_ff_tissue_track.py -- fetal-fraction tissue-of-origin track for a seqFF++ model.

Stage 1 (importance) -> Stage 3 (chromatin explanation), fetal-fraction fork:

  1. Load a per-region importance table from your FF model (|SHAP| per genomic bin,
     or any per-bin importance).  Feature keys map to the shipped 50 kb reference grid.
  2. For the top-N most important bins, ask two questions against the FF openness atlas
     (9 tissues: placenta = fetal signal; liver/endothelial + monocyte/Bcell/CD4/CD8/NK
     = maternal background; K562 = control):

       ABSOLUTE openness   -- mean combined-openness C_<tissue> over the top-N bins,
                              vs a matched random-region null.  High for *every* tissue
                              means the model rides pan-active chromatin (a confound).

       TISSUE-SPECIFIC     -- spec_<tissue> = C_<tissue> - mean over tissues, i.e.
                              cell-of-origin signal.  A genuine FF model should lead
                              with PLACENTA here; a model that only tracks total cfDNA
                              openness will not.

  3. Write results.csv (per topn x tissue z/p), proportions.csv (normalized
     tissue-of-origin share), meta.json, and a two-panel figure.

The atlas ships with the repo:
  analysis/fetal_fraction/reference/ff_openness_atlas_hg19_50kb.csv.gz
It is hg19, 57,633 autosomal 50 kb bins.  Your importance keys must therefore be
hg19 50 kb bins (chrN:start-end or chrN:start:end).  See the AWS recipe in
reference/README.md for extracting |SHAP| from a glmnet FF model in place and
mapping the mixed coverage-bin + 4-mer-motif feature space onto these keys.

Usage:
  python run_ff_tissue_track.py \
      --importance ff_shap_importance.csv \
      --outdir out_ff_tissue
  # atlas path is auto-resolved to reference/; override with --atlas
  # column autodetect: key from key/chrom+start+end/interval/bin/index,
  #                    importance from mean_abs_shap/importance/abs_shap/shap/...
  # or name them explicitly:
  python run_ff_tissue_track.py --importance imp.csv \
      --key-col feature --shap-col beta_abs --outdir out_ff_tissue
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import ff_tissue_proportion as FT

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_ATLAS = os.path.join(HERE, "reference", "ff_openness_atlas_hg19_50kb.csv.gz")


def _tissue_order(results, meta, max_tissues=None):
    """Placenta-first, then descending top-N specific z; optionally trimmed to the
    max_tissues most informative (always keeping placenta)."""
    tissues = meta["tissues"]
    topns = sorted(results["topn"].unique())
    lead = "placenta" if "placenta" in tissues else tissues[0]
    zorder = (results[results.topn == topns[0]]
              .set_index("tissue")["spec_z"].reindex(tissues))
    order = [lead] + [t for t in zorder.sort_values(ascending=False).index if t != lead]
    if max_tissues is not None and len(order) > max_tissues:
        strength = (results.assign(a=results["spec_z"].abs())
                    .groupby("tissue")["a"].max())
        keep = {lead, *[t for t in strength.sort_values(ascending=False).index
                        if t != lead][:max_tissues - 1]}
        order = [t for t in order if t in keep]
    return order, topns


def _draw_panel_a(ax, results, order, topns, show_labels, rot, ha, tick_fs):
    x = np.arange(len(order))
    w = 0.8 / len(topns)
    for i, tn in enumerate(topns):
        sub = results[results.topn == tn].set_index("tissue").reindex(order)
        ax.bar(x + i * w, sub["spec_z"].values, w, label=f"top-{tn}")
    ax.axhline(0, color="0.4", lw=0.8)
    ax.set_xticks(x + w * (len(topns) - 1) / 2)
    ax.set_xticklabels(order if show_labels else [], rotation=rot, ha=ha, fontsize=tick_fs)
    ax.set_xlim(-0.5, len(order) - 0.5 + w * len(topns))
    ax.set_ylabel("tissue-specific enrichment (z vs matched null)", fontsize=11)
    ax.set_title("A. Cell-of-origin: spec_<tissue> of top-N bins", fontsize=12)
    ax.tick_params(axis="y", labelsize=10)
    ax.legend(fontsize=10, frameon=False, ncol=len(topns))


def _draw_panel_b(ax, proportions, order, topns, show_labels, rot, ha, tick_fs):
    x = np.arange(len(order))
    w = 0.8 / len(topns)
    for i, tn in enumerate(topns):
        ax.bar(x + i * w, proportions[tn].reindex(order).values, w, label=f"top-{tn}")
    ax.set_xticks(x + w * (len(topns) - 1) / 2)
    ax.set_xticklabels(order if show_labels else [], rotation=rot, ha=ha, fontsize=tick_fs)
    ax.set_xlim(-0.5, len(order) - 0.5 + w * len(topns))
    ax.set_ylabel("tissue-of-origin share", fontsize=11)
    ax.set_title("B. Normalized proportion (positive specificity)", fontsize=12)
    ax.tick_params(axis="y", labelsize=10)
    ax.legend(fontsize=10, frameon=False, ncol=len(topns))


def _plot(results, proportions, meta, outpng, max_tissues=None, layout="auto",
          tick_fs=None):
    """Render the cell-of-origin figure.

    layout:
      "row"      -- panels A|B side by side (compact; best for few tissues)
      "col"      -- panels A over B, each panel labeled (best for many tissues)
      "separate" -- two standalone PNGs (<stem>_A.png, <stem>_B.png)
      "auto"     -- "row" for <=16 tissues, else "col"
    tick_fs:  explicit tissue-label font size; None -> size-adaptive default.
    """
    order, topns = _tissue_order(results, meta, max_tissues)
    n = len(order)
    if layout == "auto":
        layout = "row" if n <= 16 else "col"
    # label geometry scales with tissue count. Tilted (not vertical) labels let
    # the font grow: angled text spreads its footprint diagonally instead of
    # stacking straight down, so the 73-tissue atlas stays legible at a larger
    # point size than 90 deg rotation allowed.
    tick_fs = tick_fs if tick_fs is not None else (
        12 if n > 40 else (13 if n > 24 else 14))
    # tilt rather than fully vertical for large atlases; angled labels need
    # horizontal room ~ text_len * cos(angle), so widen per-tissue spacing.
    rot, ha = (55, "right") if n > 16 else (45, "right")
    # per-tissue width scales with both font size and the horizontal footprint
    # of a tilted label (bigger for shallower tilt / larger font)
    per_tissue = max(0.52, 0.045 * tick_fs)
    panel_w = max(6.0, per_tissue * n + 1.8)
    # tilted long labels drop further below the axis than vertical ones, so the
    # bottom margin needs to be deeper still to avoid clipping
    bottom = 0.52 if n > 40 else (0.40 if n > 16 else 0.24)

    if layout == "separate":
        import os
        stem, ext = os.path.splitext(outpng)
        outs = []
        for tag, drawer, src in (("A", _draw_panel_a, results),
                                 ("B", _draw_panel_b, proportions)):
            fig, ax = plt.subplots(figsize=(panel_w, 5.0))
            drawer(ax, src, order, topns, True, rot, ha, tick_fs)
            fig.tight_layout()
            fig.subplots_adjust(bottom=bottom)
            p = f"{stem}_{tag}{ext}"
            fig.savefig(p, dpi=150)
            plt.close(fig)
            outs.append(p)
        return outs

    if layout == "col":
        # stacked A over B; label BOTH panels so each is readable on its own.
        # (not sharex -> panel A prints its tissue names directly beneath it)
        # tilted labels drop further than vertical, so panel A's names need more
        # room before panel B starts: taller figure + larger inter-panel gap.
        fig_h = 13.5 if n > 40 else (11.0 if n > 16 else 8.5)
        hspace = 1.25 if n > 40 else (0.85 if n > 16 else 0.35)
        fig, axes = plt.subplots(2, 1, figsize=(panel_w, fig_h))
        _draw_panel_a(axes[0], results, order, topns, True, rot, ha, tick_fs)
        _draw_panel_b(axes[1], proportions, order, topns, True, rot, ha, tick_fs)
        fig.tight_layout()
        fig.subplots_adjust(bottom=bottom / 2 + 0.06, hspace=hspace)
        fig.savefig(outpng, dpi=150)
        plt.close(fig)
        return [outpng]

    # layout == "row"
    fig, axes = plt.subplots(1, 2, figsize=(2 * panel_w, 4.8))
    _draw_panel_a(axes[0], results, order, topns, True, rot, ha, tick_fs)
    _draw_panel_b(axes[1], proportions, order, topns, True, rot, ha, tick_fs)
    fig.tight_layout()
    fig.subplots_adjust(bottom=bottom)
    fig.savefig(outpng, dpi=150)
    plt.close(fig)
    return [outpng]


def _lollipop_order(results, meta, top_k):
    """Top-K tissues by max |spec_z| across cutoffs (placenta always kept),
    returned sorted by spec_z at the smallest cutoff ASCENDING so the strongest
    enrichment sits at the TOP of a horizontal (barh-style) axis."""
    tissues = meta["tissues"]
    topns = sorted(results["topn"].unique())
    lead = "placenta" if "placenta" in tissues else None
    strength = (results.assign(a=results["spec_z"].abs())
                .groupby("tissue")["a"].max())
    keep = list(strength.sort_values(ascending=False).index[:top_k])
    if lead is not None and lead not in keep:
        keep = keep[:top_k - 1] + [lead]
    z0 = (results[results.topn == topns[0]]
          .set_index("tissue")["spec_z"].reindex(keep))
    order = list(z0.sort_values(ascending=True).index)  # bottom->top
    n_hidden = len(tissues) - len(order)
    return order, topns, n_hidden


def _draw_lollipop(ax, series_by_topn, order, topns, xlabel, title,
                   tick_fs, ref_is_stem=True):
    """Horizontal dot-and-stem: one row per tissue (labels on y, always
    horizontal), a stem from 0 to the reference-cutoff value, and a colored
    dot per top-N cutoff so cutoff stability reads at a glance."""
    y = np.arange(len(order))
    cmap = plt.get_cmap("viridis")
    colors = [cmap(i / max(len(topns) - 1, 1)) for i in range(len(topns))]
    ref = topns[0]
    if ref_is_stem:
        ref_vals = series_by_topn[ref].reindex(order).values
        ax.hlines(y, 0, ref_vals, color="0.7", lw=1.4, zorder=1)
    for i, tn in enumerate(topns):
        vals = series_by_topn[tn].reindex(order).values
        ax.scatter(vals, y, s=42, color=colors[i], label=f"top-{tn}",
                   zorder=3, edgecolor="white", linewidth=0.5)
    ax.axvline(0, color="0.4", lw=0.8, zorder=0)
    ax.set_yticks(y)
    ax.set_yticklabels(order, fontsize=tick_fs)
    ax.set_ylim(-0.7, len(order) - 0.3)
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_title(title, fontsize=13)
    ax.tick_params(axis="x", labelsize=10)
    ax.legend(fontsize=9, frameon=False, ncol=len(topns), loc="lower right")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


def _plot_lollipop(results, proportions, meta, outpng, top_k=20, tick_fs=None):
    """Two-panel horizontal lollipop: A = tissue-specific enrichment (spec_z),
    B = tissue-of-origin share. Top-K tissues only, sorted, labels on the
    y-axis so they never need rotation. The long tail of near-zero tissues is
    dropped (noted in the caption) instead of crowding the axis."""
    order, topns, n_hidden = _lollipop_order(results, meta, top_k)
    n = len(order)
    tick_fs = tick_fs if tick_fs is not None else (11 if n > 24 else 12)
    fig_h = max(4.5, 0.34 * n + 1.6)

    specz = {tn: results[results.topn == tn].set_index("tissue")["spec_z"]
             for tn in topns}
    fig, axes = plt.subplots(1, 2, figsize=(13.0, fig_h), sharey=True)
    _draw_lollipop(axes[0], specz, order, topns,
                   "tissue-specific enrichment (z vs matched null)",
                   "A. Cell-of-origin (spec_z), top-N bins", tick_fs)
    _draw_lollipop(axes[1], {tn: proportions[tn] for tn in topns}, order, topns,
                   "tissue-of-origin share",
                   "B. Normalized proportion", tick_fs, ref_is_stem=True)
    lead = order[-1] if order else "?"
    cap = f"Top {n} of {n + n_hidden} tissues by |spec_z|"
    if n_hidden > 0:
        cap += f"  (+{n_hidden} more near zero, omitted)"
    fig.suptitle(cap, y=0.995, fontsize=10, color="0.35")
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(outpng, dpi=150)
    plt.close(fig)
    return [outpng]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--importance", required=True,
                    help="per-bin importance table (CSV); |SHAP| or any per-region score")
    ap.add_argument("--atlas", default=DEFAULT_ATLAS,
                    help="FF openness atlas csv[.gz] (default: shipped reference/)")
    ap.add_argument("--key-col", default=None,
                    help="column holding the bin key (default: autodetect)")
    ap.add_argument("--shap-col", default=None,
                    help="column holding importance (default: autodetect |SHAP|)")
    ap.add_argument("--prefix", default="C_",
                    help="atlas track prefix: C_ combined (default), D_ DNase, "
                         "H_ histone, ACC_ accessibility")
    ap.add_argument("--topns", default="500,1000,2000,5000",
                    help="comma-separated top-N cutoffs")
    ap.add_argument("--n-perm", type=int, default=2000,
                    help="matched random-region null permutations")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--outdir", default="out_ff_tissue")
    ap.add_argument("--max-tissues", type=int, default=None,
                    help="cap x-axis to the N most informative tissues (always "
                         "incl. placenta); use with the 73-tissue atlas to keep "
                         "labels legible, e.g. --max-tissues 20. Default: show all.")
    ap.add_argument("--tick-fs", type=float, default=None,
                    help="tissue-label font size (points). Default: adaptive "
                         "(10-12 by tissue count). Raise it, e.g. --tick-fs 14, "
                         "if names are still hard to read.")
    ap.add_argument("--style", choices=["bars", "lollipop"], default="bars",
                    help="'bars' = grouped bars over all/max tissues (default); "
                         "'lollipop' = horizontal dot-and-stem over the top-K "
                         "tissues sorted by |spec_z|, labels always horizontal, "
                         "near-zero tail dropped. Best for the 73-tissue atlas.")
    ap.add_argument("--top-k", type=int, default=20,
                    help="lollipop only: number of top tissues (by |spec_z|) to "
                         "show; placenta always included. Default 20.")
    ap.add_argument("--layout", choices=["auto", "row", "col", "separate"],
                    default="auto",
                    help="panel arrangement: 'row' = A|B side by side; "
                         "'col' = A over B with one shared tissue-label row "
                         "(clean for many tissues); 'separate' = two standalone "
                         "PNGs (<stem>_A.png/_B.png); 'auto' = row for <=16 "
                         "tissues else col. Default: auto.")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    topns = tuple(int(x) for x in args.topns.split(","))

    print(f"[load] atlas: {args.atlas}", file=sys.stderr)
    atlas = pd.read_csv(args.atlas)
    print(f"[load] importance: {args.importance}", file=sys.stderr)
    imp_raw = pd.read_csv(args.importance)

    results, proportions, meta = FT.tissue_proportion(
        imp_raw, atlas, topns=topns, prefix=args.prefix,
        n_perm=args.n_perm, seed=args.seed,
        shap_col=args.shap_col, key_col=args.key_col)

    # coverage sanity: how many importance keys landed on the grid
    frac = meta["n_bins_matched"] / max(meta["n_bins_in_importance"], 1)
    meta["matched_fraction"] = frac

    res_path = os.path.join(args.outdir, "ff_tissue_results.csv")
    prop_path = os.path.join(args.outdir, "ff_tissue_proportions.csv")
    meta_path = os.path.join(args.outdir, "ff_tissue_meta.json")
    fig_path = os.path.join(args.outdir, "ff_tissue_track.png")
    results.to_csv(res_path, index=False)
    proportions.to_csv(prop_path)
    with open(meta_path, "w") as fh:
        json.dump(meta, fh, indent=2, default=str)
    if args.style == "lollipop":
        fig_paths = _plot_lollipop(results, proportions, meta, fig_path,
                                   top_k=args.top_k, tick_fs=args.tick_fs)
    else:
        fig_paths = _plot(results, proportions, meta, fig_path,
                          max_tissues=args.max_tissues, layout=args.layout,
                          tick_fs=args.tick_fs)
    for _fp in fig_paths:
        print(f"[write] {_fp}")

    # console summary: leading tissue-specific enrichment at the smallest topn
    tn0 = topns[0]
    r0 = results[results.topn == tn0].sort_values("spec_z", ascending=False)
    lead = r0.iloc[0]
    print(f"\n[coverage] {meta['n_bins_matched']}/{meta['n_bins_in_importance']} "
          f"importance keys on grid ({frac:.1%}); atlas grid {meta['n_atlas_bins']} bins")
    print(f"[result] leading tissue-specific (top-{tn0}): "
          f"{lead['tissue']}  spec_z={lead['spec_z']:.2f}  p={lead['spec_p']:.3g}")
    if lead["tissue"] == "placenta":
        print("[interpret] placenta leads -> attribution carries genuine fetal "
              "tissue-of-origin signal, not just pan-openness.")
    else:
        print(f"[interpret] placenta does NOT lead (top is {lead['tissue']}); "
              "check the ABSOLUTE panel -- attribution may track total cfDNA openness.")
    print(f"\n[write] {res_path}\n[write] {prop_path}\n[write] {meta_path}")


if __name__ == "__main__":
    main()
