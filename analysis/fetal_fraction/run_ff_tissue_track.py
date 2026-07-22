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


def _plot(results, proportions, meta, outpng):
    tissues = meta["tissues"]
    topns = sorted(results["topn"].unique())
    # order tissues by placenta-first then by top-N specific z (readability)
    lead = "placenta" if "placenta" in tissues else tissues[0]
    ref_topn = topns[0]
    zorder = (results[results.topn == ref_topn]
              .set_index("tissue")["spec_z"].reindex(tissues))
    order = [lead] + [t for t in zorder.sort_values(ascending=False).index if t != lead]

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6))
    x = np.arange(len(order))
    w = 0.8 / len(topns)
    # panel A: tissue-specific z (cell-of-origin)
    ax = axes[0]
    for i, tn in enumerate(topns):
        sub = results[results.topn == tn].set_index("tissue").reindex(order)
        ax.bar(x + i * w, sub["spec_z"].values, w, label=f"top-{tn}")
    ax.axhline(0, color="0.4", lw=0.8)
    ax.set_xticks(x + w * (len(topns) - 1) / 2)
    ax.set_xticklabels(order, rotation=45, ha="right")
    ax.set_ylabel("tissue-specific enrichment (z vs matched null)")
    ax.set_title("A. Cell-of-origin: spec_<tissue> of top-N bins")
    ax.legend(fontsize=8, frameon=False)
    # panel B: proportion (share of positive spec)
    ax = axes[1]
    for i, tn in enumerate(topns):
        ax.bar(x + i * w, proportions[tn].reindex(order).values, w, label=f"top-{tn}")
    ax.set_xticks(x + w * (len(topns) - 1) / 2)
    ax.set_xticklabels(order, rotation=45, ha="right")
    ax.set_ylabel("tissue-of-origin share")
    ax.set_title("B. Normalized proportion (positive specificity)")
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(outpng, dpi=150)
    plt.close(fig)


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
    ap.add_argument("--topns", default="500,1000,2000",
                    help="comma-separated top-N cutoffs")
    ap.add_argument("--n-perm", type=int, default=2000,
                    help="matched random-region null permutations")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--outdir", default="out_ff_tissue")
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
    _plot(results, proportions, meta, fig_path)

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
    print(f"\n[write] {res_path}\n[write] {prop_path}\n[write] {meta_path}\n[write] {fig_path}")


if __name__ == "__main__":
    main()
