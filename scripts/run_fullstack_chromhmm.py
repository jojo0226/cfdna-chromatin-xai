#!/usr/bin/env python
"""
run_fullstack_chromhmm.py -- cell-type-AGNOSTIC ChromHMM enrichment of the top
importance / differential cfDNA regions, using the universal ("full-stack")
ChromHMM annotation (Vu & Ernst, Genome Biol 2022).

Unlike the per-tissue bundle (placenta/neutrophil/...), the full-stack model is a
single genome-wide 100-state segmentation integrated across >1000 datasets, so it
answers "what kind of chromatin are these regions, generally" without picking a
cell type. The 100 states collapse to ~16 functional groups (promoters, enhancers,
transcription, Polycomb, heterochromatin, quiescent, ...).

Native builds for hg19 AND hg38 -> no liftover (use the file matching your data).

Input : importance table (key/chrom-start-end + importance/direction), same
        auto-detection as run_chromhmm.py.
Output: fullstack_composition_top{N}.csv per cutoff, and a groups x cutoffs
        heatmap PNG.
"""
import argparse, csv, gzip, json, os, sys
import numpy as np
import pandas as pd


def load_state_map(csv_path):
    """state_annotations_processed.csv -> (num2group dict, ordered group list)."""
    num2group = {}
    with open(csv_path, newline="") as f:
        for r in csv.DictReader(f):
            num2group[r["state"].strip()] = r["Group"].strip()
    groups = []
    for n in sorted(num2group, key=int):
        g = num2group[n]
        if g not in groups:
            groups.append(g)
    return num2group, groups


def load_fullstack(bed_gz, num2group, group2idx, chroms=None):
    """Load full-stack segments -> {chrom: {starts, ends, gidx}} (group index per seg)."""
    from collections import defaultdict
    tmp = defaultdict(list)
    keep = set(chroms) if chroms else None
    op = gzip.open(bed_gz, "rt")
    with op as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            c = p[0]
            if keep is not None and c not in keep:
                continue
            statenum = p[3].split("_")[0]
            g = num2group[statenum]
            tmp[c].append((int(p[1]), int(p[2]), group2idx[g]))
    seg = {}
    for c, rows in tmp.items():
        rows.sort()
        arr = np.array(rows, dtype=np.int64)
        seg[c] = {"starts": arr[:, 0], "ends": arr[:, 1], "gidx": arr[:, 2]}
    return seg


def dominant_group(chrom, start, end, seg):
    """Group index covering the most bp of [start,end); None if unannotated."""
    if chrom not in seg:
        return None
    s = seg[chrom]
    lo = int(np.searchsorted(s["ends"], start, side="right"))
    hi = int(np.searchsorted(s["starts"], end, side="left"))
    bp = {}
    for i in range(lo, hi):
        ov = min(end, int(s["ends"][i])) - max(start, int(s["starts"][i]))
        if ov > 0:
            k = int(s["gidx"][i])
            bp[k] = bp.get(k, 0) + ov
    if not bp:
        return None
    return max(bp, key=bp.get)


def composition(regions, seg, n_groups):
    """regions: list of (chrom,start,end). Returns fraction vector over groups."""
    counts = np.zeros(n_groups, dtype=float)
    n = 0
    for c, s, e in regions:
        g = dominant_group(c, s, e, seg)
        if g is None:
            continue
        counts[g] += 1
        n += 1
    return (counts / n if n else counts), n


def parse_bin(k):
    chrom, rest = k.split(":")
    a, b = rest.split("-")
    return chrom, int(a), int(b)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--importance", required=True)
    ap.add_argument("--segments", required=True, help="full-stack *_segments.bed.gz (matching build)")
    ap.add_argument("--state-annot", required=True, help="state_annotations_processed.csv")
    ap.add_argument("--top-ns", default="1000,2000,5000")
    ap.add_argument("--outdir", default="out_fullstack")
    ap.add_argument("--title", default=None)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    num2group, groups = load_state_map(args.state_annot)
    group2idx = {g: i for i, g in enumerate(groups)}

    sep = "\t" if args.importance.endswith((".tsv", ".txt", ".bed")) else ","
    df = pd.read_csv(args.importance, sep=sep)
    cols = set(df.columns)
    keycol = next((c for c in ("key", "bin", "feature", "name") if c in cols), None)
    impcol = next((c for c in ("importance", "mean_abs_shap", "abs_shap") if c in cols), None)
    dircol = next((c for c in ("direction", "signed_shap", "signed") if c in cols), None)
    rankcol = impcol or dircol
    if rankcol is None:
        sys.exit("ERROR: need an importance or signed column to rank by")
    df = df.reindex(df[rankcol].abs().sort_values(ascending=False).index)

    if {"chrom", "start", "end"}.issubset(cols):
        recs_all = list(df[["chrom", "start", "end"]].itertuples(index=False, name=None))
    else:
        recs_all = [parse_bin(k) for k in df[keycol]]

    cutoffs = sorted({int(x) for x in args.top_ns.split(",") if x.strip()})
    chroms = sorted({c for c, _, _ in recs_all[:max(cutoffs)]})
    print(f"[fullstack] loading segments for {len(chroms)} chroms ...", file=sys.stderr)
    seg = load_fullstack(args.segments, num2group, group2idx, chroms=chroms)

    mat = np.zeros((len(groups), len(cutoffs)), dtype=float)
    meta_ns = {}
    for j, k in enumerate(cutoffs):
        frac, n = composition(recs_all[:k], seg, len(groups))
        mat[:, j] = frac
        meta_ns[k] = n
        comp_df = pd.DataFrame({"group": groups, "query_frac": frac})
        comp_df = comp_df.sort_values("query_frac", ascending=False)
        comp_df.to_csv(os.path.join(args.outdir, f"fullstack_composition_top{k}.csv"), index=False)
        lead = comp_df.iloc[0]
        print(f"[top {k}] n_annot={n}  lead group: {lead['group']} ({lead['query_frac']:.2%})")

    # drop all-zero groups for the plot
    nz = mat.sum(axis=1) > 0
    plot_groups = [g for g, keep in zip(groups, nz) if keep]
    pmat = mat[nz]

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    vmax = float(pmat.max()) or 1.0
    fig, ax = plt.subplots(figsize=(1.3 * len(cutoffs) + 3.2, 0.36 * len(plot_groups) + 1.4))
    im = ax.imshow(pmat, aspect="auto", cmap="viridis", vmin=0.0, vmax=vmax)
    ax.set_xticks(range(len(cutoffs)))
    ax.set_xticklabels([f"top {k}\n(n={meta_ns[k]})" for k in cutoffs], fontsize=8)
    ax.set_yticks(range(len(plot_groups)))
    ax.set_yticklabels(plot_groups, fontsize=8)
    ax.set_ylabel("full-stack ChromHMM group", fontsize=9)
    thr = 0.6 * vmax
    for i in range(len(plot_groups)):
        for j in range(len(cutoffs)):
            v = pmat[i, j]
            if v <= 0:
                continue
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7,
                    color="white" if v >= thr else "black")
    ax.set_xticks(np.arange(-.5, len(cutoffs), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(plot_groups), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.5)
    ax.tick_params(which="minor", length=0)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label("fraction of bins (dominant group)", fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    ax.set_title(args.title or "Universal (cell-type-agnostic) ChromHMM enrichment of top regions",
                 fontsize=10)
    png = os.path.join(args.outdir, "fullstack_composition_heatmap.png")
    fig.savefig(png, dpi=200, bbox_inches="tight")
    json.dump({"cutoffs": cutoffs, "n_annotated": meta_ns, "groups": groups,
               "segments": os.path.basename(args.segments)},
              open(os.path.join(args.outdir, "fullstack_meta.json"), "w"), indent=2)
    print(f"[plot] {png}")


if __name__ == "__main__":
    main()
