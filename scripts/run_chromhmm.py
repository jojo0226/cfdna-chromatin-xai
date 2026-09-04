#!/usr/bin/env python
"""
run_chromhmm.py -- ChromHMM state annotation of importance-ranked cfDNA bins.

Takes an importance table (key/chrom-start-end + optional importance/direction
columns -- e.g. the diff_importance.csv from fragbed_to_coverage.py, or a
select_signed_shap.py subset) and reports, per reference epigenome, which
ChromHMM 18-states the top bins fall in, and -- if an hg38 genome FASTA is
available -- the GC/length/N-matched state enrichment with empirical p-values.

Input
-----
  --importance FILE   CSV/TSV with a bin column ('key'/'bin'/'feature' like
                      'chr1:0-50000') OR explicit chrom,start,end columns, and
                      optionally an importance/mean_abs_shap and a
                      signed_shap/direction column.
  --top-n N           annotate only the top-N bins by importance (default 500).

Build crossing (hg19 cohorts, e.g. PRAD_14230_hg19)
---------------------------------------------------
  The bundled ChromHMM references are hg38. If your bins are hg19, pass
    --from-build hg19
  and the top-N bins are lifted hg19->hg38 (coordinates only, round-trip QC'd)
  before the reference lookup; the mapping rate is printed. To instead run a
  whole hg19 cohort natively, pre-lift the references once with
  scripts/prelift_chromhmm_hg19.py and point --data-dir at the hg19 bundle.

Reference selection
-------------------
  --panel fetal            use a named panel's tissues (default), OR
  --tissues placenta,Bcell explicit comma list. Available bundled tissues:
                           placenta, neutrophil, Bcell, K562, keratinocyte.

Enrichment (optional)
--------------------
  --fasta-dir DIR   dir of per-chrom hg38 FASTA (chrN.fa.gz). If given, the
                    script adds matched-null state enrichment (log2 fold + p_emp)
                    on top of the raw composition. Without it, only composition
                    is reported (no genome needed).

Outputs (under --outdir)
------------------------
  chromhmm_composition_<tissue>.csv   dominant-state fractions of the top bins
  chromhmm_enrichment_<tissue>.csv    matched-null enrichment (only with --fasta-dir)
  chromhmm_annotated_bins.csv         each top bin + its dominant state per tissue
  chromhmm_meta.json                  run parameters + liftover QC
"""
import argparse
import json
import os
import sys

import pandas as pd

from cfdna_chromatin import chromhmm as CH
from cfdna_chromatin import genome as G


def _read_table(path):
    sep = "\t" if path.endswith((".tsv", ".txt", ".bed")) else ","
    return pd.read_csv(path, sep=sep)


def _pick(cols, *cands):
    for c in cands:
        if c in cols:
            return c
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--importance", required=True, help="importance table (csv/tsv)")
    ap.add_argument("--data-dir", default=None,
                    help="dir with ChromHMM segmentation BEDs (default: package data/)")
    ap.add_argument("--panel", default="fetal")
    ap.add_argument("--tissues", default=None,
                    help="comma list of tissues (overrides --panel)")
    ap.add_argument("--top-n", type=int, default=500)
    ap.add_argument("--top-ns", default=None,
                    help="comma list of cutoffs for a multi-cutoff heatmap, "
                         "e.g. 1000,2000,5000 (overrides --top-n; implies --plot)")
    ap.add_argument("--plot", action="store_true",
                    help="write chromhmm_composition_heatmap.png")
    ap.add_argument("--from-build", default="hg38", choices=["hg38", "hg19"],
                    help="build of the input bins (hg19 -> lifted to hg38)")
    ap.add_argument("--chain", default=None, help="local forward UCSC chain file (from-build -> hg38)")
    ap.add_argument("--back-chain", default=None,
                    help="local reverse chain (hg38 -> from-build); enables round-trip QC offline")
    ap.add_argument("--fasta-dir", default=None,
                    help="hg38 per-chrom FASTA dir (enables matched-null enrichment)")
    ap.add_argument("--by", default="state", choices=["state", "group"])
    ap.add_argument("--n-per", type=int, default=100)
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--outdir", default="out_chromhmm")
    args = ap.parse_args()

    cutoffs = None
    if args.top_ns:
        cutoffs = sorted({int(x) for x in args.top_ns.split(",") if x.strip()})
        args.top_n = max(cutoffs)   # annotate/lift once at the largest cutoff
        args.plot = True

    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_dir = args.data_dir or os.path.join(here, "data")
    os.makedirs(args.outdir, exist_ok=True)

    df = _read_table(args.importance)
    cols = set(df.columns)
    keycol = _pick(cols, "key", "bin", "feature", "name")
    impcol = _pick(cols, "importance", "mean_abs_shap", "abs_shap")
    dircol = _pick(cols, "direction", "signed_shap", "signed")

    # rank by importance (or |signed| if that's all we have) and take top-N
    if impcol:
        df = df.reindex(df[impcol].abs().sort_values(ascending=False).index)
    elif dircol:
        df = df.reindex(df[dircol].abs().sort_values(ascending=False).index)
    top = df.head(args.top_n).copy()
    print(f"[chromhmm] {len(df)} bins in; annotating top {len(top)}"
          f" (key={keycol}, importance={impcol}, direction={dircol})")

    # build a frame the module understands: chrom/start/end + preserved cols
    if {"chrom", "start", "end"}.issubset(cols):
        bins = top
    else:
        if not keycol:
            sys.exit("ERROR: need chrom/start/end columns or a key/bin column")
        bins = top.rename(columns={keycol: "key"})

    tissues = args.tissues.split(",") if args.tissues else None

    genome = None
    if args.fasta_dir:
        chroms = sorted({CH.parse_bin(k)[0] for k in
                         (bins["key"] if "key" in bins else
                          bins.apply(lambda r: f"{r.chrom}:{r.start}-{r.end}", axis=1))})
        paths = {c: os.path.join(args.fasta_dir, f"{c}.fa.gz") for c in chroms
                 if os.path.exists(os.path.join(args.fasta_dir, f"{c}.fa.gz"))}
        if paths:
            print(f"[chromhmm] loading genome for {len(paths)} chroms (matched-null enrichment on)")
            genome = G.load_genome(paths)
        else:
            print("[chromhmm] WARNING: --fasta-dir given but no chrN.fa.gz found; "
                  "skipping enrichment", file=sys.stderr)

    res = CH.annotate_importance_bins(
        bins, data_dir, panel=args.panel, tissues=tissues,
        query_build=args.from_build, ref_build="hg38",
        chain_path=args.chain, back_chain_path=args.back_chain, by=args.by, genome=genome,
        null_kwargs=dict(seed=args.seed),
    )

    if res["liftover_qc"] is not None:
        qc = res["liftover_qc"]
        print(f"[liftover] {args.from_build}->hg38 mapping rate "
              f"{qc['mapping_rate']:.3f} ({qc['n_mapped']}/{qc['n_in']})")

    # write composition per tissue
    for t, comp in res["composition"].items():
        comp.to_csv(os.path.join(args.outdir, f"chromhmm_composition_{t}.csv"), index=False)
        lead = comp.iloc[0]
        print(f"[{t}] top state: {lead['label']} ({lead['query_frac']:.2%} of bins)")

    # write enrichment per tissue (if computed)
    for t, enr in res.get("enrichment", {}).items():
        enr.to_csv(os.path.join(args.outdir, f"chromhmm_enrichment_{t}.csv"), index=False)
        sig = enr[(enr["p_emp"] < 0.05) & (enr["log2_fold"] > 0)]
        if len(sig):
            top_e = sig.iloc[0]
            print(f"[{t}] top enriched state: {top_e['label']} "
                  f"log2fold={top_e['log2_fold']:.2f} p={top_e['p_emp']:.4f}")

    # annotated bins: dominant state per tissue, side by side
    ann = res["regions"][["chrom", "start", "end"]].copy()
    for c in ("key", "importance", "mean_abs_shap", "signed_shap", "direction"):
        if c in res["regions"].columns:
            ann[c] = res["regions"][c].values
    for t in res["tissues"]:
        a = CH.annotate(res["regions"], CH.load_reference_panel(data_dir, [t])[t])
        ann[f"state_{t}"] = a["dominant_state"].values
    ann.to_csv(os.path.join(args.outdir, "chromhmm_annotated_bins.csv"), index=False)

    # ---- optional composition heatmap (single or multi-cutoff) --------------
    if args.plot:
        seg_panel = CH.load_reference_panel(data_dir, res["tissues"])
        regions = res["regions"]           # importance-ranked, post-lift
        if cutoffs:
            # clamp cutoffs to the annotated set; drop duplicates after clamping
            eff = sorted({min(k, len(regions)) for k in cutoffs})
            comps = CH.multi_cutoff_composition(regions, seg_panel, eff, by=args.by)
            mats = {k: CH.composition_matrix(comps[k], by=args.by) for k in eff}
            for k in eff:
                for t, c in comps[k].items():
                    c.to_csv(os.path.join(
                        args.outdir, f"chromhmm_composition_top{k}_{t}.csv"), index=False)
        else:
            mats = {args.top_n: CH.composition_matrix(res["composition"], by=args.by)}
        png = os.path.join(args.outdir, "chromhmm_composition_heatmap.png")
        CH.plot_composition_heatmaps(
            mats, by=args.by, outpath=png,
            suptitle=f"ChromHMM {args.by} composition of top importance bins "
                     f"({args.panel} panel)")
        print(f"[plot] {png}")

    meta = {"n_in": int(len(df)), "n_annotated": int(len(res["regions"])),
            "top_n": args.top_n, "from_build": args.from_build,
            "panel": args.panel, "tissues": res["tissues"], "by": args.by,
            "enrichment": bool(genome is not None),
            "liftover_qc": res["liftover_qc"]}
    with open(os.path.join(args.outdir, "chromhmm_meta.json"), "w") as fh:
        json.dump(meta, fh, indent=2)
    print(f"[write] {args.outdir}/  (composition + annotated_bins"
          + (" + enrichment" if genome is not None else "") + ")")


if __name__ == "__main__":
    main()
