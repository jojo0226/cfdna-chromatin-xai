#!/usr/bin/env python
"""build_repressive_atlas.py -- H3K27me3 (Polycomb) repressive-chromatin atlas.

A companion to build_openness_atlas.py. Every track in the openness atlas is an
*active/open* signal (DNase, activating histone marks, accessibility). This
builds the complementary axis: tissue-specific *repression*, measured by
H3K27me3 (facultative Polycomb) narrowPeak coverage on the SAME hg19 50 kb grid.

WHY a separate atlas, never folded into C_:
  cfDNA coverage and openness are linked through nucleosome occupancy -- open
  chromatin is nucleosome-depleted, fragments more, and is DEPLETED in plasma
  coverage; repressed/compact chromatin is protected and ENRICHED. So a
  repressive mark points the OPPOSITE way from DNase. Averaging it into the
  composite openness score would partly cancel the active signal. Kept apart, it
  gives a two-sided cell-of-origin test: a genuine tissue-of-origin signal is
  enriched in that tissue's ACTIVE marks AND depleted in its REPRESSIVE marks; a
  pan-openness confound shows no coherent repressive-depletion.

Assay choice: H3K27me3 (facultative, cell-type-specific), NOT H3K9me3
(constitutive heterochromatin, tissue-invariant -> would add noise, not identity).

Output columns (per tissue with H3K27me3 available):
  R_<tissue>       z-mean H3K27me3 coverage over replicates (repressive signal)
  n_R_<tissue>     replicate count
  repr_spec_<tissue>   R_<tissue> - mean_over_tissues(R)   (tissue-specific repression)

Usage:
  python build_repressive_atlas.py --probe h3k27me3_probe.json \
      --outdir repressive_atlas --max-reps 3
The --probe JSON is the availability sweep (term -> accessions); if omitted the
script queries ENCODE live for the curated-72 biosample terms.
"""
import argparse
import gzip
import io
import json
import os
import sys
import time
import urllib.request

import numpy as np
import pandas as pd

# reuse the exact grid + binning + aggregation from the openness builder
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_openness_atlas import (  # noqa: E402
    aggregate_replicates,
    bin_grid,
    peaks_to_coverage,
)

ENCODE = "https://www.encodeproject.org"
HDRS = {"Accept": "application/json"}
# output-type preference: replicated call > optimal > peaks > pseudoreplicated
OUTPUT_PREF = ["replicated peaks", "optimal idr thresholded peaks",
               "peaks", "pseudoreplicated peaks"]


def _get(url, params=None, tries=3):
    for k in range(tries):
        try:
            if params:
                from urllib.parse import urlencode
                url2 = url + "?" + urlencode(params)
            else:
                url2 = url
            req = urllib.request.Request(url2, headers=HDRS)
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.load(r)
        except Exception as e:  # noqa: BLE001
            if k == tries - 1:
                raise
            time.sleep(1.5 * (k + 1))
            _ = e
    return None


def experiments_for(term):
    """H3K27me3 hg19 experiment accessions for one biosample term."""
    params = {
        "type": "Experiment", "assay_title": "Histone ChIP-seq",
        "target.label": "H3K27me3", "biosample_ontology.term_name": term,
        "assembly": "hg19", "status": "released",
        "files.file_type": "bed narrowPeak", "frame": "object", "limit": "10",
    }
    try:
        j = _get(ENCODE + "/search/", params)
    except Exception:  # noqa: BLE001
        return []
    return [g["accession"] for g in j.get("@graph", [])]


def best_peak_href(expt_acc):
    """Pick one released hg19 bed narrowPeak per experiment, by output-type pref."""
    j = _get(f"{ENCODE}/experiments/{expt_acc}/")
    cands = []
    for f in j.get("files", []):
        if (f.get("file_type") == "bed narrowPeak"
                and f.get("assembly") == "hg19"
                and f.get("status") == "released"):
            ot = (f.get("output_type") or "").lower()
            rank = OUTPUT_PREF.index(ot) if ot in OUTPUT_PREF else len(OUTPUT_PREF)
            cands.append((rank, f.get("href")))
    if not cands:
        return None
    cands.sort(key=lambda x: x[0])
    return cands[0][1]


def load_narrowpeak(href):
    """Download an ENCODE narrowPeak bed.gz -> chrom/start/end/signalValue (autosomes)."""
    url = ENCODE + href
    with urllib.request.urlopen(urllib.request.Request(url, headers={}), timeout=120) as r:
        raw = r.read()
    with gzip.open(io.BytesIO(raw), "rt") as fh:
        df = pd.read_csv(fh, sep="\t", header=None,
                         usecols=[0, 1, 2, 6],
                         names=["chrom", "start", "end", "signalValue"])
    auto = {f"chr{i}" for i in range(1, 23)}
    return df[df["chrom"].isin(auto)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--probe", default=None,
                    help="JSON from the availability sweep (term->accessions)")
    ap.add_argument("--terms", default=None,
                    help="comma-separated biosample terms (overrides --probe)")
    ap.add_argument("--outdir", default="repressive_atlas")
    ap.add_argument("--bin-size", type=int, default=50_000)
    ap.add_argument("--max-reps", type=int, default=3)
    ap.add_argument("--weight", default="signalValue",
                    help="signalValue (fold-enrichment-weighted) or none")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    weight = None if args.weight.lower() in ("none", "") else args.weight

    # resolve tissue -> experiment accessions
    if args.terms:
        terms = [t.strip() for t in args.terms.split(",")]
        tissue_expts = {t: experiments_for(t) for t in terms}
    elif args.probe:
        probe = json.load(open(args.probe))
        tissue_expts = {r["term"]: r.get("accessions", [])
                        for r in probe if r.get("n")}
    else:
        raise SystemExit("provide --probe or --terms")

    grid = bin_grid(args.bin_size, "hg19")
    out = grid.copy()
    kept, dropped = [], []

    for ti, (tissue, expts) in enumerate(sorted(tissue_expts.items()), 1):
        expts = expts[:args.max_reps]
        vecs = []
        for e in expts:
            href = best_peak_href(e)
            if not href:
                continue
            try:
                pk = load_narrowpeak(href)
            except Exception as ex:  # noqa: BLE001
                print(f"  [skip] {tissue} {e}: {ex}")
                continue
            vecs.append(peaks_to_coverage(pk, grid, args.bin_size, weight=weight))
        if not vecs:
            dropped.append(tissue)
            print(f"[{ti}/{len(tissue_expts)}] {tissue}: NO usable files -> dropped")
            continue
        track = aggregate_replicates(np.vstack(vecs), method="zmean")
        col = _clean(tissue)
        out[f"R_{col}"] = track
        out[f"n_R_{col}"] = len(vecs)
        kept.append(col)
        print(f"[{ti}/{len(tissue_expts)}] {tissue} -> R_{col}  ({len(vecs)} reps)")

    # repressive specificity: R_<t> - mean over tissues (removes pan-repression)
    R_cols = [f"R_{c}" for c in kept]
    Rmat = out[R_cols].to_numpy()
    panmean = np.nanmean(Rmat, axis=1, keepdims=True)
    for j, c in enumerate(kept):
        out[f"repr_spec_{c}"] = Rmat[:, j] - panmean[:, 0]

    atlas_path = os.path.join(args.outdir,
                              "ff_repressive_atlas_hg19_50kb_H3K27me3.csv.gz")
    out.to_csv(atlas_path, index=False, compression="gzip")

    manifest = {
        "assembly": "hg19", "bin_size": args.bin_size, "n_bins": len(grid),
        "mark": "H3K27me3", "assay": "Histone ChIP-seq (bed narrowPeak)",
        "aggregation": "zmean over replicates", "weight": args.weight,
        "n_tissues": len(kept), "tissues": kept, "dropped": dropped,
        "columns": {"R_<tissue>": "z-mean H3K27me3 coverage (repressive signal)",
                    "repr_spec_<tissue>": "R_<tissue> - mean_over_tissues(R)"},
        "note": ("SEPARATE repressive axis; NOT merged into openness C_. "
                 "Two-sided test: genuine origin = active-enriched AND "
                 "repressive-depleted."),
        "provenance": "ENCODE hg19 H3K27me3 narrowPeak, <=%d reps/biosample"
                      % args.max_reps,
    }
    json.dump(manifest, open(os.path.join(args.outdir, "manifest_repressive.json"),
                             "w"), indent=1)
    print(f"\n[atlas] {atlas_path}  ({len(grid)} bins x {len(kept)} tissues)")
    print(f"[manifest] {args.outdir}/manifest_repressive.json")
    if dropped:
        print(f"[dropped] {len(dropped)}: {dropped}")


def _clean(name):
    """ENCODE free-text biosample -> compact column token."""
    s = name.lower()
    repl = {"cd14-positive monocyte": "monocyte",
            "cd4-positive, alpha-beta t cell": "CD4",
            "cd8-positive, alpha-beta t cell": "CD8",
            "natural killer cell": "NK", "b cell": "Bcell",
            "common myeloid progenitor, cd34-positive": "CMP",
            "endothelial cell of umbilical vein": "HUVEC",
            "dermis blood vessel endothelial cell": "dermal_endo",
            "trophoblast cell": "trophoblast"}
    if s in repl:
        return repl[s]
    return (name.replace(",", "").replace(" ", "_")
            .replace("'", "").replace("(", "").replace(")", ""))


if __name__ == "__main__":
    main()
