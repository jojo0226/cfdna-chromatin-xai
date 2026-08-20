#!/usr/bin/env python3
"""fetch_ff_atlas_beds.py -- reproduce the 73-tissue FF openness atlas BEDs
from PUBLIC ENCODE (no S3, no credentials, no archive restore).

The curated73 atlas was built by querying ENCODE's public API for released
hg19 DNase-seq narrowPeak files, keeping <=3 replicates per biosample (ranked
by peak type), curating to `classification == "tissue"` plus a fixed additions
list (immune cell types, a myeloid progenitor, and cfDNA-relevant biosamples),
and finally adding one neutrophil track from LOCAL bed files. This script
reproduces the ENCODE-derived 72 tissues exactly: it re-runs the same query and
selection, downloads each replicate's narrowPeak, converts to the shipped
4-column form (chrom, start, end, signalValue), writes them into
`data/access_ff/`, and records the tissue->accession mapping into a manifest so
the repo becomes self-describing.

The 73rd track (neutrophil) is NOT on ENCODE hg19 as a tissue biosample -- it
was built from three local files under ~/Downloads (see build lineage /
add_neutrophil_dnase.py). This script reports it as a known gap; supply those
BEDs separately if you want the full 73.

Usage:
    python scripts/fetch_ff_atlas_beds.py --dry-run     # list tissues+accessions, download nothing
    python scripts/fetch_ff_atlas_beds.py               # fetch all, write BEDs + manifest
    python scripts/fetch_ff_atlas_beds.py --resume      # skip files already on disk

No third-party deps beyond the stdlib + pandas (already a package dep).
"""
from __future__ import annotations

import argparse
import collections
import gzip
import io
import json
import os
import re
import sys
import time
import urllib.request

import pandas as pd

# ---------------------------------------------------------------- constants
# hg19 autosomes only (the atlas grid); chrom filter for downloaded peaks
HG19_CHROMS = {f"chr{i}" for i in range(1, 23)}

# Same ENCODE query the atlas build used.
FIELDS = (
    "&field=accession"
    "&field=biosample_ontology.term_name"
    "&field=biosample_ontology.classification"
    "&field=files.@id"
    "&field=files.file_type"
    "&field=files.assembly"
    "&field=files.output_type"
    "&field=files.status"
    "&field=files.href"
    "&field=files.file_format"
)
QUERY = (
    "https://www.encodeproject.org/search/?type=Experiment&assay_title=DNase-seq"
    "&assembly=hg19&status=released&limit=all&format=json" + FIELDS
)

# Exactly the additions list from the atlas build lineage (non-"tissue" class
# biosamples that are deliberately kept for cfDNA relevance).
IMMUNE = [
    "CD14-positive monocyte", "B cell", "CD4-positive, alpha-beta T cell",
    "CD8-positive, alpha-beta T cell", "natural killer cell",
]
MYELOID_PROXY = ["common myeloid progenitor, CD34-positive"]
CFDNA_RELEVANT = [
    "trophoblast cell", "endothelial cell of umbilical vein",
    "dermis blood vessel endothelial cell",
]
ADDITIONS = IMMUNE + MYELOID_PROXY + CFDNA_RELEVANT

# Peak-type preference (lower = better) -- same ranking as the build.
_PEAK_RANK = {
    "optimal IDR thresholded peaks": 0,
    "stable peaks": 1,
    "pseudoreplicated peaks": 2,
    "peaks": 3,
}
_KEEP_OUTPUT_TYPES = set(_PEAK_RANK)
MAX_REPS = 3

OUT_DIR = os.path.join("data", "access_ff")
MANIFEST = os.path.join(OUT_DIR, "encode_accessions.json")


# ---------------------------------------------------------------- helpers
def _rank(output_type: str) -> int:
    return _PEAK_RANK.get(output_type, 9)


def _slug(name: str) -> str:
    """Filesystem-safe token for a biosample name (matches restore driver)."""
    return re.sub(r"[^0-9A-Za-z]+", "_", name).strip("_")


def _facc_from_href(href: str) -> str:
    """/files/ENCFF123ABC/@@download/ENCFF123ABC.bed.gz -> ENCFF123ABC"""
    m = re.search(r"(ENCFF[0-9A-Z]+)", href or "")
    return m.group(1) if m else "UNKNOWN"


def _get_json(url: str) -> dict:
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.load(r)


def _download(href: str, tries: int = 3) -> bytes:
    url = "https://www.encodeproject.org" + href
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=180) as r:
                return r.read()
        except Exception:
            if attempt == tries - 1:
                raise
            time.sleep(2 * (attempt + 1))
    raise RuntimeError("unreachable")


def _to_four_col(raw: bytes) -> pd.DataFrame:
    """narrowPeak (gz) -> DataFrame[chrom,start,end,signalValue], hg19 autosomes."""
    rows = []
    with gzip.open(io.BytesIO(raw), "rt") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            c = p[0]
            if c not in HG19_CHROMS:
                continue
            try:
                s, e = int(p[1]), int(p[2])
            except ValueError:
                continue
            sig = 1.0
            if len(p) > 6 and p[6] not in ("", "."):
                try:
                    sig = float(p[6])
                except ValueError:
                    sig = 1.0
            rows.append((c, s, e, sig))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "signalValue"])


# ---------------------------------------------------------------- selection
def build_selection() -> tuple[dict, dict]:
    """Return (selection, classif).

    selection: {tissue_term: [(file_accession, href, output_type), ...]}  (<=3)
    classif:   {tissue_term: classification}
    Reproduces the atlas build: query ENCODE, keep <=3 ranked reps per biosample,
    then restrict to class=="tissue" + ADDITIONS.
    """
    print("[query] ENCODE hg19 DNase-seq released experiments ...", flush=True)
    t = time.time()
    d = _get_json(QUERY)
    graph = d.get("@graph", [])
    print(f"[query] {len(graph)} experiments in {time.time() - t:.0f}s", flush=True)

    by_tissue = collections.defaultdict(list)
    classif = {}
    for g in graph:
        term = g.get("biosample_ontology", {}).get("term_name")
        if not term:
            continue
        classif[term] = g.get("biosample_ontology", {}).get("classification")
        for f in g.get("files", []):
            if (
                f.get("file_format") == "bed"
                and "narrowPeak" in (f.get("file_type") or "")
                and f.get("assembly") == "hg19"
                and f.get("status") == "released"
                and f.get("output_type") in _KEEP_OUTPUT_TYPES
            ):
                by_tissue[term].append(
                    (_facc_from_href(f["href"]), f["href"], f.get("output_type"))
                )

    # <=3 ranked, de-duplicated by file accession
    sel = {}
    for term, lst in by_tissue.items():
        # rank first, then accession as a stable tie-break so re-fetches are
        # reproducible when a biosample has >MAX_REPS files at the same rank
        # (all hg19 DNase narrowPeaks are output_type "peaks", i.e. all tied).
        lst = sorted(lst, key=lambda x: (_rank(x[2]), x[0]))
        seen, keep = set(), []
        for facc, href, ot in lst:
            if facc in seen:
                continue
            seen.add(facc)
            keep.append((facc, href, ot))
            if len(keep) >= MAX_REPS:
                break
        if keep:
            sel[term] = keep

    # curated panel: all class=="tissue" + the fixed additions that resolved
    tissue_only = sorted(t for t in sel if classif.get(t) == "tissue")
    additions = [t for t in ADDITIONS if t in sel]
    panel = tissue_only + additions
    selection = {t: sel[t] for t in panel}
    return selection, classif


# ---------------------------------------------------------------- main
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true",
                    help="list tissues + accessions, download nothing")
    ap.add_argument("--resume", action="store_true",
                    help="skip BEDs already present on disk")
    ap.add_argument("--out", default=OUT_DIR,
                    help=f"output directory (default {OUT_DIR})")
    args = ap.parse_args()

    selection, classif = build_selection()
    n_tissue = sum(1 for t in selection if classif.get(t) == "tissue")
    n_add = len(selection) - n_tissue
    n_files = sum(len(v) for v in selection.values())
    print(f"[panel] {len(selection)} biosamples "
          f"({n_tissue} tissue + {n_add} additions), {n_files} peak files")
    print("[note] neutrophil (73rd track) is built from LOCAL files, not ENCODE "
          "hg19 -- supply data/neutrophil_dnase/*.bed separately for the full 73.")

    if args.dry_run:
        for t in sorted(selection):
            accs = ", ".join(f for f, _, _ in selection[t])
            print(f"  {t}  [{classif.get(t)}]  ({len(selection[t])})  {accs}")
        print(f"[dry-run] would write {n_files} BEDs into {args.out}")
        return 0

    os.makedirs(args.out, exist_ok=True)
    manifest = {"assembly": "hg19", "assay": "DNase-seq", "source": "ENCODE public API",
                "query": QUERY, "max_reps": MAX_REPS, "tissues": {}}
    n_ok = n_skip = n_err = 0
    for i, term in enumerate(sorted(selection), 1):
        slug = _slug(term)
        recs = []
        for facc, href, ot in selection[term]:
            dst = os.path.join(args.out, f"{slug}_DNase_{facc}.bed.gz")
            recs.append({"file_accession": facc, "output_type": ot,
                         "path": os.path.relpath(dst, args.out)})
            if args.resume and os.path.exists(dst):
                n_skip += 1
                continue
            try:
                df = _to_four_col(_download(href))
                with gzip.open(dst, "wt") as fh:
                    df.to_csv(fh, sep="\t", header=False, index=False)
                n_ok += 1
            except Exception as e:
                n_err += 1
                print(f"  ERR {facc} ({term}): {type(e).__name__}: {str(e)[:150]}")
        manifest["tissues"][term] = {"classification": classif.get(term),
                                     "n_rep": len(selection[term]), "files": recs}
        if i % 10 == 0 or i == len(selection):
            print(f"[{i}/{len(selection)}] {term}  "
                  f"(ok={n_ok} skip={n_skip} err={n_err})", flush=True)

    manifest["n_tissues_encode"] = len(selection)
    with open(MANIFEST if args.out == OUT_DIR
              else os.path.join(args.out, "encode_accessions.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"[done] wrote {n_ok} BEDs, skipped {n_skip}, errors {n_err} -> {args.out}")
    print(f"[done] accession manifest -> {MANIFEST}")
    print("[next] un-gitignore data/access_ff/ and commit if you want the BEDs in the repo.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
