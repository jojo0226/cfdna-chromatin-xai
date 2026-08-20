#!/usr/bin/env python
"""Fetch ENCODE hg19 DNase-seq peaks for the fetal-fraction (seqFF++) openness panel.

The NIPT / fetal-fraction fork compares tissue-of-origin openness across the cell
types that contribute to cfDNA in pregnancy -- placental (fetal) signal against a
maternal background of hematopoietic + solid-tissue cell types. There is NO tumor
proxy — the fetal (placental) signal is the cell-of-origin of interest.

Accessions live in analysis/fetal_fraction/reference/ff_reference_manifest.json. Each
tissue has up to 3 DNase-seq replicates from distinct ENCODE experiments; they are
z-scored per replicate across bins then averaged (depth/threshold-robust) by
build_openness_atlas.py.

Raw 10-column narrowPeak files download to data/access_ff_raw/ (gitignored) and are
converted to the shipped 4-column (chrom,start,end,signalValue) form in
data/access_ff/. signalValue is column 7 of ENCODE narrowPeak (fold-enrichment),
NOT the saturated integer score in column 5.

neutrophil has no DNase-seq in ENCODE (hg19 or GRCh38), so it ships as a histone-only
tissue (4 active marks; C_neutrophil = z(H)). CMP (CD34+ common myeloid progenitor)
has DNase + 3 histone marks (no H3K27ac narrowPeak in hg19). Both are driven from the
same manifest — DNase accessions in `tissues`, histone in `histone_accessions`.

Usage:
    python scripts/fetch_ff_reference.py
Then rebuild the atlas:
    python -c "import build_openness_atlas as B, glob; ..."   # see reference/README
"""
import gzip
import json
import os
import urllib.request

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANIFEST = os.path.join(HERE, "analysis", "fetal_fraction", "reference",
                        "ff_reference_manifest.json")
RAW = os.path.join(HERE, "data", "access_ff_raw")
OUT = os.path.join(HERE, "data", "access_ff")


def fetch_and_convert(tissue, acc):
    """Download one ENCODE narrowPeak and convert 10-col -> 4-col (col7=signalValue)."""
    raw = os.path.join(RAW, f"{acc}.bed.gz")
    out = os.path.join(OUT, f"{tissue}_DNase_{acc}.bed.gz")
    if os.path.exists(out) and os.path.getsize(out) > 0:
        return "cached"
    if not (os.path.exists(raw) and os.path.getsize(raw) > 0):
        url = f"https://www.encodeproject.org/files/{acc}/@@download/{acc}.bed.gz"
        urllib.request.urlretrieve(url, raw)
    n = 0
    with gzip.open(raw, "rt") as fi, gzip.open(out, "wt") as fo:
        for line in fi:
            p = line.rstrip("\n").split("\t")
            if len(p) < 7:
                continue
            fo.write(f"{p[0]}\t{p[1]}\t{p[2]}\t{p[6]}\n")
            n += 1
    return f"{n} peaks"


HRAW = os.path.join(HERE, "data", "histone_ff_raw")
HOUT = os.path.join(HERE, "data", "histone_ff")


def fetch_and_convert_histone(tissue, mark, acc):
    """Download one ENCODE histone narrowPeak -> 4-col (chrom,start,end,signalValue)."""
    raw = os.path.join(HRAW, f"{acc}.bed.gz")
    out = os.path.join(HOUT, f"{tissue}_{mark}_{acc}.bed.gz")
    if os.path.exists(out) and os.path.getsize(out) > 0:
        return "cached"
    if not (os.path.exists(raw) and os.path.getsize(raw) > 0):
        url = f"https://www.encodeproject.org/files/{acc}/@@download/{acc}.bed.gz"
        urllib.request.urlretrieve(url, raw)
    n = 0
    with gzip.open(raw, "rt") as fi, gzip.open(out, "wt") as fo:
        for line in fi:
            p = line.rstrip("\n").split("\t")
            if len(p) < 7:
                continue
            fo.write(f"{p[0]}\t{p[1]}\t{p[2]}\t{p[6]}\n")
            n += 1
    return f"{n} peaks"


def main():
    man = json.load(open(MANIFEST))
    # DNase
    os.makedirs(RAW, exist_ok=True)
    os.makedirs(OUT, exist_ok=True)
    print("DNase-seq:")
    for tissue, accs in man["tissues"].items():
        for acc in accs:
            msg = fetch_and_convert(tissue, acc)
            print(f"  {tissue:12s} {acc}  {msg}")
    # Histone ChIP-seq (H3K4me3/H3K4me1/H3K27ac/H3K36me3), keyed tissue->mark->[accs]
    if man.get("histone_accessions"):
        os.makedirs(HRAW, exist_ok=True)
        os.makedirs(HOUT, exist_ok=True)
        print("Histone ChIP-seq:")
        for tissue, marks in man["histone_accessions"].items():
            for mark, accs in marks.items():
                for acc in accs:
                    msg = fetch_and_convert_histone(tissue, mark, acc)
                    print(f"  {tissue:12s} {mark:8s} {acc}  {msg}")
    print(f"done -> {OUT} , {HOUT}")


if __name__ == "__main__":
    main()
