#!/usr/bin/env python
"""Fetch genome-wide ENCODE histone-mark peak files for the reference panel.

The repository ships chr19-22 subsets of each peak file (data/histone/), which is
all the tests and the demo need. This script downloads the full-genome versions
of the same files (one ENCODE experiment per tissue x mark) for users who want to
run the histone layer genome-wide.

Files and accessions are taken from data/histone/manifest.json, so this stays in
sync with whatever is shipped. Downloads go to data/histone_raw/ (gitignored).

Usage:
    python scripts/fetch_histone.py
"""
import json
import os
import shutil
import urllib.request

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, "data")
RAW = os.path.join(DATA, "histone_raw")
CHROMS = {"chr19", "chr20", "chr21", "chr22"}


def main(subset_only=False):
    os.makedirs(RAW, exist_ok=True)
    man = json.load(open(os.path.join(DATA, "histone", "manifest.json")))
    for rec in man:
        acc = rec["file_acc"]
        url = f"https://www.encodeproject.org/files/{acc}/@@download/{acc}.bed.gz"
        dst = os.path.join(RAW, f"{acc}.bed.gz")
        if os.path.exists(dst) and os.path.getsize(dst) > 0:
            print(f"  have {acc}")
            continue
        print(f"  fetching {rec['tissue']:10s} {rec['mark']:8s} {acc} ...")
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req, timeout=180) as r, open(dst, "wb") as f:
            shutil.copyfileobj(r, f)
    print(f"done -> {RAW}")


if __name__ == "__main__":
    main()
