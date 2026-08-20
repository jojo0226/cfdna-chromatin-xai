#!/usr/bin/env python
"""Fetch genome-wide ENCODE DNase-seq peak files for the accessibility panel.

The repository ships chr19-22 subsets of each peak file (data/access/), which is
all the tests and the demo need. This script downloads the full-genome versions
of the same files (one ENCODE DNase-seq experiment per tissue) for users who want
to run the accessibility layer genome-wide.

Files and accessions are taken from data/access/manifest.json, so this stays in
sync with whatever is shipped. Downloads go to data/access_raw/ (gitignored).

Usage:
    python scripts/fetch_access.py
"""
import json
import os
import shutil
import urllib.request

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, "data")
RAW = os.path.join(DATA, "access_raw")


def main():
    os.makedirs(RAW, exist_ok=True)
    man = json.load(open(os.path.join(DATA, "access", "manifest.json")))
    for rec in man:
        acc = rec["file_acc"]
        url = f"https://www.encodeproject.org/files/{acc}/@@download/{acc}.bed.gz"
        dst = os.path.join(RAW, f"{acc}.bed.gz")
        if os.path.exists(dst) and os.path.getsize(dst) > 0:
            print(f"  have {acc}")
            continue
        proxy = f" (proxy for {rec['proxy_for']})" if rec.get("proxy_for") else ""
        print(f"  fetching {rec['tissue']:10s} DNase-seq {acc}{proxy} ...")
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req, timeout=300) as r, open(dst, "wb") as f:
            shutil.copyfileobj(r, f)
    print(f"done -> {RAW}")


if __name__ == "__main__":
    main()
