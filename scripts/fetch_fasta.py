#!/usr/bin/env python
"""Fetch hg38 per-chromosome FASTA needed for GC/mappability-matched null sampling.

The ChromHMM reference BEDs ship in the repo (data/*.bed.gz); the sequence does not
(it is large and re-downloadable). Run this once:

    python scripts/fetch_fasta.py --chroms chr19 chr20 chr21 chr22

Files land in data/fasta/ (gitignored).
"""
import argparse
import gzip
import os
import time
import urllib.request

UCSC = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/{c}.fa.gz"


def _valid_gzip(path):
    try:
        with gzip.open(path, "rt") as fh:
            for _ in fh:
                pass
        return True
    except Exception:
        return False


def fetch(chrom, out_dir, tries=5):
    dst = os.path.join(out_dir, f"{chrom}.fa.gz")
    if os.path.exists(dst) and _valid_gzip(dst):
        print(f"  {chrom}: cached")
        return dst
    for _ in range(tries):
        try:
            urllib.request.urlretrieve(UCSC.format(c=chrom), dst)
            if _valid_gzip(dst):
                print(f"  {chrom}: {os.path.getsize(dst)/1e6:.1f} MB")
                return dst
        except Exception:
            pass
        if os.path.exists(dst):
            os.remove(dst)
        time.sleep(2)
    raise RuntimeError(f"failed to fetch {chrom}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", nargs="+",
                    default=["chr19", "chr20", "chr21", "chr22"])
    ap.add_argument("--out-dir", default="data/fasta")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    print(f"Fetching {len(args.chroms)} chromosome(s) -> {args.out_dir}")
    for c in args.chroms:
        fetch(c, args.out_dir)
    print("done.")


if __name__ == "__main__":
    main()
