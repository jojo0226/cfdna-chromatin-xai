"""
genome.py -- hg38 sequence bookkeeping for GC / mappability-matched null sampling.

Loads per-chromosome FASTA into cumulative GC and N counts so that the GC content
and N-fraction (a mappability proxy) of *any* interval are O(1) lookups. This is
the backbone of the matched null: a query region is only ever compared against
background intervals of equal length drawn from positions with matched GC and
low N-fraction.
"""
from __future__ import annotations
import gzip
import numpy as np


def load_fasta_cumulative(fa_gz_path):
    """Read one gzipped chromosome FASTA -> {'n', 'cumGC', 'cumN'}.

    cumGC[i] = number of G/C bases in seq[:i]; cumN[i] = number of N bases in seq[:i]
    (length n+1 each, so interval [s,e) counts are cum[e]-cum[s]).
    """
    parts = []
    with gzip.open(fa_gz_path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            parts.append(line.strip())
    s = "".join(parts).upper()
    b = np.frombuffer(s.encode("ascii"), dtype=np.uint8)
    is_gc = ((b == ord("G")) | (b == ord("C"))).astype(np.int64)
    is_n = (b == ord("N")).astype(np.int64)
    return {
        "n": int(b.size),
        "cumGC": np.concatenate([[0], np.cumsum(is_gc)]),
        "cumN": np.concatenate([[0], np.cumsum(is_n)]),
    }


def load_genome(fasta_paths):
    """fasta_paths: {chrom: path}. Returns {chrom: {'n','cumGC','cumN'}}."""
    return {c: load_fasta_cumulative(p) for c, p in fasta_paths.items()}


def region_gc_nfrac(chrom, start, end, genome):
    """(gc_fraction over non-N bases, N-fraction over length, effective length)."""
    ga = genome[chrom]
    start = max(0, start)
    end = min(ga["n"], end)
    L = end - start
    if L <= 0:
        return np.nan, 1.0, 0
    gc = ga["cumGC"][end] - ga["cumGC"][start]
    nn = ga["cumN"][end] - ga["cumN"][start]
    eff = L - nn
    gcf = gc / eff if eff > 0 else np.nan
    return float(gcf), float(nn / L), int(L)
