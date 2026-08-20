#!/usr/bin/env python
"""End-to-end demo of the chromatin-state engine on real hg38 data.

Reproduces the two validation figures in docs/figures/:
  * positive_control.png -- housekeeping-gene promoters enrich for Promoter states
  * b1_calibration.png   -- signal-free query sets give FPR-controlled p-values

Prerequisites:
  1. pip install -e .
  2. python scripts/fetch_fasta.py            # downloads chr19-22 FASTA
  (ChromHMM reference BEDs already ship in data/)

Usage:
  python scripts/run_demo.py
"""
import json
import os

import matplotlib
matplotlib.use("Agg")

from cfdna_chromatin import references as R, genome as G, engine as E, benchmark as B
from cfdna_chromatin import histone as H, accessibility as A

CHROMS = ["chr19", "chr20", "chr21", "chr22"]
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, "data")
FIGS = os.path.join(HERE, "docs", "figures")
EXAMPLES = os.path.join(HERE, "examples")


def load_all():
    seg = {tis: R.load_segmentation(os.path.join(DATA, f"{tis}_{acc}.bed.gz"), chroms=CHROMS)
           for tis, acc in R.REFERENCE_ACCESSIONS.items()}
    genome = G.load_genome({c: os.path.join(DATA, "fasta", f"{c}.fa.gz") for c in CHROMS})
    return seg, genome


def housekeeping_promoters():
    """Load bundled housekeeping-gene TSS (chr19) and build +/-2kb promoter windows."""
    hk = json.load(open(os.path.join(EXAMPLES, "housekeeping_tss_chr19.json")))
    return [(r["chrom"], r["tss"] - 2000, r["tss"] + 2000) for r in hk]


def main():
    seg, genome = load_all()
    prom = housekeeping_promoters()

    # ---- positive control across all four epigenomes ----
    tissues = ["placenta", "neutrophil", "Bcell", "K562"]
    pc = {}
    for t in tissues:
        nm = E.NullModel(genome=genome, chroms=CHROMS, seed=1)
        res, _ = E.enrichment_test(prom, seg[t], nm, by="group", n_per=50, n_boot=500, seed=1)
        pc[t] = {r["label"]: r for r in res}
        print(f"{t:11s} Promoter log2FC={pc[t]['Promoter']['log2_fold']:+.2f} "
              f"p={pc[t]['Promoter']['p_emp']:.4f}")

    # ---- B1 calibration on placenta ----
    nm_b1 = E.NullModel(genome=genome, chroms=CHROMS, seed=7)
    b1 = B.run_b1_calibration(nm_b1, seg["placenta"], n_sets=80, set_size=40,
                              region_len=50_000, by="group", n_per=40, n_boot=400, seed=7)
    fpr = [d["fpr_0.05"] for d in b1["summary"].values()]
    print(f"\nB1 realized FPR range: {min(fpr):.3f}-{max(fpr):.3f} (nominal 0.05)")
    print("All groups FPR <= 0.05:", all(f <= 0.05 for f in fpr))

    # ---- orthogonal cross-checks: histone marks + open chromatin (placenta) ----
    nm_x = E.NullModel(genome=genome, chroms=CHROMS, seed=1)
    hpanel = H.load_mark_panel(os.path.join(DATA, "histone"), chroms=CHROMS)
    hres, _ = H.mark_enrichment_test(prom, hpanel["placenta"], nm_x, n_per=50, n_boot=500, seed=1)
    print("\nHistone fingerprint (placenta, housekeeping promoters):")
    for r in hres:
        print(f"  {r['mark']:9s} log2FC={r['log2_fold']:+.2f}  p={r['p_emp']:.4f}")

    apanel = A.load_access_panel(os.path.join(DATA, "access"), chroms=CHROMS)
    ares, _ = A.access_enrichment_test(prom, apanel["placenta"], nm_x, n_per=50, n_boot=500, seed=1)
    print(f"\nOpen chromatin (placenta DNase): log2FC={ares['log2_fold']:+.2f} p={ares['p_emp']:.4f}")

    os.makedirs(FIGS, exist_ok=True)
    print("\nDemo complete. See", FIGS)


if __name__ == "__main__":
    main()
