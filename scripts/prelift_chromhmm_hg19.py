#!/usr/bin/env python
"""
prelift_chromhmm_hg19.py -- build an hg19 ChromHMM reference bundle once.

The bundled ChromHMM segmentations are hg38. For a whole hg19 cohort it is
cleaner (and faster) to lift the FIVE static reference BEDs DOWN to hg19 a single
time than to lift every analysis's query bins UP to hg38. This script does that
with round-trip-free interval lifting (segmentations are dense; per-interval RT
QC is unnecessary and slow -- coverage-level accuracy is what matters for
dominant-state calls).

Usage
-----
  python scripts/prelift_chromhmm_hg19.py \
      --in-dir  data \
      --out-dir data_hg19 \
      [--chain hg38ToHg19.over.chain.gz]

Then annotate hg19 bins natively (no per-analysis lifting):
  python scripts/run_chromhmm.py --importance imp_hg19.csv \
      --data-dir data_hg19 --from-build hg38   # coords already match hg19 refs
"""
import argparse
import os

from cfdna_chromatin import references as R
from cfdna_chromatin import liftover as L


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in-dir", default="data", help="dir with hg38 tissue_ACC.bed.gz")
    ap.add_argument("--out-dir", default="data_hg19")
    ap.add_argument("--chain", default=None, help="hg38ToHg19 chain (optional)")
    ap.add_argument("--tissues", default=None, help="comma list (default: all bundled)")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    tissues = (args.tissues.split(",") if args.tissues
               else list(R.REFERENCE_ACCESSIONS.keys()))

    for t in tissues:
        acc = R.REFERENCE_ACCESSIONS[t]
        src = os.path.join(args.in_dir, f"{t}_{acc}.bed.gz")
        if not os.path.exists(src):
            print(f"[skip] {t}: {src} not found")
            continue
        dst = os.path.join(args.out_dir, f"{t}_{acc}.bed.gz")
        qc = L.lift_reference_bed(src, dst, from_build="hg38", to_build="hg19",
                                  chain_path=args.chain, round_trip=False)
        print(f"[{t}] hg38->hg19 mapping rate {qc['mapping_rate']:.3f} "
              f"({qc['n_mapped']}/{qc['n_in']}) -> {dst}")


if __name__ == "__main__":
    main()
