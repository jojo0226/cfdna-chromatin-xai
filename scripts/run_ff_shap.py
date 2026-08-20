#!/usr/bin/env python
"""
Fetal-fraction SHAP -> placenta cell-of-origin explanation (seqFF++).

End-to-end Stage 1 (importance) -> Stage 3 (chromatin explanation) for the
fetal-fraction application, wired to the shipped `fetal` panel (placenta =
signal / cell-of-origin; neutrophil+Bcell = maternal blood background; K562 =
control). It consumes EITHER your own fitted FF model or a precomputed SHAP
matrix, so no cfDNA data has to leave your machine to reuse this tool.

Two input modes
---------------
  (A) model + feature matrix   -- compute SHAP here:
        python scripts/run_ff_shap.py --model ff_model.pkl --matrix X.csv
      X.csv: rows = samples, columns = genomic-bin names ('chrX:start-end').
      Model is loaded with joblib/pickle; regressor (FF) or classifier.

  (B) precomputed SHAP matrix  -- skip straight to the biology:
        python scripts/run_ff_shap.py --shap shap_ff.csv
      shap_ff.csv: rows = samples, columns = genomic-bin names.

Build crossing (hg19 cohorts)
-----------------------------
  The reference layer is hg38. If your bins are hg19, pass --from-build hg19:
  the top ranked bins are lifted hg19->hg38 (coordinates only, QC'd) before the
  reference lookup, and the mapping rate is reported. The recommended
  alternative -- pre-lifting the references down to hg19 once -- is done with
  cfdna_chromatin.liftover.lift_reference_bed() and then run with --from-build hg38
  against the hg19 bundle.

Positive control
-----------------
  For a real FF model the important regions should ENRICH placenta-specific open
  chromatin (fetal cell-of-origin). The script reports the placenta(signal)-vs-
  blood(background) contrast and flags whether the expected positive fires.
"""
import argparse
import os
import pickle
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from cfdna_chromatin import (  # noqa: E402
    histone as H, selection as SEL, liftover as LO,
)

HERE = os.path.join(os.path.dirname(__file__), "..")
HIST_DIR = os.path.join(HERE, "data", "histone")
# reference tracks are shipped for chr19-22 (see docs/data_provenance.md)
REF_CHROMS = ("chr19", "chr20", "chr21", "chr22")


def _load_model(path):
    try:
        import joblib
        return joblib.load(path)
    except Exception:
        with open(path, "rb") as fh:
            return pickle.load(fh)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--shap", help="precomputed SHAP CSV (samples x bins)")
    src.add_argument("--model", help="fitted FF model (joblib/pickle)")
    ap.add_argument("--matrix", help="feature matrix CSV (samples x bins); required with --model")
    ap.add_argument("--panel", default="fetal")
    ap.add_argument("--top", type=int, default=30, help="top-N bins to display")
    ap.add_argument("--from-build", default="hg38", choices=["hg38", "hg19"],
                    help="genome build of the bin coordinates (hg19 -> lift to hg38)")
    ap.add_argument("--explainer", default="auto",
                    choices=["auto", "tree", "linear", "kernel"])
    ap.add_argument("--out", help="write ranked bins to this CSV")
    args = ap.parse_args()

    # ---- Stage 1: obtain the SHAP matrix --------------------------------
    if args.shap:
        sv = pd.read_csv(args.shap, index_col=0)
        print(f"SHAP matrix (loaded): {sv.shape[0]} samples x {sv.shape[1]} bins")
    else:
        if not args.matrix:
            ap.error("--matrix is required with --model")
        model = _load_model(args.model)
        X = pd.read_csv(args.matrix, index_col=0)
        print(f"feature matrix: {X.shape[0]} samples x {X.shape[1]} features")
        sv = SEL.shap_from_model(model, X, explainer=args.explainer)
        print(f"SHAP computed via {args.explainer} explainer: {sv.shape}")

    ranked = SEL.rank_by_shap(sv)

    # ---- optional build crossing hg19 -> hg38 ---------------------------
    if args.from_build == "hg19":
        lifted, qc = LO.lift_regions(
            ranked.reset_index().rename(columns={"index": "bin"}),
            from_build="hg19", to_build="hg38", round_trip=True)
        print(f"\nliftover hg19->hg38: mapped {qc['n_mapped']}/{qc['n_in']} "
              f"({qc['mapping_rate']:.1%}); dropped "
              f"unmapped={qc['dropped_unmapped']} diff_chrom={qc['dropped_diff_chrom']} "
              f"degenerate={qc['dropped_degenerate']} round_trip={qc['dropped_round_trip']}")
        ranked = lifted.set_index("bin") if "bin" in lifted.columns else lifted

    ranked = ranked[ranked["chrom"].isin(REF_CHROMS)]
    print(f"bins on reference chromosomes ({','.join(REF_CHROMS)}): {len(ranked)}")
    if len(ranked) == 0:
        print("no bins on reference chromosomes -- nothing to explain "
              "(ship more reference chromosomes, or check coordinates/build).")
        return

    # ---- Stage 3: placenta cell-of-origin explanation -------------------
    hpanel = H.load_mark_panel(HIST_DIR, chroms=REF_CHROMS)
    res, J = SEL.compartment_importance_test(ranked, hpanel, args.panel)
    print(f"\ncompartment_importance_test (panel={args.panel}, n={res['n']}):")
    # for the fetal panel, signal_specific == placenta-specific contrast
    sig = res["signal_specific"]
    print(f"  placenta(signal)-specific : rho={sig['rho']:+.3f}  p={sig['p']:.3f}")
    print(f"  blood(background)-specific: rho={-sig['rho']:+.3f}  p={sig['p']:.3f}")

    if sig["rho"] > 0 and sig["p"] < 0.05:
        print("\n=> POSITIVE CONTROL FIRES: FF-important regions concentrate in "
              "PLACENTA-specific chromatin (fetal cell-of-origin).")
    elif sig["rho"] < 0 and sig["p"] < 0.05:
        print("\n=> important regions concentrate in maternal BLOOD chromatin, "
              "not placenta -- the FF model may be riding a background/coverage "
              "artifact. Inspect before trusting.")
    else:
        print("\n=> no significant compartment preference (diffuse/weak signal).")

    cols = ["chrom", "start", "end", "importance", "direction"]
    print(f"\nTop {args.top} bins by |SHAP|:")
    print(ranked.head(args.top)[cols].to_string())

    if args.out:
        ranked.to_csv(args.out)
        print(f"\nwrote ranked bins -> {args.out}")


if __name__ == "__main__":
    main()
