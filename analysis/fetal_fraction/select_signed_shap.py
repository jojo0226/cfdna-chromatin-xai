#!/usr/bin/env python
"""select_signed_shap.py -- split FF SHAP bins into most-positive / most-negative sets.

Consumes the importance CSV written by extract_ff_shap.R
(columns: key, mean_abs_shap, signed_shap) and, for each cutoff N, writes the
top-N bins with the LARGEST positive signed_shap and the top-N with the most
negative signed_shap. These directed sets are the inputs to the downstream
cell-of-origin / openness analysis run separately per direction.

Signed convention (from the FF model's effective per-bin weight w_bin):
  signed_shap > 0  -> more coverage/openness in the bin raises predicted FF
                      (fetal / placental-leaning direction)
  signed_shap < 0  -> lowers predicted FF (maternal-leaning direction)

Outputs (under --outdir):
  pos_top{N}.csv / neg_top{N}.csv   one per cutoff and direction
                                    (key, mean_abs_shap, signed_shap, rank)
  signed_selection_manifest.json    cutoffs, counts, and value ranges

Usage:
  python select_signed_shap.py --importance ff_shap_importance.csv \
      --topns 5000,2000,1000,500 --outdir ff_signed_sets
"""
import argparse
import json
import os

import pandas as pd


def load_signed(path: str, key_col=None, signed_col=None, abs_col=None) -> pd.DataFrame:
    df = pd.read_csv(path)
    key = key_col or ("key" if "key" in df.columns else df.columns[0])
    if signed_col is None:
        for c in ("signed_shap", "signed", "w_bin", "weight", "mean_shap"):
            if c in df.columns:
                signed_col = c
                break
    if signed_col is None:
        raise KeyError(
            "no signed column found (looked for signed_shap/signed/w_bin/weight/"
            "mean_shap); re-run extract_ff_shap.R which now emits signed_shap, "
            "or pass --signed-col"
        )
    if abs_col is None:
        abs_col = "mean_abs_shap" if "mean_abs_shap" in df.columns else None
    out = pd.DataFrame({"key": df[key].astype(str), "signed_shap": df[signed_col].astype(float)})
    out["mean_abs_shap"] = df[abs_col].astype(float) if abs_col else out["signed_shap"].abs()
    # collapse any duplicate keys, keeping the largest-magnitude signed value
    out = (out.reindex(out["mean_abs_shap"].abs().sort_values(ascending=False).index)
              .drop_duplicates("key"))
    return out.reset_index(drop=True)


def select(df: pd.DataFrame, topns, outdir: str):
    os.makedirs(outdir, exist_ok=True)
    pos = df[df["signed_shap"] > 0].sort_values("signed_shap", ascending=False)
    neg = df[df["signed_shap"] < 0].sort_values("signed_shap", ascending=True)
    manifest = {"n_bins_total": int(len(df)),
                "n_positive": int(len(pos)), "n_negative": int(len(neg)),
                "cutoffs": {}}
    for n in topns:
        for label, src in (("pos", pos), ("neg", neg)):
            sub = src.head(n).copy()
            sub.insert(0, "rank", range(1, len(sub) + 1))
            fn = os.path.join(outdir, f"{label}_top{n}.csv")
            sub[["rank", "key", "mean_abs_shap", "signed_shap"]].to_csv(fn, index=False)
        manifest["cutoffs"][str(n)] = {
            "pos_n": int(min(n, len(pos))),
            "neg_n": int(min(n, len(neg))),
            "pos_signed_range": [float(pos.head(n)["signed_shap"].min()),
                                 float(pos.head(n)["signed_shap"].max())] if len(pos) else None,
            "neg_signed_range": [float(neg.head(n)["signed_shap"].min()),
                                 float(neg.head(n)["signed_shap"].max())] if len(neg) else None,
        }
    with open(os.path.join(outdir, "signed_selection_manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)
    return manifest


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--importance", required=True,
                    help="CSV from extract_ff_shap.R (key,mean_abs_shap,signed_shap)")
    ap.add_argument("--topns", default="5000,2000,1000,500",
                    help="comma-separated cutoffs (default 5000,2000,1000,500)")
    ap.add_argument("--outdir", default="ff_signed_sets")
    ap.add_argument("--key-col", default=None)
    ap.add_argument("--signed-col", default=None)
    ap.add_argument("--abs-col", default=None)
    args = ap.parse_args()

    topns = [int(x) for x in args.topns.split(",")]
    df = load_signed(args.importance, key_col=args.key_col,
                     signed_col=args.signed_col, abs_col=args.abs_col)
    manifest = select(df, topns, args.outdir)
    print(f"[select] {manifest['n_bins_total']} bins "
          f"({manifest['n_positive']} positive, {manifest['n_negative']} negative)")
    for n, c in manifest["cutoffs"].items():
        print(f"  top{n}: {c['pos_n']} pos, {c['neg_n']} neg")
    print(f"[write] {args.outdir}/  (pos_top*/neg_top* + signed_selection_manifest.json)")


if __name__ == "__main__":
    main()
