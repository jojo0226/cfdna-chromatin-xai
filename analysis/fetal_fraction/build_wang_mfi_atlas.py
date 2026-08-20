#!/usr/bin/env python3
"""build_wang_mfi_atlas.py -- maternal-fetal interface snATAC openness layer, 50 kb hg19.

Turns the Wang et al. 2026 maternal-fetal interface (MFI) pseudobulk snATAC atlas
into a per-cell-type open-chromatin track on the SAME 50 kb hg19 autosomal grid +
`key` format as the DNase openness atlas (ff_openness_atlas_hg19_50kb*.csv.gz). It is
a SEPARATE layer (W_<celltype> columns), kept side by side with the DNase C_ columns
and the Wang-2024 placenta A_ columns -- different assays, auditable independently.

Unlike the bulk C_ atlas, every Wang MFI cell type carries a Fetal / Maternal origin
label, so this layer also emits two ORIGIN COMPOSITE tracks (openness_fetal,
openness_maternal) built from clean-origin cell types only (origin_purity >= 0.9,
Unknown excluded). Those drive the per-PC origin score in pc_origin_map.R.

Source (wang_mfi_atlas_hg19.tar.gz, hg19; Wang et al. Nature 653, 2026;
doi:10.1038/s41586-026-10316-x; cell.ucsf.edu/snPlacenta):
    wang_mfi_win5kb_hg19.h5ad    33 cell types x 567,239 5 kb windows (uint32 counts)
    wang_mfi_win50kb_hg19.h5ad   33 cell types x  56,734 50 kb windows (10:1 sum of 5 kb)
    celltype_metadata.csv        origin, origin_purity, total_fragments, QC flags

GRID PROBLEM (measured): the Wang windows are hg19 only on their START edge (lifted
from hg38), so their starts are phase-OFFSET from the clean 50 kb FF grid. Exact
(chrom,start) key match recovers 0.7% of bins; re-binning by window MIDPOINT recovers
92.8%. We therefore aggregate from the 5 kb matrix (the 50 kb matrix is an exact 10:1
sum of it, so nothing is lost) onto the FF grid by midpoint:
    FF_bin(window) = floor((start + 2500) / 50000) * 50000     # 5 kb window midpoint
    raw[b,t] = sum over 5 kb windows w mapped to bin b of counts[w, t]

Normalization (raw counts -> comparable openness), per cell type t:
    1. CPM  = raw / total_fragments(t) * 1e6      (depth spans 1000x; must normalize)
    2. log1p(CPM)                                  (variance stabilize)
    3. z-score across bins (per cell type)         (same footing as C_ / A_ columns,
                                                    so Pearson corr & PC alignment match)

Origin composites are the per-bin MEAN of the z-scored clean-origin columns (fetal /
maternal), then themselves left as-is (already comparable). Ambiguous types
(origin_purity < 0.75) and Unknown are excluded from composites but KEPT as W_ columns.

Usage:
  python build_wang_mfi_atlas.py \
      --tar  ~/Downloads/wang_mfi_atlas_hg19.tar.gz \
      --grid analysis/fetal_fraction/reference/ff_openness_atlas_hg19_50kb_curated73.csv.gz \
      --out  analysis/fetal_fraction/reference/wang_mfi_openness_hg19_50kb.csv.gz

Outputs (next to --out):
  wang_mfi_openness_hg19_50kb.csv.gz    chrom,start,end,key, W_<ct> x33, openness_fetal, openness_maternal
  manifest_wang_mfi.json                provenance, ct->origin/purity/QC map, grid stats, recipe
"""
from __future__ import annotations
import argparse, json, os, sys, tarfile, tempfile
import numpy as np
import pandas as pd

BIN = 50000
CLEAN_PURITY = 0.90     # origin composites use only cell types this pure
AMBIG_PURITY = 0.75     # below this = origin-ambiguous, flagged
MEMBER = "wang_mfi_atlas_hg19"


def _extract(tar_path, member, dstdir):
    with tarfile.open(tar_path) as t:
        t.extract(f"{MEMBER}/{member}", path=dstdir)
    return os.path.join(dstdir, MEMBER, member)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tar", required=True, help="wang_mfi_atlas_hg19.tar.gz")
    ap.add_argument("--grid", required=True, help="FF openness atlas csv.gz (defines the 50 kb grid + key format)")
    ap.add_argument("--out", required=True, help="output csv.gz for the W_ layer")
    ap.add_argument("--res", choices=["5kb", "50kb"], default="5kb",
                    help="which h5ad to aggregate from (default 5kb -- phase-correct)")
    args = ap.parse_args()

    import anndata as ad

    # ---- FF grid (target) ----------------------------------------------------
    grid = pd.read_csv(args.grid, usecols=["chrom", "start", "end", "key"])
    grid["start"] = grid["start"].astype(int)
    assert (grid["start"] % BIN == 0).all(), "FF grid not on clean 50 kb phase"
    grid_index = pd.MultiIndex.from_arrays([grid["chrom"].values, grid["start"].values])
    grid_chroms = set(grid["chrom"].unique())          # autosomes only
    print(f"[grid] {len(grid)} bins, {len(grid_chroms)} chroms, key e.g. {grid['key'].iloc[0]}")

    # ---- Wang h5ad + metadata ------------------------------------------------
    with tempfile.TemporaryDirectory() as td:
        h5name = f"wang_mfi_win{args.res}_hg19.h5ad"
        h5 = _extract(args.tar, h5name, td)
        meta_p = _extract(args.tar, "celltype_metadata.csv", td)
        A = ad.read_h5ad(h5)
        meta = pd.read_csv(meta_p)

    win_res = 5000 if args.res == "5kb" else BIN
    cts = list(A.obs_names)
    X = np.asarray(A.X, dtype=np.float64)              # cell types x windows
    var = A.var.copy(); var["start"] = var["start"].astype(int)
    print(f"[wang] {A.shape[0]} cell types x {A.shape[1]} {args.res} windows")

    # ---- midpoint -> FF grid start ------------------------------------------
    mid = var["start"].values + win_res // 2
    gstart = (mid // BIN) * BIN
    wchrom = var["chrom"].values
    keep = np.array([c in grid_chroms for c in wchrom])
    print(f"[map] {keep.sum()}/{len(keep)} windows on autosomal FF chroms")

    # aggregate counts per (chrom, gridstart) via pandas groupby on the window axis
    # build long frame: one row per window, columns = counts per cell type -> sum
    wf = pd.DataFrame(X.T[keep], columns=cts)          # windows x cell types
    wf["chrom"] = wchrom[keep]; wf["start"] = gstart[keep]
    agg = wf.groupby(["chrom", "start"], observed=True)[cts].sum()
    print(f"[agg] {agg.shape[0]} unique grid bins carry Wang signal")

    # reindex onto the FF grid (bins with no Wang signal -> NaN, not 0)
    agg = agg.reindex(grid_index)
    covered = agg[cts[0]].notna().sum()  # any ct col; all-or-nothing per bin
    n_any = agg.notna().any(axis=1).sum()
    print(f"[cover] {n_any}/{len(grid)} FF bins mapped ({n_any/len(grid)*100:.1f}%)")

    raw = agg.values                                   # bins x cell types (NaN where unmapped)

    # ---- normalize: CPM -> log1p -> z-score per cell type --------------------
    tot = meta.set_index("cell_type")["total_fragments"].reindex(cts).values.astype(float)
    cpm = raw / tot[None, :] * 1e6
    logv = np.log1p(cpm)
    mu = np.nanmean(logv, axis=0); sd = np.nanstd(logv, axis=0) + 1e-12
    z = (logv - mu[None, :]) / sd[None, :]             # bins x cell types, NaN preserved

    out = grid.copy()
    for j, ct in enumerate(cts):
        out[f"W_{ct}"] = z[:, j]

    # ---- origin composites (clean types only) --------------------------------
    # Built from row-centered SPECIFICITY, not raw z-score: at 50 kb every cell type
    # shares a large "open-everywhere" component (README dilution caveat), so raw
    # z-score composites correlate ~+0.91 and cannot separate origin. Row-centering
    # (spec = z - per-bin mean across cell types) removes the shared component and
    # yields a genuine fetal-vs-maternal axis (composite corr ~ -0.83; trophoblast
    # +0.6..0.76 fetal / -0.5..-0.72 maternal). The W_ columns stay z-scored to match
    # the C_/A_ layers (pc_origin_map.R re-centers at analysis time, idempotently).
    m = meta.set_index("cell_type")
    fetal_cts = [c for c in cts if m.loc[c, "origin"] == "Fetal" and m.loc[c, "origin_purity"] >= CLEAN_PURITY]
    mat_cts   = [c for c in cts if m.loc[c, "origin"] == "Maternal" and m.loc[c, "origin_purity"] >= CLEAN_PURITY]
    spec = z - np.nanmean(z, axis=1, keepdims=True)     # bins x cell types, row-centered
    out["openness_fetal"] = np.nanmean(spec[:, [cts.index(c) for c in fetal_cts]], axis=1)
    out["openness_maternal"] = np.nanmean(spec[:, [cts.index(c) for c in mat_cts]], axis=1)
    print(f"[origin] fetal composite: {len(fetal_cts)} cts; maternal: {len(mat_cts)} cts (from specificity)")

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    out.to_csv(args.out, index=False, compression="gzip")
    print(f"[write] {args.out}  ({out.shape[0]} bins x {out.shape[1]} cols)")

    # ---- manifest ------------------------------------------------------------
    ct_map = {}
    for c in cts:
        r = m.loc[c]
        ct_map[c] = dict(origin=str(r["origin"]), origin_purity=float(r["origin_purity"]),
                         n_cells=int(r["n_cells_total"]), total_fragments=int(r["total_fragments"]),
                         ambiguous=bool(r["origin_purity"] < AMBIG_PURITY),
                         low_n=bool(r["low_n_flag"]))
    manifest = dict(
        atlas="wang_mfi_openness_hg19_50kb",
        source="Wang et al. Nature 653:167-179 (2026), doi:10.1038/s41586-026-10316-x; cell.ucsf.edu/snPlacenta",
        assembly="hg19", bin_size_bp=BIN, prefix="W_",
        aggregated_from=f"wang_mfi_win{args.res}_hg19.h5ad",
        n_bins=int(len(grid)), n_bins_mapped=int(n_any),
        grid_coverage_frac=round(float(n_any / len(grid)), 4),
        grid_join="window midpoint -> floor(mid/50000)*50000 (starts are hg38-lifted, phase-offset)",
        normalization="CPM (per cell-type total_fragments) -> log1p -> z-score across bins per cell type",
        n_celltypes=len(cts),
        origin_composites=dict(clean_purity_threshold=CLEAN_PURITY,
                               built_from="row-centered specificity (z - per-bin mean across cell types)",
                               fetal_celltypes=fetal_cts, maternal_celltypes=mat_cts),
        ambiguous_celltypes=[c for c in cts if m.loc[c, "origin_purity"] < AMBIG_PURITY],
        low_n_celltypes=[c for c in cts if bool(m.loc[c, "low_n_flag"])],
        celltypes=ct_map,
    )
    mpath = os.path.join(os.path.dirname(args.out), "manifest_wang_mfi.json")
    with open(mpath, "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"[write] {mpath}")


if __name__ == "__main__":
    main()
