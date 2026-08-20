#!/usr/bin/env python
"""
add_neutrophil_dnase.py  --  add a REAL neutrophil accessibility column to the
fetal-fraction openness atlases from user-supplied hg19 DNase/accessibility BEDs.

Fills the documented gap: ENCODE has no hg19 neutrophil DNase-seq, so the
11-tissue atlas carried neutrophil as histone-only and the curated-72 atlas
substituted CMP (CD34+ myeloid progenitor) as a proxy. Three hg19 neutrophil
accessibility peak files (4-col BED: chrom/start/end/signalValue) replace those
proxies with real data.

Aggregation matches the atlas exactly: signalValue-weighted covered-fraction per
50 kb bin, z-scored per replicate, mean across the 3 replicates (zmean).
  - curated-72:  C_neutrophil = signalValue-weighted zmean coverage  (73rd tissue)
  - 11-tissue :  D_neutrophil = same coverage; ACC = z(D);
                 C_neutrophil = mean(z(H_neutrophil), z(ACC_neutrophil))
                 (upgrades the prior histone-only C = z(H) to the full composite)
"""
import argparse
import json
import os
import numpy as np
import pandas as pd

BIN = 50000
NEUT_FILES = [
    "Neutrophil_fetalung.bed.gz",
    "GMPNeutrophil_CD34_bonemarrowPBMC.bed.gz",
    "GMPNeutrophil_bonemarrowPBMC.bed.gz",
]

def peaks_to_coverage(peaks, grid, bin_size, weight="signalValue"):
    idx = {(c, s): i for i, (c, s) in enumerate(zip(grid.chrom.values, grid.start.values))}
    cov = np.zeros(len(grid))
    wsum = np.zeros(len(grid))
    wcnt = np.zeros(len(grid))
    sig = peaks.get("signalValue", pd.Series(np.ones(len(peaks)))).values
    for chrom, s, e, sg in zip(peaks.chrom.values, peaks.start.values, peaks.end.values, sig):
        s = int(s)
        e = int(e)
        b0 = (s // bin_size) * bin_size
        b1 = (e // bin_size) * bin_size
        for bs in range(b0, b1 + bin_size, bin_size):
            i = idx.get((chrom, bs))
            if i is None:
                continue
            ov = min(e, bs + bin_size) - max(s, bs)
            if ov <= 0:
                continue
            cov[i] += ov / bin_size
            wsum[i] += sg * ov
            wcnt[i] += ov
    cov = np.clip(cov, 0, 1)
    if weight == "signalValue":
        mean_sig = np.divide(wsum, wcnt, out=np.zeros_like(wsum), where=wcnt > 0)
        return cov * mean_sig
    return cov

def z(a):
    a = a.astype(float)
    sd = np.nanstd(a)
    return (a - np.nanmean(a)) / (sd if sd else 1.0)

def zmean_reps(mat):
    mat = np.atleast_2d(mat).astype(float)
    if mat.shape[0] == 1:
        return z(mat[0])
    zz = np.vstack([z(r) for r in mat])
    return np.nanmean(zz, axis=0)

def build_neut_coverage(datadir, grid):
    mats = []
    for f in NEUT_FILES:
        p = os.path.join(datadir, f)
        pk = pd.read_csv(p, sep="\t", header=None,
                         names=["chrom", "start", "end", "signalValue"])
        pk = pk[pk.chrom.isin(set(grid.chrom.unique()))]
        mats.append(peaks_to_coverage(pk, grid, BIN, weight="signalValue"))
    return zmean_reps(np.vstack(mats)), len(NEUT_FILES)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--refdir", required=True, help="reference/ dir with atlases + manifests")
    ap.add_argument("--datadir", required=True, help="dir with the 3 neutrophil bed.gz files")
    args = ap.parse_args()

    a11p = os.path.join(args.refdir, "ff_openness_atlas_hg19_50kb.csv.gz")
    a72p = os.path.join(args.refdir, "ff_openness_atlas_hg19_50kb_curated73.csv.gz")
    m11p = os.path.join(args.refdir, "ff_reference_manifest.json")
    m72p = os.path.join(args.refdir, "manifest_curated73.json")

    a11 = pd.read_csv(a11p)
    a72 = pd.read_csv(a72p)
    grid = a72[["chrom", "start", "end"]].copy()
    D_neut, n_rep = build_neut_coverage(args.datadir, grid)
    ACC = z(D_neut)

    # ---- curated-72: add C_neutrophil (73rd tissue) ----
    a72["C_neutrophil"] = D_neut          # signalValue-weighted zmean coverage
    a72.to_csv(a72p, index=False, compression="gzip")

    # ---- 11-tissue: add D_/n_D_/ACC_/upgrade C_neutrophil ----
    a11["D_neutrophil"] = D_neut
    a11["n_D_neutrophil"] = n_rep
    a11["ACC_neutrophil"] = ACC
    if "H_neutrophil" in a11.columns:
        a11["C_neutrophil"] = np.nanmean(np.vstack([z(a11["H_neutrophil"].values), z(ACC)]), axis=0)
    else:
        a11["C_neutrophil"] = ACC
    a11.to_csv(a11p, index=False, compression="gzip")

    # ---- manifests ----
    m11 = json.load(open(m11p))
    m11.setdefault("tissues", {})
    m11["tissues"]["neutrophil"] = [f.replace(".bed.gz", "") for f in NEUT_FILES]
    m11["n_tissues"] = len({k for k in m11["tissues"]})
    m11.setdefault("tissue_notes", {})["neutrophil"] = (
        "granulocyte lineage, often the single largest cfDNA contributor. "
        "REAL hg19 accessibility now included (3 user-supplied neutrophil peak "
        "sets: fetal-lung neutrophil + two GMP/bone-marrow neutrophil), "
        "signalValue-weighted zmean; full D/ACC/H/C composite like every other "
        "tissue. Previously histone-only (no ENCODE neutrophil DNase in hg19).")
    json.dump(m11, open(m11p, "w"), indent=2)

    m72 = json.load(open(m72p))
    if isinstance(m72.get("tissues"), dict):
        m72["tissues"]["neutrophil"] = {"n_rep": n_rep, "class": "primary cell"}
    m72["n_tissues"] = len(m72["tissues"]) if isinstance(m72.get("tissues"), dict) else m72.get("n_tissues", 72) + 1
    m72["neutrophil_note"] = (
        "RESOLVED: real hg19 neutrophil accessibility added as C_neutrophil "
        "(3 user-supplied peak sets: fetal-lung neutrophil + two GMP/bone-marrow "
        "neutrophil, signalValue-weighted zmean). The CD34+ CMP proxy remains in "
        "the panel as the myeloid-progenitor axis; neutrophil is now a distinct "
        "mature-granulocyte column.")
    ca = m72.setdefault("curated_additions", {})
    ca["real_neutrophil_accessibility"] = [f.replace(".bed.gz", "") for f in NEUT_FILES]
    json.dump(m72, open(m72p, "w"), indent=2)

    print(f"D_neutrophil: mean={D_neut.mean():.3f} sd={D_neut.std():.3f} "
          f"nonzero={int((D_neut!=0).sum())}")
    print(f"11-tissue now {a11.shape[1]} cols, n_tissues={m11['n_tissues']}")
    print(f"curated now {a72.shape[1]} cols, n_tissues={m72['n_tissues']}")

if __name__ == "__main__":
    main()
