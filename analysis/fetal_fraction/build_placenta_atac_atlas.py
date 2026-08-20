#!/usr/bin/env python3
"""build_placenta_atac_atlas.py -- placental snATAC openness layer, 50 kb hg19.

Turns the Wang et al. 2024 pregnancy snATAC pseudobulk atlas into a per-cell-type
open-chromatin track on the SAME 50 kb hg19 autosomal grid + `key` format as the
DNase openness atlas (ff_openness_atlas_hg19_50kb*.csv.gz). It is a SEPARATE layer
(A_<celltype> columns), not merged into the C_<tissue> DNase columns -- the two are
different assays and are kept side by side so cross-modal agreement stays auditable.

Source h5 (produced upstream, hg38->hg19 lifted; see LIFTOVER note in the atlas
README): datasets counts(peaks x profiles), peaks/{chr,start,end},
profiles/{id,cell_type,stage,gw,n_cells,n_donors}. Values are per-profile CPM-like
counts over 283,406 peaks; 11 pseudobulk profiles = 6 cell types x {early GW8, term
GW38-39} (EVT early only).

Per bin b, per cell type t:
    raw[b,t] = sum over peaks p overlapping bin b of  CPM[p, t]
where CPM[p,t] pools that cell type's stage profiles weighted by n_cells (a term and
early profile of the same cell type are averaged by their nucleus counts, so a rare
early population does not dominate). A peak is assigned to every 50 kb bin it overlaps
(peaks are short << 50 kb, so this is almost always one bin). Columns are then
z-scored across bins (per cell type), matching the openness-atlas convention so that
Pearson correlations and PC alignments are on the same footing as the C_ columns.

Usage:
  python build_placenta_atac_atlas.py \
      --h5 ~/Downloads/scATACseq/Wang2024_pregnancy_snATAC_atlas_hg19.h5 \
      --grid analysis/fetal_fraction/reference/ff_openness_atlas_hg19_50kb_curated73.csv.gz \
      --out  analysis/fetal_fraction/reference/ff_placenta_atac_hg19_50kb.csv.gz

--grid supplies the canonical chrom/start/end/key bins (any openness atlas works;
its C_ columns are ignored). If omitted, the 50 kb autosomal grid is rebuilt from
hardcoded UCSC hg19 chromosome lengths (identical bins).
"""
import argparse, json, os, sys
import numpy as np
import pandas as pd
import h5py

# UCSC hg19 autosome lengths (fallback grid; identical to the openness atlas)
HG19 = {
 "chr1":249250621,"chr2":243199373,"chr3":198022430,"chr4":191154276,
 "chr5":180915260,"chr6":171115067,"chr7":159138663,"chr8":146364022,
 "chr9":141213431,"chr10":135534747,"chr11":135006516,"chr12":133851895,
 "chr13":115169878,"chr14":107349540,"chr15":102531392,"chr16":90354753,
 "chr17":81195210,"chr18":78077248,"chr19":59128983,"chr20":63025520,
 "chr21":48129895,"chr22":51304566}
BIN = 50000

def build_grid_from_lengths():
    rows=[]
    for c in [f"chr{i}" for i in range(1,23)]:
        n=HG19[c]
        for s in range(0,n,BIN):
            e=min(s+BIN,n)
            rows.append((c,s,e,f"{c}:{s}-{e}"))
    return pd.DataFrame(rows,columns=["chrom","start","end","key"])

def load_grid(path):
    if path and os.path.exists(os.path.expanduser(path)):
        g=pd.read_csv(os.path.expanduser(path),usecols=["chrom","start","end","key"])
        return g.reset_index(drop=True)
    return build_grid_from_lengths()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--h5",required=True)
    ap.add_argument("--grid",default="analysis/fetal_fraction/reference/ff_openness_atlas_hg19_50kb_curated73.csv.gz")
    ap.add_argument("--out",default="analysis/fetal_fraction/reference/ff_placenta_atac_hg19_50kb.csv.gz")
    ap.add_argument("--manifest",default=None)
    args=ap.parse_args()

    grid=load_grid(args.grid)
    autoso=set(f"chr{i}" for i in range(1,23))
    grid=grid[grid.chrom.isin(autoso)].reset_index(drop=True)
    # bin index lookup: (chrom,start)->row
    grid["bin_idx"]=(grid.start//BIN).astype(int)
    chrom_offset={}; off=0; order=[f"chr{i}" for i in range(1,23)]
    for c in order:
        chrom_offset[c]=off; off+=int(np.ceil(HG19[c]/BIN))
    def gidx(chrom,start):
        return chrom_offset.get(chrom,-1)+(start//BIN)
    grid["gidx"]=[gidx(c,s) for c,s in zip(grid.chrom,grid.start)]
    total_bins=off

    with h5py.File(os.path.expanduser(args.h5),"r") as f:
        counts=f["counts"][:].astype(np.float64)              # peaks x profiles
        pchr=[x.decode() for x in f["peaks/chr"][:]]
        pstart=f["peaks/start"][:].astype(int)
        pend=f["peaks/end"][:].astype(int)
        prof_ct=[x.decode() for x in f["profiles/cell_type"][:]]
        prof_ncell=f["profiles/n_cells"][:].astype(float)
    n_peaks,n_prof=counts.shape

    # pool profiles -> cell type, n_cells-weighted
    cell_types=sorted(set(prof_ct))
    ct_cpm=np.zeros((n_peaks,len(cell_types)))
    for j,ct in enumerate(cell_types):
        cols=[i for i,c in enumerate(prof_ct) if c==ct]
        w=prof_ncell[cols]; w=w/w.sum()
        ct_cpm[:,j]=counts[:,cols]@w

    # peak -> bin(s); accumulate
    peak_g=np.array([gidx(c,s) for c,s in zip(pchr,pstart)])
    peak_g_end=np.array([gidx(c,max(s,e-1)) for c,s,e in zip(pchr,pstart,pend)])
    raw=np.zeros((total_bins,len(cell_types)))
    for p in range(n_peaks):
        g0=peak_g[p]
        if g0<0: continue
        g1=peak_g_end[p]
        for g in range(g0,g1+1):
            raw[g,:]+=ct_cpm[p,:]

    # map onto the grid rows, z-score per cell type across grid bins
    out=grid[["chrom","start","end","key"]].copy()
    sub=raw[grid.gidx.values,:]
    for j,ct in enumerate(cell_types):
        v=sub[:,j].astype(float)
        mu,sd=v.mean(),v.std()
        z=(v-mu)/sd if sd>0 else v*0.0
        out[f"A_{ct}"]=z

    # placenta_atac_spec: max A_ minus mean of the rest (a placenta-specificity score
    # analogous to the *_specific contrasts in the DNase atlas)
    A=out[[c for c in out.columns if c.startswith("A_")]].values
    out["placenta_atac_spec"]=A.max(1)-(A.sum(1)-A.max(1))/(A.shape[1]-1)

    os.makedirs(os.path.dirname(os.path.expanduser(args.out)) or ".",exist_ok=True)
    out.to_csv(os.path.expanduser(args.out),index=False,compression="gzip")
    n_touch=int((sub.sum(1)>0).sum())
    print(f"[atac] {len(out)} bins x {len(cell_types)} cell types; "
          f"{n_touch} bins with >=1 ATAC peak ({100*n_touch/len(out):.1f}%)")
    print(f"[write] {args.out}")

    manifest={
      "atlas":"placental snATAC openness (Wang et al. 2024)",
      "genome":"hg19 (lifted from hg38; UCSC hg38ToHg19)",
      "bin_size":BIN,"n_bins":len(out),"cell_types":cell_types,
      "columns":["chrom","start","end","key"]+[f"A_{c}" for c in cell_types]+["placenta_atac_spec"],
      "value":"per-cell-type CPM summed per 50kb bin (stages pooled by n_cells), z-scored across bins",
      "source_h5":os.path.basename(args.h5),
      "study":"Wang M et al., Nat Genet 56:294-305 (2024), DOI 10.1038/s41588-023-01647-w, GSE247036",
      "n_peaks_source":int(n_peaks),
      "note":"SEPARATE layer -- do NOT merge into C_ DNase columns; use for cross-modal validation and PC alignment (pc_tissue_map.R --prefix A_)."
    }
    mpath=args.manifest or os.path.join(os.path.dirname(os.path.expanduser(args.out)),"manifest_placenta_atac.json")
    with open(mpath,"w") as fh: json.dump(manifest,fh,indent=2)
    print(f"[write] {mpath}")

if __name__=="__main__":
    main()
