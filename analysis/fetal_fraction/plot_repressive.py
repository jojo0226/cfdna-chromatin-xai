#!/usr/bin/env python
"""Two SEPARATE repressive-axis figures for the FF cell-of-origin track.

The repressive atlas (build_repressive_atlas.py -> H3K27me3 z-mean coverage
per tissue, plus pan-subtracted repr_spec_<tissue>) is a fully SEPARATE axis
from the DNase openness atlas. It is never merged into the C_ openness
columns. These two figures characterise that axis on its own and test the
two-sided cell-of-origin premise against the openness atlas.

FIG 1  repr_fig1_clustering.png
    Tissue x tissue Spearman correlation of the genome-wide H3K27me3
    *specificity* profile (repr_spec_<tissue>), hierarchically ordered.
    A genuine repressive axis clusters tissues by lineage (hematopoietic
    block, placenta family, etc.).

FIG 2  repr_fig2_twosided.png
    Per-tissue Spearman rho between active-specificity (C_ pan-subtracted,
    DNase openness atlas) and repressive-specificity (repr_spec_).  rho < 0
    means: where a tissue is *specifically* open it is *specifically*
    depleted of its own H3K27me3 -- the signature a true cell-of-origin
    signal must show, and the justification for a two-sided enrichment test
    (active-enriched AND repressive-depleted).

NOTE on why *specificity* and not raw z: at 50 kb, raw DNase and raw
H3K27me3 are both concentrated in gene-rich euchromatin, so raw-vs-raw is
positively correlated (a gene-density baseline, ~+0.2).  Pan-tissue
subtraction removes that shared baseline; the residual specificity is the
biologically interpretable quantity.

Usage:
    python plot_repressive.py \
        --repressive ff_repressive_atlas_hg19_50kb_H3K27me3.csv.gz \
        --openness   ff_openness_atlas_hg19_50kb_curated73.csv.gz \
        --outdir     figures/
"""
import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform
from scipy.stats import rankdata

CFDNA = {"placenta", "trophoblast", "chorion", "monocyte", "Bcell",
         "CD4", "CD8", "NK", "liver", "CMP", "HUVEC"}


def _clean(name):
    """ENCODE free-text biosample -> compact token (mirror of builder)."""
    s = name.lower()
    table = {
        "cd14-positive monocyte": "monocyte", "b cell": "Bcell",
        "cd4-positive, alpha-beta t cell": "CD4",
        "cd8-positive, alpha-beta t cell": "CD8",
        "natural killer cell": "NK", "common myeloid progenitor, cd34-positive": "CMP",
        "trophoblast cell": "trophoblast",
        "endothelial cell of umbilical vein": "HUVEC",
        "peyer's patch": "Peyers_patch",
    }
    if s in table:
        return table[s]
    return (s.replace("'", "").replace(",", "").replace("-", "_")
            .replace(" ", "_").replace("/", "_"))


def spearman_matrix(mat):
    ranked = np.apply_along_axis(rankdata, 0, mat)
    return np.corrcoef(ranked, rowvar=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repressive", required=True)
    ap.add_argument("--openness", required=True)
    ap.add_argument("--outdir", default="figures")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    rep = pd.read_csv(args.repressive)
    op = pd.read_csv(args.openness)
    assert (rep["key"].values == op["key"].values).all(), "grid mismatch"

    spec_cols = [c for c in rep.columns if c.startswith("repr_spec_")]
    toks = [c[len("repr_spec_"):] for c in spec_cols]

    # ---- openness specificity (pan-tissue subtracted) ----
    Cc = [c for c in op.columns if c.startswith("C_")]
    op_terms = [c[2:] for c in Cc]
    op_clean = {_clean(t): t for t in op_terms}
    Cmat = op[Cc].to_numpy()
    op_pan = np.nanmean(Cmat, axis=1, keepdims=True)
    op_spec = {op_terms[j]: Cmat[:, j] - op_pan[:, 0] for j in range(len(op_terms))}

    # ================= FIGURE 1 =================
    S = rep[spec_cols].to_numpy()
    A = spearman_matrix(S)
    D = 1.0 - A
    np.fill_diagonal(D, 0.0)
    order = leaves_list(linkage(squareform(D, checks=False), method="average"))
    Ao = A[np.ix_(order, order)]
    labs = [toks[i].replace("_", " ") for i in order]
    critset = {t.replace("_", " ") for t in CFDNA}

    fig1, ax = plt.subplots(figsize=(8.2, 7.6))
    im = ax.imshow(Ao, cmap="RdBu_r", vmin=-1, vmax=1, aspect="equal")
    ax.set_xticks(range(len(labs)))
    ax.set_xticklabels(labs, rotation=90, fontsize=5.2)
    ax.set_yticks(range(len(labs)))
    ax.set_yticklabels(labs, fontsize=5.2)
    for i, lb in enumerate(labs):
        if lb in critset:
            for tk in (ax.get_xticklabels()[i], ax.get_yticklabels()[i]):
                tk.set_color("#b8272e")
                tk.set_fontweight("bold")
    cb = fig1.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label("Spearman \u03c1  (H3K27me3 specificity profile)", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    ax.set_title("Repressive (H3K27me3) axis clusters tissues by lineage\n"
                 "(cfDNA-relevant tissues in red)", fontsize=8.5, loc="left")
    p1 = os.path.join(args.outdir, "repr_fig1_clustering.png")
    fig1.savefig(p1, dpi=200, bbox_inches="tight")

    # ================= FIGURE 2 =================
    rows = []
    for tok in toks:
        raw = op_clean.get(tok)
        if raw is None:
            continue
        a = op_spec[raw]
        r = rep[f"repr_spec_{tok}"].to_numpy()
        m = np.isfinite(a) & np.isfinite(r)
        rho = np.corrcoef(rankdata(a[m]), rankdata(r[m]))[0, 1]
        rows.append((tok, rho))
    pB = (pd.DataFrame(rows, columns=["tissue", "rho"])
          .sort_values("rho").reset_index(drop=True))
    y = np.arange(len(pB))
    colors = ["#b8272e" if t in CFDNA else "#9bb0c1" for t in pB.tissue]

    fig2, ax = plt.subplots(figsize=(6.4, 8.0))
    ax.barh(y, pB.rho, color=colors, height=0.72)
    ax.axvline(0, color=".3", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([t.replace("_", " ") for t in pB.tissue], fontsize=5.6)
    for i, t in enumerate(pB.tissue):
        if t in CFDNA:
            ax.get_yticklabels()[i].set_color("#b8272e")
            ax.get_yticklabels()[i].set_fontweight("bold")
    ax.set_xlabel("Spearman \u03c1:  active-specificity  vs  repressive-specificity",
                  fontsize=7.5)
    ax.set_ylim(-0.7, len(pB) - 0.3)
    ax.set_title("Where a tissue is specifically OPEN, its H3K27me3 is "
                 "specifically DEPLETED\n"
                 "\u03c1<0 (left) = the two-sided cell-of-origin test is valid",
                 fontsize=8.5, loc="left")
    ax.annotate(f"mean \u03c1 = {pB.rho.mean():+.3f}\n"
                f"{int((pB.rho < 0).mean() * 100)}% of tissues \u03c1<0",
                xy=(0.03, 0.06), xycoords="axes fraction", fontsize=6.5,
                ha="left", va="bottom", color=".25")
    ax.legend(handles=[Patch(color="#b8272e", label="cfDNA-relevant"),
                       Patch(color="#9bb0c1", label="other tissue")],
              fontsize=6.5, frameon=False, loc="lower right")
    p2 = os.path.join(args.outdir, "repr_fig2_twosided.png")
    fig2.savefig(p2, dpi=200, bbox_inches="tight")

    pB.to_csv(os.path.join(args.outdir, "repr_twosided_rho.csv"), index=False)
    print(f"[fig1] {p1}  ({A.shape[0]} tissues)")
    print(f"[fig2] {p2}  (mean rho {pB.rho.mean():+.3f}, "
          f"{int((pB.rho < 0).mean() * 100)}% negative)")


if __name__ == "__main__":
    main()
