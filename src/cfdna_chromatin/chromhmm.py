"""
chromhmm.py -- user-facing ChromHMM state-annotation layer for importance bins.

This is a thin façade over :mod:`references` (18-state vocabulary + segmentation
loading) and :mod:`engine` (per-region annotation + matched-null enrichment). It
gives one entry point that takes a set of importance-ranked genomic bins and
returns, per reference epigenome, (a) the ChromHMM dominant-state composition of
those bins and (b) -- optionally, when an hg38 genome is supplied -- the
GC/length/N-matched state enrichment with empirical p-values.

------------------------------------------------------------------------------
GENOME BUILD -- read this if your cohort is hg19 (e.g. PRAD_14230_hg19)
------------------------------------------------------------------------------
The bundled ChromHMM reference segmentations are **hg38 only** (ENCODE
Roadmap/EpiMap 18-state tracks: placenta, neutrophil, Bcell, K562,
keratinocyte). ChromHMM annotations are just labelled genome intervals, so the
*engine* is build-agnostic -- the only build-specific thing is which reference
BEDs exist. Two supported ways to annotate an hg19 query set:

  1. LIFT THE QUERY (per-analysis, default here). Pass ``query_build="hg19"`` and
     the top bins are lifted hg19->hg38 with round-trip QC (liftover.lift_regions)
     before annotation. Cheap because only the handful of ranked bins move; the
     importance/direction columns ride along. A ``liftover_qc`` dict is returned
     so you can check the mapping rate.

  2. PRE-LIFT THE REFERENCE (recommended for a whole hg19 cohort). Lift the five
     static hg38 segmentation BEDs DOWN to hg19 once with
     ``liftover.lift_reference_bed(..., from_build="hg38", to_build="hg19")``,
     build an hg19 reference bundle, then annotate natively with
     ``query_build="hg38"`` semantics against that hg19 bundle (no per-analysis
     lifting). See ``scripts/prelift_chromhmm_hg19.py``.

If you would rather annotate against **native hg19 ChromHMM tracks**, Roadmap
Epigenomics distributes 15/18-state hg19 segmentations per epigenome; drop them
in a data dir with the same ``tissue_ACCESSION.bed.gz`` naming and point
``data_dir`` at it (they share the same 18-state vocabulary as loaded here).

All engine coordinates are 0-based half-open.
"""
from __future__ import annotations
import os

import numpy as np
import pandas as pd

from . import references as R
from . import engine as E
from . import liftover as L
from .selection import parse_bin


# --------------------------------------------------------------------------- #
#  Reference loading                                                           #
# --------------------------------------------------------------------------- #
def reference_bed_path(data_dir, tissue):
    """Resolve the ChromHMM segmentation BED for one tissue in data_dir.

    Accepts either ``<tissue>_<ACCESSION>.bed.gz`` (bundled naming) or a bare
    ``<tissue>.bed.gz``. Raises FileNotFoundError with a helpful message.
    """
    acc = R.REFERENCE_ACCESSIONS.get(tissue)
    cands = []
    if acc:
        cands.append(os.path.join(data_dir, f"{tissue}_{acc}.bed.gz"))
    cands.append(os.path.join(data_dir, f"{tissue}.bed.gz"))
    for p in cands:
        if os.path.exists(p):
            return p
    raise FileNotFoundError(
        f"No ChromHMM segmentation for tissue {tissue!r} in {data_dir!r}; "
        f"looked for {', '.join(os.path.basename(c) for c in cands)}"
    )


def load_reference_panel(data_dir, tissues, chroms=None):
    """Load ChromHMM segmentations for a list of tissues.

    Returns {tissue: segmentation} where segmentation is the per-chrom sorted
    arrays produced by references.load_segmentation.
    """
    return {t: R.load_segmentation(reference_bed_path(data_dir, t), chroms=chroms)
            for t in tissues}


def panel_tissues(panel_name):
    """Convenience re-export: ordered tissue list for a named panel."""
    return R.panel_tissues(panel_name)


# --------------------------------------------------------------------------- #
#  Query preparation (build crossing)                                         #
# --------------------------------------------------------------------------- #
def _regions_frame(bins):
    """bins: iterable of bin names OR (chrom,start,end) OR DataFrame -> DataFrame.

    Output DataFrame has columns chrom/start/end and preserves any importance/
    direction columns when a DataFrame is passed in.
    """
    if isinstance(bins, pd.DataFrame):
        df = bins.copy()
        if not {"chrom", "start", "end"}.issubset(df.columns):
            # derive coords from the index (bin names) or a 'key'/'bin' column
            src = df.index if df.index.name or df.index.dtype == object else None
            keycol = next((c for c in ("key", "bin", "feature", "name") if c in df.columns), None)
            names = df[keycol] if keycol else src
            coords = [parse_bin(n) for n in names]
            cf = pd.DataFrame(coords, columns=["chrom", "start", "end"], index=df.index)
            df = pd.concat([cf, df], axis=1)
        return df
    rows, idx = [], []
    for i, b in enumerate(bins):
        if isinstance(b, (tuple, list)) and len(b) >= 3:
            rows.append((b[0], int(b[1]), int(b[2]))); idx.append(i)
        else:
            rows.append(parse_bin(b)); idx.append(b)
    return pd.DataFrame(rows, columns=["chrom", "start", "end"], index=idx)


def prepare_query(bins, query_build="hg38", ref_build="hg38",
                  chain_path=None, back_chain_path=None, round_trip=True):
    """Return (regions_df, liftover_qc) ready to annotate against ref_build.

    If query_build == ref_build, returns the frame unchanged and qc=None. Else
    lifts chrom/start/end query_build->ref_build with QC, preserving all extra
    columns (importance, direction, ...).
    """
    df = _regions_frame(bins)
    if query_build == ref_build:
        return df, None
    lifted, qc = L.lift_regions(df, from_build=query_build, to_build=ref_build,
                                chain_path=chain_path, back_chain_path=back_chain_path,
                                round_trip=round_trip)
    return lifted, qc


# --------------------------------------------------------------------------- #
#  Annotation (no genome needed)                                              #
# --------------------------------------------------------------------------- #
def annotate(regions_df, segmentation):
    """Per-region ChromHMM annotation in one reference tissue.

    regions_df : DataFrame with chrom/start/end.
    Returns regions_df with added columns dominant_state, active_fraction,
    covered_bp (one row per input region).
    """
    recs = list(regions_df[["chrom", "start", "end"]].itertuples(index=False, name=None))
    ann = E.annotate_regions(recs, segmentation)
    out = regions_df.copy()
    out["dominant_state"] = [a["dominant_state"] for a in ann]
    out["active_fraction"] = [a["active_fraction"] for a in ann]
    out["covered_bp"] = [a["covered_bp"] for a in ann]
    return out


def state_composition(regions_df, segmentation, by="state"):
    """Fraction of query bins whose dominant label is each state (or group).

    Returns a DataFrame indexed by state/group with a 'query_frac' column and a
    'group' column (when by='state'), sorted by query_frac descending.
    """
    recs = list(regions_df[["chrom", "start", "end"]].itertuples(index=False, name=None))
    _, frac, n = E._dominant_counts(recs, segmentation, by=by)
    df = pd.DataFrame({"label": list(frac.keys()),
                       "query_frac": list(frac.values())})
    if by == "state":
        df["group"] = df["label"].map(R.GROUP)
    df["n_bins"] = n
    return df.sort_values("query_frac", ascending=False).reset_index(drop=True)


# --------------------------------------------------------------------------- #
#  Enrichment (needs an hg38 genome for the matched null)                      #
# --------------------------------------------------------------------------- #
def enrichment(regions_df, segmentation, null_model, by="state",
               n_per=100, n_boot=1000, seed=0):
    """State/group enrichment vs. a GC/length/N-matched null (engine.enrichment_test).

    Returns (results_df, meta). results_df columns: label, query_frac, null_frac,
    log2_fold, ci_low, ci_high, p_emp, n_query, n_null (+ group when by='state').
    """
    recs = list(regions_df[["chrom", "start", "end"]].itertuples(index=False, name=None))
    res, meta = E.enrichment_test(recs, segmentation, null_model, by=by,
                                  n_per=n_per, n_boot=n_boot, seed=seed)
    df = pd.DataFrame(res)
    if by == "state" and len(df):
        df["group"] = df["label"].map(R.GROUP)
    return df.sort_values("log2_fold", ascending=False).reset_index(drop=True), meta


# --------------------------------------------------------------------------- #
#  One-call driver                                                             #
# --------------------------------------------------------------------------- #
def annotate_importance_bins(bins, data_dir, panel="fetal", tissues=None,
                             query_build="hg38", ref_build="hg38",
                             chroms=None, chain_path=None, back_chain_path=None,
                             round_trip=None, by="state",
                             genome=None, null_kwargs=None):
    """End-to-end ChromHMM readout for a set of importance bins.

    bins        : bin names / (chrom,start,end) / DataFrame (importance columns
                  preserved).
    data_dir    : dir holding tissue_ACCESSION.bed.gz ChromHMM segmentations.
    panel       : reference panel name (default 'fetal'); ignored if `tissues` given.
    tissues     : explicit tissue list (overrides panel).
    query_build : build of `bins` ('hg19' triggers a lift to ref_build).
    genome      : optional hg38 genome dict (genome.load_genome). When supplied,
                  matched-null enrichment is computed in addition to composition.
    Returns dict:
        regions        : annotated query frame (post-lift if lifted)
        liftover_qc    : QC dict or None
        composition    : {tissue: state_composition DataFrame}
        enrichment     : {tissue: enrichment DataFrame}  (only if genome given)
        tissues, panel, by
    """
    tissues = tissues or R.panel_tissues(panel)
    seg_panel = load_reference_panel(data_dir, tissues, chroms=chroms)
    # Round-trip QC needs a reverse chain; default it on only when one is available.
    if round_trip is None:
        round_trip = back_chain_path is not None
    regions, qc = prepare_query(bins, query_build=query_build, ref_build=ref_build,
                                chain_path=chain_path, back_chain_path=back_chain_path,
                                round_trip=round_trip)

    comp = {t: state_composition(regions, seg_panel[t], by=by) for t in tissues}
    out = {"regions": regions, "liftover_qc": qc, "composition": comp,
           "tissues": tissues, "panel": panel, "by": by, "enrichment": {}}

    if genome is not None:
        nk = dict(chroms=chroms or list(genome.keys()))
        nk.update(null_kwargs or {})
        null = E.NullModel(genome=genome, **nk)
        out["enrichment"] = {t: enrichment(regions, seg_panel[t], null, by=by)[0]
                             for t in tissues}
    return out


# --------------------------------------------------------------------------- #
#  Composition matrix + multi-cutoff heatmap                                   #
# --------------------------------------------------------------------------- #
def composition_matrix(comp_dict, by="state"):
    """Stack per-tissue state_composition frames into a states x tissues matrix.

    comp_dict : {tissue: DataFrame returned by state_composition}
    by        : 'state' (18 rows, canonical STATES order) or 'group' (coarse).
    Returns a DataFrame indexed by state/group (canonical order, all-zero rows
    dropped), columns = tissues, values = query_frac (0.0 where a tissue has no
    bins in that state).
    """
    if by == "state":
        order = list(R.STATES)
    else:
        order = list(dict.fromkeys(R.GROUP.values()))
    tissues = list(comp_dict.keys())
    mat = pd.DataFrame(0.0, index=order, columns=tissues, dtype=float)
    for t, df in comp_dict.items():
        s = df.set_index("label")["query_frac"]
        for lab, v in s.items():
            if lab in mat.index:
                mat.at[lab, t] = float(v)
    return mat.loc[mat.sum(axis=1) > 0]


def multi_cutoff_composition(regions, seg_panel, cutoffs, by="state"):
    """Composition at several top-N cutoffs from one importance-ranked frame.

    regions   : DataFrame already sorted by importance (row 0 = most important),
                as returned by prepare_query on a ranked input.
    seg_panel : {tissue: segmentation DataFrame} (load_reference_panel output).
    cutoffs   : iterable of ints, e.g. (1000, 2000, 5000). Cutoffs larger than
                len(regions) are clamped to len(regions).
    Returns {cutoff: {tissue: state_composition DataFrame}}.
    """
    tissues = list(seg_panel.keys())
    out = {}
    for k in cutoffs:
        sub = regions.head(int(k))
        out[int(k)] = {t: state_composition(sub, seg_panel[t], by=by) for t in tissues}
    return out


def plot_composition_heatmaps(mats_by_cutoff, by="state", outpath=None,
                              cmap="viridis", suptitle=None):
    """Draw one ChromHMM state-composition heatmap per top-N cutoff.

    mats_by_cutoff : dict {cutoff_label: states x tissues DataFrame} from
                     composition_matrix. cutoff_label may be an int (top-N) or a
                     string; panels are drawn left->right in dict order and share
                     the state (y) axis. A single DataFrame is also accepted.
    by             : 'state' or 'group' (only affects the y-axis title).
    outpath        : if given, fig.savefig(outpath) at 200 dpi (bbox tight).
    Returns (fig, axes). Cells are annotated with the fraction; the color scale
    is shared across panels (sequential, 0..max).
    """
    import matplotlib.pyplot as plt
    import numpy as np

    if isinstance(mats_by_cutoff, pd.DataFrame):
        mats_by_cutoff = {"": mats_by_cutoff}
    labels = list(mats_by_cutoff.keys())

    # shared row set: canonical order filtered to states nonzero in ANY panel
    if by == "state":
        order = list(R.STATES)
    else:
        order = list(dict.fromkeys(R.GROUP.values()))
    nz = set()
    for m in mats_by_cutoff.values():
        nz |= set(m.index[m.sum(axis=1) > 0])
    rows = [s for s in order if s in nz] or order
    tissues = list(next(iter(mats_by_cutoff.values())).columns)
    vmax = max((float(m.values.max()) for m in mats_by_cutoff.values()
                if m.size), default=1.0) or 1.0

    n = len(labels)
    fig, axes = plt.subplots(
        1, n, figsize=(max(2.6, 1.05 * len(tissues)) * n + 1.2,
                       0.34 * len(rows) + 1.4),
        squeeze=False, sharey=True)
    axes = axes[0]

    im = None
    for ax, lab in zip(axes, labels):
        m = mats_by_cutoff[lab].reindex(index=rows, columns=tissues).fillna(0.0)
        im = ax.imshow(m.values, aspect="auto", cmap=cmap, vmin=0.0, vmax=vmax)
        ax.set_xticks(range(len(tissues)))
        ax.set_xticklabels(tissues, rotation=45, ha="right", fontsize=8)
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels(rows, fontsize=7)
        title = f"top {lab}" if isinstance(lab, int) else str(lab)
        ax.set_title(title, fontsize=9)
        # annotate cells (<200 -> always readable)
        thr = 0.6 * vmax
        for i in range(len(rows)):
            for j in range(len(tissues)):
                v = m.values[i, j]
                if v <= 0:
                    continue
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        fontsize=6.2,
                        color="white" if v >= thr else "black")
        ax.set_xticks(np.arange(-.5, len(tissues), 1), minor=True)
        ax.set_yticks(np.arange(-.5, len(rows), 1), minor=True)
        ax.grid(which="minor", color="white", linewidth=0.5)
        ax.tick_params(which="minor", length=0)

    axes[0].set_ylabel("ChromHMM " + ("state" if by == "state" else "group"),
                       fontsize=9)
    cbar = fig.colorbar(im, ax=list(axes), fraction=0.025, pad=0.02)
    cbar.set_label("fraction of bins (dominant state)", fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    if suptitle:
        fig.suptitle(suptitle, fontsize=10)
    if outpath:
        fig.savefig(outpath, dpi=200, bbox_inches="tight")
    return fig, axes
