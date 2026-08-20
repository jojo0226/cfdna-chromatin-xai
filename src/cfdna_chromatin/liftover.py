"""
liftover.py -- controlled genome-build crossing for the cfDNA explainer.

The quantitative pipeline (bins, coverage, end-motifs, model, SHAP) is
build-native: an hg19 cohort stays in hg19 end-to-end. The chromatin reference
layer, however, is hg38. This module crosses that boundary in exactly ONE place
and on COORDINATES ONLY (never per-sample data), with round-trip QC so a bad
lift is caught rather than silently mislocating a region.

Recommended direction: pre-lift the (small, static) hg38 reference peak sets
DOWN to hg19 once, so an hg19 analysis runs natively against an hg19 reference
bundle. The reverse (lift the handful of top-SHAP hg19 bins UP to hg38 per
analysis) is also supported via lift_regions().

Backed by pyliftover (pure-Python; downloads/uses UCSC chain files). No binary
liftOver needed. Chain files:
    hg38ToHg19.over.chain.gz   (down)
    hg19ToHg38.over.chain.gz   (up)

Liftover hygiene enforced here:
    * drop regions that fail to map,
    * drop regions whose endpoints map to a DIFFERENT chromosome,
    * drop regions that invert or collapse (end <= start after lift),
    * optional round-trip check: lift back and require overlap within `rt_tol`.
Every call returns both the lifted frame and a QC dict with the mapping rate.
"""
from __future__ import annotations

import pandas as pd


def get_lifter(from_build, to_build, chain_path=None):
    """Return a pyliftover.LiftOver for from_build->to_build.

    chain_path : local chain file; if None, pyliftover fetches the UCSC chain
                 for the (from_build, to_build) pair on first use.
    """
    from pyliftover import LiftOver
    if chain_path is not None:
        return LiftOver(chain_path)
    return LiftOver(from_build, to_build)


def _lift_point(lifter, chrom, pos):
    """Lift a single 0-based point; return (chrom, pos) or None."""
    res = lifter.convert_coordinate(chrom, int(pos))
    if not res:
        return None
    c, p = res[0][0], res[0][1]
    return c, p


def lift_regions(regions, from_build="hg19", to_build="hg38",
                 chain_path=None, lifter=None, round_trip=True, rt_tol=2):
    """Lift a list/frame of regions between builds with QC.

    regions : list of (chrom,start,end) OR a DataFrame with chrom/start/end
              columns (any extra columns are preserved on kept rows).
    Returns (lifted, qc):
        lifted : DataFrame with lifted chrom/start/end (+ preserved columns),
                 index-aligned to the *kept* input rows.
        qc     : dict(n_in, n_mapped, mapping_rate, dropped_unmapped,
                      dropped_diff_chrom, dropped_degenerate, dropped_round_trip).
    """
    if lifter is None:
        lifter = get_lifter(from_build, to_build, chain_path)
    back = None
    if round_trip:
        back = get_lifter(to_build, from_build, None if chain_path else None)

    if isinstance(regions, pd.DataFrame):
        df = regions.copy()
        recs = list(df[["chrom", "start", "end"]].itertuples(index=True, name=None))
        extra_cols = [c for c in df.columns if c not in ("chrom", "start", "end")]
    else:
        recs = [(i, r[0], r[1], r[2]) for i, r in enumerate(regions)]
        df = None
        extra_cols = []

    kept_idx, rows = [], []
    qc = dict(n_in=len(recs), n_mapped=0, dropped_unmapped=0,
              dropped_diff_chrom=0, dropped_degenerate=0, dropped_round_trip=0)

    for rec in recs:
        idx, chrom, start, end = rec
        a = _lift_point(lifter, chrom, start)
        b = _lift_point(lifter, chrom, end)
        if a is None or b is None:
            qc["dropped_unmapped"] += 1
            continue
        (ca, pa), (cb, pb) = a, b
        if ca != cb:
            qc["dropped_diff_chrom"] += 1
            continue
        ns, ne = sorted((pa, pb))
        if ne <= ns:
            qc["dropped_degenerate"] += 1
            continue
        if round_trip:
            ra = _lift_point(back, ca, ns)
            rb = _lift_point(back, cb, ne)
            if ra is None or rb is None or ra[0] != chrom or rb[0] != chrom:
                qc["dropped_round_trip"] += 1
                continue
            if abs(ra[1] - start) > rt_tol * (end - start + 1) and \
               abs(rb[1] - end) > rt_tol * (end - start + 1):
                qc["dropped_round_trip"] += 1
                continue
        kept_idx.append(idx)
        rows.append((ca, ns, ne))
        qc["n_mapped"] += 1

    qc["mapping_rate"] = qc["n_mapped"] / qc["n_in"] if qc["n_in"] else 0.0

    lifted = pd.DataFrame(rows, columns=["chrom", "start", "end"],
                          index=kept_idx if df is not None else range(len(rows)))
    if df is not None and extra_cols:
        lifted = lifted.join(df.loc[kept_idx, extra_cols])
    return lifted, qc


def lift_reference_bed(in_bed_gz, out_bed_gz, from_build="hg38", to_build="hg19",
                       chain_path=None, round_trip=False):
    """Pre-lift a reference peak BED(.gz) between builds (the recommended path).

    Reads a (chrom,start,end,...) BED, lifts every interval, writes the kept
    intervals to out_bed_gz (gzipped, same extra columns). Returns the QC dict.
    Use this once per reference track to build an hg19 reference bundle from the
    shipped hg38 peaks.
    """
    import gzip
    op = gzip.open if str(in_bed_gz).endswith(".gz") else open
    rows = []
    with op(in_bed_gz, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            rows.append((f[0], int(f[1]), int(f[2]), "\t".join(f[3:])))
    df = pd.DataFrame(rows, columns=["chrom", "start", "end", "_rest"])
    lifted, qc = lift_regions(df, from_build, to_build, chain_path,
                              round_trip=round_trip)
    ow = gzip.open if str(out_bed_gz).endswith(".gz") else open
    with ow(out_bed_gz, "wt") as fh:
        for _, r in lifted.iterrows():
            rest = r["_rest"] if isinstance(r.get("_rest"), str) and r["_rest"] else ""
            fh.write(f"{r['chrom']}\t{int(r['start'])}\t{int(r['end'])}"
                     + (f"\t{rest}" if rest else "") + "\n")
    return qc
