"""
build_openness_atlas.py -- build a compartment-openness atlas aggregating MULTIPLE
REPLICATES per tissue x assay (the NIPT/seqFF++ fork).

TWO ATLAS SOURCES (same replicate-normalization engine, pick with `source=`)
---------------------------------------------------------------------------
  source="local"  Option 1: the here-built atlas -- the SAME ENCODE peaks the shipped
                  reference/compartment_openness_atlas.csv.gz was made from
                  (data/access/, data/histone/; 4-col bed.gz; hg38), but now you may
                  list SEVERAL replicate files per tissue x assay instead of one.
                  No S3, no archive wait; runs anywhere. Grid: hg38, default 500 kb.

  source="s3"     Option 2: the S3 atlas_of_chromatin (hg19, self-describing csv.gz).
                  Richer replicate pool, but objects are in Intelligent-Tiering
                  ARCHIVE ACCESS -- you MUST batch_restore() and wait ~3-5 h before the
                  peak reads succeed (see s3_atlas.py). Grid: hg19, default 50 kb.

Why replicates need care (both sources)
---------------------------------------
Peak files are NOT directly averageable: replicates differ in sequencing depth and in
the peak-caller threshold, so a deeply-sequenced replicate has systematically more/wider
peaks and would dominate a raw mean. This is what makes "adding replicates" tricky.

The fix: per-replicate NORMALIZATION before aggregation.
  1. bin each replicate's peaks to the grid -> per-bin coverage vector
  2. standardize EACH replicate across bins (z-score) so depth/threshold cancel
  3. aggregate the standardized replicates within tissue x assay (mean)
  4. average the four active-histone marks -> H_<tissue>; DNase -> D_<tissue>
  5. build the same *_specific contrasts as the reference atlas

For source="s3", s3_atlas.py (build_catalog / read_peaks) provides metadata-first access.
"""
from __future__ import annotations

import glob
import gzip
import os

import numpy as np
import pandas as pd

# hg19 autosome lengths (UCSC) -- the S3 atlas assembly
HG19_AUTOSOMES = {
    "chr1": 249250621, "chr2": 243199373, "chr3": 198022430, "chr4": 191154276,
    "chr5": 180915260, "chr6": 171115067, "chr7": 159138663, "chr8": 146364022,
    "chr9": 141213431, "chr10": 135534747, "chr11": 135006516, "chr12": 133851895,
    "chr13": 115169878, "chr14": 107349540, "chr15": 102531392, "chr16": 90354753,
    "chr17": 81195210, "chr18": 78077248, "chr19": 59128983, "chr20": 63025520,
    "chr21": 48129895, "chr22": 51304566,
}
# hg38 autosome lengths (UCSC) -- the local (data/access, data/histone) assembly
HG38_AUTOSOMES = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559, "chr4": 190214555,
    "chr5": 181538259, "chr6": 170805979, "chr7": 159345973, "chr8": 145138636,
    "chr9": 138394717, "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189, "chr16": 90338345,
    "chr17": 83257441, "chr18": 80373285, "chr19": 58617616, "chr20": 64444167,
    "chr21": 46709983, "chr22": 50818468,
}
ASSEMBLIES = {"hg19": HG19_AUTOSOMES, "hg38": HG38_AUTOSOMES}
ACTIVE_MARKS = ("H3K4me3", "H3K4me1", "H3K27ac", "H3K36me3")

# ── tissue-group harmonization (ENCODE free-text biosample -> clean group) ──────
# Cell-free DNA in plasma is a mixture of the tissues that shed it. In the MATERNAL
# plasma / NIPT setting the fetal fraction comes from PLACENTA (trophoblast); the
# maternal background is dominated by HEMATOPOIETIC cells, with smaller solid-tissue
# contributions (liver, endothelium, and — in disease — other organs). This map lets
# you group many ENCODE biosamples (each a free-text `Biosample term name`) into a
# fixed set of origin labels so replicates from different donors/experiments pool.
#
# Ordered: more specific patterns first (e.g. CD8 before generic "t cell", HepG2/K562
# carcinoma lines before their lineage) so the FIRST match wins. Patterns are matched
# case-insensitively as substrings of the biosample name. Verify against your own
# catalog's `Biosample term name` values and extend as needed.
TISSUE_GROUP_PATTERNS = [
    # --- fetal-derived (the fetal-fraction signal) ---
    ("placenta",      ["placenta", "trophoblast", "cytotrophoblast", "chorion", "chorionic", "amnion"]),
    # --- hematopoietic / blood (dominant maternal background) ---
    ("neutrophil",    ["neutrophil"]),
    ("monocyte",      ["monocyte", "cd14"]),
    ("macrophage",    ["macrophage"]),
    ("Bcell",         ["b cell", "b-cell", "cd19", "cd20"]),
    ("CD8",           ["cd8"]),
    ("CD4",           ["cd4"]),
    ("Tcell",         ["t cell", "t-cell", "cd3", "thymocyte"]),
    ("NK",            ["natural killer", "cd56"]),
    ("erythroblast",  ["erythroblast", "erythroid", "erythrocyte"]),
    ("megakaryocyte", ["megakaryocyte", "platelet"]),
    ("HSC",           ["hematopoietic stem", "cd34", "myeloid progenitor",
                       "lymphoid progenitor", "hematopoietic multipotent"]),
    ("K562",          ["k562"]),        # CML line, historical control
    # --- solid-organ contributors to cfDNA (minor in health, larger in disease) ---
    ("liver",         ["hepatocyte", "liver", "hepg2"]),
    ("lung",          ["lung", "pneumocyte", "bronchial", "alveolar"]),
    ("colon",         ["colon", "colonic", "sigmoid", "intestine", "intestinal", "rectal"]),
    ("kidney",        ["kidney", "renal", "nephron"]),
    ("pancreas",      ["pancreas", "pancreatic", "islet"]),
    ("heart",         ["heart", "cardiac", "cardiomyocyte", "ventricle", "atrium", "aorta"]),
    ("endothelial",   ["endothelial", "huvec", "vein", "artery", "arterial", "vascular"]),
    ("adipose",       ["adipose", "adipocyte"]),
    ("breast",        ["breast", "mammary", "mcf-7", "mcf7"]),
    ("skeletal_muscle", ["skeletal muscle", "myotube", "myoblast", "myocyte"]),
    ("brain",         ["brain", "neuron", "neural", "cortex", "astrocyte", "cerebell"]),
    ("esophagus",     ["esophagus", "esophageal"]),
    ("stomach",       ["stomach", "gastric"]),
    ("spleen",        ["spleen", "splenic"]),
    ("thyroid",       ["thyroid"]),
    ("adrenal",       ["adrenal"]),
    ("ovary",         ["ovary", "ovarian", "uterus", "uterine", "endometrium", "endometrial"]),
    ("prostate",      ["prostate", "prostatic"]),
    ("keratinocyte",  ["keratinocyte", "epidermis", "epidermal"]),  # epithelial reference
]

# Origin groups most relevant to a maternal-plasma / fetal-fraction cfDNA panel,
# in rough order of expected contribution. Use as a default panel skeleton.
CFDNA_FF_GROUPS = [
    "placenta",                                        # fetal fraction
    "neutrophil", "monocyte", "Bcell", "CD4", "CD8", "NK",  # blood (maternal bg)
    "erythroblast", "megakaryocyte", "HSC",
    "liver", "endothelial", "lung", "colon", "kidney", "pancreas", "heart",
    "breast", "adipose", "skeletal_muscle", "brain",
]


def harmonize_tissue(biosample_name: str) -> str | None:
    """Map one ENCODE free-text biosample name to a clean origin-group label.

    Returns the group label (e.g. 'placenta', 'monocyte', 'liver') or None if no
    pattern matches. First match in TISSUE_GROUP_PATTERNS wins (specific before generic).
    """
    if not isinstance(biosample_name, str):
        return None
    s = biosample_name.lower()
    for group, pats in TISSUE_GROUP_PATTERNS:
        if any(p in s for p in pats):
            return group
    return None


def add_tissue_group(catalog: pd.DataFrame, name_col: str = "Biosample term name",
                     group_col: str = "group") -> pd.DataFrame:
    """Add a harmonized `group` column to an s3_atlas catalog (see harmonize_tissue).

    Rows whose biosample name matches no pattern get group=None -- inspect those with
    catalog[catalog[group_col].isna()][name_col].value_counts() and extend the map.
    """
    col = name_col if name_col in catalog.columns else next(
        (c for c in catalog.columns if "biosample" in c.lower() and "name" in c.lower()), None)
    if col is None:
        raise KeyError(f"no biosample-name column found in {list(catalog.columns)[:8]}")
    out = catalog.copy()
    out[group_col] = out[col].map(harmonize_tissue)
    return out


# ── grid ──────────────────────────────────────────────────────────────────────
def bin_grid(bin_size: int = 50_000, assembly: str = "hg19") -> pd.DataFrame:
    """Fixed-width autosomal bin grid (chrom, start, end, key) for hg19 or hg38."""
    rows = []
    for chrom, size in ASSEMBLIES[assembly].items():
        starts = np.arange(0, size, bin_size, dtype=int)
        for s in starts:
            e = min(int(s) + bin_size, size)
            rows.append((chrom, int(s), e))
    g = pd.DataFrame(rows, columns=["chrom", "start", "end"])
    g["key"] = g.chrom + ":" + g.start.astype(str) + "-" + g.end.astype(str)
    return g


# ── local peak loader (Option 1: same format as the shipped reference) ───────
def read_peaks_local(path: str, autosomes_only=True) -> pd.DataFrame:
    """Load one local peak file into chrom/start/end/signalValue.

    Handles both shipped formats: the 4-column bed.gz in data/access & data/histone
    (chrom, start, end, signalValue) and a named-column csv.gz (S3-style) if present.
    """
    with gzip.open(path, "rt") as fh:
        head = fh.readline()
    if "," in head and "seqnames" in head:  # S3-style csv saved locally
        df = pd.read_csv(path, usecols=lambda c: c in (
            "seqnames", "start", "end", "signalValue"))
        df = df.rename(columns={"seqnames": "chrom"})
    else:  # 4-column ENCODE narrowPeak subset (chrom start end signalValue)
        df = pd.read_csv(path, sep="\t", header=None,
                         names=["chrom", "start", "end", "signalValue"],
                         usecols=[0, 1, 2, 3])
    if autosomes_only:
        auto = {f"chr{i}" for i in range(1, 23)}
        df = df[df["chrom"].isin(auto)]
    return df


# ── one replicate -> per-bin coverage ────────────────────────────────────────
def peaks_to_coverage(peaks: pd.DataFrame, grid: pd.DataFrame,
                      bin_size: int, weight: str | None = None) -> np.ndarray:
    """Per-bin coverage for one replicate, aligned to `grid` row order.

    coverage = fraction of each bin covered by peaks (0..1); if weight='signalValue',
    the covered fraction is multiplied by the mean signalValue of the overlapping peaks
    (fold-enrichment-weighted openness). Returns a vector len(grid).
    """
    # index bins by (chrom, start) for O(1) placement
    bin_index = {(c, s): i for i, (c, s) in
                 enumerate(zip(grid.chrom.values, grid.start.values))}
    cov = np.zeros(len(grid), dtype=float)
    wsum = np.zeros(len(grid), dtype=float)
    wcnt = np.zeros(len(grid), dtype=float)
    for chrom, s, e, sig in zip(peaks.chrom.values, peaks.start.values,
                                peaks.end.values, peaks.get(
                                    "signalValue", pd.Series(np.ones(len(peaks)))).values):
        s = int(s)
        e = int(e)
        b0 = (s // bin_size) * bin_size
        b1 = (e // bin_size) * bin_size
        for bs in range(b0, b1 + 1, bin_size):
            i = bin_index.get((chrom, bs))
            if i is None:
                continue
            ov = min(e, bs + bin_size) - max(s, bs)
            if ov <= 0:
                continue
            cov[i] += ov / bin_size
            wsum[i] += sig * ov
            wcnt[i] += ov
    cov = np.clip(cov, 0.0, 1.0)
    if weight == "signalValue":
        mean_sig = np.divide(wsum, wcnt, out=np.zeros_like(wsum), where=wcnt > 0)
        return cov * mean_sig
    return cov


# ── replicate aggregation (the depth-robust step) ────────────────────────────
def aggregate_replicates(mat: np.ndarray, method: str = "zmean") -> np.ndarray:
    """Combine replicate coverage vectors (rows = replicates) into one track.

    'zmean'  : z-score EACH replicate across bins, then mean -> depth/threshold-robust
               (recommended; this is what makes replicates comparable).
    'mean'   : plain mean (biased toward the deepest replicate; for comparison only).
    'median' : per-bin median across replicates (robust to one outlier replicate).
    """
    mat = np.atleast_2d(mat).astype(float)
    if mat.shape[0] == 1:
        row = mat[0]
        return (row - np.nanmean(row)) / (np.nanstd(row) or 1.0) if method == "zmean" else row
    if method == "zmean":
        mu = np.nanmean(mat, axis=1, keepdims=True)
        sd = np.nanstd(mat, axis=1, keepdims=True)
        sd[sd == 0] = 1.0
        return np.nanmean((mat - mu) / sd, axis=0)
    if method == "median":
        return np.nanmedian(mat, axis=0)
    return np.nanmean(mat, axis=0)


# ── replicate list -> one aggregated track (source-agnostic core) ────────────
def _aggregate_track(replicates: list[pd.DataFrame], grid: pd.DataFrame,
                     bin_size: int, agg: str, weight: str | None
                     ) -> tuple[np.ndarray, int]:
    """Bin + normalize + aggregate a list of already-loaded replicate peak frames."""
    vecs = [peaks_to_coverage(pk, grid, bin_size, weight=weight)
            for pk in replicates if len(pk)]
    if not vecs:
        return np.full(len(grid), np.nan), 0
    return aggregate_replicates(np.vstack(vecs), method=agg), len(vecs)


def _load_replicates(entries, source: str, s3=None) -> list[pd.DataFrame]:
    """Turn a panel's per-assay list into loaded peak frames.

    Option 1 (source="local"): entries are file PATHS or globs under the repo
        (e.g. "data/access/keratinocyte_DNase_*.bed.gz"); read with read_peaks_local.
    Option 2 (source="s3"): entries are ACCESSIONS or full s3 keys; read with
        s3_atlas.read_peaks (which handles the archive tier).
    """
    frames = []
    if source == "local":
        paths = []
        for e in entries:
            paths.extend(sorted(glob.glob(e)) if any(ch in e for ch in "*?[") else [e])
        for p in paths:
            if not os.path.exists(p):
                print(f"  skip (missing): {p}")
                continue
            try:
                frames.append(read_peaks_local(p))
            except Exception as ex:
                print(f"  skip {p}: {type(ex).__name__}: {ex}")
    elif source == "s3":
        from s3_atlas import ATLAS_PREFIX, read_peaks
        for e in entries:
            key = e if e.endswith(".csv.gz") else (
                f"{ATLAS_PREFIX}/chip/ENCODE/bed_narrowPeak_hg19/{e}_bed.csv.gz")
            try:
                frames.append(read_peaks(key, s3=s3))
            except Exception as ex:  # archived/not-restored/missing -> skip, report
                print(f"  skip {key}: {type(ex).__name__}: {ex}")
    else:
        raise ValueError(f"source must be 'local' or 's3', got {source!r}")
    return frames


def build_openness_atlas(tissue_panel: dict, source: str = "local",
                         bin_size: int | None = None, assembly: str | None = None,
                         agg: str = "zmean", weight: str | None = None,
                         s3=None) -> pd.DataFrame:
    """Build the replicate-aggregated compartment-openness atlas from EITHER source.

    source="local"  Option 1 -- entries in tissue_panel are file paths/globs under the
                    repo (same peaks the shipped reference used); assembly defaults to
                    hg38, bin_size to 500 kb. No S3, no archive wait.
    source="s3"     Option 2 -- entries are ENCODE accessions (or full s3 keys) in the
                    archived S3 atlas; assembly hg19, bin_size 50 kb. Restore first
                    (s3_atlas.batch_restore / wait_until_restored) or reads will skip.

    tissue_panel: {tissue: {"DNase": [...], "H3K4me3": [...], ...}} -- any number of
    replicate entries per assay; they are normalized then aggregated. Example (local):
        {"keratinocyte": {"DNase": ["data/access/keratinocyte_DNase_*.bed.gz"],
                          "H3K4me3": ["data/histone/keratinocyte_H3K4me3_*.bed.gz"], ...}}

    Returns chrom,start,end,key + H_<tissue>/D_<tissue> + n_H_/n_D_ replicate counts +
    the reference *_specific contrasts.
    """
    if assembly is None:
        assembly = "hg38" if source == "local" else "hg19"
    if bin_size is None:
        bin_size = 500_000 if source == "local" else 50_000
    grid = bin_grid(bin_size, assembly=assembly)
    out = grid.copy()

    for tissue, assays in tissue_panel.items():
        # accessibility assays: DNase -> D_, ATAC -> A_
        for assay, pref in (("DNase", "D"), ("ATAC", "A")):
            if assays.get(assay):
                reps = _load_replicates(assays[assay], source, s3=s3)
                vec, n = _aggregate_track(reps, grid, bin_size, agg, weight)
                out[f"{pref}_{tissue}"] = vec
                out[f"n_{pref}_{tissue}"] = n
        # combined open-chromatin composite = z-mean of whichever accessibility
        # assays are present (DNase and/or ATAC); lets tissues with only one still line up
        acc = [out[f"{p}_{tissue}"].values.astype(float)
               for p in ("D", "A") if f"{p}_{tissue}" in out]
        if acc:
            zs = [(v - np.nanmean(v)) / (np.nanstd(v) or 1.0) for v in acc]
            out[f"ACC_{tissue}"] = np.nanmean(np.vstack(zs), axis=0)
        # active histone marks -> H_
        mark_tracks, n_marks = [], 0
        for mark in ACTIVE_MARKS:
            if assays.get(mark):
                reps = _load_replicates(assays[mark], source, s3=s3)
                vec, n = _aggregate_track(reps, grid, bin_size, agg, weight)
                if n:
                    mark_tracks.append(vec)
                    n_marks += n
        if mark_tracks:
            out[f"H_{tissue}"] = np.nanmean(np.vstack(mark_tracks), axis=0)
            out[f"n_H_{tissue}"] = n_marks

    return _add_contrasts(out)


def _add_contrasts(out: pd.DataFrame) -> pd.DataFrame:
    """Add the reference atlas composites/contrasts (contrast-first z, then average)."""
    def z(col):
        x = out[col].values.astype(float)
        return (x - np.nanmean(x)) / (np.nanstd(x) or 1.0)

    for a in ("H", "D", "A", "ACC"):
        # hematopoietic background axis (resting granulocyte/B-cell)
        bg_cols = [f"{a}_neutrophil", f"{a}_Bcell"]
        if all(c in out for c in bg_cols):
            out[f"{a}_background"] = out[bg_cols].mean(axis=1)
        # activated-myeloid (inflammation) axis vs background
        if f"{a}_monocyte" in out and f"{a}_background" in out:
            out[f"{a}_infl_specific"] = z(f"{a}_monocyte") - z(f"{a}_background")
        for lym in ("Tcell", "CD8", "NK"):
            if f"{a}_{lym}" in out and f"{a}_background" in out:
                out[f"{a}_{lym}_specific"] = z(f"{a}_{lym}") - z(f"{a}_background")
    # concordance composites
    for name in ("infl_specific",):
        if f"H_{name}" in out and f"D_{name}" in out:
            out[f"C_{name}"] = 0.5 * (out[f"H_{name}"] + out[f"D_{name}"])
    # generic per-tissue combined-openness composite, comparable across a panel of
    # arbitrary tissues (e.g. the fetal-fraction panel):
    #   - histone + accessibility present -> z-mean of z(H) and z(ACC/DNase)
    #   - histone-only tissue (e.g. neutrophil: no ENCODE DNase) -> z(H) alone
    # so a tissue with only one modality still lines up on the same scale.
    # Contrast artifacts from the block above (H_background, H_*_specific) are NOT
    # real tissues, so they are excluded from the C_ composite.
    def _is_contrast(name):
        return name.endswith("_specific") or name in ("background", "signal",
                                                       "immune", "tumor")
    hist_tissues = {c[2:] for c in out.columns if c.startswith("H_")}
    acc_tissues = {c.split("_", 1)[1] for c in out.columns
                   if c.startswith("ACC_") or c.startswith("D_")}
    base_tissues = sorted(t for t in (hist_tissues | acc_tissues)
                          if not _is_contrast(t))
    for t in base_tissues:
        if f"C_{t}" in out:
            continue
        acc_col = f"ACC_{t}" if f"ACC_{t}" in out else (
            f"D_{t}" if f"D_{t}" in out else None)
        has_h = f"H_{t}" in out
        if has_h and acc_col is not None:
            out[f"C_{t}"] = 0.5 * (z(f"H_{t}") + z(acc_col))
        elif has_h:                       # histone-only tissue
            out[f"C_{t}"] = z(f"H_{t}")
        elif acc_col is not None:         # accessibility-only tissue
            out[f"C_{t}"] = z(acc_col)
    return out


def panel_from_catalog(catalog: pd.DataFrame, groups=None,
                       assay_col="assay", target_col="target",
                       acc_col="file_accession") -> dict:
    """Build a {group: {assay: [accessions]}} panel from a harmonized s3 catalog.

    Requires a `group` column (add_tissue_group). Maps the catalog's assay/target
    fields onto the builder's assay keys: DNase-seq->DNase, ATAC-seq->ATAC, and each
    histone ChIP target (H3K4me3, ...) onto its own key. Every replicate accession for
    a group x assay is collected, so tissues with many ENCODE experiments pool them all.
    """
    if "group" not in catalog.columns:
        raise KeyError("catalog needs a 'group' column -- run add_tissue_group first")
    cat = catalog[catalog["group"].notna()].copy()
    if groups is not None:
        cat = cat[cat["group"].isin(groups)]

    def assay_key(assay, target):
        a = str(assay).lower()
        if "dnase" in a:
            return "DNase"
        if "atac" in a:
            return "ATAC"
        t = str(target).strip()
        return t if t in ACTIVE_MARKS else None

    panel: dict = {}
    for _, r in cat.iterrows():
        key = assay_key(r.get(assay_col, ""), r.get(target_col, ""))
        if key is None:
            continue
        panel.setdefault(r["group"], {}).setdefault(key, []).append(r[acc_col])
    return panel


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source", choices=["local", "s3"], default="local",
                    help="local = here-built ENCODE peaks (Option 1); "
                         "s3 = archived S3 atlas (Option 2, restore first)")
    ap.add_argument("--bin-size", type=int, default=None)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    if args.source == "local":
        # Option 1: point each assay at file GLOBS -- add more replicate files (and more
        # tissues) to the data/ folders and they are picked up + depth-normalized. Name
        # files <group>_<assay>_<accession>.bed.gz so one glob collects all replicates.
        # DNase and ATAC both feed accessibility (ACC_); list either or both per tissue.
        A = os.path.join(REPO, "data", "access")
        H = os.path.join(REPO, "data", "histone")

        def local_tissue(group):
            d = {"DNase": [f"{A}/{group}_DNase_*.bed.gz"],
                 "ATAC":  [f"{A}/{group}_ATAC_*.bed.gz"]}
            for m in ACTIVE_MARKS:
                d[m] = [f"{H}/{group}_{m}_*.bed.gz"]
            return d

        # Expand to the cfDNA / fetal-fraction origin groups. Globs that match no file
        # are silently skipped, so this panel is safe even before every tissue is fetched.
        PANEL = {g: local_tissue(g) for g in CFDNA_FF_GROUPS + ["K562"]}
    else:
        # Option 2: drive the panel from the harmonized S3 catalog (restore beds first).
        # This pools EVERY hg19 replicate ENCODE has per origin group x assay.
        import s3_atlas
        cat = s3_atlas.build_catalog(("chip", "atac", "DNase"), genome="hg19")
        cat = add_tissue_group(cat)
        unmatched = cat[cat["group"].isna()]
        if len(unmatched):
            print(f"[panel] {len(unmatched)} tracks unmatched to a group; top names:")
            print(unmatched.get("tissue", unmatched.iloc[:, 0]).value_counts().head(10))
        PANEL = panel_from_catalog(cat, groups=CFDNA_FF_GROUPS + ["K562"])
        print(f"[panel] {len(PANEL)} groups: "
              + ", ".join(f"{g}({sum(len(v) for v in a.values())})" for g, a in PANEL.items()))

    atlas = build_openness_atlas(PANEL, source=args.source, bin_size=args.bin_size)
    out = args.out or f"compartment_openness_atlas_{args.source}.csv.gz"
    atlas.to_csv(out, index=False)
    print(f"wrote {out}: {atlas.shape}")
    print(atlas.filter(regex="^n_").describe())
