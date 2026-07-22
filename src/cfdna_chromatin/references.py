"""
references.py -- ChromHMM 18-state model vocabulary and reference-epigenome loading.

Reference model: Roadmap / EpiMap ChromHMM 18-state (the "core 18-state" model,
built on the 6-mark imputed panel). The four bundled reference epigenomes share
an identical state vocabulary:

    placenta      ENCFF024IDF   BSS01443  placenta male embryo (85 days)   [fetal signal]
    neutrophil    ENCFF412COT   BSS01381  neutrophil male                  [hematopoietic bg]
    Bcell         ENCFF004HFC   BSS00095  B cell female adult (27 yr)      [hematopoietic bg]
    K562          ENCFF026ZCC   BSS01039  K562                             [erythroleukemia control]
    keratinocyte  ENCFF131MHN             keratinocyte                     [squamous-epithelial reference]

The references are organised into application *panels* (see PANELS). The shipped
panel is "fetal" (placental signal on a maternal hematopoietic background, for
seqFF++ fetal fraction). The chromatin engine is panel-agnostic: any per-region
query set can be run against any panel, and additional panels can be declared
from the bundled references without touching the engine. All five segmentations
share the identical ChromHMM 18-state vocabulary, and matching histone (6-mark)
and DNase accessibility tracks are bundled per tissue.

All coordinates are 0-based half-open, hg38.
"""
from __future__ import annotations
import gzip
from collections import defaultdict
import numpy as np

# ---- 18-state vocabulary (canonical order) ----
STATES = [
    "TssA", "TssFlnk", "TssFlnkU", "TssFlnkD", "Tx", "TxWk",
    "EnhG1", "EnhG2", "EnhA1", "EnhA2", "EnhWk", "ZNF/Rpts",
    "Het", "TssBiv", "EnhBiv", "ReprPC", "ReprPCWk", "Quies",
]
S2I = {s: i for i, s in enumerate(STATES)}
I2S = {i: s for s, i in S2I.items()}

# ---- coarse functional groups (direction-aware biology in B2/B3) ----
GROUP = {
    "TssA": "Promoter", "TssFlnk": "Promoter", "TssFlnkU": "Promoter", "TssFlnkD": "Promoter",
    "TssBiv": "Bivalent", "EnhBiv": "Bivalent",
    "EnhA1": "Enhancer", "EnhA2": "Enhancer", "EnhG1": "Enhancer", "EnhG2": "Enhancer", "EnhWk": "Enhancer",
    "Tx": "Transcription", "TxWk": "Transcription",
    "ZNF/Rpts": "ZNF_Repeat", "Het": "Heterochromatin",
    "ReprPC": "Polycomb", "ReprPCWk": "Polycomb", "Quies": "Quiescent",
}

# active (open/expressed) states -- used as the positive label for active-vs-inactive AUROC
ACTIVE = {"TssA", "TssFlnk", "TssFlnkU", "TssFlnkD", "EnhA1", "EnhA2", "EnhG1", "EnhG2", "EnhWk", "Tx"}

# ENCODE file accessions for the bundled references (ChromHMM 18-state segmentations)
REFERENCE_ACCESSIONS = {
    "placenta":     "ENCFF024IDF",
    "neutrophil":   "ENCFF412COT",
    "Bcell":        "ENCFF004HFC",
    "K562":         "ENCFF026ZCC",
    "keratinocyte": "ENCFF131MHN",   # BSS keratinocyte; squamous-epithelial reference track
}

# ---- Reference panels -------------------------------------------------------
# A *panel* is a named set of reference epigenomes chosen for one cfDNA
# application, each tagged with the role it plays in interpretation:
#   "signal"     -- the tissue whose fingerprint we expect to be ENRICHED in
#                   the important/differential regions (the cell-of-origin hypothesis)
#   "background" -- tissues that dominate the plasma cfDNA pool in that context
#                   (hematopoietic in both pregnancy and cancer); their fingerprint
#                   is the null we contrast against
#   "control"    -- a cell line / stromal tissue included as a specificity check
#   "inflammation" -- an optional activated-myeloid axis (e.g. CD14+ monocyte) used
#                   to LABEL, and separate out, regions driven by a systemic immune
#                   response rather than the signal tissue. When a plasma pool carries
#                   a host inflammatory shift that moves fragmentomics genome-wide
#                   independent of the signal source, this axis is what lets
#                   attribute_signal() distinguish an "inflammation" call from a
#                   genuine "signal" (cell-of-origin) call.
#
# The chromatin engine is feature-agnostic: any per-region query set can be run
# against any panel. Panels only declare which references to load and how to read
# the resulting enrichment (which tissue is the expected positive).
PANELS = {
    # seqFF / fetal-fraction: placental signal on a maternal hematopoietic background
    "fetal": {
        "signal":     ["placenta"],
        "background": ["neutrophil", "Bcell"],
        "control":    ["K562"],
        "note": "Fetal fraction (seqFF++). Placenta is the cell-of-origin of fetal cfDNA; "
                "maternal plasma background is myeloid + lymphoid.",
    },
}

# roles that carry an expected-positive (cell-of-origin) interpretation vs. roles
# used as contrast/labels. Order matters for panel_tissues() output.
_ROLE_ORDER = ("signal", "background", "inflammation", "control")


def panel_tissues(panel_name):
    """Flat ordered list of all tissues in a panel (signal, background, inflammation, control)."""
    p = PANELS[panel_name]
    seen, out = set(), []
    for role in _ROLE_ORDER:
        for t in p.get(role, []):
            if t not in seen:
                seen.add(t)
                out.append(t)
    return out


def tissue_role(panel_name, tissue):
    """Return the role of a tissue within a panel ('signal'|'background'|'inflammation'|'control')."""
    for role in _ROLE_ORDER:
        if tissue in PANELS[panel_name].get(role, []):
            return role
    return None


def load_segmentation(bed_gz_path, chroms=None):
    """Load one ChromHMM .bed.gz into per-chromosome sorted segment arrays.

    Returns {chrom: {"starts": int64[], "ends": int64[], "code": int8[]}} with
    segments sorted by start (they are non-overlapping and gap-free by construction).
    """
    tmp = defaultdict(list)
    keep = set(chroms) if chroms is not None else None
    with gzip.open(bed_gz_path, "rt") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            c = p[0]
            if keep is not None and c not in keep:
                continue
            state = p[3]
            if state not in S2I:
                raise ValueError(f"Unknown ChromHMM state {state!r} in {bed_gz_path}")
            tmp[c].append((int(p[1]), int(p[2]), S2I[state]))
    seg = {}
    for c, rows in tmp.items():
        rows.sort()
        arr = np.array(rows, dtype=np.int64)
        seg[c] = {"starts": arr[:, 0], "ends": arr[:, 1], "code": arr[:, 2].astype(np.int8)}
    return seg
