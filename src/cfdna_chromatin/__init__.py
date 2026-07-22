"""cfdna_chromatin -- explainable chromatin-state annotation for cfDNA regions (hg38).

Public API:
    references : ChromHMM 18-state vocabulary, groups, reference loading
    genome     : GC / mappability bookkeeping for matched-null sampling
    engine     : annotation + matched-null enrichment test
    histone    : histone-mark fingerprint layer (orthogonal cross-check)
    accessibility : open-chromatin (DNase-seq) layer (orthogonal cross-check)
    benchmark  : B1 statistical-calibration harness
    selection  : importance/selection front-end (SHAP or differential -> ranked bins)
    attribute  : compartment attribution (signal / inflammation / background) of regions
    liftover   : controlled hg19<->hg38 build crossing (coordinates only, QC'd)
"""
from . import (  # noqa: F401
    references, genome, engine, histone, accessibility, benchmark, attribute, selection,
    liftover,
)

__version__ = "0.5.2"
