"""
attribute.py -- compartment attribution for differential cfDNA regions.

The chromatin/accessibility layers tell you *whether* a region's chromatin is
open in a given reference tissue. attribute_signal() turns that into the answer
the cell-of-origin question actually needs: for each region, which compartment's
active chromatin best explains it --

    "signal"       -- open in the panel's signal tissue, closed elsewhere
                      (the cell-of-origin hypothesis)
    "inflammation" -- open in the activated-myeloid axis (e.g. monocyte),
                      not signal-specific
    "background"   -- open in the resting hematopoietic background (neutrophil/Bcell);
                      i.e. ordinary blood-derived cfDNA, no cell-of-origin information
    "ambiguous"    -- open in several compartments at once (shared/housekeeping)
    "unexplained"  -- not clearly open in any reference (quiescent / novel)

This is the step that lets you state *how much* of a region set concentrates in
the cell-of-origin (signal) compartment versus a background/host-response
compartment, rather than reporting a single confounded aggregate score. It is
deliberately conservative: a region is only called "signal" when the signal
compartment is both active AND clearly more active than every non-signal
compartment.

Uses the histone active-mark fingerprint (mean peak-coverage over the 4 active
marks) as the per-compartment activity score. Accessibility is available as a
tie-breaker/annotation but is not used for the primary call when a panel's
inflammation-axis DNase shares an ENCODE experiment with the background proxy;
in that case the inflammation axis is read from the histone marks.
"""
from __future__ import annotations
import numpy as np

from . import references as R
from . import histone as H

ACTIVE_MARKS = ("H3K4me3", "H3K4me1", "H3K27ac", "H3K36me3")


def _activity(chrom, start, end, tissue_panel):
    """Mean active-mark peak coverage fraction for one region in one tissue (0..1)."""
    fp = H.region_mark_fingerprint(chrom, start, end, tissue_panel, marks=H.MARKS)
    cov = fp["coverage"]
    return float(np.mean([cov[m] for m in ACTIVE_MARKS]))


def _compartment_map(panel_name):
    """tissue -> compartment label, collapsing background tissues into one 'background'."""
    p = R.PANELS[panel_name]
    cmap = {}
    for t in p.get("signal", []):
        cmap[t] = "signal"
    for t in p.get("inflammation", []):
        cmap[t] = "inflammation"
    for t in p.get("background", []):
        cmap[t] = "background"
    for t in p.get("control", []):
        cmap[t] = "control"
    return cmap


def attribute_region(chrom, start, end, hpanel, panel_name,
                     min_active=0.05, margin=0.5):
    """Attribute one region to a compartment.

    Parameters
    ----------
    hpanel : dict tissue -> mark panel (from histone.load_mark_panel)
    panel_name : key into references.PANELS (e.g. 'fetal')
    min_active : minimum active-mark coverage for a compartment to count as "open"
                 at all; below this everywhere -> "unexplained".
    margin : how much higher (log2) the top compartment's activity must be over the
             best NON-signal compartment for a confident "signal" call; and over the
             runner-up generally. Expressed as a log2 ratio with a small pseudocount.

    Returns a dict: {chrom,start,end, compartment, scores{...}, top, runner_up,
                     log2_margin, called_open}.
    """
    cmap = _compartment_map(panel_name)
    # per-compartment activity = max over that compartment's tissues
    comp_scores = {}
    for tissue, comp in cmap.items():
        if tissue not in hpanel:
            continue
        a = _activity(chrom, start, end, hpanel[tissue])
        comp_scores[comp] = max(comp_scores.get(comp, 0.0), a)

    eps = 1e-3
    # rank the interpretive compartments (exclude 'control' from the call itself;
    # control is a specificity annotation, not a plausible cfDNA source)
    interp = {c: s for c, s in comp_scores.items() if c != "control"}
    if not interp or max(interp.values()) < min_active:
        return {"chrom": chrom, "start": start, "end": end,
                "compartment": "unexplained", "scores": comp_scores,
                "top": None, "runner_up": None, "log2_margin": 0.0,
                "called_open": False}

    ordered = sorted(interp.items(), key=lambda kv: kv[1], reverse=True)
    top, top_s = ordered[0]
    run, run_s = ordered[1] if len(ordered) > 1 else (None, 0.0)
    log2_margin = float(np.log2((top_s + eps) / (run_s + eps)))

    # signal call is held to the stricter bar: must beat the best non-signal compartment
    if top == "signal":
        best_nonsignal = max([s for c, s in interp.items() if c != "signal"], default=0.0)
        margin_vs_nonsignal = float(np.log2((top_s + eps) / (best_nonsignal + eps)))
        comp = "signal" if margin_vs_nonsignal >= margin else "ambiguous"
        log2_margin = margin_vs_nonsignal
    else:
        comp = top if log2_margin >= margin else "ambiguous"

    return {"chrom": chrom, "start": start, "end": end,
            "compartment": comp, "scores": comp_scores,
            "top": top, "runner_up": run, "log2_margin": log2_margin,
            "called_open": True}


def attribute_signal(regions, hpanel, panel_name, min_active=0.05, margin=0.5):
    """Attribute a set of differential regions and summarise the compartment mix.

    regions : iterable of (chrom, start, end)
    Returns (per_region list, summary dict). The summary's 'fractions' is the
    headline result: what fraction of the region set is signal (cell-of-origin)
    vs inflammation vs background vs ambiguous vs unexplained.
    """
    per = [attribute_region(c, s, e, hpanel, panel_name,
                            min_active=min_active, margin=margin)
           for (c, s, e) in regions]
    n = len(per)
    counts = {}
    for r in per:
        counts[r["compartment"]] = counts.get(r["compartment"], 0) + 1
    fractions = {k: v / n for k, v in counts.items()} if n else {}
    summary = {"n": n, "counts": counts, "fractions": fractions,
               "panel": panel_name, "min_active": min_active, "margin": margin,
               "signal_fraction": fractions.get("signal", 0.0),
               "background_fraction": fractions.get("inflammation", 0.0) + fractions.get("background", 0.0)}
    return per, summary
