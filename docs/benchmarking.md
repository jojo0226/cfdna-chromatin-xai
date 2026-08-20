# Benchmarking strategy

The module is two separable things, benchmarked differently:
**(A)** the enrichment *engine* — is it statistically calibrated? and
**(B)** the *biological claim* — do the state calls match reality? Four tiers,
increasing in stringency. Implementation status is marked on each.

---

## B1 — Statistical calibration &nbsp; `implemented` (`benchmark.run_b1_calibration`)

**Test.** Feed the engine random region sets drawn from the mappable genome — sets
with *no* real chromatin signal. The distribution of enrichment p-values should be
**uniform**, and the realized false-positive rate at nominal α should be ≤ α.

**Metric.** p-value QQ-plot vs. uniform (`docs/figures/b1_calibration.png`, left);
realized FPR vs. nominal α (right). *Run this first — an uncalibrated engine
invalidates every biological result.*

**Current result & known limitation.** The realized FPR is controlled at every state
group (0.000–0.037, all ≤ 0.05). However, the empirical p-values are **conservative,
not perfectly uniform** — the QQ curves sit above the diagonal and a KS test rejects
uniformity.

*Diagnosis:* for rare states (e.g. Promoter, which occupies <1 % of the genome), a
signal-free query set of ~40 regions of 50 kb contains very few dominant-Promoter
regions, and so do the bootstrap null draws. The two-sided empirical p is built from
`|log2 fold| ≥ observed`, and with many draws tied at zero fold-change the mass piles
up near p=1. This is a discreteness artifact, not a bias toward false positives —
the FPR result confirms the engine is *safe* (conservative). Conservative-but-safe is
the correct failure mode for a discovery gate.

*Planned refinement:* (i) test on larger query sets / smaller bins so rare states are
better populated; (ii) add a mid-p correction for the discrete tie mass; (iii) report
calibration per *frequency stratum* so common vs. rare states are judged separately.
Making B1 a CI gate is deferred until the discreteness correction lands, so a
legitimately conservative engine is not flagged as a regression.

---

## B2 — Positive / negative controls &nbsp; `positive control implemented`

**Positive controls.** Curated region sets with a *known* state. Implemented:
±2 kb promoter windows of 22 housekeeping genes → should call active-promoter states.
**Result:** Promoter enrichment log2FC +2.5 to +3.5 (p=0.002) across all four
epigenomes, with Quiescent/Polycomb/Transcription depleted
(`docs/figures/positive_control.png`). To add: FANTOM/VISTA validated enhancers
(→ Enhancer/H3K27ac), CTCF/insulator sites.

**Negative controls** *(planned).* Gene deserts / constitutive heterochromatin →
should call Quies/Het.

**Metric.** Per-state sensitivity/specificity and confusion matrix on the labeled
set; overall AUROC for "active vs. inactive" calls (the `ACTIVE` state set in
`references.py` defines the positive label).

**Orthogonal histone-mark cross-check.** The histone layer (`histone.py`,
`mark_enrichment_test`) tests the *same* housekeeping-promoter set against the
matched null on six ChIP-seq marks — a readout that does not use the ChromHMM
segmentation at all. Because the two annotations are independent, agreement is
genuine corroboration rather than circularity. **Result** (placenta shown;
consistent across all four epigenomes, `docs/figures/histone_fingerprint.png`):

| Mark | log2FC | direction |
|---|---|---|
| H3K4me3  | +3.0 | active promoter ↑ |
| H3K27ac  | +2.2 | active enhancer/promoter ↑ |
| H3K4me1  | +1.4 | enhancer ↑ |
| H3K36me3 | −0.5 | gene body (not promoter) ↓ |
| H3K27me3 | −3.1 | Polycomb ↓ |
| H3K9me3  | −3.8 | heterochromatin ↓ |

The active-promoter marks are enriched (p<0.05) and both repressive marks depleted —
the textbook active-promoter signature, recovered independently of the ChromHMM call
that B2 already validated.

**Orthogonal open-chromatin cross-check.** The accessibility layer
(`accessibility.py`, `access_enrichment_test`) applies the same test to ENCODE
DNase-seq peaks — a third axis, again independent of the ChromHMM segmentation.
Housekeeping promoters are nucleosome-depleted, so they should sit in open
chromatin. **Result** (`docs/figures/accessibility_enrichment.png`):

| Track | log2FC | p | note |
|---|---|---|---|
| placenta (fetal)          | +1.23 | 0.012 | significant ↑ |
| neutrophil (monocyte proxy) | +0.58 | 0.187 | positive, n.s. |
| Bcell                     | +0.77 | 0.098 | positive, n.s. |
| K562                      | +0.31 | 0.144 | positive, n.s. |

All four are positive (promoters more open than the matched null) with the fetal
track significant. The effect is smaller than the histone fingerprint by design:
DNase peaks cover much more of the genome than sharp histone marks, so the
matched-null background is already fairly open (especially for the monocyte and
K562 tracks) and the contrast compresses. This is the honest expected behaviour of
a broad-coverage assay on a small (n=22, chr19) promoter set, not a failure — the
direction is correct on every track and the layer's overlap math is unit-tested
against synthetic peaks.

---

## B3 — Recovery of known cfDNA biology &nbsp; `planned`

**Test.** Run on region sets with an expected tissue signature. Placental /
trophoblast-derived cfDNA regions should enrich for active chromatin in **placental**
epigenomes and be depleted in unrelated tissues; maternal-hematopoietic regions the
reverse. Fetal-enriched (short-fragment) vs. maternal-enriched regions should give
**opposite** signatures (direction-aware test). The four bundled references
(placenta vs. neutrophil/B cell) are exactly the contrast this tier needs.

**Metric.** Tissue-matched fold-enrichment with CI; effect-size separation between
up- and down-region sets.

**Circularity caveat.** If the same tissue reference defines both the state calls and
the "expected" signature, the test is trivially passed. Use independent reference
panels for annotation vs. evaluation.

---

## B4 — Orthogonal ground truth &nbsp; `planned`

**Held-out tissue concordance.** Annotate using one set of reference epigenomes, test
agreement on a **held-out** ENCODE tissue — guards against over-fitting one panel.

**Cross-modality.** Where matched WGBS / enzymatic-methyl data exist, compare the
chromatin-state cell-of-origin inference to **methylation-based tissue-of-origin**
deconvolution on the same samples (independent method, same truth).

**Metric.** Concordance / correlation between the two independent readouts.

---

## Metrics summary

| Benchmark | Question | Primary metric | Status |
|---|---|---|---|
| B1 calibration | Is the null well-behaved? | p-value uniformity (QQ), FPR vs. α | ✅ implemented (conservative; refinement noted) |
| B2 controls | Are state calls correct on known regions? | sensitivity/specificity, AUROC, confusion matrix | ◑ positive control done |
| B3 cfDNA biology | Does it recover expected tissue signal? | tissue-matched fold-enrichment + CI, up/down separation | ○ planned |
| B4 orthogonal | Does it agree with an independent method? | held-out concordance, WGBS correlation | ○ planned |

Report **effect sizes with CIs alongside p-values** everywhere — with large
annotation catalogs, significance is cheap and fold-enrichment is what matters.

---

## The matched-null background (non-negotiable)

Enrichment is meaningless without the right null. Every query region is compared
against background regions matched on **GC content, length, and mappability
(N-fraction)**, drawn from the mappable genome (`engine.NullModel`). Statistics are
computed against this null, **not** against the whole genome — otherwise GC and
mappability biases produce false "active chromatin" enrichment. This null is shared
by every tier and every future XAI layer, so getting it right here pays off
downstream. Coverage-stratum matching (for fragmentomics inputs) is a planned
extension.

---

## Build order

1. ✅ **Engine core** — region intersection + matched-null enrichment.
2. ✅ **B1 calibration harness** — permutation/QQ test (CI gate deferred pending
   discreteness correction).
3. ◑ **B2 control panel** — positive control done; add negative controls + AUROC.
4. ✅ **Histone + accessibility layers** — histone-mark fingerprint (H3K4me3/me1,
   H3K27ac, H3K36me3, H3K27me3, H3K9me3) implemented (`histone.py`) and open-chromatin
   DNase-seq layer implemented (`accessibility.py`); both cross-checked in B2.
5. ✅ **Application panels + compartment attribution** — `references.PANELS`
   (`fetal`) groups references by role (signal / background / inflammation /
   control); `attribute.attribute_signal()` labels each region
   signal / inflammation / background / ambiguous / unexplained. On housekeeping
   promoters (active in all cell types) it returns **0%** confident signal calls
   (100% ambiguous), as expected — the conservative behavior that lets a genuine
   cell-of-origin enrichment stand out. This is the step that separates
   cell-of-origin signal from background/host-response chromatin.
6. ○ **B3 cfDNA biology** — placental vs. hematopoietic, direction-aware.
7. ○ **B4 orthogonal** — held-out tissue + methylation cross-check.

---

## Risks & notes

- **Reference-tissue sparsity** — placental/trophoblast epigenomes are thinner than
  blood; state coverage must be reported and a nearest-tissue fallback flagged.
- **Circularity in B3** — independent reference panels for annotation vs. evaluation
  (see B3).
- **Build/coordinate errors** — all data is hg38-native here (no liftOver). The B1
  harness also surfaces systematic coordinate bugs as p-value inflation.

## References

- ENCODE / Roadmap ChromHMM chromatin-state models; ENCODE histone ChIP & DNase.
- Snyder et al. *Cell* 2016 — cfDNA tissues-of-origin (motivates B3).
- Loyfer et al. / Moss et al. — methylation atlases (orthogonal truth for B4).
- Kim et al. 2015 (PMID 25967380) — seqFF (the fetal-fraction model this interprets).
