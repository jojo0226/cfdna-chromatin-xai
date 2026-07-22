# Data provenance

## ChromHMM reference segmentations (bundled in `data/`)

All four are the **Roadmap / EpiMap ChromHMM 18-state model** on **GRCh38 (hg38)**,
downloaded from the ENCODE portal (`https://www.encodeproject.org`). Output type:
`semi-automated genome annotation`. State vocabulary is identical across all four
(verified programmatically): TssA, TssFlnk, TssFlnkU, TssFlnkD, Tx, TxWk, EnhG1,
EnhG2, EnhA1, EnhA2, EnhWk, ZNF/Rpts, Het, TssBiv, EnhBiv, ReprPC, ReprPCWk, Quies.

| File | Accession | Biosample ID | Description |
|------|-----------|--------------|-------------|
| `placenta_ENCFF024IDF.bed.gz`   | ENCFF024IDF | BSS01443 | placenta male embryo (85 days) |
| `neutrophil_ENCFF412COT.bed.gz` | ENCFF412COT | BSS01381 | neutrophil male |
| `Bcell_ENCFF004HFC.bed.gz`      | ENCFF004HFC | BSS00095 | B cell female adult (27 years) |
| `K562_ENCFF026ZCC.bed.gz`       | ENCFF026ZCC | BSS01039 | K562 |

Format: BED9 (chrom, start, end, state, score, strand, thickStart, thickEnd, itemRgb),
0-based half-open. Only column 4 (state label) is used by the engine.

**Biosample selection rationale.** cfDNA in maternal plasma is a mixture of a
placental (fetal) component and a maternal hematopoietic background. The four
references bracket that mixture: placenta for the fetal side; neutrophil and B cell
for the two main maternal myeloid/lymphoid lineages; K562 as an erythroid reference
and cross-tissue control. Additional biosamples can be added by dropping their
ChromHMM 18-state BED into `data/` and extending `REFERENCE_ACCESSIONS`.

## Histone-mark ChIP-seq peaks (bundled in `data/histone/`)

The histone-mark fingerprint layer (`histone.py`) uses ENCODE **Histone ChIP-seq**
peak calls on **GRCh38**, one experiment per (tissue, mark) across the same four
biosamples, for six marks:

| Mark | Chromatin association |
|------|-----------------------|
| H3K4me3  | active promoters (sharp) |
| H3K4me1  | enhancers (poised + active) |
| H3K27ac  | active enhancers / promoters |
| H3K36me3 | transcribed gene bodies |
| H3K27me3 | Polycomb repression |
| H3K9me3  | constitutive heterochromatin |

The full 6 × 4 file matrix (accessions, parent experiments, output type, peak counts)
is recorded in `data/histone/manifest.json`. Peak calls are preferred as
`replicated peaks` where available, else `pseudoreplicated peaks`. Each shipped file
is a **chr19–22 subset** reduced to 4 columns (chrom, start, end, signalValue),
0-based half-open — ~2 MB total. The genome-wide originals (~34 MB) are
re-downloadable with `scripts/fetch_histone.py` into `data/histone_raw/` (gitignored).

This layer is a **state-independent** readout: because the marks are called
independently of the ChromHMM segmentation, agreement between the mark fingerprint
and the ChromHMM call is genuine corroboration, not circularity.

## Open-chromatin DNase-seq peaks (bundled in `data/access/`)

The accessibility layer (`accessibility.py`) uses ENCODE **DNase-seq** narrowPeak
calls on **GRCh38**, one experiment per tissue. Open chromatin is the most direct
footprint of regulatory activity and is what makes cfDNA cell-of-origin inference
possible, so it forms a third orthogonal axis alongside ChromHMM state and the
histone fingerprint.

| Tissue (key) | Biosample assayed | File | Experiment | chr19–22 peaks |
|--------------|-------------------|------|------------|----------------|
| placenta   | placenta | ENCFF316XWU | ENCSR035QHH | 32,136 |
| neutrophil | CD14-positive monocyte *(proxy)* | ENCFF136KXE | ENCSR000ELE | 21,885 |
| Bcell      | B cell | ENCFF251ZFA | ENCSR000EII | 21,065 |
| K562       | K562 | ENCFF218UJY | ENCSR000EKN | 67,991 |

**Neutrophil substitution.** Primary neutrophils have highly condensed chromatin
and are essentially absent from ENCODE accessibility assays (no DNase-seq or
ATAC-seq on GRCh38). The `neutrophil` track is therefore a **CD14-positive
monocyte** DNase-seq experiment used as a labeled myeloid proxy for the maternal
background; the substitution is recorded in the manifest (`proxy_for: neutrophil`)
so it is never silently conflated with the true-neutrophil ChromHMM and histone
references. The ChromHMM and histone layers still use genuine neutrophil data.

The manifest (`data/access/manifest.json`) records accessions, parent experiments,
output type, and peak counts. DNase experiments of this vintage emit a single
`peaks` narrowPeak output (the replicated-peaks convention postdates them). Each
shipped file is a **chr19–22 subset** reduced to 4 columns (chrom, start, end,
signalValue), 0-based half-open — ~1.3 MB total. Genome-wide originals are
re-downloadable with `scripts/fetch_access.py` into `data/access_raw/` (gitignored).

## Genome sequence (fetched, not bundled)

hg38 per-chromosome FASTA from UCSC goldenPath
(`https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/`), fetched by
`scripts/fetch_fasta.py` into `data/fasta/` (gitignored). Used only to compute GC
content and N-fraction (mappability proxy) for the matched-null background.

The demo works on **chr19–22** (gene-dense chr19 plus the small autosomes) to keep
the footprint small; the engine itself is chromosome-agnostic — fetch the full
genome to run genome-wide.

## Housekeeping-gene TSS (`examples/housekeeping_tss_chr19.json`)

Canonical TSS coordinates for 22 housekeeping / constitutively-expressed genes on
chr19, retrieved via Ensembl REST (GRCh38). Strand-aware TSS = gene start on +
strand, gene end on − strand. Used as the B2 positive control (their ±2 kb promoter
windows should annotate as active-promoter chromatin).
## Fetal-fraction tissue-of-origin atlas (`analysis/fetal_fraction/reference/`)

The seqFF++ fork ships a precomputed **tissue-of-origin openness atlas**
(`ff_openness_atlas_hg19_50kb.csv.gz`): 57,633 autosomal 50 kb bins (hg19) × per-tissue
openness tracks for the nine cell types that contribute to maternal-plasma cfDNA.
There is **no tumor proxy** — the placenta (fetal) signal is the cell-of-origin of
interest, scored against a maternal background of hematopoietic and solid-tissue cells.

| role | tissues | assays |
|------|---------|--------|
| fetal signal (cell-of-origin) | placenta | DNase + 4 active histone marks |
| maternal solid | liver, endothelial | DNase + histone |
| maternal hematopoietic | monocyte, Bcell, CD4, CD8, NK | DNase + histone |
| control | K562 | DNase + histone |

For each tissue the atlas carries a DNase openness track (`D_`), an active-histone
fingerprint (`H_`, mean over H3K4me3/H3K4me1/H3K27ac/H3K36me3), an accessibility
composite (`ACC_`), and a cross-assay combined-openness track (`C_`, z-mean of histone
and accessibility) that is comparable across tissues and is the default the
tissue-proportion analysis reads. Peaks are z-scored per replicate (depth-robust) then
averaged; binning weights each bin by the peak's signalValue (ENCODE narrowPeak col 7),
not the saturated integer score. Full accessions and per-tissue replicate lists are in
`ff_reference_manifest.json`.

**No ATAC-seq layer.** ENCODE has zero ATAC narrowPeak (hg19 or GRCh38) for placenta,
endothelial, and monocyte. Because placenta is the critical fetal signal, an ATAC layer
would be missing exactly where it matters and non-comparable across the panel, so the
atlas is DNase + histone only.

**Build-native rule.** The atlas is hg19 to match legacy seqFF cohorts; importance keys
fed to the tissue-of-origin track must be hg19 50 kb bins. Only reference peaks ever
cross genome builds — cohort data is never lifted over.

## SHAP selection front-end (v0.4.1)

`selection.py` is the importance/selection stage (Stage 1) of the three-stage
architecture. It is model-agnostic: `rank_by_shap()` consumes any SHAP matrix
(samples × genomic bins), `rank_by_differential()` computes a per-bin Mann-Whitney
ranking with no model, and `shap_from_model()` computes a region-anchored SHAP
attribution from a fitted model locally (tree → TreeExplainer, linear →
LinearExplainer, else KernelExplainer), so no cfDNA data has to leave the machine.
`compartment_importance_test()` correlates per-bin importance with signal-specific vs
background-specific chromatin openness (Spearman), reporting **both** the raw-openness
and the specific-openness correlations so the genome-wide openness confound is explicit:
importance often tracks total openness (any active chromatin), and only the
compartment-*specific* contrast reveals cell-of-origin.
