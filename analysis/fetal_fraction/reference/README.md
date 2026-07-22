# Fetal-fraction tissue-of-origin openness atlas

`ff_openness_atlas_hg19_50kb.csv.gz` — precomputed per-bin chromatin openness for
every cell type in the fetal-fraction (seqFF++) panel, so the tissue-of-origin
analysis runs without re-deriving openness from raw ENCODE peaks each time.

- **Rows**: 57,633 bins, all 22 autosomes, 50 kb (hg19).
- **`key`**: `chrN:start-end` join key (also `chrom`,`start`,`end`).
- **Assembly**: hg19. The FF cohort is build-native — importance keys must be hg19
  50 kb bins; no liftover of cohort data (only reference peaks ever cross builds).

## Tissues (11)

Placental (fetal) **signal** against a maternal background of hematopoietic and
solid-tissue cell types. There is **no tumor proxy** — this is the seqFF++ fork,
not the cancer application.

| Role | Tissues |
|------|---------|
| Fetal signal (cell-of-origin) | `placenta` |
| Maternal solid | `liver`, `endothelial` |
| Maternal hematopoietic (mature) | `monocyte`, `Bcell`, `CD4`, `CD8`, `NK`, `neutrophil` |
| Maternal hematopoietic (progenitor) | `CMP` (CD34+ common myeloid progenitor) |
| Control | `K562` |

`neutrophil` (granulocyte lineage, often the single largest cfDNA contributor) and
`CMP` (the shared erythroid/megakaryocyte/granulocyte progenitor axis) close the two
largest gaps in a mature-cell-only panel. Two coverage caveats, both handled by the
`C_` composite:

- **neutrophil is histone-only** — ENCODE has no neutrophil DNase-seq in hg19 or
  GRCh38, so `C_neutrophil = z(H_neutrophil)` (all 4 active marks; no `D_`/`ACC_`).
- **CMP has no H3K27ac** narrowPeak in hg19 — `H_CMP` averages the other 3 active
  marks, and `C_CMP` uses its DNase-seq (3 reps) for the accessibility half.

## Columns

For each tissue `T`, four openness tracks:

- `D_T` — **DNase-seq** peak coverage (openness), z-mean over replicates, signalValue-weighted.
- `H_T` — mean coverage by the four **active histone marks** (H3K4me3, H3K4me1, H3K27ac, H3K36me3).
- `ACC_T` — accessibility composite (currently `z(D_T)`; an ATAC layer would go here but is excluded, see note).
- `C_T` — **combined openness** = mean of z(`H_T`) and z(`ACC_T`), comparable across tissues.
  Histone-only tissues (`neutrophil`) use z(`H_T`) alone; accessibility-only tissues use
  z(`ACC_T`) alone — so every tissue lands on the same scale regardless of which assays
  ENCODE has. This is the default track the tissue-proportion analysis uses (`prefix="C_"`).

Per-track replicate counts are in the `n_*` columns.

## Why no ATAC-seq

ENCODE has **zero** ATAC narrowPeak (hg19 or GRCh38) for placenta, endothelial, and
monocyte. Because placenta is the critical fetal signal, an ATAC layer would be
missing exactly where it matters and would be non-comparable across the panel, so
the atlas is DNase + histone only. (Recorded in `ff_reference_manifest.json`.)

## Provenance

Accessions and per-tissue replicate lists are in `ff_reference_manifest.json`
(DNase accessions per tissue + `histone_accessions`). Raw ENCODE narrowPeak files
are pulled and converted to 4-column (chrom,start,end,signalValue) form by
`scripts/fetch_ff_reference.py` (DNase) and `scripts/fetch_histone.py`; the raw and
intermediate peak dirs are gitignored (regenerable). Only this derived atlas +
manifest ship.

## Running the tissue-of-origin track

See `../FF_IMPLEMENTATION.md` for the full AWS recipe: extract per-feature |SHAP|
from a glmnet FF model **in place** (`extract_ff_shap.R`, no per-patient data
leaves the instance), map the mixed coverage-bin + 4-mer-motif feature space onto
these hg19 50 kb keys, then run `run_ff_tissue_track.py` against this atlas.
