# Fetal-fraction tissue-of-origin openness atlas

Precomputed per-bin chromatin openness for every cell type in the fetal-fraction
(seqFF++) panel, so the tissue-of-origin analysis runs without re-deriving openness
from raw ENCODE peaks each time. **Two openness atlases ship here** — pick by panel
breadth:

| file | panel | use |
|------|-------|-----|
| **`ff_openness_atlas_hg19_50kb_curated73.csv.gz`** | **73 tracks** (72 ENCODE + neutrophil) | **default** — the broad tissue-of-origin panel |
| `ff_openness_atlas_hg19_50kb.csv.gz` | 11 tracks | the original focused seqFF++ panel |

Both share the same grid:

- **Rows**: 57,633 bins, all 22 autosomes, 50 kb (hg19).
- **`key`**: `chrN:start-end` join key (also `chrom`,`start`,`end`).
- **Assembly**: hg19. The FF cohort is build-native — importance keys must be hg19
  50 kb bins; no liftover of cohort data (only reference peaks ever cross builds).

**Full build log + every accession:** see
[`ATLAS_BUILD_LOG.md`](ATLAS_BUILD_LOG.md) and
[`ff_atlas_accessions.csv`](ff_atlas_accessions.csv) (164 ENCODE files, one row
per replicate). Provenance in one line: **100% public ENCODE DNase-seq** (hg19
narrowPeak query → ≤3 reps/biosample → 63 tissues + 9 cfDNA-relevant additions),
plus **1 local neutrophil track** ENCODE can't supply. No S3, no AWS.

## The 73-tissue panel (default)

Placental (fetal) **signal** against a maternal background of solid tissue,
vasculature, and hematopoietic cell types — 63 curated ENCODE `tissue` biosamples
plus a fixed set of 9 cfDNA-relevant additions and 1 local neutrophil track:

| Role | Tracks |
|------|--------|
| Fetal / placental signal | `placenta`, `trophoblast cell`, `chorion` |
| Vascular (cfDNA-relevant) | `endothelial cell of umbilical vein`, `dermis blood vessel endothelial cell`, `coronary artery`, `tibial artery` |
| Maternal hematopoietic (mature) | `CD14-positive monocyte`, `B cell`, `CD4-positive, alpha-beta T cell`, `CD8-positive, alpha-beta T cell`, `natural killer cell`, `neutrophil` (local) |
| Maternal hematopoietic (progenitor) | `common myeloid progenitor, CD34-positive` |
| Solid tissue | 56 ENCODE tissues — brain/cerebellum/cortex, heart chambers, kidney/renal, lung lobes, GI tract, liver, muscle groups, endocrine, reproductive, etc. |

`placenta` ranks #1 of 73 by placental specificity, `trophoblast cell` #2. Placenta
is 3 reps; the full per-tissue accession + replicate list is in `ATLAS_BUILD_LOG.md`
§2. **neutrophil** (often the single largest cfDNA contributor) is the one non-ENCODE
track — ENCODE has no hg19 neutrophil DNase-seq, so it is built from 3 local hg19
accessibility BEDs (`data/neutrophil_dnase/`, `add_neutrophil_dnase.py`).

## The 11-tissue panel (focused)

The original seqFF++ panel — placenta signal against a compact maternal background.
There is **no tumor proxy** — this is the seqFF++ fork, not the cancer application.

| Role | Tissues |
|------|---------|
| Fetal signal (cell-of-origin) | `placenta` |
| Maternal solid | `liver`, `endothelial` |
| Maternal hematopoietic (mature) | `monocyte`, `Bcell`, `CD4`, `CD8`, `NK`, `neutrophil` |
| Maternal hematopoietic (progenitor) | `CMP` (CD34+ common myeloid progenitor) |
| Control | `K562` |

## Choosing the openness modality (DNase vs histone marks)

`build_openness_atlas.py` lets you pick what the combined-openness track
`C_<tissue>` is built from — you are **not** locked to a fixed histone+DNase blend:

| flag | `C_<tissue>` built from | when to use |
|------|-------------------------|-------------|
| `--composite dnase` | accessibility only (DNase) | **cleanest open-chromatin signal**; the default for the curated-73 atlas |
| `--composite histone` | selected active marks only | when you only trust histone ChIP for a panel |
| `--composite histone+dnase` | z-mean of both (default) | maximal coverage when DNase is sparse |
| `--marks H3K4me3[,H3K27ac,…]` | restrict which marks feed `H_<tissue>` | drop marks that don't mark *open* chromatin |

**On H3K4me1** — your concern is right. Of the four active marks
(H3K4me3, H3K4me1, H3K27ac, H3K36me3), **H3K4me1 and H3K36me3 do not cleanly mark
open/nucleosome-depleted chromatin**: H3K4me1 marks *primed/poised* enhancers
(often closed), and H3K36me3 marks transcribed gene bodies. If you want openness in
the strict accessibility sense, prefer **`--composite dnase`** (what the shipped
curated-73 atlas uses), or if you must use histone, restrict to the two
active-promoter/enhancer marks with **`--marks H3K4me3,H3K27ac`**. Avoid feeding
H3K4me1 into an "openness" track unless you specifically want poised-enhancer signal.

## Per-tissue coverage caveats (11-tissue atlas)

- **neutrophil is histone-only** — ENCODE has no neutrophil DNase-seq in hg19 or
  GRCh38, so `C_neutrophil = z(H_neutrophil)` (all 4 active marks; no `D_`/`ACC_`).
  *(The curated-73 atlas instead uses the 3 local neutrophil DNase BEDs.)*
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

## ATAC-seq: no ENCODE panel, but a dedicated placental layer

ENCODE has **zero** ATAC narrowPeak (hg19 or GRCh38) for placenta, endothelial, and
monocyte, so ATAC cannot enter the pan-tissue `C_` panel comparably — the panel stays
DNase + histone (recorded in `ff_reference_manifest.json`). The gap that mattered most
was placenta, and that is now covered by a **separate placental snATAC layer** built
from the Wang et al. 2024 pregnancy atlas (`ff_placenta_atac_hg19_50kb.csv.gz`; see
below and `PLACENTA_ATAC_LAYER.md`). It is kept as its own file — never merged into
`C_` — so it serves as an independent cross-modal check on the placental signal, not a
non-comparable panel column.

## Provenance

Accessions and per-tissue replicate lists are in `ff_reference_manifest.json`
(DNase accessions per tissue + `histone_accessions`). Raw ENCODE narrowPeak files
are pulled and converted to 4-column (chrom,start,end,signalValue) form by
`scripts/fetch_ff_reference.py` (DNase) and `scripts/fetch_histone.py`; the raw and
intermediate peak dirs are gitignored (regenerable). Only this derived atlas +
manifest ship.

## Additional atlases shipped here

- **`ff_openness_atlas_hg19_50kb_curated73.csv.gz`** (+ `manifest_curated73.json`) —
  the expanded openness panel: 72 curated **tissues** (primary cell lines, in-vitro
  differentiated, and cancer lines removed), with the blood lineages re-added
  (monocyte, B/CD4/CD8/NK) and the placenta family broadened (trophoblast ×3,
  chorion, HTR-8/SVneo). Same 57,633-bin hg19 grid and `C_<tissue>` openness
  columns. Placenta ranks #1/72, trophoblast #2. This is the panel Step 5 pairs
  against the repressive axis.
- **`ff_repressive_atlas_hg19_50kb_H3K27me3.csv.gz`** (+ `manifest_repressive.json`) —
  the **separate repressive axis** (see FF guide Step 5): z-mean ENCODE hg19
  H3K27me3 narrowPeak coverage `R_<tissue>` and pan-tissue-subtracted specificity
  `repr_spec_<tissue>` for 39 tissues (every cfDNA-critical one covered). Kept as
  its own file — never merged into the `C_` openness columns — so a genuine
  cell-of-origin can be required to be **active-enriched and repressive-depleted**.
- **`ff_placenta_atac_hg19_50kb.csv.gz`** (+ `manifest_placenta_atac.json`,
  `PLACENTA_ATAC_LAYER.md`) — the **placental snATAC openness layer** from Wang et al.
  2024 (Nat Genet, GSE247036; hg38→hg19 lifted, 283,406 peaks). Same 57,633-bin hg19
  grid; six cell-type columns `A_STB / A_vCTB / A_EVT / A_Endothelial / A_Erythroid /
  A_Fibroblast` (z across bins) plus `placenta_atac_spec`. Built by
  `analysis/fetal_fraction/build_placenta_atac_atlas.py`. Kept separate from the `C_`
  DNase columns; cross-modal agreement with `C_placenta` is Spearman 0.91–0.93 per cell
  type. Feed it to `pc_tissue_map.R --prefix A_` to correlate the model's PCs against
  placental chromatin accessibility (see FF guide / `PLACENTA_ATAC_LAYER.md`).

## Re-fetching the source BEDs (public ENCODE, no S3)

The shipped atlas CSVs are self-contained — the tissue-of-origin track reads them
directly and never needs the raw peaks. You only need the source BEDs to **rebuild**
an atlas (e.g. with different `--marks`/`--composite`), and they are **entirely
public ENCODE** — no S3 bucket, no AWS credentials, no archive restore:

```bash
python scripts/fetch_ff_atlas_beds.py --dry-run   # preview the 72-biosample panel
python scripts/fetch_ff_atlas_beds.py             # download BEDs + manifest -> data/access_ff/
```

It re-runs the same hg19 DNase-seq narrowPeak query the atlas was built from,
applies the same ≤3-replicate-per-biosample ranking and the same curation
(classification==`tissue` + the fixed ADDITIONS list — 5 immune, 1 myeloid
progenitor, 3 cfDNA-relevant biosamples), converts each narrowPeak to the shipped
4-column form (`chrom,start,end,signalValue`, hg19 autosomes), and writes
`encode_accessions.json` keyed on the **raw biosample name** so it joins directly
to the `C_<tissue>` columns here.

The panel is **72 ENCODE biosamples → 164 peak files**. The 73rd atlas column,
`neutrophil`, is the one gap: ENCODE has no hg19 neutrophil DNase-seq, so it is
supplied from local BEDs (the fetch script reports this explicitly).

> `scripts/restore_ff_atlas_beds.py` is a **deprecated tombstone** — an earlier
> release wrongly assumed the atlas accessions survived only in a private
> archived-S3 index and needed an IAM-authenticated restore. They don't; use the
> fetch script above.

## Running the tissue-of-origin track

See `../FF_IMPLEMENTATION.md` for the full AWS recipe: extract per-bin SHAP from
the PCA + linear-regression FF model **in place** (`extract_ff_shap.R`,
back-projecting `w_bin = Σ_k L[bin,k]·β_k`; no per-patient data leaves the
instance), map the model's bin names onto these hg19 50 kb keys, then run
`run_ff_tissue_track.py` against this atlas.
