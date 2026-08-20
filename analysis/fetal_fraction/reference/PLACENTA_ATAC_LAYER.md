# Placental snATAC openness layer

A **separate** open-chromatin track derived from single-nucleus ATAC-seq of the
human placenta, kept side by side with the DNase openness atlas rather than merged
into it. The two are different assays (snATAC vs bulk DNase-seq); holding them apart
keeps cross-modal agreement auditable and lets the tissue-of-origin step corroborate
a placental signal with two independent measurements.

## File
`ff_placenta_atac_hg19_50kb.csv.gz` — 57,633 rows (same 50 kb hg19 autosomal grid and
`key` format as `ff_openness_atlas_hg19_50kb*.csv.gz`), columns:

| column | meaning |
|--------|---------|
| `chrom,start,end,key` | 50 kb bin, identical to the DNase atlas (`chr1:0-50000`) |
| `A_STB` | syncytiotrophoblast openness (z across bins) |
| `A_vCTB` | villous cytotrophoblast |
| `A_EVT` | extravillous trophoblast |
| `A_Endothelial` | placental endothelial |
| `A_Erythroid` | fetal erythroid |
| `A_Fibroblast` | placental fibroblast |
| `placenta_atac_spec` | max(A_*) − mean(rest): a placenta-cell-type specificity score |

`A_` columns are z-scored across bins per cell type (same convention as the DNase
`C_` columns), so correlations and PC alignments are on the same footing.

## Source
Wang M, Zhu H, et al. *"A single-cell atlas of chromatin accessibility in the human
placenta."* Nat Genet 56:294–305 (2024). DOI 10.1038/s41588-023-01647-w; PubMed
38267607; GEO GSE247036 / BioProject PRJNA1035954. 12 donors; gestational stages
early (GW8) and term (GW38–39); 6 cell types; 11 pseudobulk profiles (EVT early only).

## hg38 → hg19 liftover
The published atlas is **hg38** (`provenance.json`: `genome = GRCh38/hg38`,
`n_peaks = 283847`), with prevalence ≥10 and the ENCODE **hg38** blacklist v2 already
applied at source. It was lifted to hg19 to match the FF model grid:

1. **Chain:** UCSC `hg38ToHg19.over.chain`, tool UCSC `liftOver` at default settings.
2. **Loss:** 283,847 hg38 peaks → **283,406** hg19 peaks retained (441 dropped, 0.16%
   — peaks with no unique hg19 image or that split across mappings).

The resulting hg19 peak BED + per-profile count matrix are what
`build_placenta_atac_atlas.py` consumes (packaged in the source `*_hg19.h5`, which
carries the 283,406 lifted peaks × 11 pseudobulk profiles).

> Note: the source `provenance.json` documents the hg38 atlas; the liftover step
> itself was performed during atlas prep and is not re-recorded in that file. The
> retained-peak count (283,406) is verified directly against the shipped h5.

## How the atlas is built
`analysis/fetal_fraction/build_placenta_atac_atlas.py`:
- pools each cell type's early+term profiles, weighted by nucleus count (`n_cells`),
  into one CPM vector per cell type;
- sums peak CPM into each overlapping 50 kb bin (peaks ≪ 50 kb, almost always one bin);
- z-scores each column across bins;
- writes the csv.gz + `manifest_placenta_atac.json`.

```
python analysis/fetal_fraction/build_placenta_atac_atlas.py \
    --h5 <Wang2024..._hg19.h5> \
    --grid analysis/fetal_fraction/reference/ff_openness_atlas_hg19_50kb_curated73.csv.gz \
    --out  analysis/fetal_fraction/reference/ff_placenta_atac_hg19_50kb.csv.gz
```

## Cross-modal validation (measured)
Each ATAC cell type vs DNase `C_placenta`, Spearman over 57,633 bins:
`A_vCTB` 0.93, `A_Erythroid` 0.93, `A_STB` 0.92, `A_Fibroblast` 0.92,
`A_Endothelial` 0.92, `A_EVT` 0.91. High agreement → the ATAC layer measures open
chromatin in the same regions as the DNase placenta column, at cell-type resolution.

## Using it in the PC / tissue-of-origin step
`pc_tissue_map.R` is atlas-agnostic — point it at this file with `--prefix A_`:

```
Rscript analysis/fetal_fraction/pc_tissue_map.R \
    --model ff_model.rds \
    --atlas analysis/fetal_fraction/reference/ff_placenta_atac_hg19_50kb.csv.gz \
    --prefix A_ --topn 2000 --outdir pc_atac_out
```

Output `pc_tissue_corr.csv` then has `corr_STB`, `corr_vCTB`, … columns: the
β-oriented Pearson correlation of each of the model's 100 PC loading vectors with
each placental cell type's ATAC openness. Run it once with `--prefix C_` (DNase) and
once with `--prefix A_` (ATAC); PCs that align with placenta in BOTH are the
cross-modally supported fetal-origin components.
