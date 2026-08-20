# FF openness atlas — build log & provenance

**File:** `ff_openness_atlas_hg19_50kb_curated73.csv.gz`
**Assembly:** hg19 (GRCh37) · **Bins:** 57,633 fixed-width 50 kb, autosomes only (chr1-chr22)
**Shape:** 57,633 rows x 77 cols = `chrom,start,end,key` + 73 `C_<tissue>` openness columns
**Source:** 100% public ENCODE portal (DNase-seq narrowPeak). No S3, no AWS, no credentials.
**Rebuild:** `python scripts/fetch_ff_atlas_beds.py` (fetch BEDs) then `build_openness_atlas.py` (bin + z-score)

---

## 1. Data source

Every source track is a **released ENCODE DNase-seq narrowPeak** file on the hg19
assembly, pulled from the public portal:

```
https://www.encodeproject.org/search/?type=Experiment&assay_title=DNase-seq&assembly=hg19&status=released&limit=all&format=json
    &field=accession
    &field=biosample_ontology.term_name
    &field=biosample_ontology.classification
    &field=files.href, files.file_type, files.file_format,
     files.assembly, files.status, files.output_type
```

Individual files are downloaded from
`https://www.encodeproject.org/files/<ENCFF...>/@@download/<ENCFF...>.bed.gz`.

The query returned **579 hg19 DNase-seq experiments**. File-level filter kept only:
`file_format == "bed"`, `"narrowPeak" in file_type`, `assembly == "hg19"`,
`status == "released"`, and `output_type` in {optimal IDR thresholded peaks, stable
peaks, pseudoreplicated peaks, peaks}. In practice **all** kept hg19 DNase
narrowPeaks are `output_type == "peaks"`.

## 2. Panel curation (72 ENCODE biosamples + 1 local = 73 tracks)

1. **<=3 replicates per biosample.** Files are ranked by peak type -- optimal IDR
   thresholded (0) < stable (1) < pseudoreplicated (2) < peaks (3) -- then, because
   all hg19 DNase files tie at "peaks", **broken deterministically by accession**
   (ascending). Top <=3 kept. *(The accession tie-break was added in v0.6.1; before
   it, re-fetches could pick different reps for biosamples with >3 files.)*
2. **Keep every `classification == "tissue"` biosample -> 63 tissues.**
3. **Plus a fixed ADDITIONS list of 9 non-"tissue" biosamples** deliberately kept
   for cfDNA relevance (immune lineages, a myeloid progenitor, and placenta/vascular
   cell types the tissue set misses).
4. **+1 neutrophil track from LOCAL files** (ENCODE has no hg19 neutrophil DNase-seq).

Total = 63 + 9 + 1 = **73 openness columns**. Panel = **164 ENCODE peak files**.

### 2a. The 63 curated tissues

| biosample | n_rep | file accessions |
|-----------|-------|-----------------|
| Peyer's patch | 1 | ENCFF055ZSF |
| adrenal gland | 3 | ENCFF094AUJ, ENCFF179WUI, ENCFF215UXZ |
| brain | 3 | ENCFF103UPI, ENCFF205KSG, ENCFF207IYJ |
| cerebellum | 2 | ENCFF455ZGC, ENCFF952ETG |
| chorion | 1 | ENCFF772FRH |
| coronary artery | 2 | ENCFF021ZDC, ENCFF739ASH |
| esophagus muscularis mucosa | 1 | ENCFF747AGO |
| esophagus squamous epithelium | 1 | ENCFF197ZNE |
| eye | 3 | ENCFF600YFF, ENCFF720MLF, ENCFF984RUT |
| femur | 1 | ENCFF760NDX |
| forelimb muscle | 1 | ENCFF850STV |
| frontal cortex | 3 | ENCFF101NVD, ENCFF144EGT, ENCFF281ASO |
| germinal center | 1 | ENCFF025VVW |
| heart | 3 | ENCFF007RSQ, ENCFF053SED, ENCFF086KUU |
| heart left ventricle | 3 | ENCFF224GQJ, ENCFF409QXF, ENCFF520HKD |
| heart right ventricle | 2 | ENCFF250FRA, ENCFF687FIQ |
| hindlimb muscle | 1 | ENCFF903ULR |
| kidney | 3 | ENCFF050XEY, ENCFF056MTI, ENCFF136JFN |
| large intestine | 3 | ENCFF058ACJ, ENCFF196VQQ, ENCFF276FKY |
| left cardiac atrium | 1 | ENCFF068ZPM |
| left forelimb | 1 | ENCFF459ENG |
| left hindlimb | 1 | ENCFF793VFP |
| left kidney | 3 | ENCFF311ACF, ENCFF360YZY, ENCFF476LDM |
| left lung | 3 | ENCFF059HII, ENCFF130KKE, ENCFF185MBE |
| left renal cortex interstitium | 3 | ENCFF206SZY, ENCFF307XDO, ENCFF531SWI |
| left renal pelvis | 3 | ENCFF132BIV, ENCFF276XHO, ENCFF352EVQ |
| liver | 3 | ENCFF286LYP, ENCFF547IOY, ENCFF572VZQ |
| lung | 3 | ENCFF169QJZ, ENCFF230CHA, ENCFF243LFH |
| muscle of arm | 3 | ENCFF099MBG, ENCFF172LKJ, ENCFF282BQK |
| muscle of back | 3 | ENCFF081UUE, ENCFF171JLD, ENCFF191JLB |
| muscle of leg | 3 | ENCFF020LII, ENCFF230WVU, ENCFF319CPJ |
| muscle of trunk | 3 | ENCFF448DSH, ENCFF562NMD, ENCFF836ETB |
| occipital lobe | 1 | ENCFF121NZP |
| omental fat pad | 1 | ENCFF925YKV |
| ovary | 2 | ENCFF498QOJ, ENCFF883WWT |
| pancreas | 2 | ENCFF735JRV, ENCFF897PRD |
| placenta | 3 | ENCFF005DMK, ENCFF029JFN, ENCFF128BFZ |
| psoas muscle | 3 | ENCFF445UFF, ENCFF523BCF, ENCFF654CVF |
| renal cortex interstitium | 3 | ENCFF027SVR, ENCFF051UGM, ENCFF131FBL |
| renal pelvis | 3 | ENCFF012SAA, ENCFF041XBL, ENCFF293PGI |
| retina | 3 | ENCFF009WLN, ENCFF015ULH, ENCFF402DMO |
| right atrium auricular region | 1 | ENCFF918YIN |
| right forelimb | 1 | ENCFF463KFI |
| right hindlimb | 1 | ENCFF369GLM |
| right kidney | 3 | ENCFF119HGJ, ENCFF186VEE, ENCFF216UBC |
| right lung | 3 | ENCFF036PGY, ENCFF046QJY, ENCFF104JWY |
| right renal cortex interstitium | 3 | ENCFF108GEI, ENCFF251NTZ, ENCFF585XKV |
| right renal pelvis | 3 | ENCFF033GXT, ENCFF360OCL, ENCFF374BCD |
| sigmoid colon | 3 | ENCFF006EIZ, ENCFF862GVN, ENCFF882CST |
| skin of body | 1 | ENCFF203VLF |
| small intestine | 3 | ENCFF007XOV, ENCFF063OCF, ENCFF321LTV |
| spinal cord | 3 | ENCFF038VGE, ENCFF162HDT, ENCFF753RBC |
| spleen | 1 | ENCFF587YNA |
| stomach | 3 | ENCFF087QAV, ENCFF104PLB, ENCFF122TCQ |
| testis | 3 | ENCFF012QTD, ENCFF833PXP, ENCFF843ZSC |
| thymus | 3 | ENCFF061UBG, ENCFF189LKP, ENCFF601PAY |
| tibial artery | 1 | ENCFF378BUF |
| tibial nerve | 1 | ENCFF424BYY |
| tongue | 3 | ENCFF173CEG, ENCFF259KDP, ENCFF814UZW |
| transverse colon | 1 | ENCFF142WAA |
| upper lobe of left lung | 1 | ENCFF271JAF |
| urinary bladder | 1 | ENCFF391VLP |
| vagina | 2 | ENCFF538ASZ, ENCFF600WRF |

### 2b. The 9 ADDITIONS (non-"tissue", kept for cfDNA relevance)

| biosample | ENCODE classification | n_rep | file accessions |
|-----------|-----------------------|-------|-----------------|
| CD14-positive monocyte | primary cell | 3 | ENCFF063IUG, ENCFF289XSM, ENCFF363XGK |
| B cell | primary cell | 3 | ENCFF159HJV, ENCFF263YRO, ENCFF403PJS |
| CD4-positive, alpha-beta T cell | primary cell | 3 | ENCFF023ZMS, ENCFF109RSW, ENCFF191WFG |
| CD8-positive, alpha-beta T cell | primary cell | 3 | ENCFF109VIJ, ENCFF512IML, ENCFF566SXQ |
| natural killer cell | primary cell | 3 | ENCFF133YNT, ENCFF675NXI, ENCFF840ZTK |
| common myeloid progenitor, CD34-positive | primary cell | 3 | ENCFF024DAR, ENCFF059VUS, ENCFF110CYC |
| trophoblast cell | in vitro differentiated cells | 3 | ENCFF520MND, ENCFF533MFM, ENCFF855CMV |
| endothelial cell of umbilical vein | primary cell | 3 | ENCFF097KBE, ENCFF330HEK, ENCFF382IRS |
| dermis blood vessel endothelial cell | primary cell | 3 | ENCFF045UNP, ENCFF084LGH, ENCFF419NYS |

Grouping: **immune** = CD14+ monocyte, B cell, CD4+ ab T cell, CD8+ ab T cell,
natural killer cell; **myeloid progenitor** = common myeloid progenitor CD34+;
**cfDNA-relevant** = trophoblast cell, endothelial cell of umbilical vein, dermis
blood vessel endothelial cell.

### 2c. The 73rd track -- neutrophil (local, not ENCODE)

ENCODE has **no hg19 (or GRCh38) neutrophil DNase-seq**, so neutrophil -- often the
single largest cfDNA contributor -- is built from 3 user-supplied hg19 accessibility
BEDs in `data/neutrophil_dnase/` (`add_neutrophil_dnase.py`):
`Neutrophil_fetalung.bed.gz`, `GMPNeutrophil_CD34_bonemarrowPBMC.bed.gz`,
`GMPNeutrophil_bonemarrowPBMC.bed.gz`. Same aggregation as the ENCODE tracks.
This is the one part of the atlas not reproducible from the public portal alone.

## 3. Binning: narrowPeak -> 50 kb openness

For every biosample track:

1. **Grid.** Fixed 50 kb autosomal bins from UCSC hg19 chromosome lengths
   (`chr1` = 249,250,621 ... `chr22` = 51,304,566), `start = k*50000`,
   `end = min(start+50000, chrom_len)`. -> **57,633 bins**, key `chrN:start-end`.
2. **Autosomes only.** chrX/chrY/chrM peaks dropped on load.
3. **signalValue-weighted covered fraction.** Each peak is spread across the bins it
   overlaps; a bin accumulates `sum(overlap_bp) / 50000` (clipped to [0,1]),
   signalValue-weighted so stronger peaks contribute more (fold-enrichment-weighted
   openness). This is the per-replicate `D_<tissue>` accessibility coverage.
4. **z-score per replicate, then mean across replicates** (`zmean`):
   `z = (x - mean(x)) / std(x)` over the 57,633 bins, averaged over the <=3 reps.
5. **Composite `C_<tissue>`** (the column the tissue-of-origin track reads,
   `prefix="C_"`): default `composite="histone+dnase"` = z-mean of z(accessibility)
   and z(histone). For the **DNase-only** curated-73 atlas, `C_<tissue>` is the
   z-scored accessibility track directly, so every tissue lands on one comparable
   z-scale regardless of assay availability.

## 4. Configurable openness (v0.6.0+)

`build_openness_atlas.py` exposes the modality choice so you are not locked to a
fixed blend (see the reference README section "Choosing the openness modality"):

- `--composite dnase` -- accessibility only (DNase). **Default for this atlas.**
- `--composite histone` -- selected active histone marks only.
- `--composite histone+dnase` -- z-mean of both.
- `--marks H3K4me3[,H3K27ac,...]` -- which of the 4 active marks
  (H3K4me3, H3K4me1, H3K27ac, H3K36me3) feed the histone track `H_<tissue>`.

## 5. Provenance summary

| property | value |
|----------|-------|
| assembly | hg19 / GRCh37 |
| assay | ENCODE DNase-seq narrowPeak (public portal) |
| experiments queried | 579 (hg19 DNase, released) |
| biosamples in panel | 72 ENCODE (63 tissue + 9 additions) + 1 local neutrophil = 73 |
| peak files | 164 ENCODE + 3 local neutrophil |
| max reps / biosample | 3 (rank by peak type, tie-break by accession) |
| bins | 57,633 x 50 kb, autosomes only |
| openness | signalValue-weighted covered fraction -> z-mean per replicate |
| S3 / AWS required | **No** -- fully public + one local track |

The complete accession list is `analysis/fetal_fraction/reference/ff_atlas_accessions.csv`
and is regenerated by `python scripts/fetch_ff_atlas_beds.py --dry-run`.
