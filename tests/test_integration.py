"""Integration test against the real bundled references.

Skipped automatically when the FASTA files are absent (they are fetched separately
by scripts/fetch_fasta.py and gitignored), so CI without the sequence still passes
on the unit tests.
"""
import json
import os

import pytest

from cfdna_chromatin import references as R, genome as G, engine as E, histone as H
from cfdna_chromatin import accessibility as A

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, "data")
FASTA = os.path.join(DATA, "fasta")
CHROMS = ["chr19"]

_have_fasta = all(os.path.exists(os.path.join(FASTA, f"{c}.fa.gz")) for c in CHROMS)
_have_beds = all(os.path.exists(os.path.join(DATA, f"{t}_{a}.bed.gz"))
                 for t, a in R.REFERENCE_ACCESSIONS.items())
_have_histone = os.path.exists(os.path.join(DATA, "histone", "manifest.json"))
_have_access = os.path.exists(os.path.join(DATA, "access", "manifest.json"))

pytestmark = pytest.mark.skipif(
    not (_have_fasta and _have_beds),
    reason="reference FASTA/BEDs not present; run scripts/fetch_fasta.py",
)


def test_housekeeping_promoters_enrich():
    seg = R.load_segmentation(
        os.path.join(DATA, f"placenta_{R.REFERENCE_ACCESSIONS['placenta']}.bed.gz"),
        chroms=CHROMS,
    )
    seg = {"chr19": seg["chr19"]}
    genome = G.load_genome({"chr19": os.path.join(FASTA, "chr19.fa.gz")})
    hk = json.load(open(os.path.join(HERE, "examples", "housekeeping_tss_chr19.json")))
    prom = [(r["chrom"], r["tss"] - 2000, r["tss"] + 2000) for r in hk if r["chrom"] == "chr19"]

    nm = E.NullModel(genome=genome, chroms=["chr19"], seed=1)
    res, _ = E.enrichment_test(prom, seg, nm, by="group", n_per=50, n_boot=300, seed=1)
    prom_res = [r for r in res if r["label"] == "Promoter"][0]
    # housekeeping promoters must be strongly Promoter-enriched vs matched null
    assert prom_res["log2_fold"] > 1.5
    assert prom_res["p_emp"] < 0.05


def test_housekeeping_promoters_histone_fingerprint():
    """Orthogonal cross-check: mark layer must give the active-promoter signature."""
    genome = G.load_genome({"chr19": os.path.join(FASTA, "chr19.fa.gz")})
    panel = H.load_mark_panel(os.path.join(DATA, "histone"), chroms=["chr19"])
    hk = json.load(open(os.path.join(HERE, "examples", "housekeeping_tss_chr19.json")))
    prom = [(r["chrom"], r["tss"] - 2000, r["tss"] + 2000) for r in hk if r["chrom"] == "chr19"]

    if not _have_histone:
        pytest.skip("histone panel not present")
    nm = E.NullModel(genome=genome, chroms=["chr19"], seed=1)
    res, _ = H.mark_enrichment_test(prom, panel["placenta"], nm, n_per=50, n_boot=300, seed=1)
    by = {r["mark"]: r for r in res}
    # active promoter marks up, repressive down
    assert by["H3K4me3"]["log2_fold"] > 1.0
    assert by["H3K27ac"]["log2_fold"] > 1.0
    assert by["H3K27me3"]["log2_fold"] < 0
    assert by["H3K9me3"]["log2_fold"] < 0


def test_housekeeping_promoters_accessible():
    """Orthogonal cross-check: housekeeping promoters must sit in open chromatin."""
    genome = G.load_genome({"chr19": os.path.join(FASTA, "chr19.fa.gz")})
    hk = json.load(open(os.path.join(HERE, "examples", "housekeeping_tss_chr19.json")))
    prom = [(r["chrom"], r["tss"] - 2000, r["tss"] + 2000) for r in hk if r["chrom"] == "chr19"]

    if not _have_access:
        pytest.skip("accessibility panel not present")
    panel = A.load_access_panel(os.path.join(DATA, "access"), chroms=["chr19"])
    nm = E.NullModel(genome=genome, chroms=["chr19"], seed=1)
    res, _ = A.access_enrichment_test(prom, panel["placenta"], nm, n_per=50, n_boot=300, seed=1)
    # housekeeping promoters are nucleosome-depleted -> strongly DNase-enriched
    assert res["log2_fold"] > 1.0
    assert res["p_emp"] < 0.05
