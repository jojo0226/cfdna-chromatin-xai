"""Unit tests for liftover.py (no network: a fake lifter stands in for pyliftover)."""
import pandas as pd

from cfdna_chromatin import liftover as LO


class FakeLifter:
    """Deterministic point lifter: shifts by +1000 on same chrom, unless the
    coordinate is in `unmapped` (returns []) or in `crosschrom` (maps to chrZ)."""
    def __init__(self, shift=1000, unmapped=(), crosschrom=()):
        self.shift = shift
        self.unmapped = set(unmapped)
        self.crosschrom = set(crosschrom)

    def convert_coordinate(self, chrom, pos):
        if pos in self.unmapped:
            return []
        if pos in self.crosschrom:
            return [("chrZ", pos + self.shift, "+", 0)]
        return [(chrom, pos + self.shift, "+", 0)]


def test_lift_regions_shift_and_preserve_columns():
    df = pd.DataFrame({"chrom": ["chr1", "chr2"], "start": [100, 200],
                       "end": [600, 700], "importance": [0.9, 0.4]})
    fwd = FakeLifter(shift=1000)
    lifted, qc = LO.lift_regions(df, lifter=fwd, round_trip=False)
    assert qc["n_mapped"] == 2 and qc["mapping_rate"] == 1.0
    assert list(lifted["start"]) == [1100, 1200]
    assert list(lifted["end"]) == [1600, 1700]
    assert "importance" in lifted.columns          # extra column preserved


def test_lift_regions_drops_unmapped_and_diff_chrom():
    df = pd.DataFrame({"chrom": ["chr1", "chr1", "chr1"],
                       "start": [100, 300, 500], "end": [600, 800, 1000]})
    # region0 start unmapped; region1 start crosses to chrZ; region2 clean
    fwd = FakeLifter(unmapped=(100,), crosschrom=(300,))
    lifted, qc = LO.lift_regions(df, lifter=fwd, round_trip=False)
    assert qc["dropped_unmapped"] == 1
    assert qc["dropped_diff_chrom"] == 1
    assert qc["n_mapped"] == 1
    assert list(lifted["start"]) == [1500]


def test_lift_regions_list_input():
    regions = [("chr3", 100, 600), ("chr3", 700, 1200)]
    lifted, qc = LO.lift_regions(regions, lifter=FakeLifter(), round_trip=False)
    assert qc["n_mapped"] == 2
    assert lifted.iloc[0]["chrom"] == "chr3"
