from __future__ import annotations

import numpy as np
import pandas as pd

import pytest

from gatk_sv_ploidy.preprocess import (
    filter_low_quality_bins,
    filter_poor_region_bins,
    normalise_depth,
    read_bed_intervals,
    read_known_sites,
)


def test_read_bed_intervals_rejects_invalid_rows(tmp_path) -> None:
    bed = tmp_path / "poor_regions.bed"
    bed.write_text("chr21\t100\t50\n")

    with pytest.raises(ValueError, match="end <= start"):
        read_bed_intervals(str(bed))


def test_filter_poor_region_bins_removes_overlapped_bins() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chr21", "chr21", "chrX"],
            "Start": [0, 100, 0],
            "End": [100, 200, 100],
            "S1": [2.0, 2.1, 1.0],
        }
    )
    poor_regions = pd.DataFrame(
        {
            "Chr": ["chr21"],
            "Start": [25],
            "End": [175],
        }
    )

    filtered = filter_poor_region_bins(
        df,
        poor_regions,
        min_poor_region_coverage=0.5,
    )

    assert filtered[["Chr", "Start"]].values.tolist() == [["chrX", 0]]


def test_normalise_depth_scales_autosome_medians_to_two() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chr18", "chr21", "chrX"],
            "Start": [0, 100, 200],
            "End": [100, 200, 300],
            "S1": [8.0, 12.0, 4.0],
            "S2": [20.0, 20.0, 20.0],
        }
    )

    normalized = normalise_depth(df)

    np.testing.assert_allclose(
        np.median(normalized.loc[normalized["Chr"] != "chrX", ["S1", "S2"]], axis=0),
        np.array([2.0, 2.0]),
    )


def test_filter_low_quality_bins_removes_cohort_shifted_bin() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chr21"] * 4,
            "Start": [0, 100, 200, 300],
            "End": [100, 200, 300, 400],
            "S1": [2.0, 2.1, 2.0, 4.5],
            "S2": [2.0, 2.2, 2.1, 4.6],
            "S3": [2.1, 2.0, 1.9, 4.4],
        }
    )

    filtered = filter_low_quality_bins(
        df,
        cohort_deviation_threshold=0.3,
        cohort_deviation_fraction_max=0.5,
        min_bins_per_chr=1,
    )

    assert filtered["Start"].tolist() == [0, 100, 200]


def test_read_known_sites_reads_vcf_and_converts_to_zero_based(tmp_path) -> None:
    vcf = tmp_path / "sites.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr21\t501\t.\tA\tG\t.\tPASS\tAF=0.3\n"
    )

    sites = read_known_sites(str(vcf))

    assert sites.loc[0, "position"] == 500
    assert sites.loc[0, "ref"] == "A"
    assert sites.loc[0, "pop_af"] == pytest.approx(0.3)
