from __future__ import annotations

import numpy as np
import pandas as pd

import pytest

from gatk_sv_ploidy.preprocess import (
    _merge_intervals_by_chr,
    _process_sd_file,
    _estimate_site_pop_af,
    _infer_cohort_major_base,
    _site_nonmajor_total,
    build_per_site_data,
    collapse_bins_per_contig,
    filter_low_quality_bins,
    filter_poor_region_bins,
    main,
    normalise_depth,
    parse_args,
    read_bed_intervals,
    read_site_depth_tsv,
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


def test_filter_low_quality_bins_surgically_removes_high_median_chry_bin() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chrY", "chrY", "chrY"],
            "Start": [0, 100, 200],
            "End": [100, 200, 300],
            "S1": [0.10, 0.86, 0.18],
            "S2": [0.12, 0.89, 0.22],
            "S3": [0.08, 0.91, 0.26],
            "S4": [0.14, 0.94, 0.24],
        }
    )

    filtered = filter_low_quality_bins(
        df,
        min_bins_per_chr=1,
    )

    assert filtered["Start"].tolist() == [0, 200]


def test_filter_low_quality_bins_keeps_all_male_chry_bins() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chrY", "chrY", "chrY"],
            "Start": [0, 100, 200],
            "End": [100, 200, 300],
            "S1": [0.95, 1.02, 0.98],
            "S2": [1.00, 1.05, 0.97],
            "S3": [0.92, 0.99, 1.01],
        }
    )

    filtered = filter_low_quality_bins(
        df,
        min_bins_per_chr=1,
    )

    assert filtered["Start"].tolist() == [0, 100, 200]


def test_filter_low_quality_bins_stratifies_chrX_and_chrY() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chrX", "chrX", "chrY", "chrY"],
            "Start": [0, 100, 0, 100],
            "End": [100, 200, 100, 200],
            "XX1": [2.00, 2.55, 0.05, 0.45],
            "XX2": [2.10, 2.60, 0.03, 0.50],
            "XY1": [1.00, 1.55, 1.02, 1.45],
            "XY2": [0.95, 1.60, 1.05, 1.50],
        }
    )

    filtered = filter_low_quality_bins(
        df,
        cohort_deviation_threshold=0.3,
        cohort_deviation_fraction_max=0.5,
        min_bins_per_chr=1,
    )

    assert filtered[["Chr", "Start"]].values.tolist() == [["chrX", 0], ["chrY", 0]]


def test_merge_intervals_and_empty_poor_region_cases() -> None:
    intervals = pd.DataFrame(
        {
            "Chr": ["chr1", "chr1", "chr1", "chrX"],
            "Start": [0, 50, 200, 10],
            "End": [50, 150, 250, 20],
        }
    )
    merged = _merge_intervals_by_chr(intervals)
    assert merged == {"chr1": [(0, 150), (200, 250)], "chrX": [(10, 20)]}

    empty_df = pd.DataFrame(columns=["Chr", "Start", "End", "S1"])
    kept = filter_poor_region_bins(empty_df, intervals)
    assert kept.empty


def test_filter_poor_region_bins_rejects_invalid_threshold() -> None:
    df = pd.DataFrame({"Chr": ["chr1"], "Start": [0], "End": [100], "S1": [2.0]})
    poor = pd.DataFrame({"Chr": ["chr1"], "Start": [0], "End": [10]})
    with pytest.raises(ValueError, match="between 0 and 1"):
        filter_poor_region_bins(df, poor, min_poor_region_coverage=1.5)


def test_filter_low_quality_bins_raises_on_insufficient_bins() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chr21"],
            "Start": [0],
            "End": [100],
            "S1": [2.0],
            "S2": [2.1],
        }
    )

    with pytest.raises(ValueError, match="Insufficient bins"):
        filter_low_quality_bins(df, min_bins_per_chr=2)


def test_collapse_bins_per_contig_uses_weighted_depth_and_limit() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chr21"] * 10 + ["chrX"] * 2,
            "Start": [0, 10, 30, 60, 100, 150, 210, 280, 360, 450, 0, 100],
            "End": [10, 30, 60, 100, 150, 210, 280, 360, 450, 550, 100, 200],
            "source_file": ["test"] * 12,
            "S1": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 0.8, 1.2],
            "S2": [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 1.0, 1.4],
        },
        index=[f"bin{i}" for i in range(12)],
    )

    collapsed = collapse_bins_per_contig(df, bins_per_contig=4)

    assert collapsed[collapsed["Chr"] == "chr21"].shape[0] == 4
    assert collapsed[collapsed["Chr"] == "chrX"].shape[0] == 2

    chr21 = collapsed[collapsed["Chr"] == "chr21"].reset_index(drop=True)
    assert chr21.loc[0, "Start"] == 0
    assert chr21.loc[0, "End"] == 60
    assert chr21.loc[0, "BinLengthBp"] == 60
    expected_s1 = (10 * 1.0 + 20 * 2.0 + 30 * 3.0) / 60.0
    expected_s2 = (10 * 2.0 + 20 * 3.0 + 30 * 4.0) / 60.0
    assert chr21.loc[0, "S1"] == pytest.approx(expected_s1)
    assert chr21.loc[0, "S2"] == pytest.approx(expected_s2)


def test_collapse_bins_per_contig_can_sum_raw_counts() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chr21"] * 6,
            "Start": [0, 100, 200, 300, 400, 500],
            "End": [100, 200, 300, 400, 500, 600],
            "source_file": ["test"] * 6,
            "S1": [10, 20, 30, 40, 50, 60],
            "S2": [1, 2, 3, 4, 5, 6],
        },
        index=[f"bin{i}" for i in range(6)],
    )

    collapsed = collapse_bins_per_contig(
        df,
        bins_per_contig=3,
        aggregation="sum",
    )

    assert collapsed.shape[0] == 3
    assert collapsed.iloc[0]["S1"] == pytest.approx(30.0)
    assert collapsed.iloc[0]["S2"] == pytest.approx(3.0)
    assert collapsed.iloc[0]["BinLengthBp"] == 200


def test_collapse_bins_per_contig_preserves_covered_length_across_gaps() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chr21"] * 6,
            "Start": [0, 1000, 2000, 3000, 4000, 5000],
            "End": [100, 1100, 2100, 3100, 4100, 5100],
            "source_file": ["test"] * 6,
            "S1": [10, 20, 30, 40, 50, 60],
            "S2": [1, 2, 3, 4, 5, 6],
        },
        index=[f"bin{i}" for i in range(6)],
    )

    collapsed = collapse_bins_per_contig(
        df,
        bins_per_contig=3,
        aggregation="sum",
    ).reset_index(drop=True)

    assert collapsed.loc[0, "Start"] == 0
    assert collapsed.loc[0, "End"] == 1100
    assert collapsed.loc[0, "BinLengthBp"] == 200
    assert collapsed.loc[0, "S1"] == pytest.approx(30.0)


def test_read_site_depth_tsv_plain_and_position_stride(tmp_path) -> None:
    sd_path = tmp_path / "test.sd.txt"
    sd_path.write_text(
        "chr21\t0\tSAMPLE_A\t10\t0\t0\t20\n"
        "chr21\t50\tSAMPLE_A\t11\t0\t0\t19\n"
        "chr21\t100\tSAMPLE_A\t12\t0\t0\t18\n"
    )

    stride_df = read_site_depth_tsv(str(sd_path), stride=2)
    pos_df = read_site_depth_tsv(str(sd_path), position_stride=25)
    empty_df = read_site_depth_tsv(str(sd_path), position_stride=13)

    assert stride_df["position"].tolist() == [0, 100]
    assert pos_df["position"].tolist() == [0, 50, 100]
    assert empty_df["position"].tolist() == [0]


def test_site_cohort_major_helpers() -> None:
    a = np.array([10, 4])
    c = np.array([0, 0])
    g = np.array([0, 0])
    t = np.array([5, 6])

    major_idx = _infer_cohort_major_base(a, c, g, t)
    alt, total = _site_nonmajor_total(a, c, g, t, major_idx)

    assert major_idx == 0
    np.testing.assert_array_equal(total, np.array([15, 10]))
    np.testing.assert_array_equal(alt, np.array([5, 6]))
    assert _estimate_site_pop_af(alt, total) == pytest.approx((11.0 + 0.5) / 26.0)


def test_process_sd_file_collects_cohort_counts(tmp_path) -> None:
    sd_path = tmp_path / "sample.sd.txt"
    sd_path.write_text(
        "chr21\t20\tS1\t10\t0\t0\t5\n"
        "chr21\t75\tS1\t8\t0\t0\t4\n"
        "chr21\t20\tUNKNOWN\t9\t0\t0\t3\n"
        "chr21\t150\tS1\t7\t0\t0\t3\n"
    )
    sample_to_idx = {"S1": 0}
    bin_lookup = {"chr21": (np.array([0]), np.array([0]), np.array([100]))}

    _, total_sites, site_data = _process_sd_file(
        str(sd_path),
        sample_to_idx,
        bin_lookup,
        position_stride=1,
        min_site_depth=10,
    )
    assert total_sites == 2
    assert sorted(site_data[0]) == [20, 75]
    assert site_data[0][20]["values"][0] == (0, 10, 0, 0, 5)


def test_build_per_site_data_uses_shared_cohort_identity_and_site_limit(tmp_path) -> None:
    bins_df = pd.DataFrame(
        {
            "Chr": ["chr21"],
            "Start": [0],
            "End": [100],
            "S1": [2.0],
            "S2": [2.1],
        }
    )
    sd_path = tmp_path / "sample.sd.txt"
    sd_path.write_text(
        "chr21\t10\tS1\t8\t0\t0\t2\n"
        "chr21\t20\tS1\t2\t0\t0\t8\n"
        "chr21\t30\tS1\t8\t0\t0\t2\n"
        "chr21\t10\tS2\t7\t0\t0\t3\n"
        "chr21\t20\tS2\t7\t0\t0\t3\n"
    )

    result = build_per_site_data(
        [str(sd_path)],
        bins_df,
        max_sites_per_bin=2,
    )

    assert result["site_alt"].shape == (1, 2, 2)
    assert result["site_mask"].sum() == 4
    np.testing.assert_array_equal(result["site_alt"][0, 0, :], np.array([2, 3]))
    np.testing.assert_array_equal(result["site_alt"][0, 1, :], np.array([2, 7]))
    np.testing.assert_allclose(
        result["site_pop_af"][0, :],
        np.array([
            (5.0 + 0.5) / 21.0,
            (9.0 + 0.5) / 21.0,
        ]),
    )


def test_preprocess_parse_args_and_main_paths(tmp_path, monkeypatch) -> None:
    raw_path = tmp_path / "raw.tsv"
    raw_path.write_text(
        "#Chr\tStart\tEnd\tS1\tS2\n"
        "chr13\t0\t100\t10\t10\n"
        "chr18\t0\t100\t10\t10\n"
        "chr21\t0\t100\t10\t10\n"
        "chrX\t0\t100\t5\t10\n"
        "chrY\t0\t100\t5\t0\n"
        "chr1\t0\t100\t10\t10\n"
    )
    poor_regions = tmp_path / "poor.bed"
    poor_regions.write_text("chr13\t0\t10\n")
    sd_path = tmp_path / "sample.sd.txt"
    sd_path.write_text(
        "chr13\t20\tS1\t8\t0\t0\t2\n"
        "chr13\t20\tS2\t7\t0\t0\t3\n"
    )
    sd_list = tmp_path / "sd_list.txt"
    sd_list.write_text(f"{sd_path}\n")
    output_dir = tmp_path / "out"

    monkeypatch.setattr(
        "sys.argv",
        [
            "gatk-sv-ploidy preprocess",
            "--input",
            str(raw_path),
            "--output-dir",
            str(output_dir),
            "--viable-only",
            "--skip-bin-filter",
            "--output-space",
            "normalized",
            "--poor-regions",
            str(poor_regions),
            "--site-depth-list",
            str(sd_list),
            "--bins-per-contig",
            "1",
            "--min-bins-per-chr",
            "1",
            "--max-sites-per-bin",
            "2",
        ],
    )

    args = parse_args()
    assert args.viable_only is True
    assert args.skip_bin_filter is True
    assert args.site_depth_list == str(sd_list)
    assert args.bins_per_contig == 1
    assert args.min_bins_per_chr == 1

    main()

    preprocessed = pd.read_csv(output_dir / "preprocessed_depth.tsv", sep="\t", index_col=0)
    assert sorted(preprocessed["Chr"].unique().tolist()) == ["chr13", "chr18", "chr21", "chrX", "chrY"]
    assert len(preprocessed) == 5
    assert (output_dir / "observation_type.txt").read_text().strip() == "normalized"
    assert (output_dir / "site_data.npz").exists()
    site_npz = np.load(output_dir / "site_data.npz")
    observed_af = site_npz["site_pop_af"][np.any(site_npz["site_mask"], axis=2)]
    assert observed_af.tolist() == [pytest.approx((5.0 + 0.5) / 21.0)]


def test_preprocess_parse_args_defaults(monkeypatch) -> None:
    monkeypatch.setattr(
        "sys.argv",
        [
            "gatk-sv-ploidy preprocess",
            "--input",
            "in.tsv",
            "--output-dir",
            "out",
        ],
    )

    args = parse_args()

    assert args.bins_per_contig == 30
    assert args.output_space == "raw"
    assert args.chrY_median_max == pytest.approx(0.85)
    assert args.min_bins_per_chr == 10


def test_preprocess_can_write_raw_counts(tmp_path, monkeypatch) -> None:
    raw_path = tmp_path / "raw.tsv"
    raw_path.write_text(
        "#Chr\tStart\tEnd\tS1\tS2\n"
        "chr13\t0\t100\t10\t20\n"
        "chr18\t0\t100\t12\t24\n"
        "chr21\t0\t100\t14\t28\n"
        "chrX\t0\t100\t6\t12\n"
        "chrY\t0\t100\t4\t0\n"
    )
    output_dir = tmp_path / "out"

    monkeypatch.setattr(
        "sys.argv",
        [
            "gatk-sv-ploidy preprocess",
            "--input",
            str(raw_path),
            "--output-dir",
            str(output_dir),
            "--skip-bin-filter",
            "--output-space",
            "raw",
            "--min-bins-per-chr",
            "1",
            "--bins-per-contig",
            "0",
        ],
    )

    main()

    preprocessed = pd.read_csv(
        output_dir / "preprocessed_depth.tsv",
        sep="\t",
        index_col=0,
    )
    assert preprocessed.loc["chr13:0-100", "S1"] == pytest.approx(10.0)
    assert preprocessed.loc["chr18:0-100", "S2"] == pytest.approx(24.0)
    assert (output_dir / "observation_type.txt").read_text().strip() == "raw"
