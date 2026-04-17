from __future__ import annotations

import numpy as np
import pandas as pd

import pytest

from gatk_sv_ploidy.preprocess import (
    _merge_intervals_by_chr,
    _process_sd_file,
    _site_alt_total,
    _site_minor_total,
    build_per_site_data,
    filter_low_quality_bins,
    filter_poor_region_bins,
    main,
    normalise_depth,
    parse_args,
    read_bed_intervals,
    read_known_sites,
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


def test_site_alt_and_minor_total_helpers() -> None:
    a = np.array([10, 0])
    c = np.array([0, 5])
    g = np.array([0, 0])
    t = np.array([5, 5])
    ref = np.array(["A", "C"])

    alt, total = _site_alt_total(a, c, g, t, ref)
    minor, total_minor = _site_minor_total(a, c, g, t)

    np.testing.assert_array_equal(total, np.array([15, 10]))
    np.testing.assert_array_equal(alt, np.array([5, 5]))
    np.testing.assert_array_equal(total_minor, np.array([15, 10]))
    np.testing.assert_array_equal(minor, np.array([5, 5]))


def test_read_known_sites_reads_tsv(tmp_path) -> None:
    path = tmp_path / "sites.tsv"
    path.write_text("chr21\t500\tA\t0.3\nchrX\t100\tC\t0.4\n")
    sites = read_known_sites(str(path))
    assert sites["contig"].tolist() == ["chr21", "chrX"]
    assert sites["pop_af"].tolist() == [0.3, 0.4]


def test_process_sd_file_handles_known_and_fallback_paths(tmp_path) -> None:
    sd_path = tmp_path / "sample.sd.txt"
    sd_path.write_text(
        "chr21\t20\tS1\t10\t0\t0\t5\n"
        "chr21\t75\tS1\t8\t0\t0\t4\n"
        "chr21\t20\tUNKNOWN\t9\t0\t0\t3\n"
        "chr21\t150\tS1\t7\t0\t0\t3\n"
    )
    sample_to_idx = {"S1": 0}
    bin_lookup = {"chr21": (np.array([0]), np.array([0]), np.array([100]))}

    known_arrays = {
        "chr21": (
            np.array([20], dtype=np.int64),
            np.array(["A"]),
            np.array([0.2], dtype=np.float32),
        )
    }
    _, known_total, known_data = _process_sd_file(
        str(sd_path),
        sample_to_idx,
        bin_lookup,
        known_arrays,
        have_known=True,
        effective_stride=1,
        min_site_depth=5,
    )
    assert known_total == 1
    assert known_data[0][20]["values"][0][:3] == (0, 5, 15)

    _, fallback_total, fallback_data = _process_sd_file(
        str(sd_path),
        sample_to_idx,
        bin_lookup,
        {},
        have_known=False,
        effective_stride=1,
        min_site_depth=10,
    )
    assert fallback_total == 2
    assert sorted(fallback_data[0]) == [20, 75]


def test_build_per_site_data_known_sites_and_site_limit(tmp_path) -> None:
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
        "chr21\t20\tS1\t8\t0\t0\t2\n"
        "chr21\t30\tS1\t8\t0\t0\t2\n"
        "chr21\t10\tS2\t7\t0\t0\t3\n"
        "chr21\t20\tS2\t7\t0\t0\t3\n"
    )
    known_sites = pd.DataFrame(
        {
            "contig": ["chr21", "chr21", "chr21"],
            "position": [10, 20, 30],
            "ref": ["A", "A", "A"],
            "pop_af": [0.1, 0.2, 0.3],
        }
    )

    result = build_per_site_data(
        [str(sd_path)],
        bins_df,
        known_sites_df=known_sites,
        max_sites_per_bin=2,
    )

    assert result["site_alt"].shape == (1, 2, 2)
    assert result["site_mask"].sum() == 4
    np.testing.assert_allclose(sorted(result["site_pop_af"][0].tolist()), [0.1, 0.2])


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
    known_sites = tmp_path / "sites.tsv"
    known_sites.write_text("chr13\t20\tA\t0.2\n")
    sd_path = tmp_path / "sample.sd.txt"
    sd_path.write_text("chr13\t20\tS1\t8\t0\t0\t2\n")
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
            "--poor-regions",
            str(poor_regions),
            "--site-depth-list",
            str(sd_list),
            "--known-sites",
            str(known_sites),
            "--max-sites-per-bin",
            "2",
        ],
    )

    args = parse_args()
    assert args.viable_only is True
    assert args.skip_bin_filter is True
    assert args.site_depth_list == str(sd_list)

    main()

    preprocessed = pd.read_csv(output_dir / "preprocessed_depth.tsv", sep="\t", index_col=0)
    assert sorted(preprocessed["Chr"].unique().tolist()) == ["chr13", "chr18", "chr21", "chrX", "chrY"]
    assert (output_dir / "site_data.npz").exists()
