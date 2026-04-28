from __future__ import annotations

import gzip

import pandas as pd
import pytest

from gatk_sv_ploidy.call import (
    _normalize_global_cn_artifact,
    _classify_aneuploidy,
    _classify_sex,
    apply_binq_filter_to_chromosome_stats,
    assign_sex_and_aneuploidy_types,
    export_ignored_bins,
    export_aneuploid_data,
    save_sex_assignments,
    summarize_binq_filter_stats,
)


def test_classify_sex_and_aneuploidy_tables() -> None:
    assert _classify_sex(1, 1) == "MALE"
    assert _classify_sex(2, 0) == "FEMALE"
    assert _classify_sex(3, 0) == "TRIPLE_X"
    assert _classify_sex(2, 2) == "OTHER"
    assert _classify_sex(8, 1) == "OTHER"

    assert _classify_aneuploidy({}, {}, 2, 0) == "NORMAL"
    assert _classify_aneuploidy({"chr21": True}, {"chr21": 3}, 2, 0) == "TRISOMY_21"
    assert _classify_aneuploidy(
        {"chr13": True, "chr18": True, "chr21": True},
        {"chr13": 3, "chr18": 3, "chr21": 3, "chrX": 1, "chrY": 2},
        1,
        2,
        sample_depth_ratio=1.5,
    ) == "MULTIPLE"
    assert _classify_aneuploidy(
        {},
        {"chr13": 2, "chr18": 2, "chr21": 2, "chrX": 1, "chrY": 1},
        1,
        1,
        sample_depth_ratio=2.0,
    ) == "NORMAL"
    assert _classify_aneuploidy({"chrX": True}, {"chrX": 1, "chrY": 0}, 1, 0) == "TURNER"
    assert _classify_aneuploidy(
        {"chr13": True, "chr18": True},
        {"chr13": 3, "chr18": 3},
        2,
        0,
    ) == "MULTIPLE"


def test_normalize_global_cn_artifact_is_no_op() -> None:
    normalized_cn_map, scale_factor = _normalize_global_cn_artifact(
        {"chr13": 4, "chr18": 4, "chr21": 4, "chrX": 2, "chrY": 2},
        sample_depth_ratio=1.0,
    )
    assert scale_factor == 1.0
    assert normalized_cn_map == {
        "chr13": 4,
        "chr18": 4,
        "chr21": 4,
        "chrX": 2,
        "chrY": 2,
    }

    unchanged_cn_map, unchanged_scale = _normalize_global_cn_artifact(
        {"chr13": 4, "chr18": 4, "chr21": 4, "chrX": 2, "chrY": 2},
        sample_depth_ratio=2.0,
    )
    assert unchanged_scale == 1.0
    assert unchanged_cn_map["chr13"] == 4

    trisomy_cn_map, trisomy_scale = _normalize_global_cn_artifact(
        {"chr13": 4, "chr18": 4, "chr21": 5, "chrX": 2, "chrY": 2},
        sample_depth_ratio=1.0,
    )
    assert trisomy_scale == 1.0
    assert trisomy_cn_map == {
        "chr13": 4,
        "chr18": 4,
        "chr21": 5,
        "chrX": 2,
        "chrY": 2,
    }

    preserved_cn_map, preserved_scale = _normalize_global_cn_artifact(
        {"chr13": 4, "chr18": 4, "chr21": 4, "chrX": 2, "chrY": 2},
        sample_depth_ratio=1.0,
        normalize_doubled_labels=False,
    )
    assert preserved_scale == 1.0
    assert preserved_cn_map["chr13"] == 4


def test_assign_and_export_outputs(tmp_path) -> None:
    df = pd.DataFrame(
        [
            {
                "sample": "S1",
                "chromosome": "chrX",
                "copy_number": 2,
                "is_aneuploid": False,
                "mean_cn_probability": 0.95,
                "median_depth": 2.0,
                "sample_var_map": 0.04,
                "chr_type": 1,
            },
            {
                "sample": "S1",
                "chromosome": "chrY",
                "copy_number": 0,
                "is_aneuploid": False,
                "mean_cn_probability": 0.93,
                "median_depth": 0.0,
                "sample_var_map": 0.04,
                "chr_type": 2,
            },
            {
                "sample": "S1",
                "chromosome": "chr21",
                "copy_number": 2,
                "is_aneuploid": False,
                "mean_cn_probability": 0.97,
                "median_depth": 2.0,
                "sample_var_map": 0.04,
                "chr_type": 0,
            },
            {
                "sample": "S2",
                "chromosome": "chrX",
                "copy_number": 1,
                "is_aneuploid": True,
                "mean_cn_probability": 0.91,
                "median_depth": 1.0,
                "sample_var_map": 0.09,
                "chr_type": 1,
            },
            {
                "sample": "S2",
                "chromosome": "chrY",
                "copy_number": 0,
                "is_aneuploid": True,
                "mean_cn_probability": 0.92,
                "median_depth": 0.0,
                "sample_var_map": 0.09,
                "chr_type": 2,
            },
            {
                "sample": "S2",
                "chromosome": "chr21",
                "copy_number": 3,
                "is_aneuploid": True,
                "mean_cn_probability": 0.90,
                "median_depth": 3.0,
                "sample_var_map": 0.09,
                "chr_type": 0,
            },
        ]
    )

    pred_df = assign_sex_and_aneuploidy_types(df, {"S2": "MULTIPLE"})

    assert pred_df.set_index("sample").loc["S1", "predicted_aneuploidy_type"] == "NORMAL"
    assert pred_df.set_index("sample").loc["S2", "predicted_aneuploidy_type"] == "MULTIPLE"

    save_sex_assignments(pred_df, str(tmp_path))
    with gzip.open(tmp_path / "sex_assignments.txt.gz", "rt") as handle:
        written = handle.read()
    assert "sample_id" in written
    assert "Assignment" in written

    export_aneuploid_data(df, str(tmp_path / "aneuploid.tsv"))
    exported = pd.read_csv(tmp_path / "aneuploid.tsv", sep="\t")
    assert "sample_overdispersion_map" in exported.columns
    assert "sample_var_map" not in exported.columns
    assert exported["sample_overdispersion_map"].iloc[0] == pytest.approx(0.09)
    assert "chr_type" not in exported.columns


def test_assign_preserves_absolute_cn_labels() -> None:
    rows = [
        {
            "sample": "reference",
            "chromosome": chrom,
            "copy_number": copy_number,
            "is_aneuploid": False,
            "mean_cn_probability": 0.95,
            "median_depth": 1.0,
            "sample_depth_map": 110.0,
        }
        for chrom, copy_number in (
            ("chr13", 2),
            ("chr18", 2),
            ("chr21", 2),
            ("chrX", 1),
            ("chrY", 1),
        )
    ]
    rows.extend(
        {
            "sample": "artifact",
            "chromosome": chrom,
            "copy_number": copy_number,
            "is_aneuploid": True,
            "mean_cn_probability": 0.95,
            "median_depth": 1.0,
            "sample_depth_map": 100.0,
        }
        for chrom, copy_number in (
            ("chr13", 4),
            ("chr18", 4),
            ("chr21", 4),
            ("chrX", 2),
            ("chrY", 2),
        )
    )
    rows.extend(
        {
            "sample": "tetraploid",
            "chromosome": chrom,
            "copy_number": copy_number,
            "is_aneuploid": True,
            "mean_cn_probability": 0.95,
            "median_depth": 1.0,
            "sample_depth_map": 240.0,
        }
        for chrom, copy_number in (
            ("chr13", 4),
            ("chr18", 4),
            ("chr21", 4),
            ("chrX", 2),
            ("chrY", 2),
        )
    )
    df = pd.DataFrame(rows)

    pred_df = assign_sex_and_aneuploidy_types(df)
    indexed_pred_df = pred_df.set_index("sample")

    assert indexed_pred_df.loc["artifact", "sex"] == "OTHER"
    assert indexed_pred_df.loc["artifact", "predicted_aneuploidy_type"] == "MULTIPLE"
    assert indexed_pred_df.loc["artifact", "chrX_CN"] == 2
    assert indexed_pred_df.loc["artifact", "chrY_CN"] == 2
    assert indexed_pred_df.loc["artifact", "global_cn_scale_factor"] == 1.0

    assert indexed_pred_df.loc["tetraploid", "sex"] == "OTHER"
    assert indexed_pred_df.loc["tetraploid", "predicted_aneuploidy_type"] == "MULTIPLE"
    assert indexed_pred_df.loc["tetraploid", "chrX_CN"] == 2
    assert indexed_pred_df.loc["tetraploid", "chrY_CN"] == 2
    assert indexed_pred_df.loc["tetraploid", "global_cn_scale_factor"] == 1.0


def test_assign_normalize_doubled_labels_flag_is_ignored() -> None:
    df = pd.DataFrame(
        [
            {
                "sample": "tetraploid",
                "chromosome": chrom,
                "copy_number": copy_number,
                "is_aneuploid": True,
                "mean_cn_probability": 0.95,
                "median_depth": 1.0,
                "sample_depth_map": 110.0,
            }
            for chrom, copy_number in (
                ("chr13", 4),
                ("chr18", 4),
                ("chr21", 4),
                ("chrX", 2),
                ("chrY", 2),
            )
        ]
    )

    pred_default = assign_sex_and_aneuploidy_types(df).set_index("sample")
    pred_df = assign_sex_and_aneuploidy_types(
        df,
        normalize_doubled_labels=False,
    ).set_index("sample")

    assert pred_default.loc["tetraploid", "sex"] == "OTHER"
    assert pred_default.loc["tetraploid", "predicted_aneuploidy_type"] == "MULTIPLE"
    assert pred_df.loc["tetraploid", "sex"] == "OTHER"
    assert pred_df.loc["tetraploid", "predicted_aneuploidy_type"] == "MULTIPLE"
    assert pred_df.loc["tetraploid", "global_cn_scale_factor"] == 1.0


def test_apply_binq_filter_to_chromosome_stats_rebuilds_autosome_calls() -> None:
    chrom_df = pd.DataFrame(
        [
            {
                "sample": "S1",
                "chromosome": "chrX",
                "copy_number": 2,
                "is_aneuploid": False,
                "mean_cn_probability": 0.95,
                "median_depth": 2.0,
                "n_bins": 1,
            },
            {
                "sample": "S1",
                "chromosome": "chrY",
                "copy_number": 0,
                "is_aneuploid": False,
                "mean_cn_probability": 0.93,
                "median_depth": 0.0,
                "n_bins": 1,
            },
            {
                "sample": "S1",
                "chromosome": "chr21",
                "copy_number": 3,
                "is_aneuploid": True,
                "mean_cn_probability": 0.90,
                "median_depth": 3.0,
                "n_bins": 2,
            },
        ]
    )
    bin_stats_df = pd.DataFrame(
        [
            {
                "chr": "chr21",
                "start": 0,
                "end": 100,
                "sample": "S1",
                "observed_depth": 3.0,
                "cn_map": 3,
                "cn_prob_0": 0.0,
                "cn_prob_1": 0.0,
                "cn_prob_2": 0.05,
                "cn_prob_3": 0.95,
                "cn_prob_4": 0.0,
                "cn_prob_5": 0.0,
                "cnq": 5,
                "sample_var": 0.04,
                "plot_depth": 3.0,
            },
            {
                "chr": "chr21",
                "start": 100,
                "end": 200,
                "sample": "S1",
                "observed_depth": 2.0,
                "cn_map": 2,
                "cn_prob_0": 0.0,
                "cn_prob_1": 0.0,
                "cn_prob_2": 0.99,
                "cn_prob_3": 0.01,
                "cn_prob_4": 0.0,
                "cn_prob_5": 0.0,
                "cnq": 20,
                "sample_var": 0.04,
                "plot_depth": 2.0,
            },
        ]
    )
    bin_quality_df = pd.DataFrame(
        [
            {"chr": "chr21", "start": 0, "end": 100, "BINQ20": 5.0, "BINQ15": 3.0},
            {"chr": "chr21", "start": 100, "end": 200, "BINQ20": 40.0, "BINQ15": 30.0},
        ]
    )

    filtered = apply_binq_filter_to_chromosome_stats(
        chrom_df,
        bin_stats_df,
        bin_quality_df,
        min_binq=20.0,
        binq_field="BINQ20",
        prob_threshold=0.5,
    )

    chr21 = filtered[filtered["chromosome"] == "chr21"].iloc[0]
    chrY = filtered[filtered["chromosome"] == "chrY"].iloc[0]

    assert int(chr21["copy_number"]) == 2
    assert bool(chr21["is_aneuploid"]) is False
    assert int(chr21["n_bins_retained"]) == 1
    assert float(chr21["frac_bins_retained"]) == 0.5
    assert int(chr21["plq"]) == 20
    assert float(chr21["ploidy_prob_2"]) == pytest.approx(1.0)
    assert int(chrY["copy_number"]) == 0


def test_apply_binq_filter_handles_fully_filtered_run() -> None:
    chrom_df = pd.DataFrame(
        [
            {
                "sample": "S1",
                "chromosome": "chrX",
                "copy_number": 2,
                "is_aneuploid": False,
                "mean_cn_probability": 0.95,
                "median_depth": 2.0,
                "n_bins": 1,
            },
            {
                "sample": "S1",
                "chromosome": "chr21",
                "copy_number": 3,
                "is_aneuploid": True,
                "mean_cn_probability": 0.90,
                "median_depth": 3.0,
                "n_bins": 2,
            },
        ]
    )
    bin_stats_df = pd.DataFrame(
        [
            {
                "chr": "chr21",
                "start": 0,
                "end": 100,
                "sample": "S1",
                "observed_depth": 3.0,
                "cn_map": 3,
                "cn_prob_0": 0.0,
                "cn_prob_1": 0.0,
                "cn_prob_2": 0.05,
                "cn_prob_3": 0.95,
                "cn_prob_4": 0.0,
                "cn_prob_5": 0.0,
                "cnq": 20,
            },
            {
                "chr": "chr21",
                "start": 100,
                "end": 200,
                "sample": "S1",
                "observed_depth": 3.0,
                "cn_map": 3,
                "cn_prob_0": 0.0,
                "cn_prob_1": 0.0,
                "cn_prob_2": 0.05,
                "cn_prob_3": 0.95,
                "cn_prob_4": 0.0,
                "cn_prob_5": 0.0,
                "cnq": 20,
            },
        ]
    )
    bin_quality_df = pd.DataFrame(
        [
            {"chr": "chr21", "start": 0, "end": 100, "BINQ20": 5.0},
            {"chr": "chr21", "start": 100, "end": 200, "BINQ20": 5.0},
        ]
    )

    filtered = apply_binq_filter_to_chromosome_stats(
        chrom_df,
        bin_stats_df,
        bin_quality_df,
        min_binq=20.0,
        binq_field="BINQ20",
        prob_threshold=0.5,
    )

    chr21 = filtered[filtered["chromosome"] == "chr21"].iloc[0]
    chrX = filtered[filtered["chromosome"] == "chrX"].iloc[0]

    assert int(chr21["copy_number"]) == 2
    assert pd.isna(chr21["mean_cn_probability"])
    assert int(chr21["n_bins_retained"]) == 0
    assert float(chr21["frac_bins_retained"]) == 0.0
    assert bool(chr21["is_aneuploid"]) is False
    assert int(chrX["copy_number"]) == 2
    assert int(chrX["n_bins_retained"]) == 0


def test_apply_binq_filter_uses_majority_vote_and_average_cnq() -> None:
    chrom_df = pd.DataFrame(
        [
            {
                "sample": "S1",
                "chromosome": "chr21",
                "copy_number": 2,
                "is_aneuploid": False,
                "mean_cn_probability": 0.75,
                "median_depth": 2.5,
                "n_bins": 3,
            },
        ]
    )
    bin_stats_df = pd.DataFrame(
        [
            {
                "chr": "chr21",
                "start": 0,
                "end": 100,
                "sample": "S1",
                "observed_depth": 2.1,
                "cn_map": 3,
                "cn_prob_0": 0.0,
                "cn_prob_1": 0.0,
                "cn_prob_2": 0.05,
                "cn_prob_3": 0.95,
                "cn_prob_4": 0.0,
                "cn_prob_5": 0.0,
                "cnq": 20,
            },
            {
                "chr": "chr21",
                "start": 100,
                "end": 200,
                "sample": "S1",
                "observed_depth": 3.1,
                "cn_map": 3,
                "cn_prob_0": 0.0,
                "cn_prob_1": 0.0,
                "cn_prob_2": 0.10,
                "cn_prob_3": 0.90,
                "cn_prob_4": 0.0,
                "cn_prob_5": 0.0,
                "cnq": 30,
            },
            {
                "chr": "chr21",
                "start": 200,
                "end": 300,
                "sample": "S1",
                "observed_depth": 1.9,
                "cn_map": 2,
                "cn_prob_0": 0.0,
                "cn_prob_1": 0.0,
                "cn_prob_2": 0.95,
                "cn_prob_3": 0.05,
                "cn_prob_4": 0.0,
                "cn_prob_5": 0.0,
                "cnq": 40,
            },
        ]
    )
    bin_quality_df = pd.DataFrame(
        [
            {"chr": "chr21", "start": 0, "end": 100, "BINQ20": 40.0, "BINQ15": 30.0},
            {"chr": "chr21", "start": 100, "end": 200, "BINQ20": 40.0, "BINQ15": 30.0},
            {"chr": "chr21", "start": 200, "end": 300, "BINQ20": 40.0, "BINQ15": 30.0},
        ]
    )

    filtered = apply_binq_filter_to_chromosome_stats(
        chrom_df,
        bin_stats_df,
        bin_quality_df,
        min_binq=20.0,
        binq_field="BINQ20",
        prob_threshold=0.5,
    )

    row = filtered.iloc[0]
    assert int(row["copy_number"]) == 3
    assert float(row["mean_cn_probability"]) == pytest.approx(2.0 / 3.0)
    assert int(row["plq"]) == 30
    assert bool(row["is_aneuploid"]) is True


def test_apply_binq_filter_preserves_original_sex_call_on_conflicting_subset() -> None:
    chrom_df = pd.DataFrame(
        [
            {
                "sample": "S1",
                "chromosome": "chrX",
                "copy_number": 1,
                "is_aneuploid": False,
                "mean_cn_probability": 0.96,
                "median_depth": 1.0,
                "n_bins": 3,
            },
            {
                "sample": "S1",
                "chromosome": "chrY",
                "copy_number": 1,
                "is_aneuploid": False,
                "mean_cn_probability": 0.53,
                "median_depth": 1.0,
                "n_bins": 10,
            },
        ]
    )
    bin_stats_rows = []
    bin_quality_rows = []
    for idx, cn_map in enumerate([1, 1, 1, 1, 2, 2, 2, 2, 2, 1]):
        start = idx * 100
        bin_stats_rows.append(
            {
                "chr": "chrY",
                "start": start,
                "end": start + 100,
                "sample": "S1",
                "observed_depth": 1.0 if cn_map == 1 else 2.0,
                "cn_map": cn_map,
                "cn_prob_0": 0.0,
                "cn_prob_1": 0.9 if cn_map == 1 else 0.1,
                "cn_prob_2": 0.1 if cn_map == 1 else 0.9,
                "cn_prob_3": 0.0,
                "cn_prob_4": 0.0,
                "cn_prob_5": 0.0,
                "cnq": 30,
            }
        )
        bin_quality_rows.append(
            {
                "chr": "chrY",
                "start": start,
                "end": start + 100,
                "BINQ20": 40.0 if idx < 9 else 5.0,
            }
        )

    filtered = apply_binq_filter_to_chromosome_stats(
        chrom_df,
        pd.DataFrame(bin_stats_rows),
        pd.DataFrame(bin_quality_rows),
        min_binq=20.0,
        binq_field="BINQ20",
        prob_threshold=0.5,
    )

    chr_y = filtered[filtered["chromosome"] == "chrY"].iloc[0]

    assert int(chr_y["copy_number"]) == 1
    assert float(chr_y["mean_cn_probability"]) == pytest.approx(0.53)
    assert int(chr_y["n_bins_retained"]) == 9
    assert float(chr_y["frac_bins_retained"]) == pytest.approx(0.9)
    assert bool(chr_y["is_aneuploid"]) is False


def test_apply_binq_filter_auto_prefers_binq_when_available() -> None:
    chrom_df = pd.DataFrame(
        [
            {
                "sample": "S1",
                "chromosome": "chr21",
                "copy_number": 3,
                "is_aneuploid": True,
                "mean_cn_probability": 0.90,
                "median_depth": 3.0,
                "n_bins": 2,
            },
        ]
    )
    bin_stats_df = pd.DataFrame(
        [
            {
                "chr": "chr21",
                "start": 0,
                "end": 100,
                "sample": "S1",
                "observed_depth": 3.0,
                "cn_map": 3,
                "cn_prob_0": 0.0,
                "cn_prob_1": 0.0,
                "cn_prob_2": 0.05,
                "cn_prob_3": 0.95,
                "cn_prob_4": 0.0,
                "cn_prob_5": 0.0,
            },
            {
                "chr": "chr21",
                "start": 100,
                "end": 200,
                "sample": "S1",
                "observed_depth": 2.0,
                "cn_map": 2,
                "cn_prob_0": 0.0,
                "cn_prob_1": 0.0,
                "cn_prob_2": 0.99,
                "cn_prob_3": 0.01,
                "cn_prob_4": 0.0,
                "cn_prob_5": 0.0,
            },
        ]
    )
    bin_quality_df = pd.DataFrame(
        [
            {"chr": "chr21", "start": 0, "end": 100, "BINQ20": 40.0, "CALLQ20": 5.0},
            {"chr": "chr21", "start": 100, "end": 200, "BINQ20": 5.0, "CALLQ20": 40.0},
        ]
    )

    filtered = apply_binq_filter_to_chromosome_stats(
        chrom_df,
        bin_stats_df,
        bin_quality_df,
        min_binq=20.0,
        binq_field="auto",
        prob_threshold=0.5,
    )

    row = filtered.iloc[0]
    assert int(row["copy_number"]) == 3
    assert bool(row["is_aneuploid"]) is True


def test_export_ignored_bins_writes_gzip_table(tmp_path) -> None:
    ignored_df = pd.DataFrame(
        [
            {
                "sample": "S1",
                "chr": "chr21",
                "start": 0,
                "end": 100,
                "ignored_in_call": True,
                "binq_field": "BINQ20",
                "binq_value": 5.0,
                "min_binq": 20.0,
                "cn_map": 3,
            },
            {
                "sample": "S1",
                "chr": "chr21",
                "start": 100,
                "end": 200,
                "ignored_in_call": False,
                "binq_field": "BINQ20",
                "binq_value": 40.0,
                "min_binq": 20.0,
                "cn_map": 2,
            },
        ]
    )

    output_path = tmp_path / "ignored_bins.tsv.gz"
    export_ignored_bins(ignored_df, str(output_path))

    exported = pd.read_csv(output_path, sep="\t", compression="gzip")
    assert len(exported) == 1
    assert exported.loc[0, "sample"] == "S1"
    assert exported.loc[0, "binq_field"] == "BINQ20"


def test_summarize_binq_filter_stats_reports_unique_bin_coverage() -> None:
    annotated_bin_df = pd.DataFrame(
        [
            {
                "sample": "S1",
                "chr": "chr1",
                "start": 0,
                "end": 100,
                "ignored_in_call": True,
            },
            {
                "sample": "S2",
                "chr": "chr1",
                "start": 0,
                "end": 100,
                "ignored_in_call": True,
            },
            {
                "sample": "S1",
                "chr": "chr1",
                "start": 100,
                "end": 300,
                "ignored_in_call": False,
            },
            {
                "sample": "S2",
                "chr": "chr1",
                "start": 100,
                "end": 300,
                "ignored_in_call": False,
            },
            {
                "sample": "S1",
                "chr": "chr2",
                "start": 0,
                "end": 50,
                "ignored_in_call": False,
            },
        ]
    )

    summary = summarize_binq_filter_stats(annotated_bin_df)

    chr1 = summary[summary["chromosome"] == "chr1"].iloc[0]
    chr2 = summary[summary["chromosome"] == "chr2"].iloc[0]

    assert int(chr1["n_bins"]) == 2
    assert int(chr1["n_bins_filtered"]) == 1
    assert float(chr1["pct_coverage_remaining"]) == pytest.approx(200.0 / 300.0 * 100.0)
    assert float(chr1["coverage_remaining_mb"]) == pytest.approx(200.0 / 1e6)

    assert int(chr2["n_bins"]) == 1
    assert int(chr2["n_bins_filtered"]) == 0
    assert float(chr2["pct_coverage_remaining"]) == pytest.approx(100.0)
    assert float(chr2["coverage_remaining_mb"]) == pytest.approx(50.0 / 1e6)
