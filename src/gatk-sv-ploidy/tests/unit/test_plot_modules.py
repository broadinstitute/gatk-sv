from __future__ import annotations

import sys
import logging

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from gatk_sv_ploidy._plot_detail import (
    _draw_af_reference_lines,
    _contiguous_bar_geometry,
    _depth_axis_label_and_limit,
    _draw_sample_depth_panel,
    _draw_site_af_scatter,
    _equal_width_x,
    _expected_cn_range,
    _samples_by_chrX,
    _samples_by_chrX_bin,
    plot_cn_per_bin_chromosome,
    plot_cn_per_contig_boxplot,
    plot_sample_with_variance,
    plot_sex_assignments,
)
from gatk_sv_ploidy._plot_ppd import run_ppd_plots
from gatk_sv_ploidy._plot_report import write_plot_report
from gatk_sv_ploidy._plot_style import (
    _suppress_chatty_plot_dependency_loggers,
    apply_theme,
    get_plot_output_format,
    set_plot_output_format,
)
from gatk_sv_ploidy.plot import (
    _annotate_binq_values,
    _annotate_ignored_bins,
    _apply_plot_depth_bin_columns,
    _apply_plot_depth_columns,
    _apply_plot_depth_sex_columns,
    _sample_baseline_ploidy_metadata,
    _run_aneuploidy_plots,
    _run_ploidy_plots,
    parse_args,
    plot_background_factor_diagnostics,
    plot_binq_genome_profile,
    plot_chromosome_cn_heatmap,
    plot_chromosome_plq_heatmap,
    plot_histograms_by_chr_type,
    plot_median_depth_distributions,
    plot_site_af_estimates,
    plot_training_loss_with_gradient,
    main,
)


def _small_chrom_df() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"sample": "S1", "chromosome": "chr21", "copy_number": 2, "mean_cn_probability": 0.9, "plq": 40, "is_aneuploid": False, "mean_depth": 2.0, "median_depth": 2.0, "std_depth": 0.1, "mad_depth": 0.1, "chrX_depth": 1.0, "chrY_depth": 1.0, "sex": "MALE"},
            {"sample": "S1", "chromosome": "chrX", "copy_number": 1, "mean_cn_probability": 0.9, "plq": 35, "is_aneuploid": False, "mean_depth": 1.0, "median_depth": 1.0, "std_depth": 0.1, "mad_depth": 0.1, "chrX_depth": 1.0, "chrY_depth": 1.0, "sex": "MALE"},
            {"sample": "S1", "chromosome": "chrY", "copy_number": 1, "mean_cn_probability": 0.9, "plq": 30, "is_aneuploid": False, "mean_depth": 1.0, "median_depth": 1.0, "std_depth": 0.1, "mad_depth": 0.1, "chrX_depth": 1.0, "chrY_depth": 1.0, "sex": "MALE"},
            {"sample": "S2", "chromosome": "chr21", "copy_number": 3, "mean_cn_probability": 0.95, "plq": 60, "is_aneuploid": True, "mean_depth": 3.0, "median_depth": 3.0, "std_depth": 0.2, "mad_depth": 0.1, "chrX_depth": 2.0, "chrY_depth": 0.0, "sex": "FEMALE"},
            {"sample": "S2", "chromosome": "chrX", "copy_number": 2, "mean_cn_probability": 0.95, "plq": 55, "is_aneuploid": False, "mean_depth": 2.0, "median_depth": 2.0, "std_depth": 0.1, "mad_depth": 0.1, "chrX_depth": 2.0, "chrY_depth": 0.0, "sex": "FEMALE"},
            {"sample": "S2", "chromosome": "chrY", "copy_number": 0, "mean_cn_probability": 0.95, "plq": 50, "is_aneuploid": False, "mean_depth": 0.0, "median_depth": 0.0, "std_depth": 0.1, "mad_depth": 0.1, "chrX_depth": 2.0, "chrY_depth": 0.0, "sex": "FEMALE"},
        ]
    )


def _small_bin_df() -> pd.DataFrame:
    rows = []
    for sample, x_cn, y_cn, chr21_cn in [("S1", 1, 1, 2), ("S2", 2, 0, 3)]:
        sample_depth = 1.0 if sample == "S1" else 1.5
        for chrom, cn, depth in [("chr21", chr21_cn, float(chr21_cn)), ("chrX", x_cn, float(x_cn)), ("chrY", y_cn, float(y_cn))]:
            for start in [0, 100]:
                probs = [0.01] * 6
                probs[cn] = 0.95
                rows.append(
                    {
                        "chr": chrom,
                        "start": start,
                        "end": start + 100,
                        "sample": sample,
                        "observed_depth": depth,
                        "cn_map": cn,
                        "cnq": 10,
                        "binq_field": "BINQ20",
                        "binq_value": 25.0 if start == 0 else 40.0,
                        "bin_bias": 1.0,
                        "bin_var": 0.04,
                        "bin_epsilon": 0.0,
                        "bin_length_kb": 1.0,
                        "sample_depth": sample_depth,
                        "sample_var": 0.09 if sample == "S1" else 0.16,
                        "mean_het_af": 0.5,
                        "n_het_sites": 2,
                        **{f"cn_prob_{i}": probs[i] for i in range(6)},
                    }
                )
    return pd.DataFrame(rows)


def _small_loss_df() -> pd.DataFrame:
    return pd.DataFrame({"epoch": [0, 1, 2], "elbo": [100.0, 90.0, 85.0]})


def _small_ignored_bins_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample": ["S2"],
            "chr": ["chr21"],
            "start": [0],
            "end": [100],
            "ignored_in_call": [True],
            "binq_field": ["BINQ20"],
            "binq_value": [5.0],
            "min_binq": [20.0],
        }
    )


def _small_ppd_frames() -> tuple[pd.DataFrame, pd.DataFrame]:
    ppd_bin_df = pd.DataFrame(
        {
            "chr": ["chr21", "chr21", "chrX", "chrY"],
            "start": [0, 100, 0, 0],
            "end": [100, 200, 100, 100],
            "sample": ["S1", "S2", "S1", "S2"],
            "observed_depth": [2.0, 3.0, 1.0, 0.0],
            "ppd_mean": [2.1, 2.8, 1.1, 0.1],
            "cn_map": [2, 3, 1, 0],
            "ppd_q05": [1.5, 2.0, 0.5, -0.1],
            "ppd_q25": [1.8, 2.5, 0.8, 0.0],
            "ppd_q50": [2.0, 2.8, 1.0, 0.1],
            "ppd_q75": [2.2, 3.1, 1.2, 0.2],
            "ppd_q95": [2.5, 3.5, 1.5, 0.4],
            "residual": [-0.1, 0.2, -0.1, -0.1],
            "z_score": [-0.2, 0.3, -0.1, -0.2],
            "tail_prob": [0.4, 0.6, 0.5, 0.7],
            "two_tail_prob": [0.5, 0.5, 0.6, 0.8],
        }
    )
    ppd_chr_df = pd.DataFrame(
        {
            "sample": ["S1", "S1", "S2", "S2"],
            "chromosome": ["chr21", "chrX", "chr21", "chrY"],
            "bayesian_pvalue": [0.5, 0.6, 0.4, 0.7],
            "rmse": [0.1, 0.2, 0.3, 0.4],
        }
    )
    return ppd_bin_df, ppd_chr_df


def _small_site_af_df() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "chr": ["chr21", "chr21", "chrX", "chrY"],
            "bin_start": [0, 0, 0, 0],
            "bin_end": [100, 100, 100, 100],
            "site_index": [0, 1, 0, 0],
            "n_observed_samples": [2, 2, 2, 1],
            "sum_alt": [2, 8, 6, 1],
            "sum_total": [20, 20, 12, 8],
            "pooled_observed_af": [0.10, 0.40, 0.50, 0.125],
            "input_site_pop_af": [0.5, 0.5, 0.5, 0.5],
            "naive_bayes_site_pop_af": [0.136, 0.409, 0.500, 0.200],
            "effective_site_pop_af": [0.136, 0.409, 0.500, 0.200],
            "posterior_sd": [0.07, 0.10, 0.11, 0.12],
            "site_pop_af_changed": [True, True, False, True],
        }
    )


def test_plot_detail_helpers_and_outputs(tmp_path, tiny_site_data) -> None:
    chrom_df = _small_chrom_df()
    bin_df = _small_bin_df()

    assert _expected_cn_range("chrX") == (1, 2)
    assert _expected_cn_range("chrY") == (0, 1)
    assert _expected_cn_range("chr21") == (2, 2)
    raw_label, raw_limit = _depth_axis_label_and_limit(
        np.array([12.0, 18.0]),
        "observed_depth",
    )
    assert raw_label == "Normalized Depth"
    assert raw_limit == pytest.approx(18.9)

    cn_label, cn_limit = _depth_axis_label_and_limit(np.array([2.0, 3.0]), "cn_map")
    assert cn_label == "Estimated Copy Number"
    assert cn_limit == pytest.approx(5.0)
    assert _samples_by_chrX(chrom_df, lambda cn: cn < 2) == ["S1"]
    assert _samples_by_chrX_bin(bin_df, lambda cn: cn >= 2) == ["S2"]
    x = _equal_width_x(np.array(["chr1", "chr1", "chr2", "chr2"]), 4)
    assert list(x) == [0.0, 0.5, 1.0, 1.5]
    left_edges, widths = _contiguous_bar_geometry(x)
    np.testing.assert_allclose(left_edges, [-0.25, 0.25, 0.75, 1.25])
    np.testing.assert_allclose(widths, [0.5, 0.5, 0.5, 0.5])
    np.testing.assert_allclose(left_edges[1:], left_edges[:-1] + widths[:-1])

    site_data = dict(tiny_site_data)
    site_data["sample_ids"] = np.array(["S1", "S2"], dtype=object)
    sample_idx_map = {"S1": 0, "S2": 1}

    plot_sample_with_variance(
        bin_df[bin_df["sample"] == "S1"],
        np.array([0.09, 0.16]),
        str(tmp_path),
        site_data=site_data,
        sample_idx_map=sample_idx_map,
    )
    plot_cn_per_contig_boxplot(chrom_df, str(tmp_path), sex_subset="male", connect_samples=True, highlight_sample="S1")
    plot_cn_per_contig_boxplot(chrom_df, str(tmp_path), sex_subset="female", connect_samples=False, highlight_sample="S2")
    plot_cn_per_bin_chromosome(bin_df, str(tmp_path), "chr21", sex_subset="female", highlight_sample="S2")
    plot_sex_assignments(
        pd.DataFrame(
            {
                "sample": ["S1", "S2"],
                "sex": ["MALE", "FEMALE"],
                "chrX_depth": [1.0, 2.0],
                "chrY_depth": [1.0, 0.0],
            }
        ),
        str(tmp_path),
        highlight_sample="S1",
    )

    assert any((tmp_path / "sample_plots").iterdir())
    assert any((tmp_path / "raw_depth_contig").iterdir())
    assert any((tmp_path / "raw_depth_binned").iterdir())
    assert (tmp_path / "sex_assignments.png").exists()


def test_af_reference_lines_are_solid() -> None:
    fig, ax = plt.subplots()
    try:
        _draw_af_reference_lines(ax)

        y_values = sorted(float(line.get_ydata()[0]) for line in ax.lines)
        np.testing.assert_allclose(
            y_values,
            [0.25, 1.0 / 3.0, 0.5, 2.0 / 3.0, 0.75],
        )
        assert {line.get_linestyle() for line in ax.lines} == {"-"}
        labels = [
            line.get_label()
            for line in ax.lines
            if not line.get_label().startswith("_")
        ]
        assert labels == ["AF grid (1/4, 1/3, 1/2, 2/3, 3/4)"]
    finally:
        plt.close(fig)


def test_plot_sample_ignores_aggregate_af_without_site_data(
    tmp_path,
    monkeypatch,
) -> None:
    captured: dict[str, int] = {}
    original_subplots = plt.subplots

    def capture_subplots(nrows, *args, **kwargs):
        captured["nrows"] = nrows
        return original_subplots(nrows, *args, **kwargs)

    monkeypatch.setattr(plt, "subplots", capture_subplots)

    bin_df = _small_bin_df()
    plot_sample_with_variance(
        bin_df[bin_df["sample"] == "S1"],
        np.array([0.09, 0.16]),
        str(tmp_path),
    )

    assert captured["nrows"] == 4
    assert any((tmp_path / "sample_plots").iterdir())


def test_draw_sample_depth_panel_omits_ignored_bins() -> None:
    fig, ax = plt.subplots()
    x = np.array([0.0, 0.5, 1.0], dtype=float)
    obs = np.array([2.0, 4.0, 3.0], dtype=float)
    cn_map = np.array([2.0, 4.0, 3.0], dtype=float)
    retained_mask = np.array([True, False, True], dtype=bool)

    _draw_sample_depth_panel(
        ax,
        x,
        obs,
        cn_map,
        x_min=-0.25,
        x_max=1.25,
        retained_mask=retained_mask,
        panel_label="Retained after BINQ filter",
    )

    cn_line, obs_line = ax.lines
    np.testing.assert_allclose(cn_line.get_xdata(), [0.0, 1.0])
    np.testing.assert_allclose(cn_line.get_ydata(), [2.0, 3.0])
    np.testing.assert_allclose(obs_line.get_xdata(), [0.0, 0.5, 1.0])
    assert np.isnan(obs_line.get_ydata()[1])
    assert tuple(ax.get_xlim()) == pytest.approx((-0.25, 1.25))
    plt.close(fig)


def test_plot_sample_with_variance_omits_filtered_panel_and_filtered_af_cnq(
    tmp_path,
    tiny_site_data,
    monkeypatch,
) -> None:
    sample_data = pd.DataFrame(
        {
            "sample": ["S1", "S1"],
            "chr": ["chr21", "chr18"],
            "start": [200, 100],
            "end": [300, 200],
            "observed_depth": [2.0, 1.9],
            "cn_map": [2, 2],
            "cnq": [20.0, 40.0],
            "binq_field": ["BINQ20", "BINQ20"],
            "binq_value": [35.0, 5.0],
            "sample_var": [0.09, 0.09],
            "ignored_in_call": [False, True],
        }
    )
    site_data = dict(tiny_site_data)
    site_data["sample_ids"] = np.array(["S1", "S2"], dtype=object)
    captured: dict[str, plt.Figure] = {}

    def fake_save_publication_figure(fig, path, dpi=None, **kwargs):
        captured["fig"] = fig

    monkeypatch.setattr(
        "gatk_sv_ploidy._plot_detail.save_publication_figure",
        fake_save_publication_figure,
    )

    plot_sample_with_variance(
        sample_data,
        np.array([0.09, 0.16]),
        str(tmp_path),
        site_data=site_data,
        sample_idx_map={"S1": 0, "S2": 1},
    )

    fig = captured["fig"]
    assert len(fig.axes) == 5

    ax_depth = fig.axes[0]
    ax_af = fig.axes[1]
    ax_cnq = fig.axes[2]
    ax_binq = fig.axes[3]
    assert ax_af.get_ylabel() == "Allele fraction"
    assert all(text.get_text() != "Retained after BINQ filter" for text in ax_depth.texts)
    assert tuple(ax_cnq.get_ylim()) == pytest.approx((0.0, 105.0))
    assert tuple(ax_binq.get_ylim()) == pytest.approx((0.0, 105.0))

    af_point_counts = [
        len(collection.get_offsets())
        for collection in ax_af.collections
        if hasattr(collection, "get_offsets") and len(collection.get_offsets()) > 0
    ]
    assert sum(af_point_counts) == 2
    cnq_score_lines = [
        line for line in ax_cnq.lines if line.get_drawstyle() == "steps-mid"
    ]
    assert len(cnq_score_lines) == 1
    assert np.isfinite(cnq_score_lines[0].get_ydata()).sum() == 1


def test_plot_helpers_prefer_plot_normalized_depth_columns() -> None:
    chrom_df = _small_chrom_df()
    chrom_df["plot_mean_depth"] = chrom_df["mean_depth"] + 10.0
    chrom_df["plot_std_depth"] = chrom_df["std_depth"] + 1.0
    chrom_df["plot_median_depth"] = chrom_df["median_depth"] + 10.0
    chrom_df["plot_mad_depth"] = chrom_df["mad_depth"] + 1.0
    bin_df = _small_bin_df()
    bin_df["plot_depth"] = bin_df["observed_depth"] + 5.0
    sex_df = pd.DataFrame(
        {
            "sample": ["S1", "S2"],
            "sex": ["MALE", "FEMALE"],
            "chrX_depth": [99.0, 99.0],
            "chrY_depth": [99.0, 99.0],
        }
    )

    plot_chrom_df = _apply_plot_depth_columns(chrom_df)
    plot_bin_df = _apply_plot_depth_bin_columns(bin_df)
    plot_sex_df = _apply_plot_depth_sex_columns(sex_df, plot_chrom_df)

    assert plot_chrom_df["median_depth"].iloc[0] == 12.0
    assert plot_bin_df["observed_depth"].iloc[0] == 7.0
    assert plot_bin_df["raw_observed_depth"].iloc[0] == 2.0
    assert plot_sex_df.loc[plot_sex_df["sample"] == "S1", "chrX_depth"].iloc[0] == 11.0
    assert plot_sex_df.loc[plot_sex_df["sample"] == "S2", "chrY_depth"].iloc[0] == 10.0


def test_plot_baseline_ploidy_metadata_prefers_call_columns() -> None:
    chrom_df = _small_chrom_df()
    chrom_df["autosomal_baseline_cn"] = [2, 2, 2, 3, 3, 3]
    sex_df = pd.DataFrame(
        {
            "sample": ["S2"],
            "baseline_ploidy_type": ["TRIPLOID"],
            "autosomal_baseline_cn": [3],
        }
    )

    metadata = _sample_baseline_ploidy_metadata(chrom_df, sex_df)

    assert metadata["S1"] == {
        "baseline_ploidy_type": "DIPLOID",
        "autosomal_baseline_cn": 2,
    }
    assert metadata["S2"] == {
        "baseline_ploidy_type": "TRIPLOID",
        "autosomal_baseline_cn": 3,
    }


def test_plot_annotate_ignored_bins_and_overlay_branches(tmp_path) -> None:
    bin_df = _small_bin_df()
    annotated = _annotate_ignored_bins(bin_df, _small_ignored_bins_df())

    ignored_row = annotated[
        (annotated["sample"] == "S2") &
        (annotated["chr"] == "chr21") &
        (annotated["start"] == 0)
    ].iloc[0]
    retained_row = annotated[
        (annotated["sample"] == "S1") &
        (annotated["chr"] == "chr21") &
        (annotated["start"] == 0)
    ].iloc[0]

    assert bool(ignored_row["ignored_in_call"]) is True
    assert bool(retained_row["ignored_in_call"]) is False
    assert float(ignored_row["ignored_fraction_in_call"]) == pytest.approx(0.5)

    plot_sample_with_variance(
        annotated[annotated["sample"] == "S2"],
        np.array([0.09, 0.16]),
        str(tmp_path),
        aneuploid_chrs=[("chr21", 3, 0.95)],
    )
    plot_cn_per_bin_chromosome(
        annotated,
        str(tmp_path),
        "chr21",
        highlight_sample="S2",
    )

    assert any((tmp_path / "sample_plots").iterdir())
    assert any((tmp_path / "raw_depth_binned").iterdir())


def test_plot_annotate_binq_values_merges_quality_table() -> None:
    bin_df = _small_bin_df().drop(columns=["binq_field", "binq_value"])
    bin_quality_df = pd.DataFrame(
        {
            "chr": ["chr21", "chr21", "chrX", "chrX", "chrY", "chrY"],
            "start": [0, 100, 0, 100, 0, 100],
            "end": [100, 200, 100, 200, 100, 200],
            "BINQ20": [15.0, 35.0, 25.0, 45.0, 55.0, 65.0],
            "CALLQ20": [99.0, 98.0, 97.0, 96.0, 95.0, 94.0],
        }
    )

    annotated = _annotate_binq_values(bin_df, bin_quality_df, "BINQ20")

    assert "binq_value" in annotated.columns
    assert "binq_field" in annotated.columns
    assert annotated.loc[
        (annotated["chr"] == "chr21") & (annotated["start"] == 100),
        "binq_value",
    ].iloc[0] == pytest.approx(35.0)

    auto_annotated = _annotate_binq_values(bin_df, bin_quality_df, "auto")
    assert auto_annotated.loc[
        (auto_annotated["chr"] == "chr21") & (auto_annotated["start"] == 100),
        "binq_field",
    ].iloc[0] == "BINQ20"

    callq_annotated = _annotate_binq_values(bin_df, bin_quality_df, "CALLQ20")
    assert callq_annotated.loc[
        (callq_annotated["chr"] == "chr21") & (callq_annotated["start"] == 100),
        "binq_field",
    ].iloc[0] == "CALLQ20"


def test_plot_detail_private_branches(tmp_path) -> None:
    fig, ax = plt.subplots()
    sample_data = pd.DataFrame(
        {
            "sample": ["S1", "S1"],
            "chr": ["chr21", "chr21"],
            "start": [0, 100],
            "end": [100, 200],
        }
    )
    site_data = {
        "bin_chr": np.array(["chr21"], dtype=object),
        "bin_start": np.array([0]),
        "bin_end": np.array([100]),
        "site_alt": np.array([[[1], [0]]], dtype=np.int32),
        "site_total": np.array([[[10], [0]]], dtype=np.int32),
        "site_mask": np.array([[[True], [False]]]),
    }
    _draw_site_af_scatter(ax, sample_data, np.array([0.1, 1.1]), site_data, 0, min_het_alt=3)
    plt.close(fig)

    no_af_sample = pd.DataFrame(
        {
            "sample": ["S1", "S1"],
            "chr": ["chr21", "chr21"],
            "start": [0, 100],
            "end": [100, 200],
            "observed_depth": [3.0, 3.1],
            "cn_map": [3, 3],
            "cnq": [12, 18],
            "sample_var": [0.09, 0.09],
        }
    )
    plot_sample_with_variance(
        no_af_sample,
        np.array([0.09, 0.16]),
        str(tmp_path),
        aneuploid_chrs=[("chr21", 3, 0.95)],
        baseline_ploidy_type="TRIPLOID",
        autosomal_baseline_cn=3,
    )
    assert any((tmp_path / "sample_plots").iterdir())

    female_only = pd.DataFrame(
        [{"sample": "S2", "chromosome": "chrX", "copy_number": 2, "median_depth": 2.0}]
    )
    plot_cn_per_contig_boxplot(female_only, str(tmp_path), sex_subset="male")

    no_chr_x_bin = pd.DataFrame(
        {
            "chr": ["chr21", "chr21"],
            "start": [0, 100],
            "end": [100, 200],
            "sample": ["S1", "S2"],
            "observed_depth": [2.0, 2.1],
            "cn_map": [2, 2],
        }
    )
    assert _samples_by_chrX_bin(no_chr_x_bin, lambda cn: cn >= 2) == ["S1", "S2"]
    plot_cn_per_bin_chromosome(no_chr_x_bin, str(tmp_path), "chr1")


def test_plot_module_orchestrators_and_skip_branches(tmp_path) -> None:
    chrom_df = _small_chrom_df()
    bin_df = _annotate_ignored_bins(_small_bin_df(), _small_ignored_bins_df())
    loss_df = _small_loss_df()
    ppd_bin_df, ppd_chr_df = _small_ppd_frames()

    plot_median_depth_distributions(chrom_df.drop(columns=["median_depth", "mad_depth"]), str(tmp_path))
    _run_ploidy_plots(
        chrom_df,
        bin_df,
        pd.DataFrame({"sample": ["S1", "S2"], "sex": ["MALE", "FEMALE"], "chrX_depth": [1.0, 2.0], "chrY_depth": [1.0, 0.0]}),
        str(tmp_path),
        "S1",
    )
    _run_aneuploidy_plots(chrom_df, bin_df, loss_df, str(tmp_path), skip_per_sample=True)
    run_ppd_plots(ppd_bin_df, ppd_chr_df, str(tmp_path), highlight_sample="S1")

    assert (tmp_path / "ppd" / "ppd_obs_vs_predicted.png").exists()
    assert (tmp_path / "diagnostics" / "chromosome_plq_heatmap.png").exists()
    assert (tmp_path / "diagnostics" / "binq_genome_profile.png").exists()


def test_plot_histograms_by_chr_type_uses_curated_unique_specs(monkeypatch, tmp_path) -> None:
    chrom_df = _small_chrom_df().assign(
        sample_var_map=[0.11, 0.11, 0.11, 0.23, 0.23, 0.23],
        sample_overdispersion_map=[0.12, 0.12, 0.12, 0.24, 0.24, 0.24],
        sample_depth_map=[1.0, 1.0, 1.0, 1.5, 1.5, 1.5],
        n_bins_retained=[8, 8, 8, 6, 6, 6],
        frac_bins_retained=[0.8, 0.8, 0.8, 0.6, 0.6, 0.6],
        autosomal_baseline_cn=[2, 2, 2, 3, 3, 3],
        sample_depth_ratio=[1.0, 1.0, 1.0, 1.45, 1.45, 1.45],
        global_cn_scale_factor=[1.0, 1.0, 1.0, 1.5, 1.5, 1.5],
    )

    calls: list[dict[str, object]] = []

    def _capture_hist(*args, **kwargs):
        calls.append(
            {
                "x_col": args[1],
                "filename": args[6],
                "bins": kwargs.get("bins"),
            }
        )

    monkeypatch.setattr("gatk_sv_ploidy.plot._hist_by_hue", _capture_hist)

    plot_histograms_by_chr_type(chrom_df, str(tmp_path), highlight_sample="S1")

    emitted_files = {str(call["filename"]) for call in calls}
    assert emitted_files == {
        "hist_sample_overdispersion.png",
        "hist_sample_depth_map.png",
        "hist_copy_number.png",
        "hist_plq.png",
        "hist_n_bins_retained.png",
        "hist_frac_bins_retained.png",
        "hist_autosomal_baseline_cn.png",
        "hist_sample_depth_ratio.png",
        "hist_global_cn_scale_factor.png",
    }
    assert "hist_sample_var_map.png" not in emitted_files
    assert all(int(call["bins"]) == 30 for call in calls)

    overdispersion_call = next(
        call for call in calls if call["filename"] == "hist_sample_overdispersion.png"
    )
    assert overdispersion_call["x_col"] == "sample_overdispersion_map"


def test_plot_report_generates_manifest_and_html(tmp_path) -> None:
    diagnostics = tmp_path / "diagnostics"
    ppd_dir = tmp_path / "ppd"
    sample_dir = tmp_path / "sample_plots"
    stale_figures_dir = tmp_path / "figures"
    stale_tables_dir = tmp_path / "tables"
    diagnostics.mkdir()
    ppd_dir.mkdir()
    sample_dir.mkdir()
    stale_figures_dir.mkdir()
    stale_tables_dir.mkdir()
    (stale_figures_dir / "stale.txt").write_text("stale", encoding="utf-8")
    (stale_tables_dir / "stale.txt").write_text("stale", encoding="utf-8")

    (tmp_path / "sex_assignments.png").write_bytes(b"plot")
    (diagnostics / "training_loss_gradient.png").write_bytes(b"plot")
    (ppd_dir / "ppd_obs_vs_predicted.png").write_bytes(b"plot")
    (sample_dir / "S1.png").write_bytes(b"plot")
    (tmp_path / "median_depth_distribution_stats.tsv").write_text(
        "sample\tchromosome\nS1\tchr21\n",
        encoding="utf-8",
    )

    manifest_df = write_plot_report(
        str(tmp_path),
        chrom_df=_small_chrom_df(),
        bin_df=_small_bin_df(),
        sex_df=pd.DataFrame({"sample": ["S1"], "sex": ["MALE"]}),
        highlight_sample="S1",
    )

    assert (tmp_path / "report" / "index.html").exists()
    assert (tmp_path / "plot_manifest.tsv").exists()
    assert not (tmp_path / "figures").exists()
    assert not (tmp_path / "tables").exists()

    assert set(manifest_df["category"]) >= {"overview", "model", "ppd", "per_sample", "tables"}
    training_loss_row = manifest_df.loc[
        manifest_df["source_path"] == "diagnostics/training_loss_gradient.png"
    ].iloc[0]
    assert training_loss_row["report_path"] == "../diagnostics/training_loss_gradient.png"
    html = (tmp_path / "report" / "index.html").read_text(encoding="utf-8")
    assert "GATK-SV Ploidy Report" in html
    assert "Highlighted sample" in html
    assert "../diagnostics/training_loss_gradient.png" in html


def test_apply_theme_uses_nature_publication_defaults() -> None:
    apply_theme()
    assert plt.rcParams["font.size"] == 8
    assert plt.rcParams["axes.labelsize"] == 8
    assert plt.rcParams["xtick.labelsize"] == 8
    assert plt.rcParams["legend.fontsize"] == 8
    assert plt.rcParams["axes.grid"] is False
    assert plt.rcParams["pdf.fonttype"] == 42
    assert plt.rcParams["svg.fonttype"] == "none"


def test_plot_style_suppresses_fonttools_info_logging() -> None:
    fonttools_logger = logging.getLogger("fontTools")
    fonttools_subset_logger = logging.getLogger("fontTools.subset")
    original_levels = (fonttools_logger.level, fonttools_subset_logger.level)

    try:
        fonttools_logger.setLevel(logging.INFO)
        fonttools_subset_logger.setLevel(logging.NOTSET)

        _suppress_chatty_plot_dependency_loggers()

        assert fonttools_logger.level == logging.WARNING
        assert fonttools_subset_logger.level == logging.WARNING
    finally:
        fonttools_logger.setLevel(original_levels[0])
        fonttools_subset_logger.setLevel(original_levels[1])


def test_run_aneuploidy_plots_passes_non_diploid_baseline_metadata(tmp_path, monkeypatch) -> None:
    chrom_df = pd.concat(
        [
            _small_chrom_df(),
            pd.DataFrame(
                [
                    {"sample": "S3", "chromosome": "chr21", "copy_number": 3, "mean_cn_probability": 0.95, "plq": 60, "is_aneuploid": False, "mean_depth": 3.0, "median_depth": 3.0, "std_depth": 0.1, "mad_depth": 0.1, "sex": "TRIPLOID_MALE", "autosomal_baseline_cn": 3},
                    {"sample": "S3", "chromosome": "chrX", "copy_number": 2, "mean_cn_probability": 0.95, "plq": 55, "is_aneuploid": False, "mean_depth": 2.0, "median_depth": 2.0, "std_depth": 0.1, "mad_depth": 0.1, "sex": "TRIPLOID_MALE", "autosomal_baseline_cn": 3},
                    {"sample": "S3", "chromosome": "chrY", "copy_number": 1, "mean_cn_probability": 0.95, "plq": 50, "is_aneuploid": False, "mean_depth": 1.0, "median_depth": 1.0, "std_depth": 0.1, "mad_depth": 0.1, "sex": "TRIPLOID_MALE", "autosomal_baseline_cn": 3},
                ]
            ),
        ],
        ignore_index=True,
    )
    chrom_df["autosomal_baseline_cn"] = chrom_df["autosomal_baseline_cn"].fillna(2)
    bin_df = pd.concat(
        [
            _small_bin_df(),
            pd.DataFrame(
                [
                    {
                        "chr": chrom,
                        "start": start,
                        "end": start + 100,
                        "sample": "S3",
                        "observed_depth": float(cn),
                        "cn_map": cn,
                        "cnq": 20,
                        "binq_value": 40.0,
                        "sample_var": 0.2,
                    }
                    for chrom, cn in [("chr21", 3), ("chrX", 2), ("chrY", 1)]
                    for start in [0, 100]
                ]
            ),
        ],
        ignore_index=True,
    )
    sex_df = pd.DataFrame(
        {
            "sample": ["S3"],
            "baseline_ploidy_type": ["TRIPLOID"],
            "autosomal_baseline_cn": [3],
        }
    )
    captured: dict[str, dict[str, object]] = {}

    for fn_name in (
        "plot_training_loss_with_gradient",
        "plot_bin_variance_bias",
        "plot_chromosome_cn_heatmap",
        "plot_chromosome_plq_heatmap",
        "plot_binq_genome_profile",
    ):
        monkeypatch.setattr(
            f"gatk_sv_ploidy.plot.{fn_name}",
            lambda *args, **kwargs: None,
        )

    def fake_plot_sample_with_variance(sample_data, *args, **kwargs):
        captured[str(sample_data["sample"].iloc[0])] = kwargs

    monkeypatch.setattr(
        "gatk_sv_ploidy.plot.plot_sample_with_variance",
        fake_plot_sample_with_variance,
    )

    _run_aneuploidy_plots(
        chrom_df,
        bin_df,
        _small_loss_df(),
        str(tmp_path),
        sex_df=sex_df,
    )

    assert captured["S3"]["baseline_ploidy_type"] == "TRIPLOID"
    assert captured["S3"]["autosomal_baseline_cn"] == 3
    assert captured["S3"]["aneuploid_chrs"] == []


def test_plot_site_af_estimates_outputs(tmp_path) -> None:
    plot_site_af_estimates(_small_site_af_df(), str(tmp_path))

    assert (tmp_path / "diagnostics" / "site_af_effective_distribution.png").exists()
    assert (tmp_path / "diagnostics" / "site_af_pooled_vs_naive_bayes.png").exists()
    assert (tmp_path / "diagnostics" / "site_af_input_vs_effective.png").exists()


def test_plot_public_helpers_cover_remaining_branches(tmp_path) -> None:
    chrom_df = _small_chrom_df()
    plot_chromosome_plq_heatmap(chrom_df, str(tmp_path))
    assert (tmp_path / "diagnostics" / "chromosome_plq_heatmap.png").exists()

    plot_binq_genome_profile(_small_bin_df(), str(tmp_path))
    assert (tmp_path / "diagnostics" / "binq_genome_profile.png").exists()

    short_loss = pd.DataFrame({"epoch": [0, 1], "elbo": [100.0, 90.0]})
    plot_training_loss_with_gradient(short_loss, str(tmp_path), window=50)
    assert (tmp_path / "diagnostics" / "training_loss_gradient.png").exists()

    big_df = pd.DataFrame(
        {
            "sample": [f"S{i}" for i in range(61)],
            "chromosome": ["chr21"] * 61,
            "copy_number": [2] * 61,
            "mean_cn_probability": [0.9] * 61,
        }
    )
    plot_chromosome_cn_heatmap(big_df, str(tmp_path))
    assert (tmp_path / "diagnostics" / "chromosome_cn_heatmap.png").exists()


def test_plot_parse_args_and_warning_only_paths(tmp_path, monkeypatch) -> None:
    chrom_stats = tmp_path / "chrom.tsv"
    _small_chrom_df().to_csv(chrom_stats, sep="\t", index=False)
    bin_stats = tmp_path / "bin.tsv.gz"
    _small_bin_df().to_csv(bin_stats, sep="\t", index=False, compression="gzip")
    ppd_bin_quality = tmp_path / "ppd_bin_quality.tsv"
    _small_bin_df()[["chr", "start", "end", "binq_value"]].drop_duplicates().rename(
        columns={"binq_value": "BINQ20"}
    ).to_csv(ppd_bin_quality, sep="\t", index=False)
    ignored_bins = tmp_path / "ignored_bins.tsv.gz"
    _small_ignored_bins_df().to_csv(ignored_bins, sep="\t", index=False, compression="gzip")
    site_af_estimates = tmp_path / "site_af_estimates.tsv.gz"
    _small_site_af_df().to_csv(site_af_estimates, sep="\t", index=False, compression="gzip")
    ppd_bin_df, _ = _small_ppd_frames()
    ppd_bin = tmp_path / "ppd_bin.tsv.gz"
    ppd_bin_df.to_csv(ppd_bin, sep="\t", index=False, compression="gzip")

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "gatk-sv-ploidy plot",
            "--chrom-stats",
            str(chrom_stats),
            "--output-dir",
            str(tmp_path / "out"),
            "--bin-stats",
            str(bin_stats),
            "--ppd-bin-quality",
            str(ppd_bin_quality),
            "--ignored-bins",
            str(ignored_bins),
            "--site-af-estimates",
            str(site_af_estimates),
            "--ppd-bin-summary",
            str(ppd_bin),
        ],
    )
    args = parse_args()
    assert args.bin_stats.endswith("bin.tsv.gz")
    assert args.ppd_bin_quality.endswith("ppd_bin_quality.tsv")
    assert args.site_af_estimates.endswith("site_af_estimates.tsv.gz")

    main()
    assert (tmp_path / "out" / "sex_assignments.png").exists()
    assert (tmp_path / "out" / "diagnostics" / "site_af_input_vs_effective.png").exists()
    assert (tmp_path / "out" / "report" / "index.html").exists()
    assert (tmp_path / "out" / "plot_manifest.tsv").exists()


def test_plot_main_writes_pdf_outputs_when_requested(tmp_path, monkeypatch) -> None:
    chrom_stats = tmp_path / "chrom.tsv"
    _small_chrom_df().to_csv(chrom_stats, sep="\t", index=False)

    original_format = get_plot_output_format()
    try:
        monkeypatch.setattr(
            sys,
            "argv",
            [
                "gatk-sv-ploidy plot",
                "--chrom-stats",
                str(chrom_stats),
                "--output-dir",
                str(tmp_path / "out_pdf"),
                "--pdf",
            ],
        )
        main()
    finally:
        set_plot_output_format(original_format)

    assert (tmp_path / "out_pdf" / "sex_assignments.pdf").exists()
    assert not (tmp_path / "out_pdf" / "sex_assignments.png").exists()
    assert (tmp_path / "out_pdf" / "report" / "index.html").exists()
    assert not (tmp_path / "out_pdf" / "figures").exists()
    html = (tmp_path / "out_pdf" / "report" / "index.html").read_text(encoding="utf-8")
    assert "../sex_assignments.pdf" in html


def test_plot_main_runs_ppd_before_sample_plot_branches(tmp_path, monkeypatch) -> None:
    chrom_stats = tmp_path / "chrom.tsv"
    _small_chrom_df().to_csv(chrom_stats, sep="\t", index=False)

    bin_stats = tmp_path / "bin.tsv.gz"
    _small_bin_df().to_csv(bin_stats, sep="\t", index=False, compression="gzip")

    training_loss = tmp_path / "training_loss.tsv"
    _small_loss_df().to_csv(training_loss, sep="\t", index=False)

    sex_assignments = tmp_path / "sex.tsv"
    pd.DataFrame(
        {
            "sample": ["S1", "S2"],
            "sex": ["MALE", "FEMALE"],
            "chrX_depth": [1.0, 2.0],
            "chrY_depth": [1.0, 0.0],
        }
    ).to_csv(sex_assignments, sep="\t", index=False)

    ppd_bin_df, ppd_chr_df = _small_ppd_frames()
    ppd_bin = tmp_path / "ppd_bin.tsv.gz"
    ppd_chr = tmp_path / "ppd_chr.tsv"
    ppd_bin_df.to_csv(ppd_bin, sep="\t", index=False, compression="gzip")
    ppd_chr_df.to_csv(ppd_chr, sep="\t", index=False)

    call_order: list[str] = []

    monkeypatch.setattr(
        "gatk_sv_ploidy.plot.plot_histograms_by_chr_type",
        lambda *args, **kwargs: call_order.append("hist"),
    )
    monkeypatch.setattr(
        "gatk_sv_ploidy.plot.plot_median_depth_distributions",
        lambda *args, **kwargs: call_order.append("median_depth"),
    )
    monkeypatch.setattr(
        "gatk_sv_ploidy.plot.run_ppd_plots",
        lambda *args, **kwargs: call_order.append("ppd"),
    )
    monkeypatch.setattr(
        "gatk_sv_ploidy.plot._run_ploidy_plots",
        lambda *args, **kwargs: call_order.append("ploidy"),
    )
    monkeypatch.setattr(
        "gatk_sv_ploidy.plot._run_aneuploidy_plots",
        lambda *args, **kwargs: call_order.append("sample_plots"),
    )

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "gatk-sv-ploidy plot",
            "--chrom-stats",
            str(chrom_stats),
            "--bin-stats",
            str(bin_stats),
            "--training-loss",
            str(training_loss),
            "--sex-assignments",
            str(sex_assignments),
            "--ppd-bin-summary",
            str(ppd_bin),
            "--ppd-chr-summary",
            str(ppd_chr),
            "--output-dir",
            str(tmp_path / "out"),
        ],
    )

    main()

    assert "ppd" in call_order
    assert "ploidy" in call_order
    assert "sample_plots" in call_order
    assert call_order.index("ppd") < call_order.index("ploidy")
    assert call_order.index("ppd") < call_order.index("sample_plots")


def test_plot_background_factor_diagnostics_outputs(tmp_path) -> None:
    bin_df = _small_bin_df()
    map_estimates = {
        "background_bin_factors": np.array(
            [
                [0.8, 0.3],
                [0.9, 0.2],
                [1.0, 0.4],
                [1.1, 0.5],
                [2.0, 1.6],
                [2.2, 1.8],
            ],
            dtype=np.float32,
        ),
        "background_sample_factors": np.array(
            [
                [0.02, 0.06],
                [0.01, 0.04],
            ],
            dtype=np.float32,
        ),
        "multiplicative_bin_factors": np.array(
            [
                [-0.4, 0.2],
                [-0.2, 0.1],
                [0.0, -0.1],
                [0.1, -0.2],
                [0.3, 0.4],
                [0.5, 0.6],
            ],
            dtype=np.float32,
        ),
        "multiplicative_sample_factors": np.array(
            [
                [-0.05, 0.08],
                [0.02, -0.03],
            ],
            dtype=np.float32,
        ),
    }

    plot_background_factor_diagnostics(bin_df, map_estimates, str(tmp_path))

    assert (tmp_path / "diagnostics" / "background_factor_diagnostics.png").exists()
    assert (tmp_path / "diagnostics" / "factor_mode_weight_diagnostics.png").exists()
