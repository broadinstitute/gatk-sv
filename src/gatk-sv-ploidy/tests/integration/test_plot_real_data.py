from __future__ import annotations

import sys

import pandas as pd
import pytest

from gatk_sv_ploidy import call, infer, plot, ppd

pytestmark = pytest.mark.integration


def test_plot_pipeline_runs_on_medium_fixture(
    medium_fixture_root,
    tmp_path,
    monkeypatch,
) -> None:
    infer_out = tmp_path / "infer"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "gatk-sv-ploidy infer",
            "--input",
            str(medium_fixture_root / "preprocessed_depth.tsv"),
            "--site-data",
            str(medium_fixture_root / "site_data.npz"),
            "--output-dir",
            str(infer_out),
            "--device",
            "cpu",
            "--max-iter",
            "10",
            "--log-freq",
            "5",
            "--no-early-stopping",
        ],
    )
    infer.main()

    call_out = tmp_path / "call"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "gatk-sv-ploidy call",
            "--chrom-stats",
            str(infer_out / "chromosome_stats.tsv"),
            "--truth-json",
            str(medium_fixture_root / "truth.json"),
            "--output-dir",
            str(call_out),
        ],
    )
    call.main()

    ppd_out = tmp_path / "ppd"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "gatk-sv-ploidy ppd",
            "--input",
            str(medium_fixture_root / "preprocessed_depth.tsv"),
            "--site-data",
            str(medium_fixture_root / "site_data.npz"),
            "--artifacts",
            str(infer_out / "inference_artifacts.npz"),
            "--output-dir",
            str(ppd_out),
            "--draws",
            "4",
        ],
    )
    ppd.main()

    plot_out = tmp_path / "plot"
    highlight_sample = "HG00150"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "gatk-sv-ploidy plot",
            "--chrom-stats",
            str(infer_out / "chromosome_stats.tsv"),
            "--bin-stats",
            str(infer_out / "bin_stats.tsv.gz"),
            "--training-loss",
            str(infer_out / "training_loss.tsv"),
            "--sex-assignments",
            str(call_out / "aneuploidy_type_predictions.tsv"),
            "--site-data",
            str(medium_fixture_root / "site_data.npz"),
            "--ppd-bin-summary",
            str(ppd_out / "ppd_bin_summary.tsv.gz"),
            "--ppd-chr-summary",
            str(ppd_out / "ppd_chromosome_summary.tsv"),
            "--highlight-sample",
            highlight_sample,
            "--output-dir",
            str(plot_out),
        ],
    )
    plot.main()

    expected_files = [
        plot_out / "sex_assignments.png",
        plot_out / "median_depth_distributions.pdf",
        plot_out / "diagnostics" / "training_loss.png",
        plot_out / "diagnostics" / "training_loss_gradient.png",
        plot_out / "diagnostics" / "bin_posteriors.png",
        plot_out / "diagnostics" / "cn_posterior_entropy.png",
        plot_out / "diagnostics" / "chromosome_cn_heatmap.png",
        plot_out / "diagnostics" / "parameter_diagnostics.png",
        plot_out / "ppd" / "ppd_obs_vs_predicted.png",
        plot_out / "ppd" / "ppd_residual_genome.png",
        plot_out / "ppd" / f"ppd_residual_{highlight_sample}.png",
        plot_out / "ppd" / "ppd_calibration_histogram.png",
        plot_out / "ppd" / "ppd_qq_plot.png",
        plot_out / "ppd" / "ppd_interval_coverage.png",
        plot_out / "ppd" / "ppd_bayesian_pvalue_heatmap.png",
        plot_out / "ppd" / "ppd_chromosome_rmse.png",
    ]
    for path in expected_files:
        assert path.exists(), path

    assert any((plot_out / "sample_plots").iterdir())
    assert any((plot_out / "raw_depth_contig").iterdir())
    assert any((plot_out / "raw_depth_binned").iterdir())

    stats_df = pd.read_csv(plot_out / "median_depth_distribution_stats.tsv", sep="\t")
    assert stats_df["sample"].nunique() == 6
