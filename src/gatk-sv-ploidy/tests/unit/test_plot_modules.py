from __future__ import annotations

import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from gatk_sv_ploidy._plot_detail import (
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
from gatk_sv_ploidy.plot import (
    _run_aneuploidy_plots,
    _run_ploidy_plots,
    parse_args,
    plot_aneuploid_histograms,
    plot_chromosome_cn_heatmap,
    plot_median_depth_distributions,
    plot_training_loss_with_gradient,
    main,
)


def _small_chrom_df() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"sample": "S1", "chromosome": "chr21", "copy_number": 2, "mean_cn_probability": 0.9, "is_aneuploid": False, "mean_depth": 2.0, "median_depth": 2.0, "std_depth": 0.1, "mad_depth": 0.1, "chrX_depth": 1.0, "chrY_depth": 1.0, "sex": "MALE"},
            {"sample": "S1", "chromosome": "chrX", "copy_number": 1, "mean_cn_probability": 0.9, "is_aneuploid": False, "mean_depth": 1.0, "median_depth": 1.0, "std_depth": 0.1, "mad_depth": 0.1, "chrX_depth": 1.0, "chrY_depth": 1.0, "sex": "MALE"},
            {"sample": "S1", "chromosome": "chrY", "copy_number": 1, "mean_cn_probability": 0.9, "is_aneuploid": False, "mean_depth": 1.0, "median_depth": 1.0, "std_depth": 0.1, "mad_depth": 0.1, "chrX_depth": 1.0, "chrY_depth": 1.0, "sex": "MALE"},
            {"sample": "S2", "chromosome": "chr21", "copy_number": 3, "mean_cn_probability": 0.95, "is_aneuploid": True, "mean_depth": 3.0, "median_depth": 3.0, "std_depth": 0.2, "mad_depth": 0.1, "chrX_depth": 2.0, "chrY_depth": 0.0, "sex": "FEMALE"},
            {"sample": "S2", "chromosome": "chrX", "copy_number": 2, "mean_cn_probability": 0.95, "is_aneuploid": False, "mean_depth": 2.0, "median_depth": 2.0, "std_depth": 0.1, "mad_depth": 0.1, "chrX_depth": 2.0, "chrY_depth": 0.0, "sex": "FEMALE"},
            {"sample": "S2", "chromosome": "chrY", "copy_number": 0, "mean_cn_probability": 0.95, "is_aneuploid": False, "mean_depth": 0.0, "median_depth": 0.0, "std_depth": 0.1, "mad_depth": 0.1, "chrX_depth": 2.0, "chrY_depth": 0.0, "sex": "FEMALE"},
        ]
    )


def _small_bin_df() -> pd.DataFrame:
    rows = []
    for sample, x_cn, y_cn, chr21_cn in [("S1", 1, 1, 2), ("S2", 2, 0, 3)]:
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
                        "bin_bias": 1.0,
                        "bin_var": 0.04,
                        "sample_var": 0.09 if sample == "S1" else 0.16,
                        "mean_het_af": 0.5,
                        "n_het_sites": 2,
                        **{f"cn_prob_{i}": probs[i] for i in range(6)},
                    }
                )
    return pd.DataFrame(rows)


def _small_loss_df() -> pd.DataFrame:
    return pd.DataFrame({"epoch": [0, 1, 2], "elbo": [100.0, 90.0, 85.0]})


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


def test_plot_detail_helpers_and_outputs(tmp_path, tiny_site_data) -> None:
    chrom_df = _small_chrom_df()
    bin_df = _small_bin_df()

    assert _expected_cn_range("chrX") == (1, 2)
    assert _expected_cn_range("chrY") == (0, 1)
    assert _expected_cn_range("chr21") == (2, 2)
    assert _samples_by_chrX(chrom_df, lambda cn: cn < 2) == ["S1"]
    assert _samples_by_chrX_bin(bin_df, lambda cn: cn >= 2) == ["S2"]
    x = _equal_width_x(np.array(["chr1", "chr1", "chr2", "chr2"]), 4)
    assert list(x) == [0.0, 0.5, 1.0, 1.5]

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
    bin_df = _small_bin_df()
    loss_df = _small_loss_df()
    ppd_bin_df, ppd_chr_df = _small_ppd_frames()

    plot_aneuploid_histograms(chrom_df[~chrom_df["is_aneuploid"]], str(tmp_path))
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
    assert (tmp_path / "diagnostics" / "parameter_diagnostics.png").exists()


def test_plot_public_helpers_cover_remaining_branches(tmp_path) -> None:
    chrom_df = _small_chrom_df()
    plot_aneuploid_histograms(chrom_df, str(tmp_path), highlight_sample="S2")
    assert (tmp_path / "diagnostics" / "hist_aneuploid_chromosome.png").exists()

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
            "--ppd-bin-summary",
            str(ppd_bin),
        ],
    )
    args = parse_args()
    assert args.bin_stats.endswith("bin.tsv.gz")

    main()
    assert (tmp_path / "out" / "sex_assignments.png").exists()
