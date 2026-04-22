from __future__ import annotations
import numpy as np
import pandas as pd
import pytest
import torch

from gatk_sv_ploidy._util import compute_cnq_from_probabilities
from gatk_sv_ploidy.data import DepthData
from gatk_sv_ploidy.infer import (
    _inference_tensor_dtype,
    _resolve_observation_model,
    build_site_af_estimates,
    build_bin_stats,
    build_chromosome_stats,
    estimate_site_pop_af_naive_bayes,
    is_fallback_minor_site_data,
    parse_args,
    resolve_site_af_estimator_application,
    should_apply_site_af_estimator,
)


def test_infer_parse_args_defaults(monkeypatch) -> None:
    monkeypatch.setattr(
        "sys.argv",
        [
            "infer",
            "--input",
            "depth.tsv",
            "--output-dir",
            "outdir",
        ],
    )

    args = parse_args()

    assert args.autosome_prior_mode == "shrinkage"
    assert args.var_bias_bin == 0.05
    assert args.alpha_ref == 50.0
    assert args.epsilon_mean == 1e-2
    assert args.depth_space == "auto"
    assert args.obs_likelihood == "auto"
    assert args.sample_depth_max == 10000.0
    assert args.af_evidence_mode == "absolute"
    assert args.cn_inference_method == "multi-draw"
    assert args.cn_inference_draws == 100
    assert args.site_af_estimator == "auto"
    assert args.site_af_prior_alpha == 1.0
    assert args.site_af_prior_beta == 1.0


def test_infer_parse_args_accepts_relative_af_evidence_mode(monkeypatch) -> None:
    monkeypatch.setattr(
        "sys.argv",
        [
            "infer",
            "--input",
            "depth.tsv",
            "--output-dir",
            "outdir",
            "--af-evidence-mode",
            "relative",
        ],
    )

    args = parse_args()

    assert args.af_evidence_mode == "relative"


def test_inference_tensor_dtype_always_uses_float64() -> None:
    assert _inference_tensor_dtype() == torch.float64


def test_compute_cnq_from_probabilities_phred_scales_top_two_gap() -> None:
    probs = np.array(
        [
            [[0.90, 0.10, 0.0], [0.55, 0.45, 0.0]],
            [[1.00, 0.00, 0.0], [0.34, 0.33, 0.33]],
        ],
        dtype=np.float64,
    )

    cnq = compute_cnq_from_probabilities(probs)

    expected = np.array(
        [
            [7, 0],
            [99, 0],
        ],
        dtype=np.int16,
    )
    np.testing.assert_array_equal(cnq, expected)


def test_build_bin_stats_includes_cnq(tiny_depth_df) -> None:
    data = DepthData(tiny_depth_df, device="cpu")
    cn_posterior = np.zeros((data.n_bins, data.n_samples, 6), dtype=np.float32)
    cn_posterior[:, :, 2] = 0.85
    cn_posterior[:, :, 1] = 0.15

    bin_df = build_bin_stats(
        data,
        map_estimates={
            "bin_bias": np.ones(data.n_bins, dtype=np.float32),
            "bin_epsilon": np.full(data.n_bins, 0.01, dtype=np.float32),
            "bin_var": np.full(data.n_bins, 0.04, dtype=np.float32),
            "sample_var": np.full(data.n_samples, 0.09, dtype=np.float32),
        },
        cn_posterior={
            "cn_posterior": cn_posterior,
            "cnq": compute_cnq_from_probabilities(cn_posterior),
            "cn_map_stability": np.full((data.n_bins, data.n_samples), 0.8, dtype=np.float32),
        },
    )

    assert "cnq" in bin_df.columns
    assert "cn_map_stability" in bin_df.columns
    assert "bin_epsilon" in bin_df.columns
    assert bin_df["cnq"].between(0, 99).all()
    assert set(bin_df["cnq"].unique()) == {5}
    assert bin_df["cn_map_stability"].unique() == pytest.approx([0.8])


def test_estimate_site_pop_af_naive_bayes_matches_beta_posterior_mean() -> None:
    site_alt = np.array([[[1, 0], [2, 0]]], dtype=np.int32)
    site_total = np.array([[[10, 10], [4, 0]]], dtype=np.int32)
    site_mask = np.array([[[True, True], [True, False]]], dtype=bool)

    result = estimate_site_pop_af_naive_bayes(
        site_alt,
        site_total,
        site_mask,
        prior_alpha=1.0,
        prior_beta=1.0,
    )

    np.testing.assert_allclose(
        result["posterior_mean"],
        np.array([[2.0 / 22.0, 3.0 / 6.0]], dtype=np.float32),
        rtol=1e-6,
        atol=1e-6,
    )
    np.testing.assert_array_equal(
        result["n_observed_samples"],
        np.array([[2, 1]], dtype=np.int32),
    )
    np.testing.assert_array_equal(
        result["sum_total"],
        np.array([[20, 4]], dtype=np.int64),
    )


def test_should_apply_site_af_estimator_auto_is_conservative() -> None:
    site_mask = np.array([[[True, True], [True, False]]], dtype=bool)
    fallback = np.array([[0.5, 0.5]], dtype=np.float32)
    known = np.array([[0.25, 0.5]], dtype=np.float32)

    assert should_apply_site_af_estimator("auto", fallback, site_mask) is False
    assert should_apply_site_af_estimator("auto", known, site_mask) is False
    assert should_apply_site_af_estimator("naive-bayes", known, site_mask) is True
    assert should_apply_site_af_estimator("off", fallback, site_mask) is False


def test_is_fallback_minor_site_data_detects_uniform_half_af_inputs() -> None:
    site_mask = np.array([[[True, True], [True, False]]], dtype=bool)
    fallback = np.array([[0.5, 0.5]], dtype=np.float32)
    known = np.array([[0.25, 0.5]], dtype=np.float32)

    assert is_fallback_minor_site_data(fallback, site_mask) is True
    assert is_fallback_minor_site_data(known, site_mask) is False


def test_resolve_site_af_estimator_application_rejects_fallback_naive_bayes() -> None:
    site_mask = np.array([[[True, True], [True, False]]], dtype=bool)
    fallback = np.array([[0.5, 0.5]], dtype=np.float32)

    assert resolve_site_af_estimator_application("auto", fallback, site_mask) is False
    assert resolve_site_af_estimator_application("off", fallback, site_mask) is False
    with pytest.raises(ValueError, match="legacy fallback site_data"):
        resolve_site_af_estimator_application("naive-bayes", fallback, site_mask)


def test_build_site_af_estimates_records_effective_afs(
    tiny_depth_df: pd.DataFrame,
    tiny_site_data: dict[str, np.ndarray],
) -> None:
    data = DepthData(tiny_depth_df, device="cpu", site_data=tiny_site_data)
    input_af = tiny_site_data["site_pop_af"].astype(np.float32)
    naive_af = np.full_like(input_af, 0.125, dtype=np.float32)
    effective_af = naive_af.copy()
    posterior_sd = np.full_like(input_af, 0.05, dtype=np.float32)
    pooled_observed_af = np.full_like(input_af, 0.1, dtype=np.float32)
    sum_alt = np.ones_like(input_af, dtype=np.int64)
    sum_total = np.full_like(input_af, 8, dtype=np.int64)
    n_observed_samples = data.site_mask.cpu().numpy().sum(axis=2).astype(np.int32)

    af_df = build_site_af_estimates(
        data,
        input_site_pop_af=input_af,
        naive_bayes_site_pop_af=naive_af,
        effective_site_pop_af=effective_af,
        posterior_sd=posterior_sd,
        pooled_observed_af=pooled_observed_af,
        sum_alt=sum_alt,
        sum_total=sum_total,
        n_observed_samples=n_observed_samples,
    )

    assert not af_df.empty
    assert "effective_site_pop_af" in af_df.columns
    assert af_df["site_pop_af_changed"].all()
    assert af_df["effective_site_pop_af"].iloc[0] == pytest.approx(0.125)


def test_build_bin_stats_includes_sample_depth_for_raw_model() -> None:
    depth_df = DepthData(
        pd.DataFrame(
            {
                "Chr": ["chr21"],
                "Start": [0],
                "End": [1000],
                "S1": np.array([20], dtype=np.float32),
            },
            index=["chr21:0-1000"],
        ),
        depth_space="raw",
        clamp_threshold=None,
        device="cpu",
    )
    cn_posterior = np.zeros((depth_df.n_bins, depth_df.n_samples, 6), dtype=np.float32)
    cn_posterior[:, :, 2] = 1.0

    bin_df = build_bin_stats(
        depth_df,
        map_estimates={
            "bin_bias": np.ones(depth_df.n_bins, dtype=np.float32),
            "bin_var": np.full(depth_df.n_bins, 0.02, dtype=np.float32),
            "sample_var": np.full(depth_df.n_samples, 0.03, dtype=np.float32),
            "sample_depth": np.full(depth_df.n_samples, 10.0, dtype=np.float32),
        },
        cn_posterior={"cn_posterior": cn_posterior},
    )

    assert "sample_depth" in bin_df.columns
    assert "bin_length_kb" in bin_df.columns
    assert "plot_depth" in bin_df.columns
    assert bin_df["plot_depth"].iloc[0] == 2.0
    assert "sample_overdispersion" in bin_df.columns


def test_build_bin_stats_includes_mean_site_pop_af(
    tiny_depth_df: pd.DataFrame,
    tiny_site_data: dict[str, np.ndarray],
) -> None:
    data = DepthData(tiny_depth_df, device="cpu", site_data=tiny_site_data)
    cn_posterior = np.zeros((data.n_bins, data.n_samples, 6), dtype=np.float32)
    cn_posterior[:, :, 2] = 1.0

    bin_df = build_bin_stats(
        data,
        map_estimates={
            "bin_bias": np.ones(data.n_bins, dtype=np.float32),
            "bin_var": np.full(data.n_bins, 0.02, dtype=np.float32),
            "sample_var": np.full(data.n_samples, 0.03, dtype=np.float32),
        },
        cn_posterior={"cn_posterior": cn_posterior},
    )

    assert "mean_site_pop_af" in bin_df.columns
    assert bin_df["mean_site_pop_af"].notna().any()


def test_build_chromosome_stats_includes_plot_normalized_depth_for_raw_model() -> None:
    depth_df = DepthData(
        pd.DataFrame(
            {
                "Chr": ["chr21", "chr21"],
                "Start": [0, 1000],
                "End": [1000, 2000],
                "S1": np.array([10, 10], dtype=np.float32),
            },
            index=["chr21:0-1000", "chr21:1000-2000"],
        ),
        depth_space="raw",
        clamp_threshold=None,
        device="cpu",
    )
    cn_posterior = np.zeros((depth_df.n_bins, depth_df.n_samples, 6), dtype=np.float32)
    cn_posterior[:, :, 2] = 1.0

    chr_df = build_chromosome_stats(
        depth_df,
        map_estimates={
            "bin_bias": np.ones(depth_df.n_bins, dtype=np.float32),
            "bin_var": np.full(depth_df.n_bins, 0.02, dtype=np.float32),
            "sample_var": np.full(depth_df.n_samples, 0.03, dtype=np.float32),
            "sample_depth": np.full(depth_df.n_samples, 10.0, dtype=np.float32),
        },
        cn_posterior={"cn_posterior": cn_posterior},
        aneuploid_map={0: []},
    )

    assert "plot_mean_depth" in chr_df.columns
    assert "plot_median_depth" in chr_df.columns
    assert "plq" in chr_df.columns
    assert "ploidy_prob_2" in chr_df.columns
    assert chr_df["plot_mean_depth"].iloc[0] == 2.0
    assert chr_df["plot_median_depth"].iloc[0] == 2.0
    assert int(chr_df["plq"].iloc[0]) == 99
    assert float(chr_df["ploidy_prob_2"].iloc[0]) == 1.0


def test_build_chromosome_stats_uses_majority_vote_and_average_cnq() -> None:
    depth_df = DepthData(
        pd.DataFrame(
            {
                "Chr": ["chr21", "chr21", "chr21"],
                "Start": [0, 1000, 2000],
                "End": [1000, 2000, 3000],
                "S1": np.array([3.0, 3.0, 2.0], dtype=np.float32),
            },
            index=["chr21:0-1000", "chr21:1000-2000", "chr21:2000-3000"],
        ),
        device="cpu",
    )
    cn_posterior = np.zeros((depth_df.n_bins, depth_df.n_samples, 6), dtype=np.float32)
    cn_posterior[0, 0, 2] = 0.05
    cn_posterior[0, 0, 3] = 0.95
    cn_posterior[1, 0, 2] = 0.10
    cn_posterior[1, 0, 3] = 0.90
    cn_posterior[2, 0, 2] = 0.95
    cn_posterior[2, 0, 3] = 0.05
    cnq = np.array([[20], [30], [40]], dtype=np.int16)

    chr_df = build_chromosome_stats(
        depth_df,
        map_estimates={
            "bin_bias": np.ones(depth_df.n_bins, dtype=np.float32),
            "bin_var": np.full(depth_df.n_bins, 0.02, dtype=np.float32),
            "sample_var": np.full(depth_df.n_samples, 0.03, dtype=np.float32),
        },
        cn_posterior={"cn_posterior": cn_posterior, "cnq": cnq},
        aneuploid_map={0: [("chr21", 3, 2.0 / 3.0)]},
    )

    row = chr_df.iloc[0]
    assert int(row["copy_number"]) == 3
    assert float(row["mean_cn_probability"]) == pytest.approx(2.0 / 3.0)
    assert int(row["plq"]) == 30
    assert float(row["ploidy_prob_3"]) == pytest.approx(2.0 / 3.0)


def test_resolve_observation_model_prefers_preprocess_marker(tmp_path) -> None:
    depth_path = tmp_path / "preprocessed_depth.tsv"
    depth_path.write_text("Bin\tChr\tStart\tEnd\tS1\n", encoding="ascii")
    (tmp_path / "observation_type.txt").write_text("raw\n", encoding="ascii")

    depth_space, obs_likelihood = _resolve_observation_model(
        "auto",
        "auto",
        str(depth_path),
    )

    assert depth_space == "raw"
    assert obs_likelihood == "negative_binomial"


def test_resolve_observation_model_defaults_to_normalized_normal_without_marker(
    tmp_path,
) -> None:
    depth_path = tmp_path / "preprocessed_depth.tsv"
    depth_path.write_text("Bin\tChr\tStart\tEnd\tS1\n", encoding="ascii")

    depth_space, obs_likelihood = _resolve_observation_model(
        "auto",
        "auto",
        str(depth_path),
    )

    assert depth_space == "normalized"
    assert obs_likelihood == "normal"
