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
    build_safe_inference_diagnostic_messages,
    build_site_af_estimates,
    build_bin_stats,
    build_chromosome_stats,
    estimate_site_pop_af_naive_bayes,
    is_fallback_minor_site_data,
    parse_args,
    resolve_site_af_estimator_application,
    should_apply_site_af_estimator,
)
from gatk_sv_ploidy.models import DEFAULT_EPSILON_MEAN
from gatk_sv_ploidy.models import DEFAULT_EPSILON_CONCENTRATION
from gatk_sv_ploidy.models import DEFAULT_BACKGROUND_FACTORS


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

    assert args.autosome_prior_mode == "dirichlet"
    assert args.var_bias_bin == pytest.approx(1e-9)
    assert args.var_sample == 0.01
    assert args.var_bin == pytest.approx(1e-9)
    assert args.alpha_ref == 50.0
    assert args.epsilon_mean == pytest.approx(DEFAULT_EPSILON_MEAN)
    assert args.epsilon_concentration == pytest.approx(
        DEFAULT_EPSILON_CONCENTRATION
    )
    assert DEFAULT_BACKGROUND_FACTORS == 2
    assert args.background_factors == DEFAULT_BACKGROUND_FACTORS
    assert args.depth_space == "auto"
    assert args.obs_likelihood == "auto"
    assert args.svi_init_restarts == 100
    assert args.grad_clip_norm == pytest.approx(10.0)
    assert args.sample_depth_max == 10000.0
    assert args.freeze_bin_bias is False
    assert args.freeze_cn_prior is False
    assert args.af_evidence_mode == "relative"
    assert args.cn_inference_method == "multi-draw"
    assert args.cn_inference_draws == 100
    assert args.site_af_estimator == "auto"
    assert args.learn_af_temperature is True
    assert args.learn_site_af is False
    assert args.site_af_prior_strength == 20.0
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


def test_infer_parse_args_accepts_freeze_cn_prior(monkeypatch) -> None:
    monkeypatch.setattr(
        "sys.argv",
        [
            "infer",
            "--input",
            "depth.tsv",
            "--output-dir",
            "outdir",
            "--freeze-cn-prior",
        ],
    )

    args = parse_args()

    assert args.freeze_cn_prior is True


def test_infer_parse_args_accepts_learn_site_af(monkeypatch) -> None:
    monkeypatch.setattr(
        "sys.argv",
        [
            "infer",
            "--input",
            "depth.tsv",
            "--output-dir",
            "outdir",
            "--learn-site-af",
            "--site-af-prior-strength",
            "12.5",
        ],
    )

    args = parse_args()

    assert args.learn_site_af is True
    assert args.site_af_prior_strength == pytest.approx(12.5)


def test_infer_parse_args_accepts_fixed_af_temperature(monkeypatch) -> None:
    monkeypatch.setattr(
        "sys.argv",
        [
            "infer",
            "--input",
            "depth.tsv",
            "--output-dir",
            "outdir",
            "--fixed-af-temperature",
        ],
    )

    args = parse_args()

    assert args.learn_af_temperature is False


def test_infer_parse_args_accepts_svi_init_restarts(monkeypatch) -> None:
    monkeypatch.setattr(
        "sys.argv",
        [
            "infer",
            "--input",
            "depth.tsv",
            "--output-dir",
            "outdir",
            "--svi-init-restarts",
            "7",
        ],
    )

    args = parse_args()

    assert args.svi_init_restarts == 7


def test_infer_parse_args_accepts_gradient_clip_norm(monkeypatch) -> None:
    monkeypatch.setattr(
        "sys.argv",
        [
            "infer",
            "--input",
            "depth.tsv",
            "--output-dir",
            "outdir",
            "--grad-clip-norm",
            "2.5",
        ],
    )

    args = parse_args()

    assert args.grad_clip_norm == pytest.approx(2.5)


def test_infer_parse_args_defaults_include_elbo_convergence_controls(
    monkeypatch,
) -> None:
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

    assert args.elbo_window == 50
    assert args.elbo_rtol == pytest.approx(1e-4)


def test_infer_parse_args_accepts_elbo_convergence_controls(monkeypatch) -> None:
    monkeypatch.setattr(
        "sys.argv",
        [
            "infer",
            "--input",
            "depth.tsv",
            "--output-dir",
            "outdir",
            "--elbo-window",
            "75",
            "--elbo-rtol",
            "1e-5",
        ],
    )

    args = parse_args()

    assert args.elbo_window == 75
    assert args.elbo_rtol == pytest.approx(1e-5)


def test_infer_parse_args_rejects_removed_sample_baseline_flag(
    monkeypatch,
) -> None:
    monkeypatch.setattr(
        "sys.argv",
        [
            "infer",
            "--input",
            "depth.tsv",
            "--output-dir",
            "outdir",
            "--learn-sample-baseline-autosome-ploidy",
        ],
    )

    with pytest.raises(SystemExit):
        parse_args()


def test_infer_parse_args_rejects_removed_sample_baseline_rebasing_flag(
    monkeypatch,
) -> None:
    monkeypatch.setattr(
        "sys.argv",
        [
            "infer",
            "--input",
            "depth.tsv",
            "--output-dir",
            "outdir",
            "--disable-sample-baseline-autosome-prior-rebasing",
        ],
    )

    with pytest.raises(SystemExit):
        parse_args()


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
    bin_epsilon = np.array(
        [
            [0.01, 0.02],
            [0.03, 0.04],
            [0.05, 0.06],
            [0.07, 0.08],
        ],
        dtype=np.float32,
    )

    bin_df = build_bin_stats(
        data,
        map_estimates={
            "bin_bias": np.ones(data.n_bins, dtype=np.float32),
            "bin_epsilon": bin_epsilon,
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
    epsilon_lookup = {
        (row["chr"], row["sample"]): row["bin_epsilon"]
        for row in bin_df[["chr", "sample", "bin_epsilon"]].to_dict("records")
    }
    sorted_chr_to_idx = {chrom: idx for idx, chrom in enumerate(data.chr)}
    chr21_idx = sorted_chr_to_idx["chr21"]
    assert epsilon_lookup[("chr21", "SAMPLE_B")] == pytest.approx(
        bin_epsilon[chr21_idx, 0]
    )
    assert epsilon_lookup[("chr21", "SAMPLE_A")] == pytest.approx(
        bin_epsilon[chr21_idx, 1]
    )


def test_build_bin_stats_expands_contig_shared_bin_epsilon() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chr21", "chr21"],
            "Start": [0, 1000],
            "End": [1000, 2000],
            "SAMPLE_A": [2.0, 2.1],
            "SAMPLE_B": [1.9, 2.0],
        }
    )
    df["Bin"] = df["Chr"].astype(str) + ":" + df["Start"].astype(str) + "-" + df["End"].astype(str)
    data = DepthData(df.set_index("Bin"), device="cpu")
    cn_posterior = np.zeros((data.n_bins, data.n_samples, 6), dtype=np.float32)
    cn_posterior[:, :, 2] = 0.85
    cn_posterior[:, :, 1] = 0.15
    contig_epsilon = np.array([[0.01, 0.02]], dtype=np.float32)

    bin_df = build_bin_stats(
        data,
        map_estimates={
            "bin_bias": np.ones(data.n_bins, dtype=np.float32),
            "bin_epsilon": contig_epsilon,
            "bin_var": np.full(data.n_bins, 0.04, dtype=np.float32),
            "sample_var": np.full(data.n_samples, 0.09, dtype=np.float32),
        },
        cn_posterior={
            "cn_posterior": cn_posterior,
            "cnq": compute_cnq_from_probabilities(cn_posterior),
        },
    )

    sample_a = bin_df.loc[bin_df["sample"] == "SAMPLE_A", "bin_epsilon"]
    sample_b = bin_df.loc[bin_df["sample"] == "SAMPLE_B", "bin_epsilon"]
    sample_a_expected = contig_epsilon[0, data.sample_ids.index("SAMPLE_A")]
    sample_b_expected = contig_epsilon[0, data.sample_ids.index("SAMPLE_B")]

    assert sample_a.to_list() == pytest.approx([sample_a_expected, sample_a_expected])
    assert sample_b.to_list() == pytest.approx([sample_b_expected, sample_b_expected])


def test_build_bin_stats_uses_background_factors() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chrY", "chrY"],
            "Start": [25, 125],
            "End": [125, 225],
            "SAMPLE_A": [0.0, 0.0],
            "SAMPLE_B": [0.0, 0.0],
        }
    )
    df["Bin"] = df["Chr"].astype(str) + ":" + df["Start"].astype(str) + "-" + df["End"].astype(str)
    data = DepthData(df.set_index("Bin"), device="cpu")
    cn_posterior = np.zeros((data.n_bins, data.n_samples, 6), dtype=np.float32)
    cn_posterior[:, :, 0] = 1.0

    bin_df = build_bin_stats(
        data,
        {
            "bin_bias": np.ones(data.n_bins, dtype=np.float32),
            "bin_var": np.full(data.n_bins, 0.01, dtype=np.float32),
            "sample_var": np.full(data.n_samples, 0.01, dtype=np.float32),
            "background_bin_factors": np.array([[1.0], [2.0]], dtype=np.float32),
            "background_sample_factors": np.array([[0.1, 0.2]], dtype=np.float32),
        },
        {
            "cn_posterior": cn_posterior,
            "cnq": compute_cnq_from_probabilities(cn_posterior),
        },
    )

    observed = {
        (row["chr"], row["sample"], int(row["start"])): row["bin_epsilon"]
        for row in bin_df[["chr", "sample", "start", "bin_epsilon"]].to_dict("records")
    }
    np.testing.assert_allclose(
        observed[("chrY", "SAMPLE_A", 25)],
        2.0 / 3.0 * 0.1,
        rtol=1e-6,
        atol=1e-6,
    )
    np.testing.assert_allclose(
        observed[("chrY", "SAMPLE_A", 125)],
        4.0 / 3.0 * 0.1,
        rtol=1e-6,
        atol=1e-6,
    )
    np.testing.assert_allclose(
        observed[("chrY", "SAMPLE_B", 25)],
        2.0 / 3.0 * 0.2,
        rtol=1e-6,
        atol=1e-6,
    )
    np.testing.assert_allclose(
        observed[("chrY", "SAMPLE_B", 125)],
        4.0 / 3.0 * 0.2,
        rtol=1e-6,
        atol=1e-6,
    )


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


def test_build_safe_inference_diagnostic_messages_excludes_sample_ids(
    tiny_depth_df: pd.DataFrame,
) -> None:
    data = DepthData(tiny_depth_df, device="cpu")
    cn_posterior = np.zeros((data.n_bins, data.n_samples, 6), dtype=np.float32)
    cn_posterior[:, :, 2] = 0.85
    cn_posterior[:, :, 3] = 0.15

    messages = build_safe_inference_diagnostic_messages(
        data,
        map_estimates={
            "sample_var": np.array([0.1, 0.2], dtype=np.float32),
            "bin_bias": np.ones(data.n_bins, dtype=np.float32),
            "bin_var": np.full(data.n_bins, 0.01, dtype=np.float32),
        },
        cn_posterior={
            "cn_posterior": cn_posterior,
            "cnq": compute_cnq_from_probabilities(cn_posterior),
            "cn_map_stability": np.full(
                (data.n_bins, data.n_samples),
                0.9,
                dtype=np.float32,
            ),
        },
    )

    joined = "\n".join(messages)
    assert "S1" not in joined
    assert "S2" not in joined
    assert "Autosomal dominant CN across bin-sample pairs" in joined
    assert "Sample variance latent across samples" in joined


def test_build_safe_inference_diagnostic_messages_reports_allosomal_signals(
    tiny_depth_df: pd.DataFrame,
) -> None:
    data = DepthData(tiny_depth_df, device="cpu")
    cn_posterior = np.zeros((data.n_bins, data.n_samples, 6), dtype=np.float32)
    cn_posterior[:, :, 2] = 1.0

    chr_x_idx = int(np.flatnonzero(data.chr == "chrX")[0])
    chr_y_idx = int(np.flatnonzero(data.chr == "chrY")[0])
    # SAMPLE_B is XX-like in the tiny fixture and SAMPLE_A is XY-like.  Use
    # the sample order from DepthData, but do not emit identifiers in logs.
    cn_posterior[chr_x_idx, 0, :] = 0.0
    cn_posterior[chr_x_idx, 0, 2] = 1.0
    cn_posterior[chr_x_idx, 1, :] = 0.0
    cn_posterior[chr_x_idx, 1, 1] = 1.0
    cn_posterior[chr_y_idx, 0, :] = 0.0
    cn_posterior[chr_y_idx, 0, 1] = 1.0
    cn_posterior[chr_y_idx, 1, :] = 0.0
    cn_posterior[chr_y_idx, 1, 1] = 1.0

    cn_probs = np.full((data.n_bins, 6), 1e-4, dtype=np.float32)
    cn_probs[:, 2] = 0.9995
    cn_probs[chr_x_idx, :] = np.array(
        [1e-4, 0.49, 0.51, 1e-4, 1e-4, 1e-4],
        dtype=np.float32,
    )
    cn_probs[chr_y_idx, :] = np.array(
        [0.40, 0.60, 1e-4, 1e-4, 1e-4, 1e-4],
        dtype=np.float32,
    )

    messages = build_safe_inference_diagnostic_messages(
        data,
        map_estimates={
            "sample_var": np.array([0.1, 0.2], dtype=np.float32),
            "bin_bias": np.ones(data.n_bins, dtype=np.float32),
            "bin_var": np.array([0.02, 0.03, 0.04, 0.8], dtype=np.float32),
            "bin_epsilon": np.zeros((data.n_bins, data.n_samples), dtype=np.float32),
            "cn_probs": cn_probs,
        },
        cn_posterior={
            "cn_posterior": cn_posterior,
            "sex_posterior": np.array(
                [[0.99, 0.01], [0.02, 0.98]],
                dtype=np.float32,
            ),
        },
    )

    joined = "\n".join(messages)
    assert "SAMPLE_A" not in joined
    assert "SAMPLE_B" not in joined
    assert "Bin variance latent across chrY bins" in joined
    assert "Bin epsilon latent across chrY bin-sample pairs" in joined
    assert "Per-bin CN0 prior mass across chrY bins" in joined
    assert "Per-bin CN1 prior mass across chrY bins" in joined
    assert "chrY dominant CN across XX-assigned bin-sample pairs" in joined
    assert "chrY unexpected dominant CN among XX-assigned bin-sample pairs" in joined
    assert "chrY bins with unexpected dominant CN sample fraction among XX-assigned samples thresholds" in joined
    assert "chrY observed-minus-expected plot depth for expected XX copy number" in joined


def test_build_safe_inference_diagnostic_messages_reports_background_factor_chrY_residuals() -> None:
    depth_df = pd.DataFrame(
        {
            "Chr": ["chr21", "chrY", "chrY", "chrY"],
            "Start": [0, 0, 100, 200],
            "End": [100, 100, 200, 300],
            "SAMPLE_XX1": [2.0, 0.02, 0.05, 0.08],
            "SAMPLE_XX2": [2.1, 0.03, 0.05, 0.09],
            "SAMPLE_XY1": [2.0, 0.85, 1.00, 1.15],
            "SAMPLE_XY2": [1.9, 0.88, 1.02, 1.12],
        }
    )
    depth_df["Bin"] = (
        depth_df["Chr"].astype(str) + ":" +
        depth_df["Start"].astype(str) + "-" +
        depth_df["End"].astype(str)
    )
    data = DepthData(depth_df.set_index("Bin"), device="cpu")

    cn_posterior = np.zeros((data.n_bins, data.n_samples, 6), dtype=np.float32)
    cn_posterior[:, :, 2] = 1.0
    chr_y_mask = data.chr == "chrY"
    cn_posterior[chr_y_mask, 0, :] = 0.0
    cn_posterior[chr_y_mask, 0, 0] = 1.0
    cn_posterior[chr_y_mask, 1, :] = 0.0
    cn_posterior[chr_y_mask, 1, 0] = 1.0
    cn_posterior[chr_y_mask, 2, :] = 0.0
    cn_posterior[chr_y_mask, 2, 1] = 1.0
    cn_posterior[chr_y_mask, 3, :] = 0.0
    cn_posterior[chr_y_mask, 3, 1] = 1.0

    messages = build_safe_inference_diagnostic_messages(
        data,
        map_estimates={
            "bin_bias": np.array([1.0, 0.989, 0.992, 0.995], dtype=np.float32),
            "bin_var": np.array([0.01, 0.18, 0.21, 0.24], dtype=np.float32),
            "bin_epsilon": np.zeros((data.n_bins, data.n_samples), dtype=np.float32),
            "background_bin_factors": np.array([[1.0], [0.2], [1.0], [2.4]], dtype=np.float32),
            "background_sample_factors": np.array([[0.02, 0.03, 0.05, 0.06]], dtype=np.float32),
        },
        cn_posterior={
            "cn_posterior": cn_posterior,
            "sex_posterior": np.array(
                [
                    [0.99, 0.01],
                    [0.98, 0.02],
                    [0.02, 0.98],
                    [0.01, 0.99],
                ],
                dtype=np.float32,
            ),
        },
    )

    joined = "\n".join(messages)
    assert "Background sample-factor amplitudes across factor-sample pairs" in joined
    assert "Background factor 1 normalized loadings across chrY bins" in joined
    assert "Correlation(background_factor_1_loading, chrY mean plot depth among XX-assigned samples by bin)" in joined
    assert "Correlation(bin_bias, chrY mean absolute observed-minus-expected plot depth among XX-assigned samples by bin)" in joined
    assert "Correlation(bin_var, chrY mean absolute observed-minus-expected plot depth among XX-assigned samples by bin)" in joined


def test_build_safe_inference_diagnostic_messages_reports_cn2_cn4_ridge_signals(
    tiny_depth_df: pd.DataFrame,
) -> None:
    data = DepthData(tiny_depth_df, device="cpu")
    cn_posterior = np.zeros((data.n_bins, data.n_samples, 6), dtype=np.float32)
    cn_posterior[:, :, 1] = 0.05
    cn_posterior[:, :, 2] = 0.55
    cn_posterior[:, :, 4] = 0.40
    cn_posterior[0, 0, 2] = 0.30
    cn_posterior[0, 0, 4] = 0.65

    messages = build_safe_inference_diagnostic_messages(
        data,
        map_estimates={
            "sample_var": np.array([0.1, 0.2], dtype=np.float32),
            "bin_bias": np.array([0.8, 1.1, 1.0, 1.0], dtype=np.float32),
            "bin_var": np.full(data.n_bins, 0.01, dtype=np.float32),
        },
        cn_posterior={
            "cn_posterior": cn_posterior,
            "cnq": compute_cnq_from_probabilities(cn_posterior),
        },
    )

    joined = "\n".join(messages)
    assert "Autosomal CN2+CN4 posterior mass across bin-sample pairs" in joined
    assert "Autosomal bin-sample pairs with CN2+CN4 posterior mass thresholds" in joined
    assert "Autosomal top-two posterior pair is CN2/CN4" in joined
    assert "Autosomal CN2 posterior on CN4-MAP bin-sample pairs" in joined
    assert "Autosomal CN4 posterior on CN2-MAP bin-sample pairs" in joined


def test_build_safe_inference_diagnostic_messages_reports_recurrent_bin_signals() -> None:
    data = DepthData(
        pd.DataFrame(
            {
                "Chr": ["chr21", "chr18", "chr13", "chrX", "chrY"],
                "Start": [200, 100, 50, 25, 10],
                "End": [300, 200, 150, 125, 110],
                "SAMPLE_B": [3.0, 2.6, 1.8, 2.0, 0.0],
                "SAMPLE_A": [2.8, 2.4, 1.7, 1.0, 1.0],
                "Bin": [
                    "chr21:200-300",
                    "chr18:100-200",
                    "chr13:50-150",
                    "chrX:25-125",
                    "chrY:10-110",
                ],
            }
        ).set_index("Bin"),
        device="cpu",
    )
    cn_posterior = np.zeros((data.n_bins, data.n_samples, 6), dtype=np.float32)
    cn_posterior[:, :, 2] = 1.0
    cn_posterior[0, :, :] = 0.0
    cn_posterior[0, :, 5] = 1.0
    cn_posterior[1, :, :] = 0.0
    cn_posterior[1, :, 4] = 1.0

    messages = build_safe_inference_diagnostic_messages(
        data,
        map_estimates={
            "sample_var": np.array([0.0, 0.0], dtype=np.float32),
            "bin_bias": np.ones(data.n_bins, dtype=np.float32),
            "bin_var": np.array([0.5, 0.3, 0.1, 0.05, 0.02], dtype=np.float32),
            "bin_epsilon": np.array([0.4, 0.2, 0.0, 0.0, 0.0], dtype=np.float32),
            "cn_probs": np.array(
                [
                    [[0.0, 0.0, 0.0, 0.0, 0.1, 0.9]],
                    [[0.0, 0.0, 0.1, 0.0, 0.8, 0.1]],
                    [[0.0, 0.0, 0.9, 0.0, 0.1, 0.0]],
                    [[0.0, 0.0, 1.0, 0.0, 0.0, 0.0]],
                    [[0.0, 0.0, 1.0, 0.0, 0.0, 0.0]],
                ],
                dtype=np.float32,
            ),
        },
        cn_posterior={
            "cn_posterior": cn_posterior,
            "cnq": compute_cnq_from_probabilities(cn_posterior),
        },
    )

    joined = "\n".join(messages)
    assert "Autosomal dominant CN>=4 sample fraction across bins" in joined
    assert "Autosomal bins with dominant CN>=4 sample fraction thresholds" in joined
    assert "Autosomal mean plot depth across bins" in joined
    assert (
        "Correlation(empirical mean autosomal plot depth by bin, autosomal dominant CN>=4 sample fraction by bin)"
        in joined
    )
    assert (
        "Autosomal mean plot depth for bins with dominant CN>=4 sample fraction >=50%"
        in joined
    )
    assert "Correlation(bin_var, autosomal dominant CN>=4 sample fraction by bin)" in joined
    assert "Correlation(bin_epsilon, autosomal dominant CN>=4 sample fraction by bin)" in joined
    assert (
        "Correlation(per-bin CN>=4 prior mass, autosomal dominant CN>=4 sample fraction by bin)"
        in joined
    )


def test_build_safe_inference_diagnostic_messages_reports_af_induced_shift_signals(
    tiny_depth_df: pd.DataFrame,
) -> None:
    data = DepthData(tiny_depth_df, device="cpu")
    depth_only_cn_posterior = np.zeros((data.n_bins, data.n_samples, 6), dtype=np.float32)
    depth_only_cn_posterior[:, :, 2] = 0.90
    depth_only_cn_posterior[:, :, 4] = 0.10

    cn_posterior = depth_only_cn_posterior.copy()
    cn_posterior[0, 0, 2] = 0.20
    cn_posterior[0, 0, 4] = 0.80
    cn_posterior[1, 1, 2] = 0.25
    cn_posterior[1, 1, 3] = 0.75

    messages = build_safe_inference_diagnostic_messages(
        data,
        map_estimates={
            "sample_var": np.array([0.1, 0.2], dtype=np.float32),
            "bin_bias": np.array([0.9, 1.1, 1.0, 1.0], dtype=np.float32),
            "bin_var": np.array([0.02, 0.03, 0.01, 0.01], dtype=np.float32),
        },
        cn_posterior={
            "cn_posterior": cn_posterior,
            "cnq": compute_cnq_from_probabilities(cn_posterior),
        },
        depth_only_cn_posterior={
            "cn_posterior": depth_only_cn_posterior,
            "cnq": compute_cnq_from_probabilities(depth_only_cn_posterior),
        },
    )

    joined = "\n".join(messages)
    assert "Autosomal depth-only to AF-enabled CN MAP transition counts across bin-sample pairs" in joined
    assert "Autosomal AF-induced CN MAP shift burden across samples" in joined
    assert "Concentration of autosomal AF-induced CN MAP shift burden across samples" in joined
    assert "Autosomal AF-induced CN MAP shift sample fraction across bins" in joined
    assert "Autosomal AF-enabled minus depth-only max posterior confidence across bin-sample pairs" in joined
    assert "Autosomal AF-induced CN MAP shift rate by chromosome" in joined


def test_build_safe_inference_diagnostic_messages_reports_af_margin_signals(
    tiny_depth_df: pd.DataFrame,
) -> None:
    data = DepthData(tiny_depth_df, device="cpu")
    cn_posterior = np.zeros((data.n_bins, data.n_samples, 6), dtype=np.float32)
    cn_posterior[:, :, 2] = 0.90
    cn_posterior[:, :, 4] = 0.10

    af_table = np.zeros((6, data.n_bins, data.n_samples), dtype=np.float32)
    af_table[4, 0, 0] = 2.5
    af_table[2, 0, 0] = 0.5
    af_table[5, 1, 1] = 1.5
    af_table[2, 1, 1] = 0.25

    messages = build_safe_inference_diagnostic_messages(
        data,
        map_estimates={
            "sample_var": np.array([0.1, 0.2], dtype=np.float32),
            "bin_bias": np.array([0.9, 1.1, 1.0, 1.0], dtype=np.float32),
            "bin_var": np.array([0.02, 0.03, 0.01, 0.01], dtype=np.float32),
            "af_temperature": np.asarray(0.0125, dtype=np.float32),
        },
        cn_posterior={
            "cn_posterior": cn_posterior,
            "cnq": compute_cnq_from_probabilities(cn_posterior),
        },
        af_table=af_table,
    )

    joined = "\n".join(messages)
    assert "Learned AF temperature scalar=0.012500" in joined
    assert "Autosomal AF CN4-CN2 log-likelihood margin across bin-sample pairs" in joined
    assert "Autosomal bin-sample pairs with AF CN4-CN2 margin thresholds" in joined
    assert "Autosomal AF CN5-CN2 log-likelihood margin across bin-sample pairs" in joined
    assert "Autosomal bin-sample pairs with AF CN5-CN2 margin thresholds" in joined
    assert "Autosomal AF preference for CN4 over CN2 across samples" in joined
    assert "Autosomal AF preference for CN5 over CN2 across bins" in joined


def test_build_safe_inference_diagnostic_messages_reports_raw_af_grid_anomaly_signals(
    tiny_depth_df: pd.DataFrame,
    tiny_site_data: dict[str, np.ndarray],
) -> None:
    data = DepthData(tiny_depth_df, device="cpu", site_data=tiny_site_data)
    depth_only_cn_posterior = np.zeros((data.n_bins, data.n_samples, 6), dtype=np.float32)
    depth_only_cn_posterior[:, :, 2] = 0.90
    depth_only_cn_posterior[:, :, 4] = 0.10

    cn_posterior = depth_only_cn_posterior.copy()
    cn_posterior[0, 0, 2] = 0.15
    cn_posterior[0, 0, 4] = 0.85
    cn_posterior[1, 1, 2] = 0.20
    cn_posterior[1, 1, 5] = 0.80

    messages = build_safe_inference_diagnostic_messages(
        data,
        map_estimates={
            "sample_var": np.array([0.1, 0.2], dtype=np.float32),
            "bin_bias": np.array([0.9, 1.1, 1.0, 1.0], dtype=np.float32),
            "bin_var": np.array([0.02, 0.03, 0.01, 0.01], dtype=np.float32),
        },
        cn_posterior={
            "cn_posterior": cn_posterior,
            "cnq": compute_cnq_from_probabilities(cn_posterior),
        },
        depth_only_cn_posterior={
            "cn_posterior": depth_only_cn_posterior,
            "cnq": compute_cnq_from_probabilities(depth_only_cn_posterior),
        },
    )

    joined = "\n".join(messages)
    assert "Autosomal raw AF CN4-vs-diploid grid advantage across observed sites" in joined
    assert "Autosomal observed sites with raw AF CN4-vs-diploid grid advantage thresholds" in joined
    assert "Autosomal raw AF CN5-vs-diploid grid advantage across observed sites" in joined
    assert "Autosomal raw AF preference for CN4 genotype grid over diploid grid across samples" in joined
    assert "Autosomal raw AF preference for CN5 genotype grid over diploid grid across bins" in joined


def test_build_safe_inference_diagnostic_messages_reports_site_pop_af_prior_signals(
    tiny_depth_df: pd.DataFrame,
    tiny_site_data: dict[str, np.ndarray],
) -> None:
    data = DepthData(tiny_depth_df, device="cpu", site_data=tiny_site_data)
    cn_posterior = np.zeros((data.n_bins, data.n_samples, 6), dtype=np.float32)
    cn_posterior[:, :, 2] = 0.90
    cn_posterior[:, :, 4] = 0.10

    af_table = np.zeros((6, data.n_bins, data.n_samples), dtype=np.float32)
    af_table[4, 0, :] = 1.5
    af_table[2, 0, :] = 0.25
    af_table[5, 1, :] = 0.75
    af_table[2, 1, :] = 0.10

    messages = build_safe_inference_diagnostic_messages(
        data,
        map_estimates={
            "sample_var": np.array([0.1, 0.2], dtype=np.float32),
            "bin_bias": np.array([0.9, 1.1, 1.0, 1.0], dtype=np.float32),
            "bin_var": np.array([0.02, 0.03, 0.01, 0.01], dtype=np.float32),
        },
        cn_posterior={
            "cn_posterior": cn_posterior,
            "cnq": compute_cnq_from_probabilities(cn_posterior),
        },
        af_table=af_table,
    )

    joined = "\n".join(messages)
    assert "Autosomal effective site_pop_af uses leave-one-out self-pooled resolution=" in joined
    assert "Autosomal effective site_pop_af across observed site-sample values" in joined
    assert "Autosomal diploid heterozygous prior mass 2p(1-p) across observed site-sample values" in joined
    assert "Autosomal observed site-sample pairs with diploid heterozygous prior mass thresholds" in joined
    assert "Mean autosomal diploid heterozygous prior mass across bins" in joined


def test_build_safe_inference_diagnostic_messages_reports_site_level_af_margin_signals(
    tiny_depth_df: pd.DataFrame,
    tiny_site_data: dict[str, np.ndarray],
) -> None:
    data = DepthData(tiny_depth_df, device="cpu", site_data=tiny_site_data)
    cn_posterior = np.zeros((data.n_bins, data.n_samples, 6), dtype=np.float32)
    cn_posterior[:, :, 2] = 0.90
    cn_posterior[:, :, 4] = 0.10

    messages = build_safe_inference_diagnostic_messages(
        data,
        map_estimates={
            "sample_var": np.array([0.1, 0.2], dtype=np.float32),
            "bin_bias": np.array([0.9, 1.1, 1.0, 1.0], dtype=np.float32),
            "bin_var": np.array([0.02, 0.03, 0.01, 0.01], dtype=np.float32),
        },
        cn_posterior={
            "cn_posterior": cn_posterior,
            "cnq": compute_cnq_from_probabilities(cn_posterior),
        },
        af_concentration=50.0,
    )

    joined = "\n".join(messages)
    assert "Autosomal site-level AF CN4-CN2 margin across observed site-sample pairs" in joined
    assert "Autosomal observed site-sample pairs with site-level AF CN4-CN2 margin thresholds" in joined
    assert "Autosomal site-level AF CN5-CN2 margin across observed site-sample pairs" in joined
    assert "Autosomal observed site_total across observed site-sample pairs" in joined
    assert "Site-total-weighted autosomal site-level AF CN4-CN2 margin=" in joined


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
