from __future__ import annotations

import numpy as np
import pandas as pd
import torch

from gatk_sv_ploidy.data import DepthData
from gatk_sv_ploidy.models import (
    CNVModel,
    _center_af_table_numpy,
    _depth_log_lik_numpy,
    _matched_residual_scale,
    _marginalized_af_log_lik,
    _marginalized_af_log_lik_numpy,
    _negative_binomial_log_lik_numpy,
    _precompute_af_table,
)


def test_marginalized_af_log_lik_matches_numpy(tiny_site_data: dict[str, np.ndarray]) -> None:
    cn = torch.full((4, 2), 2, dtype=torch.long)
    torch_result = _marginalized_af_log_lik(
        torch.tensor(tiny_site_data["site_alt"], dtype=torch.float32),
        torch.tensor(tiny_site_data["site_total"], dtype=torch.float32),
        torch.tensor(tiny_site_data["site_pop_af"], dtype=torch.float32),
        torch.tensor(tiny_site_data["site_mask"], dtype=torch.bool),
        cn=cn,
        n_states=6,
        concentration=50.0,
    ).cpu().numpy()

    numpy_result = _marginalized_af_log_lik_numpy(
        tiny_site_data["site_alt"],
        tiny_site_data["site_total"],
        tiny_site_data["site_pop_af"],
        tiny_site_data["site_mask"],
        cn_state=2,
        n_states=6,
        concentration=50.0,
    )

    np.testing.assert_allclose(torch_result, numpy_result, rtol=1e-5, atol=1e-5)


def test_precompute_af_table_contains_expected_state_slice(
    tiny_site_data: dict[str, np.ndarray],
) -> None:
    table = _precompute_af_table(
        torch.tensor(tiny_site_data["site_alt"], dtype=torch.float32),
        torch.tensor(tiny_site_data["site_total"], dtype=torch.float32),
        torch.tensor(tiny_site_data["site_pop_af"], dtype=torch.float32),
        torch.tensor(tiny_site_data["site_mask"], dtype=torch.bool),
        n_states=6,
        concentration=50.0,
    ).cpu().numpy()

    assert table.shape == (6, 4, 2)
    np.testing.assert_allclose(
        table[2],
        _marginalized_af_log_lik_numpy(
            tiny_site_data["site_alt"],
            tiny_site_data["site_total"],
            tiny_site_data["site_pop_af"],
            tiny_site_data["site_mask"],
            cn_state=2,
            n_states=6,
            concentration=50.0,
        ),
        rtol=1e-5,
        atol=1e-5,
    )


def test_relative_af_evidence_centers_against_fixed_reference_mixture(
    tiny_site_data: dict[str, np.ndarray],
) -> None:
    raw_table = _precompute_af_table(
        torch.tensor(tiny_site_data["site_alt"], dtype=torch.float32),
        torch.tensor(tiny_site_data["site_total"], dtype=torch.float32),
        torch.tensor(tiny_site_data["site_pop_af"], dtype=torch.float32),
        torch.tensor(tiny_site_data["site_mask"], dtype=torch.bool),
        n_states=6,
        concentration=50.0,
    ).cpu().numpy()
    model = CNVModel(
        autosome_prior_mode="shrinkage",
        af_evidence_mode="relative",
        sex_cn_weight=0.0,
        device="cpu",
    )
    chr_type_np = np.zeros(raw_table.shape[1], dtype=np.int64)

    transformed = model._apply_af_evidence_mode_numpy(raw_table, chr_type_np)
    reference = model._af_reference_cn_probs_numpy(chr_type_np, raw_table.shape[1])
    expected = _center_af_table_numpy(raw_table, reference)

    np.testing.assert_allclose(transformed, expected, rtol=1e-6, atol=1e-6)
    np.testing.assert_allclose(
        np.log(np.sum(reference.T[:, :, np.newaxis] * np.exp(transformed), axis=0)),
        np.zeros((raw_table.shape[1], raw_table.shape[2]), dtype=np.float64),
        rtol=1e-6,
        atol=1e-6,
    )


def test_marginalized_af_log_lik_scales_with_informative_site_count() -> None:
    single_alt = torch.tensor([[[5]]], dtype=torch.float32)
    single_total = torch.tensor([[[10]]], dtype=torch.float32)
    duplicated_alt = torch.tensor([[[5], [5]]], dtype=torch.float32)
    duplicated_total = torch.tensor([[[10], [10]]], dtype=torch.float32)
    single_pop_af = torch.tensor([[0.5]], dtype=torch.float32)
    duplicated_pop_af = torch.tensor([[0.5, 0.5]], dtype=torch.float32)
    single_mask = torch.tensor([[[True]]], dtype=torch.bool)
    duplicated_mask = torch.tensor([[[True], [True]]], dtype=torch.bool)
    cn = torch.tensor([[2]], dtype=torch.long)

    single = _marginalized_af_log_lik(
        single_alt,
        single_total,
        single_pop_af,
        single_mask,
        cn=cn,
        n_states=6,
        concentration=50.0,
    )
    duplicated = _marginalized_af_log_lik(
        duplicated_alt,
        duplicated_total,
        duplicated_pop_af,
        duplicated_mask,
        cn=cn,
        n_states=6,
        concentration=50.0,
    )

    torch.testing.assert_close(duplicated, 2.0 * single, rtol=1e-5, atol=1e-5)


def test_matched_residual_scale_matches_target_variance() -> None:
    target_variance = np.array([0.25, 1.0, 4.0], dtype=np.float64)

    normal_scale = _matched_residual_scale(target_variance, "normal", 3.5)
    laplace_scale = _matched_residual_scale(target_variance, "laplace", 3.5)
    studentt_scale = _matched_residual_scale(target_variance, "studentt", 3.5)

    np.testing.assert_allclose(normal_scale ** 2, target_variance)
    np.testing.assert_allclose(2.0 * laplace_scale ** 2, target_variance)
    np.testing.assert_allclose(
        (3.5 / (3.5 - 2.0)) * studentt_scale ** 2,
        target_variance,
    )


def test_depth_log_lik_numpy_studentt_is_heavier_tailed_than_normal() -> None:
    obs = np.array([[[8.0]]], dtype=np.float64)
    expected = np.array([[[2.0]]], dtype=np.float64)
    variance = np.array([[[1.0]]], dtype=np.float64)

    normal_ll = _depth_log_lik_numpy(obs, expected, variance, "normal", 3.5)
    studentt_ll = _depth_log_lik_numpy(obs, expected, variance, "studentt", 3.5)

    assert float(studentt_ll[0, 0, 0]) > float(normal_ll[0, 0, 0])


def test_negative_binomial_log_lik_numpy_prefers_matching_mean() -> None:
    obs = np.array([[[20.0]]], dtype=np.float64)
    overdispersion = np.array([[[0.05]]], dtype=np.float64)

    matching_ll = _negative_binomial_log_lik_numpy(
        obs,
        np.array([[[20.0]]], dtype=np.float64),
        overdispersion,
    )
    mismatched_ll = _negative_binomial_log_lik_numpy(
        obs,
        np.array([[[5.0]]], dtype=np.float64),
        overdispersion,
    )

    assert float(matching_ll[0, 0, 0]) > float(mismatched_ll[0, 0, 0])


def test_run_discrete_inference_uses_learned_af_temperature_map() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chr21"],
            "Start": [0],
            "End": [1000],
            "S1": [2.0],
        }
    )
    df["Bin"] = "chr21:0-1000"
    data = DepthData(df.set_index("Bin"), device="cpu")
    model = CNVModel(
        autosome_prior_mode="dirichlet",
        learn_af_temperature=True,
        af_weight=0.25,
        var_sample=1000.0,
        var_bin=1000.0,
        sex_cn_weight=0.0,
        epsilon_mean=0.0,
    )
    base_map = {
        "bin_bias": np.array([1.0], dtype=np.float32),
        "sample_var": np.array([1000.0], dtype=np.float32),
        "bin_var": np.array([1000.0], dtype=np.float32),
        "cn_probs": np.full((1, 6), 1.0 / 6.0, dtype=np.float32),
    }
    af_table = np.full((6, 1, 1), -5.0, dtype=np.float32)
    af_table[2, 0, 0] = 0.0

    low_temp = model.run_discrete_inference(
        data,
        map_estimates={
            **base_map,
            "af_temperature": np.asarray(0.0, dtype=np.float32),
        },
        af_table=af_table,
    )
    high_temp = model.run_discrete_inference(
        data,
        map_estimates={
            **base_map,
            "af_temperature": np.asarray(2.0, dtype=np.float32),
        },
        af_table=af_table,
    )

    assert high_temp["cn_posterior"][0, 0, 2] > low_temp["cn_posterior"][0, 0, 2]


def test_cnv_model_accepts_lowrank_guide() -> None:
    model = CNVModel(guide_type="lowrank")
    assert model.guide.__class__.__name__ == "AutoLowRankMultivariateNormal"


def test_cnv_model_negative_binomial_adds_sample_depth_latent() -> None:
    model = CNVModel(obs_likelihood="negative_binomial")
    assert "sample_depth" in model._latent_sites
    assert "bin_epsilon" in model._latent_sites


def test_estimate_sample_depth_init_uses_autosomal_counts_per_kb() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chr21", "chrX"],
            "Start": [0, 0],
            "End": [1000, 1000],
            "S1": [10, 5],
            "S2": [20, 10],
        },
        index=["chr21:0-1000", "chrX:0-1000"],
    )
    data = DepthData(
        df,
        device="cpu",
        depth_space="raw",
        clamp_threshold=None,
    )
    model = CNVModel(obs_likelihood="negative_binomial", device="cpu")

    init = model._estimate_sample_depth_init(data).cpu().numpy()

    np.testing.assert_allclose(init, np.array([10.0, 20.0], dtype=np.float32))


def test_run_discrete_inference_uses_bin_epsilon() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chrY"],
            "Start": [25],
            "End": [125],
            "S1": [0.2],
        }
    )
    df["Bin"] = "chrY:25-125"
    data = DepthData(df.set_index("Bin"), device="cpu")
    model = CNVModel(sex_cn_weight=0.0, epsilon_mean=0.1)
    base_map = {
        "bin_bias": np.array([1.0], dtype=np.float32),
        "sample_var": np.array([0.01], dtype=np.float32),
        "bin_var": np.array([0.01], dtype=np.float32),
        "cn_probs": np.full((1, 6), 1.0 / 6.0, dtype=np.float32),
    }

    no_epsilon = model.run_discrete_inference(data, map_estimates=base_map)
    with_epsilon = model.run_discrete_inference(
        data,
        map_estimates={
            **base_map,
            "bin_epsilon": np.array([0.2], dtype=np.float32),
        },
    )

    assert with_epsilon["cn_posterior"][0, 0, 0] > no_epsilon["cn_posterior"][0, 0, 0]


def test_run_discrete_inference_negative_binomial_uses_sample_depth() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chr21"],
            "Start": [0],
            "End": [1000],
            "S1": [10],
        }
    )
    df["Bin"] = "chr21:0-1000"
    data = DepthData(
        df.set_index("Bin"),
        device="cpu",
        depth_space="raw",
        clamp_threshold=None,
    )
    model = CNVModel(
        sex_cn_weight=0.0,
        epsilon_mean=0.0,
        obs_likelihood="negative_binomial",
    )
    posterior = model.run_discrete_inference(
        data,
        map_estimates={
            "bin_bias": np.array([1.0], dtype=np.float32),
            "sample_var": np.array([1e-4], dtype=np.float32),
            "bin_var": np.array([1e-4], dtype=np.float32),
            "sample_depth": np.array([10.0], dtype=np.float32),
            "cn_probs": np.full((1, 6), 1.0 / 6.0, dtype=np.float32),
        },
    )

    assert int(np.argmax(posterior["cn_posterior"][0, 0, :])) == 2
    assert float(posterior["cn_posterior"][0, 0, 2]) > float(
        posterior["cn_posterior"][0, 0, 1]
    )


def test_run_discrete_inference_supports_autosome_shrinkage_maps() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chr21"],
            "Start": [0],
            "End": [1000],
            "S1": [2.0],
        }
    )
    df["Bin"] = "chr21:0-1000"
    data = DepthData(df.set_index("Bin"), device="cpu")
    model = CNVModel(
        sex_cn_weight=0.0,
        epsilon_mean=0.0,
        autosome_prior_mode="shrinkage",
    )

    common_maps = {
        "bin_bias": np.array([1.0], dtype=np.float32),
        "sample_var": np.array([0.01], dtype=np.float32),
        "bin_var": np.array([0.01], dtype=np.float32),
    }
    shrinkage_maps = {
        **common_maps,
        "autosome_nonref_mean": np.asarray(0.05, dtype=np.float32),
        "autosome_nonref_prob": np.array([0.2], dtype=np.float32),
        "autosome_alt_cn_probs": np.array([[0.0, 0.0, 1.0, 0.0, 0.0]], dtype=np.float32),
    }
    explicit_maps = {
        **common_maps,
        "cn_probs": np.array([[0.0, 0.0, 0.8, 0.2, 0.0, 0.0]], dtype=np.float32),
    }

    shrinkage_post = model.run_discrete_inference(data, map_estimates=shrinkage_maps)
    explicit_post = model.run_discrete_inference(data, map_estimates=explicit_maps)

    np.testing.assert_allclose(
        shrinkage_post["cn_posterior"],
        explicit_post["cn_posterior"],
        rtol=1e-6,
        atol=1e-6,
    )


def test_shrinkage_prior_trains_with_autosome_and_sex_bins() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chr21", "chrX"],
            "Start": [0, 0],
            "End": [1000, 1000],
            "S1": [2.0, 1.0],
        },
        index=["chr21:0-1000", "chrX:0-1000"],
    )
    data = DepthData(df, device="cpu")
    model = CNVModel(
        sex_cn_weight=0.0,
        epsilon_mean=0.0,
        guide_type="diagonal",
        autosome_prior_mode="shrinkage",
    )

    model.train(data, max_iter=1, log_freq=1, early_stopping=False)
    estimates = model.get_map_estimates(data, estimate_method="median")

    assert estimates["cn_probs"].shape == (2, 6)
    assert "autosome_nonref_prob" in estimates
    assert "sex_cn_probs" in estimates


def test_get_map_estimates_median_uses_guide_median(monkeypatch) -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chr21"],
            "Start": [0],
            "End": [1000],
            "S1": [2.0],
        }
    )
    df["Bin"] = "chr21:0-1000"
    data = DepthData(df.set_index("Bin"), device="cpu")
    model = CNVModel(
        sex_cn_weight=0.0,
        epsilon_mean=0.0,
        guide_type="diagonal",
        autosome_prior_mode="dirichlet",
    )

    class DummyGuide:
        def median(self, *args, **kwargs):
            return {
                "bin_bias": torch.tensor([1.25], dtype=torch.float32),
                "sample_var": torch.tensor([0.05], dtype=torch.float32),
                "bin_var": torch.tensor([0.02], dtype=torch.float32),
                "cn_probs": torch.full((1, 6), 1.0 / 6.0, dtype=torch.float32),
            }

    model.guide = DummyGuide()
    monkeypatch.setattr(
        model,
        "_infer_discrete_assignments",
        lambda data, continuous_estimates, model_kw=None: {
            "cn": np.array([[2]], dtype=np.int64),
        },
    )

    estimates = model.get_map_estimates(data, estimate_method="median")

    assert float(estimates["bin_bias"][0]) == 1.25
    assert int(estimates["cn"][0, 0]) == 2


def test_run_discrete_inference_multi_draw_averages_posteriors(monkeypatch) -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chr21"],
            "Start": [0],
            "End": [1000],
            "S1": [2.0],
        }
    )
    df["Bin"] = "chr21:0-1000"
    data = DepthData(df.set_index("Bin"), device="cpu")
    model = CNVModel(sex_cn_weight=0.0, epsilon_mean=0.0, guide_type="diagonal")

    draw_biases = iter([1.0, 2.0, 3.0])

    def fake_get_continuous_estimates(data, estimate_method="current", model_kw=None):
        bias = next(draw_biases)
        return {
            "bin_bias": torch.tensor([bias], dtype=torch.float32),
            "sample_var": torch.tensor([0.05], dtype=torch.float32),
            "bin_var": torch.tensor([0.02], dtype=torch.float32),
            "cn_probs": torch.full((1, 6), 1.0 / 6.0, dtype=torch.float32),
        }

    def fake_run_fixed(data, maps, af_table=None):
        bias = float(np.asarray(maps["bin_bias"]).reshape(-1)[0])
        if bias < 1.5:
            return {
                "cn_posterior": np.array([[[1.0, 0.0, 0.0, 0.0, 0.0, 0.0]]], dtype=np.float32),
                "sex_posterior": np.array([[1.0, 0.0]], dtype=np.float32),
            }
        if bias < 2.5:
            return {
                "cn_posterior": np.array([[[0.0, 1.0, 0.0, 0.0, 0.0, 0.0]]], dtype=np.float32),
                "sex_posterior": np.array([[0.0, 1.0]], dtype=np.float32),
            }
        return {
            "cn_posterior": np.array([[[0.0, 0.0, 1.0, 0.0, 0.0, 0.0]]], dtype=np.float32),
            "sex_posterior": np.array([[0.5, 0.5]], dtype=np.float32),
        }

    monkeypatch.setattr(model, "_get_continuous_estimates", fake_get_continuous_estimates)
    monkeypatch.setattr(model, "_run_discrete_inference_fixed_latents", fake_run_fixed)

    collected_draws = {}
    posterior = model.run_discrete_inference_multi_draw(
        data,
        n_draws=3,
        draw_estimate_collector=collected_draws,
    )

    np.testing.assert_allclose(
        posterior["cn_posterior"],
        np.array([[[1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 0.0, 0.0, 0.0]]], dtype=np.float32),
    )
    np.testing.assert_allclose(
        posterior["sex_posterior"],
        np.array([[0.5, 0.5]], dtype=np.float32),
    )
    np.testing.assert_allclose(
        posterior["cn_map_stability"],
        np.array([[1.0 / 3.0]], dtype=np.float32),
    )
    np.testing.assert_allclose(
        posterior["sex_map_stability"],
        np.array([2.0 / 3.0], dtype=np.float32),
    )
    assert [
        float(np.asarray(draw).reshape(-1)[0])
        for draw in collected_draws["bin_bias"]
    ] == [1.0, 2.0, 3.0]
    assert len(collected_draws["cn_probs"]) == 3
