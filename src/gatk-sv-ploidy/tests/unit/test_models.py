from __future__ import annotations
import pytest

import numpy as np
import pandas as pd
import pyro
import pyro.poutine as poutine
import torch
from pyro.infer import TraceEnum_ELBO

from gatk_sv_ploidy.data import DepthData
from gatk_sv_ploidy.models import (
    CNVModel,
    DEFAULT_BACKGROUND_FACTORS,
    DEFAULT_EPSILON_CONCENTRATION,
    DEFAULT_EPSILON_MEAN,
    _center_af_table_numpy,
    _depth_log_lik_numpy,
    _matched_residual_scale,
    _marginalized_af_log_lik,
    _marginalized_af_log_lik_numpy,
    _negative_binomial_log_lik_numpy,
    _precompute_af_table,
    _raw_expected_depth_units,
    _resolve_fixed_site_pop_af_numpy,
    _resolve_fixed_site_pop_af_torch,
    _windowed_relative_elbo_change,
)


def test_cnv_model_defaults_match_current_preferred_configuration() -> None:
    model = CNVModel()

    assert model.autosome_prior_mode == "dirichlet"
    assert model.var_bias_bin == 0.02
    assert model.var_sample == 0.00025
    assert model.var_bin == 0.001
    assert model.epsilon_mean == pytest.approx(DEFAULT_EPSILON_MEAN)
    assert model.epsilon_concentration == pytest.approx(DEFAULT_EPSILON_CONCENTRATION)
    assert DEFAULT_BACKGROUND_FACTORS == 2
    assert model.background_factors == DEFAULT_BACKGROUND_FACTORS
    assert model.af_evidence_mode == "relative"
    assert model.learn_af_temperature is True
    assert model.freeze_bin_bias is False
    assert model.freeze_cn_prior is False
    assert model.freeze_sample_var is False
    assert "af_temperature" in model._latent_sites
    assert "background_bin_factors" in model._latent_sites
    assert "background_sample_factors" in model._latent_sites


def test_freeze_sample_var_removes_latent_and_fixes_to_zero() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chr21", "chr21", "chrX"],
            "Start": [0, 1000, 0],
            "End": [1000, 2000, 1000],
            "S1": [10, 12, 5],
            "S2": [20, 24, 10],
        },
        index=["chr21:0-1000", "chr21:1000-2000", "chrX:0-1000"],
    )
    data = DepthData(
        df,
        device="cpu",
        depth_space="raw",
        clamp_threshold=None,
        dtype=torch.float64,
    )

    model = CNVModel(
        obs_likelihood="negative_binomial",
        device="cpu",
        dtype=torch.float64,
        freeze_sample_var=True,
    )
    assert "sample_var" not in model._latent_sites
    model_kw = model._model_kwargs(data)
    assert "sample_var_fixed" in model_kw
    np.testing.assert_allclose(
        model_kw["sample_var_fixed"].cpu().numpy(),
        np.array([0.0, 0.0], dtype=np.float64),
    )

    model.train(data, max_iter=1, log_freq=1, early_stopping=False)
    estimates = model.get_map_estimates(data, estimate_method="median")

    np.testing.assert_allclose(
        estimates["sample_var"],
        np.array([0.0, 0.0], dtype=np.float64),
    )


def test_freeze_cn_prior_removes_latent_and_uses_default_prior() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chr21", "chrX", "chrY"],
            "Start": [0, 0, 0],
            "End": [1000, 1000, 1000],
            "S1": [10.0, 5.0, 1.0],
            "S2": [12.0, 6.0, 1.0],
        },
        index=["chr21:0-1000", "chrX:0-1000", "chrY:0-1000"],
    )
    data = DepthData(
        df,
        device="cpu",
        depth_space="raw",
        clamp_threshold=None,
        dtype=torch.float64,
    )

    model = CNVModel(
        obs_likelihood="negative_binomial",
        device="cpu",
        dtype=torch.float64,
        freeze_cn_prior=True,
    )
    assert "cn_probs" not in model._latent_sites

    model.train(data, max_iter=1, log_freq=1, early_stopping=False)
    estimates = model.get_map_estimates(data, estimate_method="median")

    expected = np.array(
        [
            [1.0 / 55.0, 1.0 / 55.0, 50.0 / 55.0, 1.0 / 55.0, 1.0 / 55.0, 1.0 / 55.0],
            [1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0],
            [1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0],
        ],
        dtype=np.float64,
    )
    np.testing.assert_allclose(estimates["cn_probs"], expected)


def test_freeze_cn_prior_rejects_shrinkage_mode() -> None:
    with pytest.raises(ValueError, match="freeze_cn_prior"):
        CNVModel(autosome_prior_mode="shrinkage", freeze_cn_prior=True)


def test_multi_draw_inference_backfills_frozen_latents(monkeypatch) -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chr21", "chr21", "chrX"],
            "Start": [0, 1000, 0],
            "End": [1000, 2000, 1000],
            "S1": [10, 12, 5],
            "S2": [20, 24, 10],
        },
        index=["chr21:0-1000", "chr21:1000-2000", "chrX:0-1000"],
    )
    data = DepthData(
        df,
        device="cpu",
        depth_space="raw",
        clamp_threshold=None,
        dtype=torch.float64,
    )

    model = CNVModel(
        obs_likelihood="negative_binomial",
        guide_type="diagonal",
        device="cpu",
        dtype=torch.float64,
        freeze_bin_bias=True,
        freeze_sample_var=True,
        freeze_sample_depth=True,
    )
    model.train(data, max_iter=1, log_freq=1, early_stopping=False)

    def fake_run_fixed(data, maps, af_table=None):
        np.testing.assert_allclose(
            np.asarray(maps["bin_bias"]).reshape(-1),
            np.ones(data.n_bins, dtype=np.float64),
        )
        np.testing.assert_allclose(
            np.asarray(maps["sample_var"]).reshape(-1),
            np.zeros(data.n_samples, dtype=np.float64),
        )
        assert np.asarray(maps["sample_depth"]).reshape(-1).shape == (data.n_samples,)

        cn_posterior = np.zeros(
            (data.n_bins, data.n_samples, model.n_states),
            dtype=np.float64,
        )
        cn_posterior[:, :, model.ref_state] = 1.0
        sex_posterior = np.zeros((data.n_samples, 2), dtype=np.float64)
        sex_posterior[:, 0] = 1.0
        return {
            "cn_posterior": cn_posterior,
            "sex_posterior": sex_posterior,
        }

    monkeypatch.setattr(model, "_run_discrete_inference_fixed_latents", fake_run_fixed)

    posterior = model.run_discrete_inference_multi_draw(data, n_draws=1)

    assert posterior["cn_posterior"].shape == (data.n_bins, data.n_samples, model.n_states)
    assert posterior["sex_posterior"].shape == (data.n_samples, 2)


def test_sample_var_uses_exponential_prior() -> None:
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
        var_sample=0.123,
        device="cpu",
    )

    trace = poutine.trace(model.model).get_trace(**model._model_kwargs(data))
    sample_var_fn = trace.nodes["sample_var"]["fn"]

    assert sample_var_fn.__class__.__name__ == "Exponential"
    np.testing.assert_allclose(
        np.asarray(sample_var_fn.rate.detach().cpu().numpy()).squeeze(),
        1.0 / 0.123,
        rtol=1e-6,
        atol=1e-6,
    )


def test_bin_epsilon_uses_gamma_prior_with_requested_mean_and_concentration() -> None:
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
        epsilon_mean=0.123,
        epsilon_concentration=0.5,
        device="cpu",
    )

    trace = poutine.trace(model.model).get_trace(**model._model_kwargs(data))
    bin_epsilon_fn = trace.nodes["bin_epsilon"]["fn"]
    gamma_dist = getattr(bin_epsilon_fn, "base_dist", bin_epsilon_fn)

    assert gamma_dist.__class__.__name__ == "Gamma"
    np.testing.assert_allclose(
        np.asarray(gamma_dist.concentration.detach().cpu().numpy()).squeeze(),
        0.5,
        rtol=1e-6,
        atol=1e-6,
    )
    np.testing.assert_allclose(
        np.asarray(gamma_dist.rate.detach().cpu().numpy()).squeeze(),
        0.5 / 0.123,
        rtol=1e-6,
        atol=1e-6,
    )


def test_select_svi_initialization_restores_best_candidate(monkeypatch) -> None:
    model = CNVModel()
    candidate_losses = [7.0, 2.5, 5.0]
    init_calls = []

    class DummyGuide:
        def __init__(self, candidate: int, loss: float) -> None:
            self.candidate = candidate
            self.loss = loss

    class FakeElbo:
        def loss(self, model_fn, guide, **kwargs):
            pyro.param("selected_candidate", torch.tensor(float(guide.candidate)))
            return guide.loss

    def fake_make_init_loc_fn(data, randomize: bool = False):
        candidate = len(init_calls)
        init_calls.append(randomize)
        return {"candidate": candidate}

    def fake_build_guide(init_loc_fn):
        candidate = init_loc_fn["candidate"]
        return DummyGuide(candidate, candidate_losses[candidate])

    pyro.clear_param_store()
    monkeypatch.setattr(model, "_make_init_loc_fn", fake_make_init_loc_fn)
    monkeypatch.setattr(model, "_build_guide", fake_build_guide)

    model._select_svi_initialization(
        object(),
        model_kw={},
        elbo=FakeElbo(),
        init_restarts=3,
    )

    assert init_calls == [False, True, True]
    assert model.guide.candidate == 1
    assert pyro.param("selected_candidate").item() == pytest.approx(1.0)


def test_windowed_relative_elbo_change_requires_two_complete_windows() -> None:
    assert _windowed_relative_elbo_change([100.0, 99.0, 98.0], window=2) is None


def test_windowed_relative_elbo_change_uses_window_means() -> None:
    relative_change = _windowed_relative_elbo_change(
        [100.0, 100.0, 100.01, 100.01],
        window=2,
    )

    assert relative_change == pytest.approx(1e-4)


def test_train_early_stopping_uses_windowed_relative_elbo_change(monkeypatch) -> None:
    model = CNVModel()
    losses = iter([100.0, 100.0, 100.01, 100.01, 100.01001, 100.01001, 99.0])
    model.loss_history["epoch"] = [999]
    model.loss_history["elbo"] = [123.0]

    class FakeScheduler:
        def step(self) -> None:
            return None

    class FakeSVI:
        def __init__(self, *args, **kwargs) -> None:
            pass

        def step(self, **kwargs) -> float:
            return next(losses)

    monkeypatch.setattr(model, "_model_kwargs", lambda data: {})
    monkeypatch.setattr(
        model,
        "_select_svi_initialization",
        lambda data, model_kw, elbo, init_restarts: None,
    )
    monkeypatch.setattr(
        "gatk_sv_ploidy.models.pyro.optim.LambdaLR",
        lambda config, clip_args=None: FakeScheduler(),
    )
    monkeypatch.setattr("gatk_sv_ploidy.models.TraceEnum_ELBO", lambda: object())
    monkeypatch.setattr("gatk_sv_ploidy.models.SVI", FakeSVI)

    history = model.train(
        object(),
        max_iter=20,
        log_freq=100,
        early_stopping=True,
        patience=2,
        convergence_window=2,
        convergence_rtol=2e-4,
    )

    assert history == pytest.approx([100.0, 100.0, 100.01, 100.01, 100.01001])
    assert model.loss_history["epoch"] == [0, 1, 2, 3, 4]


def test_train_passes_gradient_clip_norm_to_scheduler(monkeypatch) -> None:
    model = CNVModel()
    captured: dict[str, object] = {}

    class FakeScheduler:
        def step(self) -> None:
            return None

    class FakeSVI:
        def __init__(self, *args, **kwargs) -> None:
            pass

        def step(self, **kwargs) -> float:
            return 100.0

    def fake_lambda_lr(config, clip_args=None):
        captured["config"] = config
        captured["clip_args"] = clip_args
        return FakeScheduler()

    monkeypatch.setattr(model, "_model_kwargs", lambda data: {})
    monkeypatch.setattr(
        model,
        "_select_svi_initialization",
        lambda data, model_kw, elbo, init_restarts: None,
    )
    monkeypatch.setattr("gatk_sv_ploidy.models.pyro.optim.LambdaLR", fake_lambda_lr)
    monkeypatch.setattr("gatk_sv_ploidy.models.TraceEnum_ELBO", lambda: object())
    monkeypatch.setattr("gatk_sv_ploidy.models.SVI", FakeSVI)

    history = model.train(
        object(),
        max_iter=1,
        log_freq=100,
        early_stopping=False,
        grad_clip_norm=3.5,
    )

    assert history == pytest.approx([100.0])
    assert captured["clip_args"] == {"clip_norm": pytest.approx(3.5)}


def test_train_disables_gradient_clipping_for_non_positive_threshold(monkeypatch) -> None:
    model = CNVModel()
    captured: dict[str, object] = {}

    class FakeScheduler:
        def step(self) -> None:
            return None

    class FakeSVI:
        def __init__(self, *args, **kwargs) -> None:
            pass

        def step(self, **kwargs) -> float:
            return 100.0

    def fake_lambda_lr(config, clip_args=None):
        captured["clip_args"] = clip_args
        return FakeScheduler()

    monkeypatch.setattr(model, "_model_kwargs", lambda data: {})
    monkeypatch.setattr(
        model,
        "_select_svi_initialization",
        lambda data, model_kw, elbo, init_restarts: None,
    )
    monkeypatch.setattr("gatk_sv_ploidy.models.pyro.optim.LambdaLR", fake_lambda_lr)
    monkeypatch.setattr("gatk_sv_ploidy.models.TraceEnum_ELBO", lambda: object())
    monkeypatch.setattr("gatk_sv_ploidy.models.SVI", FakeSVI)

    history = model.train(
        object(),
        max_iter=1,
        log_freq=100,
        early_stopping=False,
        grad_clip_norm=0.0,
    )

    assert history == pytest.approx([100.0])
    assert captured["clip_args"] is None


def test_sample_depth_prior_uses_per_sample_anchor() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chr21", "chr21"],
            "Start": [0, 1000],
            "End": [1000, 2000],
            "S1": [100, 100],
            "S2": [200, 200],
        }
    )
    df["Bin"] = ["chr21:0-1000", "chr21:1000-2000"]
    data = DepthData(
        df.set_index("Bin"),
        depth_space="raw",
        clamp_threshold=None,
        device="cpu",
        dtype=torch.float64,
    )

    model = CNVModel(
        obs_likelihood="negative_binomial",
        dtype=torch.float64,
        device="cpu",
    )

    observed_center = model._estimate_sample_depth_center(data)
    sample_depth_prior = model._estimate_sample_depth_prior(data)

    np.testing.assert_allclose(
        np.array([100.0, 200.0]),
        observed_center,
        rtol=1e-6,
        atol=1e-6,
    )
    np.testing.assert_allclose(
        sample_depth_prior["center"].cpu().numpy(),
        np.array([100.0, 200.0]),
        rtol=1e-6,
        atol=1e-6,
    )
    assert sample_depth_prior["scale"][1].item() > 0.0


def test_raw_expected_depth_units_assumes_diploid_baseline() -> None:
    copy_number = np.array([0.0, 2.0, 4.0], dtype=np.float64)
    bin_bias = np.ones_like(copy_number)
    additive_background = np.zeros_like(copy_number)

    expected = _raw_expected_depth_units(copy_number, bin_bias, additive_background)

    np.testing.assert_allclose(expected, np.array([0.0, 1.0, 2.0], dtype=np.float64))


def test_model_trace_supports_negative_binomial_cn_sampling() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chr21", "chr21"],
            "Start": [0, 1000],
            "End": [1000, 2000],
            "S1": [30.0, 30.0],
        }
    )
    df["Bin"] = ["chr21:0-1000", "chr21:1000-2000"]
    data = DepthData(
        df.set_index("Bin"),
        device="cpu",
        dtype=torch.float64,
    )
    model = CNVModel(
        obs_likelihood="negative_binomial",
        sex_cn_weight=0.0,
        epsilon_mean=0.0,
        dtype=torch.float64,
        device="cpu",
    )

    trace = poutine.trace(model.model).get_trace(**model._model_kwargs(data))

    assert trace.nodes["cn_probs"]["value"].shape == (data.n_bins, 1, model.n_states)
    assert trace.nodes["cn"]["value"].shape == (data.n_bins, data.n_samples)


def test_traceenum_elbo_supports_negative_binomial_and_sex_together() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chr21", "chrX", "chrY"],
            "Start": [0, 0, 0],
            "End": [100, 100, 100],
            "S1": [100.0, 50.0, 10.0],
            "S2": [110.0, 5.0, 60.0],
        }
    )
    df["Bin"] = ["chr21:0-100", "chrX:0-100", "chrY:0-100"]
    data = DepthData(
        df.set_index("Bin"),
        depth_space="raw",
        clamp_threshold=None,
        device="cpu",
        dtype=torch.float64,
    )
    model = CNVModel(
        obs_likelihood="negative_binomial",
        dtype=torch.float64,
        device="cpu",
    )
    guide = model._build_guide(model._make_init_loc_fn(data))
    kwargs = model._model_kwargs(data)
    loss = TraceEnum_ELBO(max_plate_nesting=2).loss(model.model, guide, **kwargs)

    assert np.isfinite(loss)


def test_discrete_inference_returns_cn_posterior() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chr21", "chr21"],
            "Start": [0, 1000],
            "End": [1000, 2000],
            "S1": [30.0, 30.0],
        }
    )
    df["Bin"] = ["chr21:0-1000", "chr21:1000-2000"]
    data = DepthData(
        df.set_index("Bin"),
        device="cpu",
        dtype=torch.float64,
    )
    model = CNVModel(
        obs_likelihood="negative_binomial",
        sex_cn_weight=0.0,
        epsilon_mean=0.0,
        dtype=torch.float64,
        device="cpu",
    )
    maps = {
        "bin_bias": np.ones(data.n_bins, dtype=np.float64),
        "sample_var": np.array([1e-6], dtype=np.float64),
        "bin_var": np.full(data.n_bins, 1e-6, dtype=np.float64),
        "sample_depth": np.array([30.0], dtype=np.float64),
        "cn_probs": np.tile(
            np.array([[0.02, 0.03, 0.88, 0.04, 0.02, 0.01]], dtype=np.float64),
            (data.n_bins, 1),
        ),
    }

    posterior = model._run_discrete_inference_fixed_latents(data, maps)

    assert posterior["cn_posterior"].shape == (data.n_bins, data.n_samples, model.n_states)
    assert set(posterior) == {"cn_posterior"}


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


def test_resolve_fixed_site_pop_af_uses_leave_one_out_for_self_pooled_input() -> None:
    site_alt = torch.tensor([[[10.0, 0.0]]], dtype=torch.float64)
    site_total = torch.tensor([[[10.0, 10.0]]], dtype=torch.float64)
    site_mask = torch.tensor([[[True, True]]], dtype=torch.bool)
    pooled_site_pop_af = torch.tensor([[(10.0 + 0.5) / (20.0 + 1.0)]], dtype=torch.float64)

    resolved_torch, used_torch = _resolve_fixed_site_pop_af_torch(
        site_alt,
        site_total,
        pooled_site_pop_af,
        site_mask,
    )
    resolved_numpy, used_numpy = _resolve_fixed_site_pop_af_numpy(
        site_alt.cpu().numpy(),
        site_total.cpu().numpy(),
        pooled_site_pop_af.cpu().numpy(),
        site_mask.cpu().numpy(),
    )

    expected = np.array([[[0.5 / 11.0, 10.5 / 11.0]]], dtype=np.float64)
    assert used_torch is True
    assert used_numpy is True
    np.testing.assert_allclose(resolved_torch.cpu().numpy(), expected, rtol=1e-8, atol=1e-8)
    np.testing.assert_allclose(resolved_numpy, expected, rtol=1e-8, atol=1e-8)


def test_prepare_af_table_uses_leave_one_out_for_self_pooled_site_pop_af() -> None:
    site_alt = torch.tensor([[[10.0, 0.0]]], dtype=torch.float64)
    site_total = torch.tensor([[[10.0, 10.0]]], dtype=torch.float64)
    site_mask = torch.tensor([[[True, True]]], dtype=torch.bool)
    pooled_site_pop_af = torch.tensor([[(10.0 + 0.5) / (20.0 + 1.0)]], dtype=torch.float64)
    expected_site_pop_af = torch.tensor(
        [[[0.5 / 11.0, 10.5 / 11.0]]],
        dtype=torch.float64,
    )
    model = CNVModel(
        af_evidence_mode="absolute",
        sex_cn_weight=0.0,
        device="cpu",
        dtype=torch.float64,
    )

    auto_table = model._prepare_af_table_torch(
        site_alt,
        site_total,
        pooled_site_pop_af,
        site_mask,
        chr_type=None,
    )
    expected_table = _precompute_af_table(
        site_alt,
        site_total,
        expected_site_pop_af,
        site_mask,
        n_states=model.n_states,
        concentration=model.af_concentration,
    )

    torch.testing.assert_close(auto_table, expected_table, rtol=1e-6, atol=1e-6)


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


def test_cnv_model_uses_mixed_guide_when_learning_site_af() -> None:
    model = CNVModel(guide_type="lowrank", learn_site_pop_af=True)
    assert model.guide.__class__.__name__ == "AutoGuideList"
    assert "site_pop_af_latent" in model._latent_sites


def test_cnv_model_negative_binomial_adds_sample_depth_latent() -> None:
    model = CNVModel(obs_likelihood="negative_binomial")
    assert "sample_depth" in model._latent_sites
    assert "bin_epsilon" in model._latent_sites
    assert "background_bin_factors" in model._latent_sites
    assert "background_sample_factors" in model._latent_sites


def test_cnv_model_can_disable_background_factors() -> None:
    model = CNVModel(background_factors=0)

    assert "background_bin_factors" not in model._latent_sites
    assert "background_sample_factors" not in model._latent_sites


def test_freeze_sample_depth_removes_latent_and_fixes_anchor() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chr21", "chr21", "chrX"],
            "Start": [0, 1000, 0],
            "End": [1000, 2000, 1000],
            "S1": [10, 12, 5],
            "S2": [20, 24, 10],
        },
        index=["chr21:0-1000", "chr21:1000-2000", "chrX:0-1000"],
    )
    data = DepthData(
        df,
        device="cpu",
        depth_space="raw",
        clamp_threshold=None,
        dtype=torch.float64,
    )

    model = CNVModel(
        obs_likelihood="negative_binomial",
        device="cpu",
        dtype=torch.float64,
        freeze_sample_depth=True,
    )
    assert "sample_depth" not in model._latent_sites
    model_kw = model._model_kwargs(data)
    assert "sample_depth_fixed" in model_kw
    np.testing.assert_allclose(
        model_kw["sample_depth_fixed"].cpu().numpy(),
        np.array([11.0, 22.0], dtype=np.float64),
    )

    model.train(data, max_iter=1, init_restarts=3, log_freq=1, early_stopping=False)
    estimates = model.get_map_estimates(data, estimate_method="median")

    np.testing.assert_allclose(
        estimates["sample_depth"],
        np.array([11.0, 22.0], dtype=np.float64),
    )


def test_freeze_bin_bias_removes_latent_and_fixes_to_one() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chr21", "chr21", "chrX"],
            "Start": [0, 1000, 0],
            "End": [1000, 2000, 1000],
            "S1": [10, 12, 5],
            "S2": [20, 24, 10],
        },
        index=["chr21:0-1000", "chr21:1000-2000", "chrX:0-1000"],
    )
    data = DepthData(
        df,
        device="cpu",
        depth_space="raw",
        clamp_threshold=None,
        dtype=torch.float64,
    )

    model = CNVModel(
        obs_likelihood="negative_binomial",
        device="cpu",
        dtype=torch.float64,
        freeze_bin_bias=True,
    )
    assert "bin_bias" not in model._latent_sites
    model_kw = model._model_kwargs(data)
    assert "bin_bias_fixed" in model_kw
    np.testing.assert_allclose(
        model_kw["bin_bias_fixed"].cpu().numpy().squeeze(),
        np.ones(3, dtype=np.float64),
    )

    model.train(data, max_iter=1, log_freq=1, early_stopping=False)
    estimates = model.get_map_estimates(data, estimate_method="median")

    np.testing.assert_allclose(
        np.asarray(estimates["bin_bias"]).squeeze(),
        np.ones(3, dtype=np.float64),
    )


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


def test_estimate_sample_depth_prior_centers_on_autosomal_medians() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chr21", "chr21", "chr21", "chrX"],
            "Start": [0, 1000, 2000, 0],
            "End": [1000, 2000, 3000, 1000],
            "S1": [10, 12, 14, 5],
            "S2": [20, 24, 28, 10],
        },
        index=["chr21:0-1000", "chr21:1000-2000", "chr21:2000-3000", "chrX:0-1000"],
    )
    data = DepthData(
        df,
        device="cpu",
        depth_space="raw",
        clamp_threshold=None,
    )
    model = CNVModel(obs_likelihood="negative_binomial", device="cpu")

    prior = model._estimate_sample_depth_prior(data)
    center = prior["center"].cpu().numpy()
    raw_sd = prior["raw_sd"].cpu().numpy()
    loc = prior["loc"].cpu().numpy()
    scale = prior["scale"].cpu().numpy()
    mean = np.exp(loc + 0.5 * scale ** 2)

    np.testing.assert_allclose(center, np.array([12.0, 24.0], dtype=np.float32))
    np.testing.assert_allclose(mean, center, rtol=1e-6, atol=1e-6)
    assert raw_sd.shape == (2,)
    assert np.all(raw_sd > 0.0)
    assert raw_sd[1] > raw_sd[0]


def test_estimate_sample_depth_prior_bootstrap_sd_tracks_sample_distribution() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chr21", "chr21", "chr21", "chr21", "chr21"],
            "Start": [0, 1000, 2000, 3000, 4000],
            "End": [1000, 2000, 3000, 4000, 5000],
            "S1": [10, 10, 10, 10, 10],
            "S2": [10, 10, 30, 30, 30],
        },
        index=[
            "chr21:0-1000",
            "chr21:1000-2000",
            "chr21:2000-3000",
            "chr21:3000-4000",
            "chr21:4000-5000",
        ],
    )
    data = DepthData(
        df,
        device="cpu",
        depth_space="raw",
        clamp_threshold=None,
    )
    model = CNVModel(obs_likelihood="negative_binomial", device="cpu")

    prior = model._estimate_sample_depth_prior(data)
    center = prior["center"].cpu().numpy()
    raw_sd = prior["raw_sd"].cpu().numpy()

    np.testing.assert_allclose(center, np.array([10.0, 30.0], dtype=np.float32))
    assert raw_sd[0] >= 0.5
    assert raw_sd[1] > raw_sd[0]


def test_model_kwargs_include_sample_depth_prior_for_negative_binomial() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chr21"],
            "Start": [0],
            "End": [1000],
            "S1": [10],
        },
        index=["chr21:0-1000"],
    )
    data = DepthData(
        df,
        device="cpu",
        depth_space="raw",
        clamp_threshold=None,
    )
    model = CNVModel(obs_likelihood="negative_binomial", device="cpu")

    model_kw = model._model_kwargs(data)

    assert "sample_depth_prior_loc" in model_kw
    assert "sample_depth_prior_scale" in model_kw
    assert model_kw["sample_depth_prior_loc"].shape == (1,)
    assert model_kw["sample_depth_prior_scale"].shape == (1,)


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


def test_run_discrete_inference_supports_sample_specific_bin_epsilon() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chrY"],
            "Start": [25],
            "End": [125],
            "S1": [0.2],
            "S2": [0.2],
        }
    )
    df["Bin"] = "chrY:25-125"
    data = DepthData(df.set_index("Bin"), device="cpu")
    model = CNVModel(sex_cn_weight=0.0, epsilon_mean=0.1)

    posterior = model.run_discrete_inference(
        data,
        map_estimates={
            "bin_bias": np.array([1.0], dtype=np.float32),
            "sample_var": np.array([0.01, 0.01], dtype=np.float32),
            "bin_var": np.array([0.01], dtype=np.float32),
            "bin_epsilon": np.array([[0.0, 0.3]], dtype=np.float32),
            "cn_probs": np.full((1, 6), 1.0 / 6.0, dtype=np.float32),
        },
    )

    assert posterior["cn_posterior"][0, 1, 0] > posterior["cn_posterior"][0, 0, 0]


def test_run_discrete_inference_supports_contig_shared_bin_epsilon() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chrY", "chrY"],
            "Start": [25, 125],
            "End": [125, 225],
            "S1": [0.2, 0.2],
        }
    )
    df["Bin"] = df["Chr"].astype(str) + ":" + df["Start"].astype(str) + "-" + df["End"].astype(str)
    data = DepthData(df.set_index("Bin"), device="cpu")
    model = CNVModel(sex_cn_weight=0.0, epsilon_mean=0.1)
    base_map = {
        "bin_bias": np.ones(data.n_bins, dtype=np.float32),
        "sample_var": np.array([0.01], dtype=np.float32),
        "bin_var": np.full(data.n_bins, 0.01, dtype=np.float32),
        "cn_probs": np.full((data.n_bins, 6), 1.0 / 6.0, dtype=np.float32),
    }

    no_epsilon = model.run_discrete_inference(data, map_estimates=base_map)
    shared_epsilon = model.run_discrete_inference(
        data,
        map_estimates={
            **base_map,
            "bin_epsilon": np.array([0.2], dtype=np.float32),
        },
    )

    assert shared_epsilon["cn_posterior"][0, 0, 0] > no_epsilon["cn_posterior"][0, 0, 0]
    assert shared_epsilon["cn_posterior"][1, 0, 0] > no_epsilon["cn_posterior"][1, 0, 0]
    assert shared_epsilon["cn_posterior"][0, 0, 0] == pytest.approx(
        shared_epsilon["cn_posterior"][1, 0, 0]
    )


def test_run_discrete_inference_uses_background_factors() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chrY"],
            "Start": [25],
            "End": [125],
            "S1": [0.05],
            "S2": [0.25],
        }
    )
    df["Bin"] = "chrY:25-125"
    data = DepthData(df.set_index("Bin"), device="cpu")
    model = CNVModel(sex_cn_weight=0.0, epsilon_mean=0.0)

    posterior = model.run_discrete_inference(
        data,
        map_estimates={
            "bin_bias": np.array([1.0], dtype=np.float32),
            "sample_var": np.array([0.01, 0.01], dtype=np.float32),
            "bin_var": np.array([0.01], dtype=np.float32),
            "background_bin_factors": np.array([[2.0]], dtype=np.float32),
            "background_sample_factors": np.array([[0.0, 0.25]], dtype=np.float32),
            "cn_probs": np.full((1, 6), 1.0 / 6.0, dtype=np.float32),
        },
    )

    assert posterior["cn_posterior"][0, 1, 0] > posterior["cn_posterior"][0, 0, 0]


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


def test_negative_binomial_model_trains_one_step_with_sample_depth_prior() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chr21", "chr21"],
            "Start": [0, 1000],
            "End": [1000, 2000],
            "S1": [10, 12],
            "S2": [20, 19],
        },
        index=["chr21:0-1000", "chr21:1000-2000"],
    )
    data = DepthData(
        df,
        device="cpu",
        depth_space="raw",
        clamp_threshold=None,
    )
    model = CNVModel(
        obs_likelihood="negative_binomial",
        guide_type="diagonal",
        sex_cn_weight=0.0,
        autosome_prior_mode="dirichlet",
        device="cpu",
    )

    model.train(data, max_iter=1, log_freq=1, early_stopping=False)
    estimates = model.get_map_estimates(data, estimate_method="median")

    assert "sample_depth" in estimates
    assert estimates["sample_depth"].shape == (2,)


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


def test_model_kwargs_skip_precomputed_af_table_when_learning_site_af(
    tiny_depth_df: pd.DataFrame,
    tiny_site_data: dict[str, np.ndarray],
) -> None:
    data = DepthData(tiny_depth_df, device="cpu", site_data=tiny_site_data)
    model = CNVModel(
        sex_cn_weight=0.0,
        epsilon_mean=0.0,
        guide_type="diagonal",
        learn_site_pop_af=True,
        af_weight=0.5,
    )

    model_kw = model._model_kwargs(data)

    assert "af_table" not in model_kw
    assert model_kw["site_pop_af"].shape == data.site_pop_af.shape


def test_learn_site_pop_af_trains_with_mixed_guide(
    tiny_depth_df: pd.DataFrame,
    tiny_site_data: dict[str, np.ndarray],
) -> None:
    data = DepthData(tiny_depth_df, device="cpu", site_data=tiny_site_data)
    model = CNVModel(
        sex_cn_weight=0.0,
        epsilon_mean=0.0,
        guide_type="diagonal",
        autosome_prior_mode="dirichlet",
        learn_site_pop_af=True,
        af_weight=0.5,
    )

    model.train(data, max_iter=1, log_freq=1, early_stopping=False)
    estimates = model.get_map_estimates(data, estimate_method="median")

    assert "site_pop_af_latent" in estimates
    assert estimates["site_pop_af_latent"].shape == tiny_site_data["site_pop_af"].shape


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
        background_factors=0,
        guide_type="diagonal",
        autosome_prior_mode="dirichlet",
        learn_af_temperature=False,
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
