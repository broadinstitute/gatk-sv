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
    DEFAULT_AF_BACKGROUND_CONCENTRATION,
    DEFAULT_AF_OUTLIER_WEIGHT,
    DEFAULT_BACKGROUND_FACTORS,
    DEFAULT_EPSILON_CONCENTRATION,
    DEFAULT_EPSILON_MEAN,
    DEFAULT_GUIDE_INIT_SCALE,
    DEFAULT_MULTIPLICATIVE_FACTORS,
    DEFAULT_RAW_VARIANCE_POWER,
    _center_af_table_numpy,
    _effective_negative_binomial_overdispersion_numpy,
    _precompute_af_table_from_observed_genotype_log_lik,
    _precompute_af_genotype_log_lik,
    _precompute_af_table_from_genotype_log_lik,
    _precompute_observed_af_genotype_log_lik,
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
    assert model.var_bin == 0.0
    assert model.raw_variance_power == pytest.approx(DEFAULT_RAW_VARIANCE_POWER)
    assert model.epsilon_mean == pytest.approx(DEFAULT_EPSILON_MEAN)
    assert model.epsilon_concentration == pytest.approx(DEFAULT_EPSILON_CONCENTRATION)
    assert DEFAULT_BACKGROUND_FACTORS == 0
    assert DEFAULT_MULTIPLICATIVE_FACTORS == 0
    assert model.background_factors == DEFAULT_BACKGROUND_FACTORS
    assert model.multiplicative_factors == DEFAULT_MULTIPLICATIVE_FACTORS
    assert model.guide_init_scale == pytest.approx(DEFAULT_GUIDE_INIT_SCALE)
    assert model.lowrank_guide_rank is None
    assert model.af_evidence_mode == "relative"
    assert model.learn_af_temperature is False
    assert model.af_outlier_weight == pytest.approx(DEFAULT_AF_OUTLIER_WEIGHT)
    assert model.af_background_concentration == pytest.approx(
        DEFAULT_AF_BACKGROUND_CONCENTRATION
    )
    assert "bin_var" not in model._latent_sites
    assert "allosome_var" not in model._latent_sites
    assert "bin_bias" not in model._latent_sites
    assert "multiplicative_bin_factors" not in model._latent_sites
    assert "multiplicative_sample_factors" not in model._latent_sites
    assert "af_temperature" not in model._latent_sites
    assert "background_bin_factors" not in model._latent_sites
    assert "background_sample_factors" not in model._latent_sites

@pytest.mark.parametrize("raw_variance_power", [0.99, 2.01])
def test_cnv_model_rejects_invalid_raw_variance_power(
    raw_variance_power: float,
) -> None:
    with pytest.raises(ValueError, match="raw_variance_power"):
        CNVModel(raw_variance_power=raw_variance_power)


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
            "Chr": ["chr21", "chr21"],
            "Start": [0, 1000],
            "End": [1000, 2000],
            "S1": [2.0, 3.0],
            "S2": [2.5, 3.5],
        }
    )
    df["Bin"] = ["chr21:0-1000", "chr21:1000-2000"]
    data = DepthData(df.set_index("Bin"), device="cpu")
    model = CNVModel(
        sex_cn_weight=0.0,
        epsilon_mean=0.123,
        epsilon_concentration=0.5,
        device="cpu",
    )

    trace = poutine.trace(model.model).get_trace(**model._model_kwargs(data))
    bin_epsilon_fn = trace.nodes["bin_epsilon"]["fn"]
    bin_epsilon_value = trace.nodes["bin_epsilon"]["value"]
    gamma_dist = getattr(bin_epsilon_fn, "base_dist", bin_epsilon_fn)

    assert gamma_dist.__class__.__name__ == "Gamma"
    assert tuple(bin_epsilon_value.shape) == (data.n_bins, data.n_samples)
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

    def fake_make_init_loc_fn(
        data,
        randomize: bool = False,
        initial_values=None,
    ):
        candidate = len(init_calls)
        init_calls.append(randomize)
        return {"candidate": candidate}

    def fake_build_guide(init_loc_fn, guide_type=None):
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


def test_select_svi_initialization_uses_single_anchor_for_expressive_guides(monkeypatch) -> None:
    model = CNVModel(guide_type="diagonal")
    init_calls = []

    class DummyGuide:
        candidate = 0

    class FakeElbo:
        def loss(self, model_fn, guide, **kwargs):
            return 1.0

    def fake_make_init_loc_fn(
        data,
        randomize: bool = False,
        initial_values=None,
    ):
        init_calls.append((randomize, initial_values))
        return {"candidate": len(init_calls) - 1}

    monkeypatch.setattr(model, "_make_init_loc_fn", fake_make_init_loc_fn)
    monkeypatch.setattr(model, "_build_guide", lambda init_loc_fn, guide_type=None: DummyGuide())

    model._select_svi_initialization(
        object(),
        model_kw={},
        elbo=FakeElbo(),
        init_restarts=10,
    )

    assert init_calls == [(False, None)]


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
        lambda data, model_kw, elbo, init_restarts, **kwargs: None,
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
        lambda data, model_kw, elbo, init_restarts, **kwargs: None,
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
        lambda data, model_kw, elbo, init_restarts, **kwargs: None,
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

    model = CNVModel(dtype=torch.float64, device="cpu")

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


def test_raw_expected_depth_units_assumes_diploid_baseline_and_cn0_background() -> None:
    copy_number = np.array([0.0, 1.0, 2.0, 4.0], dtype=np.float64)
    bin_bias = np.ones_like(copy_number)
    additive_background = np.full_like(copy_number, 0.2)

    expected = _raw_expected_depth_units(copy_number, bin_bias, additive_background)

    np.testing.assert_allclose(
        expected,
        np.array([0.1, 0.5, 1.0, 2.0], dtype=np.float64),
    )

    expected_torch = _raw_expected_depth_units(
        torch.tensor(copy_number),
        torch.tensor(bin_bias),
        torch.tensor(additive_background),
    )
    np.testing.assert_allclose(expected_torch.numpy(), expected)


def test_raw_expected_depth_units_supports_sample_specific_baselines() -> None:
    copy_number = np.array(
        [[2.0, 3.0], [4.0, 0.0]],
        dtype=np.float64,
    )
    bin_bias = np.ones_like(copy_number)
    additive_background = np.full_like(copy_number, 0.3)
    baseline = np.array([2.0, 3.0], dtype=np.float64)

    expected = _raw_expected_depth_units(
        copy_number,
        bin_bias,
        additive_background,
        autosomal_baseline_copy_number=baseline,
    )

    np.testing.assert_allclose(
        expected,
        np.array([[1.0, 1.0], [2.0, 0.1]], dtype=np.float64),
    )

    expected_torch = _raw_expected_depth_units(
        torch.tensor(copy_number),
        torch.tensor(bin_bias),
        torch.tensor(additive_background),
        autosomal_baseline_copy_number=torch.tensor(baseline),
    )
    np.testing.assert_allclose(expected_torch.numpy(), expected)


def test_compose_autosome_cn_probs_supports_sample_specific_reference_states() -> None:
    model = CNVModel(dtype=torch.float64, device="cpu")
    autosome_nonref_prob = torch.tensor([0.3, 0.4], dtype=torch.float64)
    autosome_alt_cn_probs = torch.tensor(
        [
            [0.05, 0.10, 0.20, 0.25, 0.40],
            [0.05, 0.10, 0.20, 0.25, 0.40],
        ],
        dtype=torch.float64,
    )

    composed_torch = model._compose_autosome_cn_probs_torch(
        autosome_nonref_prob,
        autosome_alt_cn_probs,
        reference_state=torch.tensor([2, 3]),
    )
    expected = np.array(
        [
            [0.015, 0.03, 0.7, 0.06, 0.075, 0.12],
            [0.02, 0.04, 0.08, 0.6, 0.10, 0.16],
        ],
        dtype=np.float64,
    )
    np.testing.assert_allclose(composed_torch.cpu().numpy(), expected)

    composed_numpy = model._compose_autosome_cn_probs_numpy(
        autosome_nonref_prob.cpu().numpy(),
        autosome_alt_cn_probs.cpu().numpy(),
        reference_state=np.array([2, 3], dtype=np.int64),
    )
    np.testing.assert_allclose(composed_numpy, expected)


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
    model = CNVModel(dtype=torch.float64, device="cpu")
    guide = model._build_guide(model._make_init_loc_fn(data))
    kwargs = model._model_kwargs(data)
    loss = TraceEnum_ELBO(max_plate_nesting=2).loss(model.model, guide, **kwargs)

    assert np.isfinite(loss)


def test_traceenum_elbo_supports_joint_af_and_sex_factors() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chrX", "chrX"],
            "Start": [82462542, 86462542],
            "End": [86462542, 94462542],
            "S1": [899291.0, 893183.0],
            "S2": [450085.0, 447599.0],
        },
        index=[
            "chrX:82462542-86462542",
            "chrX:86462542-94462542",
        ],
    )
    data = DepthData(
        df,
        depth_space="raw",
        clamp_threshold=None,
        device="cpu",
        dtype=torch.float64,
    )
    model = CNVModel(
        sex_cn_weight=3.5,
        af_weight=0.8,
        learn_af_temperature=True,
        device="cpu",
        dtype=torch.float64,
    )
    kwargs = model._model_kwargs(data)
    kwargs["af_table"] = torch.tensor(
        [
            [
                [-8.731556539984444e01, -6.612377062059409e01],
                [-1.368110966644657e02, -1.814617429338733e02],
            ],
            [
                [-1.488844565971347e02, -1.656738939179022e01],
                [-3.580629440725404e02, 1.632826355252483e00],
            ],
            [
                [1.790295811632028e00, -1.692863593990445e00],
                [1.789305086134846e00, -1.259268134340488e-01],
            ],
            [
                [-4.993605053063245e00, 1.707493550236123e00],
                [-4.247519362663468e00, -1.576731826990715e01],
            ],
            [
                [-6.255349242860806e00, -1.214220162548450e00],
                [-7.804633313324430e00, -3.016906861661726e01],
            ],
            [
                [-9.508844668952335e00, -5.535106331824970e00],
                [-1.411907518276722e01, -4.268211546796122e01],
            ],
        ],
        dtype=torch.float64,
    )

    pyro.clear_param_store()
    guide = model._build_guide(model._make_init_loc_fn(data))
    loss = TraceEnum_ELBO().loss(model.model, guide, **kwargs)

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


def test_marginalized_af_log_lik_matches_numpy_with_outlier_component(
    tiny_site_data: dict[str, np.ndarray],
) -> None:
    cn = torch.full((4, 2), 2, dtype=torch.long)
    torch_result = _marginalized_af_log_lik(
        torch.tensor(tiny_site_data["site_alt"], dtype=torch.float32),
        torch.tensor(tiny_site_data["site_total"], dtype=torch.float32),
        torch.tensor(tiny_site_data["site_pop_af"], dtype=torch.float32),
        torch.tensor(tiny_site_data["site_mask"], dtype=torch.bool),
        cn=cn,
        n_states=6,
        concentration=50.0,
        outlier_weight=0.05,
    ).cpu().numpy()

    numpy_result = _marginalized_af_log_lik_numpy(
        tiny_site_data["site_alt"],
        tiny_site_data["site_total"],
        tiny_site_data["site_pop_af"],
        tiny_site_data["site_mask"],
        cn_state=2,
        n_states=6,
        concentration=50.0,
        outlier_weight=0.05,
    )

    np.testing.assert_allclose(torch_result, numpy_result, rtol=1e-5, atol=1e-5)


def test_marginalized_af_log_lik_outlier_component_softens_mismatch() -> None:
    site_alt = np.array([[[9]]], dtype=np.int32)
    site_total = np.array([[[10]]], dtype=np.int32)
    site_pop_af = np.array([[0.5]], dtype=np.float32)
    site_mask = np.array([[[True]]], dtype=bool)

    no_outlier = _marginalized_af_log_lik_numpy(
        site_alt,
        site_total,
        site_pop_af,
        site_mask,
        cn_state=2,
        n_states=6,
        concentration=400.0,
        outlier_weight=0.0,
    )
    with_outlier = _marginalized_af_log_lik_numpy(
        site_alt,
        site_total,
        site_pop_af,
        site_mask,
        cn_state=2,
        n_states=6,
        concentration=400.0,
        outlier_weight=0.10,
    )

    assert with_outlier[0, 0] > no_outlier[0, 0]


def test_marginalized_af_log_lik_matches_scalar_when_using_repeated_sample_concentration(
    tiny_site_data: dict[str, np.ndarray],
) -> None:
    cn = torch.full((4, 2), 2, dtype=torch.long)
    scalar_torch = _marginalized_af_log_lik(
        torch.tensor(tiny_site_data["site_alt"], dtype=torch.float32),
        torch.tensor(tiny_site_data["site_total"], dtype=torch.float32),
        torch.tensor(tiny_site_data["site_pop_af"], dtype=torch.float32),
        torch.tensor(tiny_site_data["site_mask"], dtype=torch.bool),
        cn=cn,
        n_states=6,
        concentration=50.0,
    ).cpu().numpy()
    vector_torch = _marginalized_af_log_lik(
        torch.tensor(tiny_site_data["site_alt"], dtype=torch.float32),
        torch.tensor(tiny_site_data["site_total"], dtype=torch.float32),
        torch.tensor(tiny_site_data["site_pop_af"], dtype=torch.float32),
        torch.tensor(tiny_site_data["site_mask"], dtype=torch.bool),
        cn=cn,
        n_states=6,
        concentration=torch.full((2,), 50.0, dtype=torch.float32),
    ).cpu().numpy()

    scalar_numpy = _marginalized_af_log_lik_numpy(
        tiny_site_data["site_alt"],
        tiny_site_data["site_total"],
        tiny_site_data["site_pop_af"],
        tiny_site_data["site_mask"],
        cn_state=2,
        n_states=6,
        concentration=50.0,
    )
    vector_numpy = _marginalized_af_log_lik_numpy(
        tiny_site_data["site_alt"],
        tiny_site_data["site_total"],
        tiny_site_data["site_pop_af"],
        tiny_site_data["site_mask"],
        cn_state=2,
        n_states=6,
        concentration=np.full(2, 50.0, dtype=np.float64),
    )

    np.testing.assert_allclose(vector_torch, scalar_torch, rtol=1e-5, atol=1e-5)
    np.testing.assert_allclose(vector_numpy, scalar_numpy, rtol=1e-5, atol=1e-5)


def test_marginalized_af_log_lik_sanitizes_invalid_site_counts() -> None:
    site_alt = np.array([[[2]]], dtype=np.int32)
    site_total = np.array([[[1]]], dtype=np.int32)
    site_pop_af = np.array([[0.5]], dtype=np.float32)
    site_mask = np.array([[[True]]], dtype=bool)

    torch_result = _marginalized_af_log_lik(
        torch.tensor(site_alt, dtype=torch.float32),
        torch.tensor(site_total, dtype=torch.float32),
        torch.tensor(site_pop_af, dtype=torch.float32),
        torch.tensor(site_mask, dtype=torch.bool),
        cn=torch.tensor([[2]], dtype=torch.long),
        n_states=6,
        concentration=50.0,
    ).cpu().numpy()

    numpy_result = _marginalized_af_log_lik_numpy(
        site_alt,
        site_total,
        site_pop_af,
        site_mask,
        cn_state=2,
        n_states=6,
        concentration=50.0,
    )

    assert np.isfinite(torch_result).all()
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


def test_precompute_af_table_matches_scalar_when_using_repeated_sample_concentration(
    tiny_site_data: dict[str, np.ndarray],
) -> None:
    scalar_table = _precompute_af_table(
        torch.tensor(tiny_site_data["site_alt"], dtype=torch.float32),
        torch.tensor(tiny_site_data["site_total"], dtype=torch.float32),
        torch.tensor(tiny_site_data["site_pop_af"], dtype=torch.float32),
        torch.tensor(tiny_site_data["site_mask"], dtype=torch.bool),
        n_states=6,
        concentration=50.0,
    )
    vector_table = _precompute_af_table(
        torch.tensor(tiny_site_data["site_alt"], dtype=torch.float32),
        torch.tensor(tiny_site_data["site_total"], dtype=torch.float32),
        torch.tensor(tiny_site_data["site_pop_af"], dtype=torch.float32),
        torch.tensor(tiny_site_data["site_mask"], dtype=torch.bool),
        n_states=6,
        concentration=torch.full((2,), 50.0, dtype=torch.float32),
    )

    torch.testing.assert_close(vector_table, scalar_table, rtol=1e-5, atol=1e-5)


def test_precompute_af_table_softens_only_target_sample_with_sample_specific_concentration() -> None:
    site_alt = torch.tensor([[[18.0, 20.0]]], dtype=torch.float32)
    site_total = torch.tensor([[[20.0, 20.0]]], dtype=torch.float32)
    site_pop_af = torch.tensor([[0.5]], dtype=torch.float32)
    site_mask = torch.tensor([[[True, True]]], dtype=torch.bool)

    scalar_table = _precompute_af_table(
        site_alt,
        site_total,
        site_pop_af,
        site_mask,
        n_states=6,
        concentration=50.0,
    )
    mixed_table = _precompute_af_table(
        site_alt,
        site_total,
        site_pop_af,
        site_mask,
        n_states=6,
        concentration=torch.tensor([5.0, 50.0], dtype=torch.float32),
    )

    assert float(mixed_table[2, 0, 0]) > float(scalar_table[2, 0, 0])
    torch.testing.assert_close(
        mixed_table[:, :, 1],
        scalar_table[:, :, 1],
        rtol=1e-5,
        atol=1e-5,
    )


def test_precomputed_af_genotype_log_lik_rebuilds_af_table(
    tiny_site_data: dict[str, np.ndarray],
) -> None:
    site_alt = torch.tensor(tiny_site_data["site_alt"], dtype=torch.float32)
    site_total = torch.tensor(tiny_site_data["site_total"], dtype=torch.float32)
    site_pop_af = torch.tensor(tiny_site_data["site_pop_af"], dtype=torch.float32)
    site_mask = torch.tensor(tiny_site_data["site_mask"], dtype=torch.bool)

    genotype_log_lik = _precompute_af_genotype_log_lik(
        site_alt,
        site_total,
        site_mask,
        n_states=6,
        concentration=50.0,
    )
    cached_table = _precompute_af_table_from_genotype_log_lik(
        site_pop_af,
        site_alt,
        site_total,
        site_mask,
        genotype_log_lik,
        n_states=6,
    )
    direct_table = _precompute_af_table(
        site_alt,
        site_total,
        site_pop_af,
        site_mask,
        n_states=6,
        concentration=50.0,
    )

    torch.testing.assert_close(cached_table, direct_table, rtol=1e-5, atol=1e-5)

    site_pop_af_grad = site_pop_af.detach().clone().requires_grad_(True)
    gradient_table = _precompute_af_table_from_genotype_log_lik(
        site_pop_af_grad,
        site_alt,
        site_total,
        site_mask,
        genotype_log_lik,
        n_states=6,
    )
    gradient_table.sum().backward()
    assert site_pop_af_grad.grad is not None
    assert torch.isfinite(site_pop_af_grad.grad).all()


def test_precompute_af_genotype_log_lik_sanitizes_invalid_site_counts() -> None:
    site_alt = torch.tensor([[[2]]], dtype=torch.float32)
    site_total = torch.tensor([[[1]]], dtype=torch.float32)
    site_mask = torch.tensor([[[True]]], dtype=torch.bool)

    genotype_log_lik = _precompute_af_genotype_log_lik(
        site_alt,
        site_total,
        site_mask,
        n_states=6,
        concentration=50.0,
    )
    observed_log_lik, *_ = _precompute_observed_af_genotype_log_lik(
        site_alt,
        site_total,
        site_mask,
        n_states=6,
        concentration=50.0,
    )

    assert torch.isfinite(genotype_log_lik).all()
    assert torch.isfinite(observed_log_lik).all()


def test_observed_af_genotype_log_lik_rebuilds_af_table(
    tiny_site_data: dict[str, np.ndarray],
) -> None:
    site_alt = torch.tensor(tiny_site_data["site_alt"], dtype=torch.float32)
    site_total = torch.tensor(tiny_site_data["site_total"], dtype=torch.float32)
    site_pop_af = torch.tensor(tiny_site_data["site_pop_af"], dtype=torch.float32)
    site_mask = torch.tensor(tiny_site_data["site_mask"], dtype=torch.bool)

    (
        observed_genotype_log_lik,
        observed_alt,
        observed_total,
        observed_bin_idx,
        observed_site_idx,
        observed_sample_idx,
        observed_site_slot_idx,
        site_slot_bin_idx,
        site_slot_site_idx,
    ) = _precompute_observed_af_genotype_log_lik(
        site_alt,
        site_total,
        site_mask,
        n_states=6,
        concentration=50.0,
    )
    observed_table = _precompute_af_table_from_observed_genotype_log_lik(
        site_pop_af,
        observed_alt,
        observed_total,
        observed_bin_idx,
        observed_site_idx,
        observed_sample_idx,
        observed_site_slot_idx,
        observed_genotype_log_lik,
        n_bins=site_alt.shape[0],
        n_samples=site_alt.shape[2],
        n_states=6,
    )
    direct_table = _precompute_af_table(
        site_alt,
        site_total,
        site_pop_af,
        site_mask,
        n_states=6,
        concentration=50.0,
    )

    torch.testing.assert_close(observed_table, direct_table, rtol=1e-5, atol=1e-5)

    observed_site_pop_af = (
        site_pop_af[site_slot_bin_idx, site_slot_site_idx]
        .detach()
        .clone()
        .requires_grad_(True)
    )
    gradient_table = _precompute_af_table_from_observed_genotype_log_lik(
        observed_site_pop_af,
        observed_alt,
        observed_total,
        observed_bin_idx,
        observed_site_idx,
        observed_sample_idx,
        observed_site_slot_idx,
        observed_genotype_log_lik,
        n_bins=site_alt.shape[0],
        n_samples=site_alt.shape[2],
        n_states=6,
    )
    gradient_table.sum().backward()
    assert observed_site_pop_af.grad is not None
    assert torch.isfinite(observed_site_pop_af.grad).all()


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
        sex_cn_weight=0.0,
        device="cpu",
    )
    chr_type_np = np.zeros(raw_table.shape[1], dtype=np.int64)

    transformed = model._apply_af_evidence_mode_numpy(raw_table, chr_type_np)
    reference = model._af_reference_cn_probs_numpy(
        chr_type_np,
        raw_table.shape[1],
        raw_table.shape[2],
    )
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
        outlier_weight=model.af_outlier_weight,
        background_concentration=model.af_background_concentration,
        informative_weight=model.af_weight,
    )
    expected_table = model._apply_af_evidence_mode_torch(expected_table, None)

    torch.testing.assert_close(auto_table, expected_table, rtol=1e-6, atol=1e-6)


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


def test_raw_variance_power_rebalances_low_and_high_count_overdispersion() -> None:
    base_overdispersion = np.array([[[0.04, 0.04]]], dtype=np.float64)
    mean = np.array([[[1.0e4, 1.0e6]]], dtype=np.float64)

    nb2 = _effective_negative_binomial_overdispersion_numpy(
        base_overdispersion,
        mean,
        raw_variance_power=2.0,
    )
    power_law = _effective_negative_binomial_overdispersion_numpy(
        base_overdispersion,
        mean,
        raw_variance_power=1.5,
    )

    np.testing.assert_allclose(nb2, base_overdispersion)
    assert float(power_law[0, 0, 0]) > float(power_law[0, 0, 1])
    assert float(power_law[0, 0, 1]) < float(nb2[0, 0, 1])


def test_negative_binomial_log_lik_rejects_large_fractional_counts() -> None:
    with pytest.raises(ValueError, match="integer-valued raw counts"):
        _negative_binomial_log_lik_numpy(
            np.array([[[1_000_000.25]]], dtype=np.float64),
            np.array([[[1_000_000.0]]], dtype=np.float64),
            np.array([[[0.04]]], dtype=np.float64),
        )


def test_run_discrete_inference_uses_learned_af_temperature_map() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chr21"],
            "Start": [0],
            "End": [1000],
            "S1": [20],
        }
    )
    df["Bin"] = "chr21:0-1000"
    site_data = {
        "site_alt": np.asarray([[[5]]], dtype=np.int32),
        "site_total": np.asarray([[[10]]], dtype=np.int32),
        "site_pop_af": np.asarray([[0.5]], dtype=np.float32),
        "site_mask": np.asarray([[[True]]], dtype=bool),
    }
    data = DepthData(
        df.set_index("Bin"),
        device="cpu",
        depth_space="raw",
        clamp_threshold=None,
        site_data=site_data,
    )
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
        "sample_depth": np.array([20.0], dtype=np.float32),
        "cn_probs": np.full((1, 6), 1.0 / 6.0, dtype=np.float32),
    }
    low_temp = model.run_discrete_inference(
        data,
        map_estimates={
            **base_map,
            "af_temperature": np.asarray(0.0, dtype=np.float32),
        },
    )
    high_temp = model.run_discrete_inference(
        data,
        map_estimates={
            **base_map,
            "af_temperature": np.asarray(1.0, dtype=np.float32),
        },
    )

    assert high_temp["cn_posterior"][0, 0, 2] > low_temp["cn_posterior"][0, 0, 2]


def test_cnv_model_accepts_lowrank_guide() -> None:
    model = CNVModel(guide_type="lowrank", lowrank_guide_rank=2)
    assert model.guide.__class__.__name__ == "AutoLowRankMultivariateNormal"
    assert model.guide._init_scale == pytest.approx(DEFAULT_GUIDE_INIT_SCALE)
    assert model.lowrank_guide_rank == 2


def test_cnv_model_passes_init_scale_to_diagonal_guide() -> None:
    model = CNVModel(guide_type="diagonal", guide_init_scale=0.005)

    assert model.guide.__class__.__name__ == "AutoDiagonalNormal"
    assert model.guide._init_scale == pytest.approx(0.005)


def test_cnv_model_uses_mixed_guide_when_learning_site_af() -> None:
    model = CNVModel(guide_type="lowrank", learn_site_pop_af=True)
    assert model.guide.__class__.__name__ == "AutoGuideList"
    assert "site_pop_af_latent_observed" in model._latent_sites


def test_cnv_model_negative_binomial_fixes_sample_depth_by_default() -> None:
    model = CNVModel()
    assert "sample_depth" not in model._latent_sites
    assert "bin_epsilon" in model._latent_sites
    assert "background_bin_factors" not in model._latent_sites
    assert "background_sample_factors" not in model._latent_sites


def test_cnv_model_can_disable_background_factors() -> None:
    model = CNVModel(background_factors=0)

    assert "background_bin_factors" not in model._latent_sites
    assert "background_sample_factors" not in model._latent_sites


def test_sample_depth_uses_autosomal_anchor() -> None:
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
        device="cpu",
        dtype=torch.float64,
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
    model = CNVModel(device="cpu")

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
    model = CNVModel(device="cpu")

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
    model = CNVModel(device="cpu")

    prior = model._estimate_sample_depth_prior(data)
    center = prior["center"].cpu().numpy()
    raw_sd = prior["raw_sd"].cpu().numpy()

    np.testing.assert_allclose(center, np.array([10.0, 30.0], dtype=np.float32))
    assert raw_sd[0] >= 0.5
    assert raw_sd[1] > raw_sd[0]


def test_model_kwargs_include_fixed_sample_depth_by_default_for_negative_binomial() -> None:
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
    model = CNVModel(device="cpu")

    model_kw = model._model_kwargs(data)

    assert "sample_depth_fixed" in model_kw
    assert "sample_depth_prior_loc" not in model_kw
    assert "sample_depth_prior_scale" not in model_kw
    assert model_kw["sample_depth_fixed"].shape == (1,)


def test_run_discrete_inference_uses_bin_epsilon() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chrY"],
            "Start": [25],
            "End": [125],
            "S1": [2],
        }
    )
    df["Bin"] = "chrY:25-125"
    data = DepthData(
        df.set_index("Bin"),
        device="cpu",
        depth_space="raw",
        clamp_threshold=None,
    )
    model = CNVModel(sex_cn_weight=0.0, epsilon_mean=0.1)
    base_map = {
        "bin_bias": np.array([1.0], dtype=np.float32),
        "sample_var": np.array([0.01], dtype=np.float32),
        "bin_var": np.array([0.01], dtype=np.float32),
        "sample_depth": np.array([20.0], dtype=np.float32),
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
            "S1": [2],
            "S2": [2],
        }
    )
    df["Bin"] = "chrY:25-125"
    data = DepthData(
        df.set_index("Bin"),
        device="cpu",
        depth_space="raw",
        clamp_threshold=None,
    )
    model = CNVModel(sex_cn_weight=0.0, epsilon_mean=0.1)

    posterior = model.run_discrete_inference(
        data,
        map_estimates={
            "bin_bias": np.array([1.0], dtype=np.float32),
            "sample_var": np.array([0.01, 0.01], dtype=np.float32),
            "bin_var": np.array([0.01], dtype=np.float32),
            "bin_epsilon": np.array([[0.0, 0.3]], dtype=np.float32),
            "sample_depth": np.array([20.0, 20.0], dtype=np.float32),
            "cn_probs": np.full((1, 6), 1.0 / 6.0, dtype=np.float32),
        },
    )

    assert posterior["cn_posterior"][0, 1, 0] > posterior["cn_posterior"][0, 0, 0]


def test_run_discrete_inference_rejects_contig_shared_bin_epsilon() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chrY", "chrY"],
            "Start": [25, 125],
            "End": [125, 225],
            "S1": [2, 2],
        }
    )
    df["Bin"] = df["Chr"].astype(str) + ":" + df["Start"].astype(str) + "-" + df["End"].astype(str)
    data = DepthData(
        df.set_index("Bin"),
        device="cpu",
        depth_space="raw",
        clamp_threshold=None,
    )
    model = CNVModel(sex_cn_weight=0.0, epsilon_mean=0.1)
    base_map = {
        "bin_bias": np.ones(data.n_bins, dtype=np.float32),
        "sample_var": np.array([0.01], dtype=np.float32),
        "bin_var": np.full(data.n_bins, 0.01, dtype=np.float32),
        "sample_depth": np.array([20.0], dtype=np.float32),
        "cn_probs": np.full((data.n_bins, 6), 1.0 / 6.0, dtype=np.float32),
    }

    with pytest.raises(ValueError, match="length n_bins"):
        model.run_discrete_inference(
            data,
            map_estimates={
                **base_map,
                "bin_epsilon": np.array([0.2], dtype=np.float32),
            },
        )


def test_run_discrete_inference_uses_background_factors() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chrY"],
            "Start": [25],
            "End": [125],
            "S1": [1],
            "S2": [5],
        }
    )
    df["Bin"] = "chrY:25-125"
    data = DepthData(
        df.set_index("Bin"),
        device="cpu",
        depth_space="raw",
        clamp_threshold=None,
    )
    model = CNVModel(sex_cn_weight=0.0, epsilon_mean=0.0)

    posterior = model.run_discrete_inference(
        data,
        map_estimates={
            "bin_bias": np.array([1.0], dtype=np.float32),
            "sample_var": np.array([0.01, 0.01], dtype=np.float32),
            "bin_var": np.array([0.01], dtype=np.float32),
            "background_bin_factors": np.array([[2.0]], dtype=np.float32),
            "background_sample_factors": np.array([[0.0, 0.25]], dtype=np.float32),
            "sample_depth": np.array([20.0, 20.0], dtype=np.float32),
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


def test_run_discrete_inference_negative_binomial_supports_triploid_baseline_cn() -> None:
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
        autosomal_baseline_cn=[3],
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

    assert int(np.argmax(posterior["cn_posterior"][0, 0, :])) == 3
    assert float(posterior["cn_posterior"][0, 0, 3]) > float(
        posterior["cn_posterior"][0, 0, 2]
    )


def test_sex_cn_score_uses_sample_baseline_cn() -> None:
    model = CNVModel(
        sex_cn_weight=3.0,
        epsilon_mean=0.0,
        autosomal_baseline_cn=[2, 3],
        dtype=torch.float64,
    )

    score = model._sex_cn_score_for_samples_numpy(2)

    assert score[1, 1, 1, 0] == 0.0
    assert score[1, 2, 1, 0] == 0.0
    assert score[1, 1, 2, 1] == 0.0
    assert score[1, 2, 1, 1] == 0.0
    assert score[1, 1, 1, 1] < 0.0
    assert score[0, 1, 3, 1] == 0.0
    assert score[0, 2, 0, 1] == 0.0


def test_baseline_cn_expansion_squeezes_singleton_sample_axis() -> None:
    model = CNVModel(
        sex_cn_weight=0.0,
        epsilon_mean=0.0,
        autosomal_baseline_cn=[2, 3],
        dtype=torch.float64,
    )
    canonical = torch.zeros((2, 1, model.n_states), dtype=torch.float64)
    canonical[..., 1] = 0.2
    canonical[..., 2] = 0.8
    chr_type = torch.zeros(2, dtype=torch.long)

    expanded = model._expand_cn_probs_to_samples_torch(
        canonical,
        chr_type,
        n_samples=2,
    )

    assert expanded.shape == (2, 2, model.n_states)
    np.testing.assert_allclose(expanded[:, 0, 2].detach().cpu().numpy(), 0.8)
    np.testing.assert_allclose(expanded[:, 1, 3].detach().cpu().numpy(), 0.8)
    np.testing.assert_allclose(expanded[:, :, 1].detach().cpu().numpy(), 0.2)


def test_elbo_loss_supports_baseline_cn_with_dirichlet_singleton_cn_probs() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chr21", "chr21"],
            "Start": [0, 1000],
            "End": [1000, 2000],
            "S1": [20, 20],
            "S2": [30, 30],
        },
        index=["chr21:0-1000", "chr21:1000-2000"],
    )
    data = DepthData(
        df,
        device="cpu",
        dtype=torch.float64,
        depth_space="raw",
        clamp_threshold=None,
    )
    model = CNVModel(
        sex_cn_weight=0.0,
        epsilon_mean=0.0,
        autosomal_baseline_cn=[2, 3],
        dtype=torch.float64,
        device="cpu",
    )
    model_kw = model._model_kwargs(data)

    pyro.clear_param_store()
    guide = model._build_guide(model._make_init_loc_fn(data))
    loss = TraceEnum_ELBO().loss(model.model, guide, **model_kw)

    assert np.isfinite(loss)


def test_elbo_loss_supports_triploid_baseline_sex_cn_prior() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chrX", "chrY"],
            "Start": [0, 0],
            "End": [1000, 1000],
            "S1": [20, 10],
        },
        index=["chrX:0-1000", "chrY:0-1000"],
    )
    data = DepthData(
        df,
        device="cpu",
        dtype=torch.float64,
        depth_space="raw",
        clamp_threshold=None,
    )
    model = CNVModel(
        sex_cn_weight=3.0,
        epsilon_mean=0.0,
        autosomal_baseline_cn=[3],
        dtype=torch.float64,
        device="cpu",
    )
    model_kw = model._model_kwargs(data)

    pyro.clear_param_store()
    guide = model._build_guide(model._make_init_loc_fn(data))
    loss = TraceEnum_ELBO().loss(model.model, guide, **model_kw)

    assert np.isfinite(loss)


def test_run_discrete_inference_negative_binomial_cn0_background() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chr21"],
            "Start": [0],
            "End": [1000],
            "S1": [2],
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
    )
    base_map = {
        "bin_bias": np.array([1.0], dtype=np.float32),
        "sample_var": np.array([1e-4], dtype=np.float32),
        "bin_var": np.array([0.0], dtype=np.float32),
        "sample_depth": np.array([20.0], dtype=np.float32),
        "cn_probs": np.full((1, 6), 1.0 / 6.0, dtype=np.float32),
    }

    no_background = model.run_discrete_inference(data, map_estimates=base_map)
    with_background = model.run_discrete_inference(
        data,
        map_estimates={
            **base_map,
            "bin_epsilon": np.array([[0.2]], dtype=np.float32),
        },
    )

    assert int(np.argmax(with_background["cn_posterior"][0, 0, :])) == 0
    assert with_background["cn_posterior"][0, 0, 0] > no_background["cn_posterior"][0, 0, 0]
    assert with_background["cn_posterior"][0, 0, 2] < no_background["cn_posterior"][0, 0, 2]


def test_run_discrete_inference_supports_autosome_shrinkage_maps() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chr21"],
            "Start": [0],
            "End": [1000],
            "S1": [20],
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
        autosome_prior_mode="shrinkage",
    )

    common_maps = {
        "bin_bias": np.array([1.0], dtype=np.float32),
        "sample_var": np.array([0.01], dtype=np.float32),
        "bin_var": np.array([0.01], dtype=np.float32),
        "sample_depth": np.array([20.0], dtype=np.float32),
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
        guide_type="diagonal",
        sex_cn_weight=0.0,
        autosome_prior_mode="dirichlet",
        device="cpu",
    )

    model.train(data, max_iter=1, log_freq=1, early_stopping=False)
    estimates = model.get_map_estimates(data, estimate_method="median")

    assert "sample_depth" in estimates
    assert estimates["sample_depth"].shape == (2,)


def test_diagonal_guide_trains_with_delta_warm_start() -> None:
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
        dtype=torch.float64,
    )
    model = CNVModel(
        guide_type="diagonal",
        sex_cn_weight=0.0,
        autosome_prior_mode="dirichlet",
        device="cpu",
        dtype=torch.float64,
    )

    history = model.train(
        data,
        max_iter=1,
        guide_warmup_iter=1,
        log_freq=1,
        early_stopping=False,
    )

    assert len(history) == 1
    assert np.isfinite(history).all()
    assert model.guide.__class__.__name__ == "AutoDiagonalNormal"


def test_shrinkage_prior_trains_with_autosome_and_sex_bins() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chr21", "chrX"],
            "Start": [0, 0],
            "End": [1000, 1000],
            "S1": [20, 10],
        },
        index=["chr21:0-1000", "chrX:0-1000"],
    )
    data = DepthData(df, device="cpu", depth_space="raw", clamp_threshold=None)
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
    tiny_raw_depth_df: pd.DataFrame,
    tiny_site_data: dict[str, np.ndarray],
) -> None:
    data = DepthData(
        tiny_raw_depth_df,
        device="cpu",
        depth_space="raw",
        clamp_threshold=None,
        site_data=tiny_site_data,
    )
    model = CNVModel(
        sex_cn_weight=0.0,
        epsilon_mean=0.0,
        guide_type="diagonal",
        learn_site_pop_af=True,
        af_weight=0.5,
    )

    model_kw = model._model_kwargs(data)

    assert "af_table" not in model_kw
    assert "af_observed_genotype_log_lik" in model_kw
    assert "af_observed_bin_idx" in model_kw
    assert "af_observed_sample_idx" in model_kw
    assert "af_observed_site_slot_idx" in model_kw
    assert "af_site_bin_idx" in model_kw
    assert model_kw["af_site_bin_idx"].numel() <= data.site_pop_af.numel()
    assert "site_alt" not in model_kw
    assert "site_total" not in model_kw
    assert model_kw["site_pop_af"].shape == data.site_pop_af.shape


def test_learn_site_pop_af_trains_with_mixed_guide(
    tiny_raw_depth_df: pd.DataFrame,
    tiny_site_data: dict[str, np.ndarray],
) -> None:
    data = DepthData(
        tiny_raw_depth_df,
        device="cpu",
        depth_space="raw",
        clamp_threshold=None,
        site_data=tiny_site_data,
    )
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

    assert "site_pop_af_latent_observed" in estimates
    assert "site_pop_af_latent" in estimates
    assert estimates["site_pop_af_latent"].shape == tiny_site_data["site_pop_af"].shape
    assert estimates["site_pop_af_latent_observed"].ndim == 1


def test_get_map_estimates_median_uses_guide_median(monkeypatch) -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chr21", "chr21"],
            "Start": [0, 1000],
            "End": [1000, 2000],
            "S1": [20, 20],
        }
    )
    df["Bin"] = ["chr21:0-1000", "chr21:1000-2000"]
    data = DepthData(
        df.set_index("Bin"),
        device="cpu",
        depth_space="raw",
        clamp_threshold=None,
    )
    model = CNVModel(
        sex_cn_weight=0.0,
        epsilon_mean=0.0,
        background_factors=0,
        multiplicative_factors=1,
        var_bin=0.0,
        guide_type="diagonal",
        autosome_prior_mode="dirichlet",
        learn_af_temperature=False,
    )

    class DummyGuide:
        def median(self, *args, **kwargs):
            return {
                "multiplicative_bin_factors": torch.tensor(
                    [[-1.0], [1.0]],
                    dtype=torch.float32,
                ),
                "multiplicative_sample_factors": torch.tensor(
                    [[0.2]],
                    dtype=torch.float32,
                ),
                "sample_var": torch.tensor([0.05], dtype=torch.float32),
                "cn_probs": torch.full((2, 6), 1.0 / 6.0, dtype=torch.float32),
                "sample_depth": torch.tensor([20.0], dtype=torch.float32),
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

    expected_bias = np.exp(np.array([-1.0, 1.0], dtype=np.float32) * 0.2)
    np.testing.assert_allclose(
        estimates["bin_bias"],
        expected_bias,
        rtol=1e-6,
        atol=1e-6,
    )
    np.testing.assert_allclose(
        estimates["bin_bias_matrix"].reshape(2, 1),
        expected_bias[:, np.newaxis],
        rtol=1e-6,
        atol=1e-6,
    )
    assert int(estimates["cn"][0, 0]) == 2


def test_run_discrete_inference_multi_draw_averages_posteriors(monkeypatch) -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chr21"],
            "Start": [0],
            "End": [1000],
            "S1": [20],
        }
    )
    df["Bin"] = "chr21:0-1000"
    data = DepthData(
        df.set_index("Bin"),
        device="cpu",
        depth_space="raw",
        clamp_threshold=None,
    )
    model = CNVModel(sex_cn_weight=0.0, epsilon_mean=0.0, guide_type="diagonal")

    draw_biases = iter([1.0, 2.0, 3.0])

    def fake_get_continuous_estimates(data, estimate_method="current", model_kw=None):
        bias = next(draw_biases)
        return {
            "bin_bias": torch.tensor([bias], dtype=torch.float32),
            "sample_var": torch.tensor([0.05], dtype=torch.float32),
            "bin_var": torch.tensor([0.02], dtype=torch.float32),
            "cn_probs": torch.full((1, 6), 1.0 / 6.0, dtype=torch.float32),
            "sample_depth": torch.tensor([20.0], dtype=torch.float32),
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
