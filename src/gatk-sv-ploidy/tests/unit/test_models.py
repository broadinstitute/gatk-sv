from __future__ import annotations

import numpy as np
import pandas as pd
import torch

from gatk_sv_ploidy.data import DepthData
from gatk_sv_ploidy.models import (
    CNVModel,
    _depth_log_lik_numpy,
    _matched_residual_scale,
    _marginalized_af_log_lik,
    _marginalized_af_log_lik_numpy,
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


def test_cnv_model_accepts_lowrank_guide() -> None:
    model = CNVModel(guide_type="lowrank")
    assert model.guide.__class__.__name__ == "AutoLowRankMultivariateNormal"


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
