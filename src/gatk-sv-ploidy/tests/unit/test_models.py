from __future__ import annotations

import numpy as np
import torch

from gatk_sv_ploidy.models import (
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
