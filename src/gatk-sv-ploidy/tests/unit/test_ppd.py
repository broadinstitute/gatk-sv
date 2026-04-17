from __future__ import annotations

import numpy as np
import pandas as pd

from gatk_sv_ploidy.data import DepthData
from gatk_sv_ploidy.ppd import (
    compute_ppd_bin_summary,
    compute_ppd_chromosome_summary,
    generate_ppd_depth,
)


def _fake_cn_posterior(n_bins: int, n_samples: int, cn_state: int) -> dict[str, np.ndarray]:
    posterior = np.zeros((n_bins, n_samples, 6), dtype=np.float32)
    posterior[:, :, cn_state] = 1.0
    return {"cn_posterior": posterior}


def test_generate_ppd_depth_is_seeded(tiny_depth_df: pd.DataFrame) -> None:
    data = DepthData(tiny_depth_df, clamp_threshold=None)
    map_est = {
        "bin_bias": np.full(data.n_bins, 1.0, dtype=np.float32),
        "sample_var": np.full(data.n_samples, 0.1, dtype=np.float32),
        "bin_var": np.full(data.n_bins, 0.05, dtype=np.float32),
    }
    cn_post = _fake_cn_posterior(data.n_bins, data.n_samples, cn_state=2)

    draws_one = generate_ppd_depth(data, map_est, cn_post, n_draws=5, seed=7)
    draws_two = generate_ppd_depth(data, map_est, cn_post, n_draws=5, seed=7)

    np.testing.assert_allclose(draws_one, draws_two)


def test_ppd_summaries_have_expected_shapes(tiny_depth_df: pd.DataFrame) -> None:
    data = DepthData(tiny_depth_df, clamp_threshold=None)
    map_est = {
        "bin_bias": np.full(data.n_bins, 1.0, dtype=np.float32),
        "sample_var": np.full(data.n_samples, 0.1, dtype=np.float32),
        "bin_var": np.full(data.n_bins, 0.05, dtype=np.float32),
    }
    cn_post = _fake_cn_posterior(data.n_bins, data.n_samples, cn_state=2)
    draws = generate_ppd_depth(data, map_est, cn_post, n_draws=8, seed=3)

    bin_summary = compute_ppd_bin_summary(data, draws, map_est, cn_post)
    chrom_summary = compute_ppd_chromosome_summary(data, draws, cn_post)

    assert len(bin_summary) == data.n_bins * data.n_samples
    assert set(["tail_prob", "two_tail_prob", "ppd_mean"]).issubset(bin_summary.columns)
    assert ((bin_summary["tail_prob"] >= 0.0) & (bin_summary["tail_prob"] <= 1.0)).all()
    assert len(chrom_summary) == data.n_samples * len(np.unique(data.chr))
