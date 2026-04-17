from __future__ import annotations

import numpy as np
import pandas as pd

from gatk_sv_ploidy.data import DepthData
from gatk_sv_ploidy.ppd import (
    compute_ppd_global_summary,
    compute_ppd_bin_summary,
    compute_ppd_chromosome_summary,
    generate_ppd_depth,
    main,
    parse_args,
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


def test_compute_ppd_global_summary_has_expected_columns() -> None:
    ppd_bin_df = pd.DataFrame(
        {
            "tail_prob": [0.2, 0.8, 0.6],
            "z_score": [0.1, -0.2, 0.3],
            "observed_depth": [1.0, 2.0, 3.0],
            "ppd_q05": [0.5, 1.0, 2.5],
            "ppd_q25": [0.8, 1.5, 2.8],
            "ppd_q75": [1.2, 2.5, 3.2],
            "ppd_q95": [1.5, 3.0, 3.5],
            "residual": [0.1, -0.2, 0.3],
        }
    )

    summary = compute_ppd_global_summary(ppd_bin_df)

    assert summary.loc[0, "n_bins_x_samples"] == 3
    assert 0.0 <= summary.loc[0, "tail_prob_ks_pval"] <= 1.0
    assert "frac_outside_90pct_interval" in summary.columns


def test_ppd_parse_args_and_main_write_outputs(
    tiny_depth_df: pd.DataFrame,
    tiny_site_data: dict[str, np.ndarray],
    tmp_path,
    monkeypatch,
) -> None:
    depth_path = tmp_path / "depth.tsv"
    tiny_depth_df.to_csv(depth_path, sep="\t")
    site_path = tmp_path / "site_data.npz"
    np.savez_compressed(site_path, **tiny_site_data)

    n_bins = len(tiny_depth_df)
    n_samples = 2
    posterior = np.zeros((n_bins, n_samples, 6), dtype=np.float32)
    posterior[:, :, 2] = 1.0
    artifact_path = tmp_path / "artifacts.npz"
    np.savez_compressed(
        artifact_path,
        bin_bias=np.full(n_bins, 1.0, dtype=np.float32),
        sample_var=np.full(n_samples, 0.1, dtype=np.float32),
        bin_var=np.full(n_bins, 0.05, dtype=np.float32),
        cn_post_cn_posterior=posterior,
    )
    output_dir = tmp_path / "out"

    monkeypatch.setattr(
        "sys.argv",
        [
            "gatk-sv-ploidy ppd",
            "--input",
            str(depth_path),
            "--artifacts",
            str(artifact_path),
            "--output-dir",
            str(output_dir),
            "--site-data",
            str(site_path),
            "--draws",
            "4",
            "--seed",
            "9",
        ],
    )

    args = parse_args()
    assert args.draws == 4
    assert args.seed == 9

    main()

    assert (output_dir / "ppd_draws.npz").exists()
    assert (output_dir / "ppd_bin_summary.tsv.gz").exists()
    assert (output_dir / "ppd_chromosome_summary.tsv").exists()
    assert (output_dir / "ppd_global_summary.tsv").exists()

    global_df = pd.read_csv(output_dir / "ppd_global_summary.tsv", sep="\t")
    assert global_df.loc[0, "n_bins_x_samples"] == n_bins * n_samples
