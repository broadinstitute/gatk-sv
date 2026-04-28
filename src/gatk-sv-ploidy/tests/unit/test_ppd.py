from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from gatk_sv_ploidy.data import DepthData
from gatk_sv_ploidy.ppd import (
    _build_model_from_artifacts,
    _compute_randomized_pit,
    _resolve_input_depth_space,
    apply_effective_site_pop_af,
    compute_call_stability_quality_summary,
    compute_ppd_bin_quality_summary,
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


def test_compute_randomized_pit_randomizes_discrete_ties_deterministically() -> None:
    draws = np.array([[[0.0]], [[1.0]], [[1.0]], [[2.0]]], dtype=np.float32)
    observed = np.array([[1.0]], dtype=np.float32)

    pit_one = _compute_randomized_pit(draws, observed, seed=0)
    pit_two = _compute_randomized_pit(draws, observed, seed=0)

    np.testing.assert_allclose(pit_one, pit_two)
    assert float(pit_one[0, 0]) >= 0.25
    assert float(pit_one[0, 0]) <= 0.75


def test_generate_ppd_depth_supports_studentt_metadata(tiny_depth_df: pd.DataFrame) -> None:
    data = DepthData(tiny_depth_df, clamp_threshold=None)
    cn_post = _fake_cn_posterior(data.n_bins, data.n_samples, cn_state=2)
    map_est_studentt = {
        "bin_bias": np.full(data.n_bins, 1.0, dtype=np.float32),
        "sample_var": np.full(data.n_samples, 0.1, dtype=np.float32),
        "bin_var": np.full(data.n_bins, 0.05, dtype=np.float32),
        "obs_likelihood": np.asarray("studentt"),
        "obs_df": np.asarray(3.5, dtype=np.float32),
    }
    map_est_normal = {
        "bin_bias": np.full(data.n_bins, 1.0, dtype=np.float32),
        "sample_var": np.full(data.n_samples, 0.1, dtype=np.float32),
        "bin_var": np.full(data.n_bins, 0.05, dtype=np.float32),
        "obs_likelihood": np.asarray("normal"),
        "obs_df": np.asarray(3.5, dtype=np.float32),
    }

    draws_studentt = generate_ppd_depth(
        data,
        map_est_studentt,
        cn_post,
        n_draws=12,
        seed=11,
    )
    draws_normal = generate_ppd_depth(
        data,
        map_est_normal,
        cn_post,
        n_draws=12,
        seed=11,
    )

    assert draws_studentt.shape == (12, data.n_bins, data.n_samples)
    assert np.isfinite(draws_studentt).all()
    assert not np.allclose(draws_studentt, draws_normal)


def test_apply_effective_site_pop_af_prefers_saved_infer_values(
    tiny_site_data: dict[str, np.ndarray],
) -> None:
    updated = apply_effective_site_pop_af(
        tiny_site_data,
        {
            "site_pop_af_effective": np.full_like(
                tiny_site_data["site_pop_af"],
                0.125,
                dtype=np.float32,
            ),
            "site_af_estimation_applied": np.asarray(True),
        },
    )

    assert updated is not None
    np.testing.assert_allclose(updated["site_pop_af"], 0.125)


def test_generate_ppd_depth_supports_negative_binomial_metadata() -> None:
    raw_df = pd.DataFrame(
        {
            "Chr": ["chr21", "chrX"],
            "Start": [0, 1000],
            "End": [1000, 2000],
            "S1": [10, 5],
            "S2": [20, 10],
        },
        index=["chr21:0-1000", "chrX:1000-2000"],
    )
    data = DepthData(raw_df, depth_space="raw", clamp_threshold=None)
    cn_post = _fake_cn_posterior(data.n_bins, data.n_samples, cn_state=2)
    map_est = {
        "bin_bias": np.full(data.n_bins, 1.0, dtype=np.float32),
        "sample_var": np.full(data.n_samples, 1e-3, dtype=np.float32),
        "bin_var": np.full(data.n_bins, 1e-3, dtype=np.float32),
        "sample_depth": np.array([10.0, 20.0], dtype=np.float32),
        "obs_likelihood": np.asarray("negative_binomial"),
        "depth_space": np.asarray("raw"),
    }

    draws = generate_ppd_depth(data, map_est, cn_post, n_draws=40, seed=17)

    assert draws.shape == (40, data.n_bins, data.n_samples)
    assert np.allclose(draws, np.rint(draws))
    assert float(draws[:, 0, 0].mean()) == pytest.approx(10.0, rel=0.3)
    assert float(draws[:, 0, 1].mean()) == pytest.approx(20.0, rel=0.3)
    assert float(draws[:, 0, 1].mean()) > float(draws[:, 0, 0].mean())


def test_generate_ppd_depth_uses_saved_posterior_draws() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chr21"],
            "Start": [0],
            "End": [1000],
            "S1": [3.0],
        },
        index=["chr21:0-1000"],
    )
    data = DepthData(df, clamp_threshold=None)
    cn_post = _fake_cn_posterior(data.n_bins, data.n_samples, cn_state=2)
    base_map = {
        "bin_bias": np.array([1.0], dtype=np.float32),
        "sample_var": np.array([1e-6], dtype=np.float32),
        "bin_var": np.array([1e-6], dtype=np.float32),
        "obs_likelihood": np.asarray("normal"),
        "obs_df": np.asarray(3.5, dtype=np.float32),
        "model_n_states": np.asarray(6),
        "model_alpha_ref": np.asarray(50.0),
        "model_alpha_non_ref": np.asarray(1.0),
        "model_var_bias_bin": np.asarray(0.05),
        "model_var_sample": np.asarray(0.001),
        "model_var_bin": np.asarray(0.001),
        "model_af_concentration": np.asarray(50.0),
        "model_af_weight": np.asarray(0.0),
        "model_alpha_sex_ref": np.asarray(1.0),
        "model_alpha_sex_non_ref": np.asarray(1.0),
        "model_sex_prior": np.asarray([0.5, 0.5], dtype=np.float32),
        "model_sex_cn_weight": np.asarray(0.0),
        "model_guide_type": np.asarray("diagonal"),
        "posterior_draws_bin_bias": np.array([[1.0], [2.0]], dtype=np.float32),
        "posterior_draws_sample_var": np.array([[1e-6], [1e-6]], dtype=np.float32),
        "posterior_draws_bin_var": np.array([[1e-6], [1e-6]], dtype=np.float32),
    }

    plugin_map = {
        key: value
        for key, value in base_map.items()
        if not key.startswith("posterior_draws_")
    }
    plugin_draws = generate_ppd_depth(data, plugin_map, cn_post, n_draws=200, seed=19)
    posterior_draws = generate_ppd_depth(
        data,
        {
            **base_map,
            "posterior_draws_bin_bias": base_map["posterior_draws_bin_bias"],
            "posterior_draws_sample_var": base_map["posterior_draws_sample_var"],
            "posterior_draws_bin_var": base_map["posterior_draws_bin_var"],
        },
        cn_post,
        n_draws=200,
        seed=19,
    )

    assert float(plugin_draws.mean()) == pytest.approx(2.0, abs=0.2)
    assert float(posterior_draws.mean()) > float(plugin_draws.mean()) + 1.0
    assert float(posterior_draws.mean()) < 4.2


def test_build_model_from_artifacts_restores_af_evidence_mode() -> None:
    model = _build_model_from_artifacts(
        {
            "model_af_weight": np.asarray(0.25),
            "model_af_evidence_mode": np.asarray("relative"),
        },
        device="cpu",
        obs_likelihood="normal",
        obs_df=3.5,
    )

    assert model.af_evidence_mode == "relative"


def test_build_model_from_artifacts_defaults_old_epsilon_prior_to_exponential() -> None:
    model = _build_model_from_artifacts(
        {
            "epsilon_mean": np.asarray(0.1),
        },
        device="cpu",
        obs_likelihood="normal",
        obs_df=3.5,
    )

    assert model.epsilon_mean == pytest.approx(0.1)
    assert model.epsilon_concentration == pytest.approx(1.0)


def test_build_model_from_artifacts_restores_sparse_gamma_epsilon_prior() -> None:
    model = _build_model_from_artifacts(
        {
            "epsilon_mean": np.asarray(0.1),
            "model_epsilon_concentration": np.asarray(0.5),
        },
        device="cpu",
        obs_likelihood="normal",
        obs_df=3.5,
    )

    assert model.epsilon_mean == pytest.approx(0.1)
    assert model.epsilon_concentration == pytest.approx(0.5)


def test_build_model_from_artifacts_defaults_old_background_factors_to_disabled() -> None:
    model = _build_model_from_artifacts(
        {
            "epsilon_mean": np.asarray(0.1),
        },
        device="cpu",
        obs_likelihood="normal",
        obs_df=3.5,
    )

    assert model.background_factors == 0


def test_generate_ppd_depth_uses_bin_epsilon(tiny_depth_df: pd.DataFrame) -> None:
    data = DepthData(tiny_depth_df.iloc[[3]], clamp_threshold=None)
    cn_post = _fake_cn_posterior(data.n_bins, data.n_samples, cn_state=0)
    base_map = {
        "bin_bias": np.full(data.n_bins, 1.0, dtype=np.float32),
        "sample_var": np.full(data.n_samples, 1e-4, dtype=np.float32),
        "bin_var": np.full(data.n_bins, 1e-4, dtype=np.float32),
    }

    draws_without_epsilon = generate_ppd_depth(
        data,
        base_map,
        cn_post,
        n_draws=200,
        seed=13,
    )
    draws_with_epsilon = generate_ppd_depth(
        data,
        {
            **base_map,
            "bin_epsilon": np.array([[0.0, 0.2]], dtype=np.float32),
        },
        cn_post,
        n_draws=200,
        seed=13,
    )

    assert float(draws_with_epsilon[:, 0, 1].mean()) > float(draws_without_epsilon[:, 0, 1].mean()) + 0.15
    assert abs(float(draws_with_epsilon[:, 0, 0].mean()) - float(draws_without_epsilon[:, 0, 0].mean())) < 0.15


def test_generate_ppd_depth_expands_contig_shared_bin_epsilon() -> None:
    df = pd.DataFrame(
        {
            "Chr": ["chrY", "chrY"],
            "Start": [25, 125],
            "End": [125, 225],
            "SAMPLE_B": [0.0, 0.0],
            "SAMPLE_A": [0.0, 0.0],
        }
    )
    df["Bin"] = df["Chr"].astype(str) + ":" + df["Start"].astype(str) + "-" + df["End"].astype(str)
    data = DepthData(df.set_index("Bin"), clamp_threshold=None)
    cn_post = _fake_cn_posterior(data.n_bins, data.n_samples, cn_state=0)
    base_map = {
        "bin_bias": np.full(data.n_bins, 1.0, dtype=np.float32),
        "sample_var": np.full(data.n_samples, 1e-4, dtype=np.float32),
        "bin_var": np.full(data.n_bins, 1e-4, dtype=np.float32),
    }
    sample_with_epsilon = data.sample_ids.index("SAMPLE_A")
    sample_without_epsilon = data.sample_ids.index("SAMPLE_B")
    contig_epsilon = np.zeros((1, data.n_samples), dtype=np.float32)
    contig_epsilon[0, sample_with_epsilon] = 0.2

    draws_without_epsilon = generate_ppd_depth(
        data,
        base_map,
        cn_post,
        n_draws=200,
        seed=17,
    )
    draws_with_epsilon = generate_ppd_depth(
        data,
        {
            **base_map,
            "bin_epsilon": contig_epsilon,
        },
        cn_post,
        n_draws=200,
        seed=17,
    )

    for bin_index in range(data.n_bins):
        assert float(draws_with_epsilon[:, bin_index, sample_with_epsilon].mean()) > float(
            draws_without_epsilon[:, bin_index, sample_with_epsilon].mean()
        ) + 0.15
        assert abs(
            float(draws_with_epsilon[:, bin_index, sample_without_epsilon].mean()) -
            float(draws_without_epsilon[:, bin_index, sample_without_epsilon].mean())
        ) < 0.15


def test_generate_ppd_depth_uses_background_factors(tiny_depth_df: pd.DataFrame) -> None:
    data = DepthData(tiny_depth_df.iloc[[3]], clamp_threshold=None)
    cn_post = _fake_cn_posterior(data.n_bins, data.n_samples, cn_state=0)
    base_map = {
        "bin_bias": np.full(data.n_bins, 1.0, dtype=np.float32),
        "sample_var": np.full(data.n_samples, 1e-4, dtype=np.float32),
        "bin_var": np.full(data.n_bins, 1e-4, dtype=np.float32),
    }

    draws_without_background = generate_ppd_depth(
        data,
        base_map,
        cn_post,
        n_draws=200,
        seed=23,
    )
    draws_with_background = generate_ppd_depth(
        data,
        {
            **base_map,
            "background_bin_factors": np.array([[3.0]], dtype=np.float32),
            "background_sample_factors": np.array([[0.0, 0.2]], dtype=np.float32),
        },
        cn_post,
        n_draws=200,
        seed=23,
    )

    assert float(draws_with_background[:, 0, 1].mean()) > float(
        draws_without_background[:, 0, 1].mean()
    ) + 0.15
    assert abs(
        float(draws_with_background[:, 0, 0].mean()) -
        float(draws_without_background[:, 0, 0].mean())
    ) < 0.15


def test_ppd_resolve_input_depth_space_prefers_preprocess_marker(tmp_path) -> None:
    depth_path = tmp_path / "preprocessed_depth.tsv"
    depth_path.write_text("Bin\tChr\tStart\tEnd\tS1\n", encoding="ascii")
    (tmp_path / "observation_type.txt").write_text("raw\n", encoding="ascii")

    depth_space = _resolve_input_depth_space(
        "auto",
        "negative_binomial",
        str(depth_path),
        {"depth_space": np.asarray("normalized")},
    )

    assert depth_space == "raw"


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
    bin_quality = compute_ppd_bin_quality_summary(bin_summary)
    chrom_summary = compute_ppd_chromosome_summary(data, draws, cn_post)

    assert len(bin_summary) == data.n_bins * data.n_samples
    assert set([
        "randomized_pit",
        "tail_prob",
        "two_tail_prob",
        "ppd_mean",
        "outside_90pct_interval",
        "outside_50pct_interval",
    ]).issubset(bin_summary.columns)
    assert ((bin_summary["randomized_pit"] >= 0.0) & (bin_summary["randomized_pit"] <= 1.0)).all()
    assert ((bin_summary["tail_prob"] >= 0.0) & (bin_summary["tail_prob"] <= 1.0)).all()
    assert len(bin_quality) == data.n_bins
    assert set([
        "frac_outside_90pct_interval",
        "frac_outside_50pct_interval",
        "BINQ15",
        "BINQ20",
    ]).issubset(bin_quality.columns)
    assert len(chrom_summary) == data.n_samples * len(np.unique(data.chr))


def test_compute_ppd_bin_quality_summary_penalizes_miscalibrated_bins() -> None:
    good_rows = [
        {
            "chr": "chr1",
            "start": 0,
            "end": 100,
            "sample": f"S{i}",
            "observed_depth": 1.0,
            "ppd_q05": 0.5,
            "ppd_q25": 0.8,
            "ppd_q75": 1.2,
            "ppd_q95": 1.5,
            "two_tail_prob": 0.9,
            "outside_90pct_interval": False,
            "outside_50pct_interval": False,
        }
        for i in range(10)
    ]
    bad_rows = [
        {
            "chr": "chr2",
            "start": 0,
            "end": 100,
            "sample": f"S{i}",
            "observed_depth": 3.0,
            "ppd_q05": 0.5,
            "ppd_q25": 0.8,
            "ppd_q75": 1.2,
            "ppd_q95": 1.5,
            "two_tail_prob": 0.01,
            "outside_90pct_interval": True,
            "outside_50pct_interval": True,
        }
        for i in range(10)
    ]
    ppd_bin_df = pd.DataFrame(good_rows + bad_rows)

    summary = compute_ppd_bin_quality_summary(ppd_bin_df)

    good = summary[summary["chr"] == "chr1"].iloc[0]
    bad = summary[summary["chr"] == "chr2"].iloc[0]

    assert float(good["BINQ20"]) > float(bad["BINQ20"])
    assert float(good["frac_outside_90pct_interval"]) < float(bad["frac_outside_90pct_interval"])


def test_compute_call_stability_quality_summary_penalizes_unstable_bins(
    tiny_depth_df: pd.DataFrame,
) -> None:
    data = DepthData(tiny_depth_df.iloc[:2], clamp_threshold=None)
    cn_post = {
        "cn_posterior": _fake_cn_posterior(data.n_bins, data.n_samples, cn_state=2)["cn_posterior"],
        "cn_map_stability": np.array(
            [
                [0.99, 0.98],
                [0.55, 0.60],
            ],
            dtype=np.float32,
        ),
    }

    summary = compute_call_stability_quality_summary(data, cn_post)

    good = summary.iloc[0]
    bad = summary.iloc[1]
    assert float(good["CALLQ20"]) > float(bad["CALLQ20"])
    assert float(good["mean_call_instability"]) < float(bad["mean_call_instability"])


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
        cn_post_cn_map_stability=np.full((n_bins, n_samples), 1.0, dtype=np.float32),
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
    assert args.depth_space == "auto"

    main()

    assert (output_dir / "ppd_draws.npz").exists()
    assert (output_dir / "ppd_bin_summary.tsv.gz").exists()
    assert (output_dir / "ppd_bin_quality.tsv").exists()
    assert (output_dir / "ppd_chromosome_summary.tsv").exists()
    assert (output_dir / "ppd_global_summary.tsv").exists()

    global_df = pd.read_csv(output_dir / "ppd_global_summary.tsv", sep="\t")
    assert global_df.loc[0, "n_bins_x_samples"] == n_bins * n_samples
    quality_df = pd.read_csv(output_dir / "ppd_bin_quality.tsv", sep="\t")
    assert "CALLQ20" in quality_df.columns
