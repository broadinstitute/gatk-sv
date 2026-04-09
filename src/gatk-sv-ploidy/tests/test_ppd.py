"""Tests for the posterior predictive check (ppd) subcommand."""

from __future__ import annotations

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Fixtures (reuse the pattern from test_allele_fraction.py)
# ---------------------------------------------------------------------------

_SAMPLE_COLS = ["SAMPLE_A", "SAMPLE_B"]
_CHROMS = ["chr21", "chr21", "chr21", "chrX", "chrX"]
_STARTS = [0, 1000, 2000, 0, 1000]
_ENDS = [1000, 2000, 3000, 1000, 2000]


def _make_depth_df() -> pd.DataFrame:
    """Return a small depth-like DataFrame for testing."""
    data = {
        "Chr": _CHROMS,
        "Start": _STARTS,
        "End": _ENDS,
        "SAMPLE_A": [2.0, 2.1, 1.9, 1.0, 0.9],
        "SAMPLE_B": [2.0, 1.8, 2.2, 2.0, 2.1],
    }
    df = pd.DataFrame(data)
    df["Bin"] = (
        df["Chr"].astype(str) + ":" + df["Start"].astype(str)
        + "-" + df["End"].astype(str)
    )
    df = df.set_index("Bin")
    return df


def _make_site_data(n_bins: int = 5, max_sites: int = 3,
                    n_samples: int = 2) -> dict:
    """Return synthetic per-site allele data arrays."""
    rng = np.random.RandomState(42)
    site_total = rng.randint(10, 40, size=(n_bins, max_sites, n_samples)
                             ).astype(np.int32)
    mask = rng.random((n_bins, max_sites, n_samples)) > 0.3
    site_total[~mask] = 0
    site_alt = (site_total * rng.uniform(0.2, 0.8, size=site_total.shape)
                ).astype(np.int32)
    site_alt[~mask] = 0
    site_pop_af = rng.uniform(0.1, 0.5, size=(n_bins, max_sites)
                              ).astype(np.float32)
    site_mask = site_total > 0
    return {
        "site_alt": site_alt,
        "site_total": site_total,
        "site_pop_af": site_pop_af,
        "site_mask": site_mask,
    }


def _train_model_and_get_artifacts(data):
    """Train a small model and return (map_estimates, cn_posterior, model)."""
    import pyro
    from gatk_sv_ploidy.models import CNVModel

    pyro.clear_param_store()
    model = CNVModel(n_states=6, guide_type="delta", af_weight=0.0,
                     sex_cn_weight=0.0)
    model.train(data, max_iter=10, log_freq=100)
    map_est = model.get_map_estimates(data)
    cn_post = model.run_discrete_inference(data, map_estimates=map_est)
    return map_est, cn_post, model


# ---------------------------------------------------------------------------
# Inference artifact I/O
# ---------------------------------------------------------------------------


class TestLoadInferenceArtifacts:
    """Tests for saving and loading inference artifacts."""

    def test_roundtrip(self, tmp_path):
        from gatk_sv_ploidy.data import DepthData
        from gatk_sv_ploidy.infer import load_inference_artifacts

        df = _make_depth_df()
        data = DepthData(df)
        map_est, cn_post, _ = _train_model_and_get_artifacts(data)

        # Save like infer.py does
        artifact_dict = {k: v for k, v in map_est.items()}
        for k, v in cn_post.items():
            artifact_dict[f"cn_post_{k}"] = v
        path = tmp_path / "inference_artifacts.npz"
        np.savez_compressed(str(path), **artifact_dict)

        # Load back
        loaded_map, loaded_cn = load_inference_artifacts(str(path))

        for key in map_est:
            np.testing.assert_array_almost_equal(
                loaded_map[key], map_est[key], decimal=5,
            )
        for key in cn_post:
            np.testing.assert_array_almost_equal(
                loaded_cn[key], cn_post[key], decimal=5,
            )


# ---------------------------------------------------------------------------
# PPD draw generation
# ---------------------------------------------------------------------------


class TestGeneratePPDDepth:
    """Tests for :func:`ppd.generate_ppd_depth`."""

    def test_shape(self):
        from gatk_sv_ploidy.data import DepthData
        from gatk_sv_ploidy.ppd import generate_ppd_depth

        df = _make_depth_df()
        data = DepthData(df)
        map_est, cn_post, _ = _train_model_and_get_artifacts(data)

        draws = generate_ppd_depth(data, map_est, cn_post, n_draws=50)
        assert draws.shape == (50, data.n_bins, data.n_samples)

    def test_finite(self):
        from gatk_sv_ploidy.data import DepthData
        from gatk_sv_ploidy.ppd import generate_ppd_depth

        df = _make_depth_df()
        data = DepthData(df)
        map_est, cn_post, _ = _train_model_and_get_artifacts(data)

        draws = generate_ppd_depth(data, map_est, cn_post, n_draws=20)
        assert np.isfinite(draws).all()

    def test_reproducible(self):
        from gatk_sv_ploidy.data import DepthData
        from gatk_sv_ploidy.ppd import generate_ppd_depth

        df = _make_depth_df()
        data = DepthData(df)
        map_est, cn_post, _ = _train_model_and_get_artifacts(data)

        draws1 = generate_ppd_depth(data, map_est, cn_post, n_draws=10, seed=99)
        draws2 = generate_ppd_depth(data, map_est, cn_post, n_draws=10, seed=99)
        np.testing.assert_array_equal(draws1, draws2)


# ---------------------------------------------------------------------------
# PPD summary statistics
# ---------------------------------------------------------------------------


class TestPPDBinSummary:
    """Tests for :func:`ppd.compute_ppd_bin_summary`."""

    def test_shape_and_columns(self):
        from gatk_sv_ploidy.data import DepthData
        from gatk_sv_ploidy.ppd import compute_ppd_bin_summary, generate_ppd_depth

        df = _make_depth_df()
        data = DepthData(df)
        map_est, cn_post, _ = _train_model_and_get_artifacts(data)
        draws = generate_ppd_depth(data, map_est, cn_post, n_draws=50)

        summary = compute_ppd_bin_summary(data, draws, map_est, cn_post)

        assert len(summary) == data.n_bins * data.n_samples
        expected_cols = {
            "chr", "start", "end", "sample", "observed_depth",
            "cn_map", "ppd_mean", "ppd_std", "residual", "z_score",
            "tail_prob", "two_tail_prob",
            "ppd_q05", "ppd_q25", "ppd_q50", "ppd_q75", "ppd_q95",
        }
        assert expected_cols.issubset(set(summary.columns))

    def test_tail_prob_range(self):
        from gatk_sv_ploidy.data import DepthData
        from gatk_sv_ploidy.ppd import compute_ppd_bin_summary, generate_ppd_depth

        df = _make_depth_df()
        data = DepthData(df)
        map_est, cn_post, _ = _train_model_and_get_artifacts(data)
        draws = generate_ppd_depth(data, map_est, cn_post, n_draws=100)

        summary = compute_ppd_bin_summary(data, draws, map_est, cn_post)

        assert (summary["tail_prob"] >= 0).all()
        assert (summary["tail_prob"] <= 1).all()
        assert (summary["two_tail_prob"] >= 0).all()
        assert (summary["two_tail_prob"] <= 1).all()


class TestPPDChromosomeSummary:
    """Tests for :func:`ppd.compute_ppd_chromosome_summary`."""

    def test_shape_and_columns(self):
        from gatk_sv_ploidy.data import DepthData
        from gatk_sv_ploidy.ppd import (
            compute_ppd_chromosome_summary,
            generate_ppd_depth,
        )

        df = _make_depth_df()
        data = DepthData(df)
        map_est, cn_post, _ = _train_model_and_get_artifacts(data)
        draws = generate_ppd_depth(data, map_est, cn_post, n_draws=50)

        summary = compute_ppd_chromosome_summary(data, draws, cn_post)

        n_chrs = len(np.unique(data.chr))
        assert len(summary) == data.n_samples * n_chrs
        assert "bayesian_pvalue" in summary.columns
        assert (summary["bayesian_pvalue"] >= 0).all()
        assert (summary["bayesian_pvalue"] <= 1).all()


class TestPPDGlobalSummary:
    """Tests for :func:`ppd.compute_ppd_global_summary`."""

    def test_single_row(self):
        from gatk_sv_ploidy.data import DepthData
        from gatk_sv_ploidy.ppd import (
            compute_ppd_bin_summary,
            compute_ppd_global_summary,
            generate_ppd_depth,
        )

        df = _make_depth_df()
        data = DepthData(df)
        map_est, cn_post, _ = _train_model_and_get_artifacts(data)
        draws = generate_ppd_depth(data, map_est, cn_post, n_draws=50)
        bin_summary = compute_ppd_bin_summary(data, draws, map_est, cn_post)

        global_summary = compute_ppd_global_summary(bin_summary)

        assert len(global_summary) == 1
        assert "rmse" in global_summary.columns
        assert "z_score_mean" in global_summary.columns
        assert "tail_prob_ks_pval" in global_summary.columns
        assert "frac_outside_90pct_interval" in global_summary.columns


# ---------------------------------------------------------------------------
# PPD plotting (smoke tests — just ensure no exceptions)
# ---------------------------------------------------------------------------


class TestPPDPlots:
    """Smoke tests for PPD plotting functions."""

    def _get_ppd_data(self):
        from gatk_sv_ploidy.data import DepthData
        from gatk_sv_ploidy.ppd import (
            compute_ppd_bin_summary,
            compute_ppd_chromosome_summary,
            generate_ppd_depth,
        )

        df = _make_depth_df()
        data = DepthData(df)
        map_est, cn_post, _ = _train_model_and_get_artifacts(data)
        draws = generate_ppd_depth(data, map_est, cn_post, n_draws=30)
        ppd_bin = compute_ppd_bin_summary(data, draws, map_est, cn_post)
        ppd_chr = compute_ppd_chromosome_summary(data, draws, cn_post)
        return ppd_bin, ppd_chr

    def test_run_ppd_plots(self, tmp_path):
        from gatk_sv_ploidy._plot_ppd import run_ppd_plots

        ppd_bin, ppd_chr = self._get_ppd_data()
        run_ppd_plots(ppd_bin, ppd_chr, str(tmp_path))

        # Check that at least some plots were created
        ppd_dir = tmp_path / "ppd"
        assert ppd_dir.exists()
        png_files = list(ppd_dir.glob("*.png"))
        assert len(png_files) >= 5

    def test_run_ppd_plots_with_highlight(self, tmp_path):
        from gatk_sv_ploidy._plot_ppd import run_ppd_plots

        ppd_bin, ppd_chr = self._get_ppd_data()
        run_ppd_plots(ppd_bin, ppd_chr, str(tmp_path),
                      highlight_sample="SAMPLE_A")

        ppd_dir = tmp_path / "ppd"
        assert ppd_dir.exists()


# ---------------------------------------------------------------------------
# New diagnostic plot smoke tests
# ---------------------------------------------------------------------------


class TestNewDiagnosticPlots:
    """Smoke tests for the new non-PPD diagnostic plots."""

    def _get_test_data(self):
        from gatk_sv_ploidy.data import DepthData
        from gatk_sv_ploidy.infer import (
            build_bin_stats,
            build_chromosome_stats,
            detect_aneuploidies,
        )

        df = _make_depth_df()
        data = DepthData(df)
        map_est, cn_post, _ = _train_model_and_get_artifacts(data)
        aneuploid_map = detect_aneuploidies(data, map_est, cn_post)
        bin_df = build_bin_stats(data, map_est, cn_post)
        chr_df = build_chromosome_stats(data, map_est, cn_post, aneuploid_map)
        loss_df = pd.DataFrame({"epoch": list(range(10)),
                                "elbo": list(range(100, 90, -1))})
        return bin_df, chr_df, loss_df

    def test_cn_posterior_entropy(self, tmp_path):
        from gatk_sv_ploidy.plot import plot_cn_posterior_entropy

        bin_df, _, _ = self._get_test_data()
        plot_cn_posterior_entropy(bin_df, str(tmp_path))
        assert (tmp_path / "diagnostics" / "cn_posterior_entropy.png").exists()

    def test_chromosome_cn_heatmap(self, tmp_path):
        from gatk_sv_ploidy.plot import plot_chromosome_cn_heatmap

        _, chr_df, _ = self._get_test_data()
        plot_chromosome_cn_heatmap(chr_df, str(tmp_path))
        assert (tmp_path / "diagnostics" / "chromosome_cn_heatmap.png").exists()

    def test_training_loss_gradient(self, tmp_path):
        from gatk_sv_ploidy.plot import plot_training_loss_with_gradient

        _, _, loss_df = self._get_test_data()
        plot_training_loss_with_gradient(loss_df, str(tmp_path))
        assert (tmp_path / "diagnostics" / "training_loss_gradient.png").exists()

    def test_parameter_diagnostics(self, tmp_path):
        from gatk_sv_ploidy.plot import plot_parameter_diagnostics

        bin_df, _, _ = self._get_test_data()
        plot_parameter_diagnostics(bin_df, str(tmp_path))
        assert (tmp_path / "diagnostics" / "parameter_diagnostics.png").exists()
