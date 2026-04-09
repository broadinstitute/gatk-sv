"""Tests for the per-site joint genotype/CN allele fraction model."""

from __future__ import annotations

import numpy as np
import pandas as pd
import torch

# ---------------------------------------------------------------------------
# Fixtures
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
        df["Chr"].astype(str) + ":" + df["Start"].astype(str) + "-" + df["End"].astype(str)
    )
    df = df.set_index("Bin")
    return df


def _make_site_data(n_bins: int = 5, max_sites: int = 3, n_samples: int = 2) -> dict:
    """Return synthetic per-site allele data arrays."""
    rng = np.random.RandomState(42)
    site_total = rng.randint(10, 40, size=(n_bins, max_sites, n_samples)).astype(np.int32)
    # Make about half the sites "active" (non-zero)
    mask = rng.random((n_bins, max_sites, n_samples)) > 0.3
    site_total[~mask] = 0
    # Alt counts = roughly half of total for het sites
    site_alt = (site_total * rng.uniform(0.2, 0.8, size=site_total.shape)).astype(np.int32)
    site_alt[~mask] = 0
    site_pop_af = rng.uniform(0.1, 0.5, size=(n_bins, max_sites)).astype(np.float32)
    site_mask = site_total > 0
    return {
        "site_alt": site_alt,
        "site_total": site_total,
        "site_pop_af": site_pop_af,
        "site_mask": site_mask,
    }


# ---------------------------------------------------------------------------
# SD file reading
# ---------------------------------------------------------------------------


class TestReadSiteDepthTsv:
    """Tests for :func:`preprocess.read_site_depth_tsv`."""

    def test_basic_read(self, tmp_path):
        from gatk_sv_ploidy.preprocess import read_site_depth_tsv

        sd_path = tmp_path / "test.sd.txt"
        sd_path.write_text(
            "chr21\t500\tSAMPLE_A\t10\t0\t0\t20\n"
            "chr21\t1500\tSAMPLE_A\t15\t0\t0\t15\n"
            "chrX\t200\tSAMPLE_A\t5\t0\t25\t0\n"
        )
        df = read_site_depth_tsv(str(sd_path))
        assert len(df) == 3
        assert list(df.columns) == ["contig", "position", "sample", "A", "C", "G", "T"]
        assert df["A"].iloc[0] == 10
        assert df["T"].iloc[0] == 20

    def test_gzipped_read(self, tmp_path):
        import gzip

        from gatk_sv_ploidy.preprocess import read_site_depth_tsv

        content = b"chr21\t500\tSAMPLE_A\t10\t0\t0\t20\n"
        sd_path = tmp_path / "test.sd.txt.gz"
        with gzip.open(str(sd_path), "wb") as fh:
            fh.write(content)
        df = read_site_depth_tsv(str(sd_path))
        assert len(df) == 1

    def test_position_stride(self, tmp_path):
        """position_stride keeps only rows where position % stride == 0."""
        from gatk_sv_ploidy.preprocess import read_site_depth_tsv

        lines = []
        for pos in [0, 50, 100, 150, 200, 250, 300]:
            lines.append(f"chr21\t{pos}\tSAMPLE_A\t10\t0\t0\t20\n")
        sd_path = tmp_path / "test.sd.txt"
        sd_path.write_text("".join(lines))

        df = read_site_depth_tsv(str(sd_path), position_stride=100)
        assert sorted(df["position"].tolist()) == [0, 100, 200, 300]

    def test_position_stride_overrides_stride(self, tmp_path):
        """When position_stride > 0, file-row stride is ignored."""
        from gatk_sv_ploidy.preprocess import read_site_depth_tsv

        lines = []
        for pos in range(0, 500, 10):
            lines.append(f"chr21\t{pos}\tSAMPLE_A\t10\t0\t0\t20\n")
        sd_path = tmp_path / "test.sd.txt"
        sd_path.write_text("".join(lines))

        # position_stride=100 should keep {0, 100, 200, 300, 400} regardless
        # of the file-row stride value.
        df = read_site_depth_tsv(str(sd_path), stride=3, position_stride=100)
        assert sorted(df["position"].tolist()) == [0, 100, 200, 300, 400]


# ---------------------------------------------------------------------------
# Minor/total computation
# ---------------------------------------------------------------------------


class TestSiteMinorTotal:
    """Tests for the internal ``_site_minor_total`` helper."""

    def test_basic(self):
        from gatk_sv_ploidy.preprocess import _site_minor_total

        a = np.array([10, 0, 5])
        c = np.array([0, 20, 5])
        g = np.array([0, 0, 5])
        t = np.array([20, 0, 5])
        minor, total = _site_minor_total(a, c, g, t)
        # Site 0: major=T(20), minor=10, total=30
        assert minor[0] == 10
        assert total[0] == 30
        # Site 1: major=C(20), minor=0, total=20
        assert minor[1] == 0
        assert total[1] == 20
        # Site 2: all equal (5 each), major=5, minor=15, total=20
        assert total[2] == 20
        assert minor[2] == 15


# ---------------------------------------------------------------------------
# Site alt/total with reference base
# ---------------------------------------------------------------------------


class TestSiteAltTotal:
    """Tests for :func:`preprocess._site_alt_total`."""

    def test_basic(self):
        from gatk_sv_ploidy.preprocess import _site_alt_total

        a = np.array([10, 0, 5])
        c = np.array([0, 20, 5])
        g = np.array([0, 0, 5])
        t = np.array([20, 0, 5])
        ref = np.array(["A", "C", "T"])

        alt, total = _site_alt_total(a, c, g, t, ref)
        # Site 0: ref=A(10), alt=total-10=20, total=30
        assert alt[0] == 20
        assert total[0] == 30
        # Site 1: ref=C(20), alt=0, total=20
        assert alt[1] == 0
        assert total[1] == 20
        # Site 2: ref=T(5), alt=15, total=20
        assert alt[2] == 15
        assert total[2] == 20


# ---------------------------------------------------------------------------
# Known-sites reading
# ---------------------------------------------------------------------------


class TestReadKnownSites:
    """Tests for :func:`preprocess.read_known_sites`."""

    def test_read_simple_tsv(self, tmp_path):
        from gatk_sv_ploidy.preprocess import read_known_sites

        path = tmp_path / "sites.tsv"
        path.write_text(
            "chr21\t500\tA\t0.3\n"
            "chr21\t1500\tG\t0.15\n"
            "chrX\t200\tC\t0.45\n"
        )
        df = read_known_sites(str(path))
        assert len(df) == 3
        assert df["pop_af"].iloc[0] == 0.3
        assert df["ref"].iloc[1] == "G"
        assert df["position"].iloc[0] == 500

    def test_read_vcf(self, tmp_path):
        from gatk_sv_ploidy.preprocess import read_known_sites

        vcf_path = tmp_path / "sites.vcf"
        vcf_path.write_text(
            "##fileformat=VCFv4.1\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr21\t501\t.\tA\tG\t.\tPASS\tAF=0.3\n"
            "chrX\t201\t.\tC\tT\t.\tPASS\tAF=0.45\n"
        )
        df = read_known_sites(str(vcf_path))
        assert len(df) == 2
        # VCF is 1-based, should be converted to 0-based
        assert df["position"].iloc[0] == 500
        assert df["pop_af"].iloc[0] == 0.3
        assert df["ref"].iloc[0] == "A"


# ---------------------------------------------------------------------------
# build_per_site_data
# ---------------------------------------------------------------------------


class TestBuildPerSiteData:
    """Tests for :func:`preprocess.build_per_site_data`."""

    def test_no_known_sites_fallback(self, tmp_path):
        from gatk_sv_ploidy.preprocess import build_per_site_data

        bins_df = _make_depth_df()

        sd_path = tmp_path / "sample_a.sd.txt"
        sd_path.write_text(
            "chr21\t500\tSAMPLE_A\t10\t0\t0\t20\n"
            "chr21\t1500\tSAMPLE_A\t8\t0\t0\t12\n"
        )

        result = build_per_site_data(
            [str(sd_path)], bins_df,
            known_sites_df=None,
            min_site_depth=5,
            max_sites_per_bin=10,
        )

        assert "site_alt" in result
        assert "site_total" in result
        assert "site_pop_af" in result
        assert "site_mask" in result
        assert result["site_alt"].shape[0] == len(bins_df)
        assert result["site_alt"].shape[1] == 10  # max_sites_per_bin
        assert result["site_alt"].shape[2] == 2   # n_samples
        # Pop AF should be 0.5 for fallback
        mask_any = result["site_mask"].any(axis=2)
        if mask_any.any():
            assert np.allclose(result["site_pop_af"][mask_any], 0.5)

    def test_with_known_sites(self, tmp_path):
        from gatk_sv_ploidy.preprocess import build_per_site_data, read_known_sites

        bins_df = _make_depth_df()

        sites_path = tmp_path / "sites.tsv"
        sites_path.write_text(
            "chr21\t500\tA\t0.3\n"
            "chr21\t1500\tA\t0.2\n"
        )
        known = read_known_sites(str(sites_path))

        sd_path = tmp_path / "sample_a.sd.txt"
        sd_path.write_text(
            "chr21\t500\tSAMPLE_A\t10\t0\t0\t20\n"   # ref=A, alt=20, total=30
            "chr21\t1500\tSAMPLE_A\t8\t0\t0\t12\n"    # ref=A, alt=12, total=20
            "chr21\t9999\tSAMPLE_A\t10\t0\t0\t20\n"   # not in known sites
        )

        result = build_per_site_data(
            [str(sd_path)], bins_df,
            known_sites_df=known,
            min_site_depth=5,
            max_sites_per_bin=10,
        )

        # Only known sites should be included — site at 9999 should be excluded
        # Check that pop_af for included sites is 0.3 or 0.2, not 0.5
        mask_any = result["site_mask"].any(axis=2)  # (n_bins, max_sites)
        if mask_any.any():
            active_af = result["site_pop_af"][mask_any]
            assert not np.allclose(active_af, 0.5)

    def test_known_sites_preserve_zero_and_low_depth(self, tmp_path):
        from gatk_sv_ploidy.preprocess import build_per_site_data, read_known_sites

        bins_df = _make_depth_df()

        sites_path = tmp_path / "sites.tsv"
        sites_path.write_text(
            "chr21\t500\tA\t0.3\n"
            "chr21\t1500\tA\t0.2\n"
        )
        known = read_known_sites(str(sites_path))

        sd_path = tmp_path / "sample_a.sd.txt"
        sd_path.write_text(
            "chr21\t500\tSAMPLE_A\t0\t0\t0\t0\n"
            "chr21\t1500\tSAMPLE_A\t1\t0\t0\t1\n"
        )

        result = build_per_site_data(
            [str(sd_path)], bins_df,
            known_sites_df=known,
            min_site_depth=5,
            max_sites_per_bin=10,
        )

        bin_sites = np.where(result["site_mask"][0, :, 0])[0]
        assert len(bin_sites) == 2

        first_slot = bin_sites[0]
        second_slot = bin_sites[1]

        assert result["site_total"][0, first_slot, 0] == 0
        assert result["site_alt"][0, first_slot, 0] == 0
        assert result["site_mask"][0, first_slot, 0]

        assert result["site_total"][0, second_slot, 0] == 2
        assert result["site_alt"][0, second_slot, 0] == 1
        assert result["site_mask"][0, second_slot, 0]

    def test_npz_roundtrip(self, tmp_path):
        """Per-site data can be saved and loaded via npz."""
        from gatk_sv_ploidy.preprocess import build_per_site_data

        bins_df = _make_depth_df()
        sd_path = tmp_path / "sample_a.sd.txt"
        sd_path.write_text(
            "chr21\t500\tSAMPLE_A\t10\t0\t0\t20\n"
        )

        result = build_per_site_data(
            [str(sd_path)], bins_df,
            min_site_depth=5, max_sites_per_bin=5,
        )

        npz_path = tmp_path / "site_data.npz"
        np.savez_compressed(str(npz_path), **result)

        from gatk_sv_ploidy.data import load_site_data
        loaded = load_site_data(str(npz_path))
        for key in ("site_alt", "site_total", "site_pop_af", "site_mask"):
            np.testing.assert_array_equal(result[key], loaded[key])


# ---------------------------------------------------------------------------
# DepthData with per-site allele data
# ---------------------------------------------------------------------------


class TestDepthDataSiteData:
    """Tests for :class:`DepthData` with per-site allele data."""

    def test_without_site_data(self):
        from gatk_sv_ploidy.data import DepthData

        df = _make_depth_df()
        data = DepthData(df)
        assert data.site_alt is None
        assert data.site_total is None
        assert data.site_pop_af is None
        assert data.site_mask is None
        assert data.max_sites == 0

    def test_with_site_data(self):
        from gatk_sv_ploidy.data import DepthData

        df = _make_depth_df()
        sd = _make_site_data(n_bins=5, max_sites=3, n_samples=2)
        data = DepthData(df, site_data=sd)

        assert data.site_alt is not None
        assert data.site_total is not None
        assert data.site_pop_af is not None
        assert data.site_mask is not None
        assert data.site_alt.shape == (5, 3, 2)
        assert data.site_total.shape == (5, 3, 2)
        assert data.site_pop_af.shape == (5, 3)
        assert data.site_mask.dtype == torch.bool
        assert data.max_sites == 3

    def test_site_data_survives_subsampling(self):
        from gatk_sv_ploidy.data import DepthData

        df = _make_depth_df()
        sd = _make_site_data(n_bins=5, max_sites=3, n_samples=2)
        data = DepthData(df, site_data=sd, subsample_bins=3)

        assert data.site_alt.shape[0] == 3
        assert data.site_alt.shape[1] == 3  # max_sites unchanged
        assert data.site_alt.shape[2] == 2  # n_samples unchanged


# ---------------------------------------------------------------------------
# Marginalized AF log-likelihood (torch)
# ---------------------------------------------------------------------------


class TestMarginalizedAFLogLik:
    """Tests for :func:`models._marginalized_af_log_lik`."""

    def test_returns_correct_shape(self):
        from gatk_sv_ploidy.models import _marginalized_af_log_lik

        n_bins, max_sites, n_samples = 3, 4, 2
        site_alt = torch.randint(0, 10, (n_bins, max_sites, n_samples)).float()
        site_total = torch.randint(10, 30, (n_bins, max_sites, n_samples)).float()
        site_pop_af = torch.rand(n_bins, max_sites) * 0.5
        site_mask = torch.ones(n_bins, max_sites, n_samples, dtype=torch.bool)
        cn = torch.full((n_bins, n_samples), 2, dtype=torch.long)

        result = _marginalized_af_log_lik(
            site_alt, site_total, site_pop_af, site_mask,
            cn=cn, n_states=6, concentration=50.0,
        )
        assert result.shape == (n_bins, n_samples)

    def test_masked_sites_contribute_zero(self):
        from gatk_sv_ploidy.models import _marginalized_af_log_lik

        site_alt = torch.tensor([[[5], [3], [7]], [[2], [8], [1]]], dtype=torch.float32)
        site_total = torch.tensor([[[20], [15], [25]], [[10], [30], [5]]], dtype=torch.float32)
        site_pop_af = torch.full((2, 3), 0.3)
        cn = torch.full((2, 1), 2, dtype=torch.long)

        # All sites active
        mask_all = torch.ones(2, 3, 1, dtype=torch.bool)
        ll_all = _marginalized_af_log_lik(
            site_alt, site_total, site_pop_af, mask_all,
            cn=cn, n_states=6, concentration=50.0,
        )

        # Mask out middle site in each bin
        mask_partial = mask_all.clone()
        mask_partial[:, 1, :] = False
        ll_partial = _marginalized_af_log_lik(
            site_alt, site_total, site_pop_af, mask_partial,
            cn=cn, n_states=6, concentration=50.0,
        )

        # Partial should have *higher* (less negative) log-lik since fewer sites
        assert (ll_partial >= ll_all).all()

    def test_cn0_gives_finite_likelihood(self):
        """At CN=0, likelihood should be finite (not -inf or nan)."""
        from gatk_sv_ploidy.models import _marginalized_af_log_lik

        site_alt = torch.tensor([[[5], [3]]], dtype=torch.float32)
        site_total = torch.tensor([[[20], [15]]], dtype=torch.float32)
        site_pop_af = torch.tensor([[0.3, 0.3]])
        site_mask = torch.ones(1, 2, 1, dtype=torch.bool)

        cn0 = torch.zeros(1, 1, dtype=torch.long)
        ll0 = _marginalized_af_log_lik(
            site_alt, site_total, site_pop_af, site_mask,
            cn=cn0, n_states=6, concentration=50.0,
        )
        assert torch.isfinite(ll0).all()


# ---------------------------------------------------------------------------
# Marginalized AF log-likelihood (numpy)
# ---------------------------------------------------------------------------


class TestMarginalizedAFLogLikNumpy:
    """Tests for :func:`models._marginalized_af_log_lik_numpy`."""

    def test_returns_correct_shape(self):
        from gatk_sv_ploidy.models import _marginalized_af_log_lik_numpy

        n_bins, max_sites, n_samples = 3, 4, 2
        rng = np.random.RandomState(0)
        site_alt = rng.randint(0, 10, (n_bins, max_sites, n_samples))
        site_total = rng.randint(10, 30, (n_bins, max_sites, n_samples))
        site_pop_af = rng.rand(n_bins, max_sites).astype(np.float32) * 0.5
        site_mask = np.ones((n_bins, max_sites, n_samples), dtype=bool)

        result = _marginalized_af_log_lik_numpy(
            site_alt, site_total, site_pop_af, site_mask,
            cn_state=2, n_states=6, concentration=50.0,
        )
        assert result.shape == (n_bins, n_samples)

    def test_agrees_with_torch(self):
        """NumPy and torch versions should give close results."""
        from gatk_sv_ploidy.models import (
            _marginalized_af_log_lik,
            _marginalized_af_log_lik_numpy,
        )

        rng = np.random.RandomState(123)
        n_bins, max_sites, n_samples = 2, 3, 2
        site_alt_np = rng.randint(0, 10, (n_bins, max_sites, n_samples)).astype(np.int32)
        site_total_np = rng.randint(15, 30, (n_bins, max_sites, n_samples)).astype(np.int32)
        site_pop_af_np = rng.uniform(0.1, 0.4, (n_bins, max_sites)).astype(np.float32)
        site_mask_np = np.ones((n_bins, max_sites, n_samples), dtype=bool)

        for cn_val in [1, 2, 3]:
            np_result = _marginalized_af_log_lik_numpy(
                site_alt_np, site_total_np, site_pop_af_np, site_mask_np,
                cn_state=cn_val, n_states=6, concentration=50.0,
            )

            cn_t = torch.full((n_bins, n_samples), cn_val, dtype=torch.long)
            torch_result = _marginalized_af_log_lik(
                torch.tensor(site_alt_np, dtype=torch.float32),
                torch.tensor(site_total_np, dtype=torch.float32),
                torch.tensor(site_pop_af_np),
                torch.tensor(site_mask_np),
                cn=cn_t, n_states=6, concentration=50.0,
            ).numpy()

            np.testing.assert_allclose(np_result, torch_result, atol=0.1, rtol=0.05)


# ---------------------------------------------------------------------------
# Model with per-site allele data
# ---------------------------------------------------------------------------


class TestCNVModelSiteData:
    """Tests for :class:`CNVModel` with per-site allele data."""

    def test_model_without_site_data(self):
        """Model should work exactly as before when no site data is present."""
        from gatk_sv_ploidy.data import DepthData
        from gatk_sv_ploidy.models import CNVModel
        import pyro

        pyro.clear_param_store()
        df = _make_depth_df()
        data = DepthData(df)
        model = CNVModel(n_states=6, guide_type="delta", af_weight=0.0)

        model.train(data, max_iter=5, log_freq=100)
        estimates = model.get_map_estimates(data)
        assert "cn" in estimates
        assert estimates["cn"].shape == (data.n_bins, data.n_samples)

    def test_model_with_site_data(self):
        """Model should train and produce MAP estimates with per-site data."""
        from gatk_sv_ploidy.data import DepthData
        from gatk_sv_ploidy.models import CNVModel
        import pyro

        pyro.clear_param_store()
        df = _make_depth_df()
        sd = _make_site_data(n_bins=5, max_sites=3, n_samples=2)
        data = DepthData(df, site_data=sd)
        model = CNVModel(
            n_states=6, guide_type="delta",
            af_concentration=50.0, af_weight=1.0,
        )

        model.train(data, max_iter=5, log_freq=100)
        estimates = model.get_map_estimates(data)
        assert "cn" in estimates

    def test_discrete_inference_with_site_data(self):
        """Analytical discrete inference should incorporate per-site AF likelihood."""
        from gatk_sv_ploidy.data import DepthData
        from gatk_sv_ploidy.models import CNVModel
        import pyro

        pyro.clear_param_store()
        df = _make_depth_df()
        sd = _make_site_data(n_bins=5, max_sites=3, n_samples=2)
        data = DepthData(df, site_data=sd)
        model = CNVModel(
            n_states=6, guide_type="delta",
            af_concentration=50.0, af_weight=1.0,
        )
        model.train(data, max_iter=5, log_freq=100)
        estimates = model.get_map_estimates(data)
        posterior = model.run_discrete_inference(data, map_estimates=estimates)

        probs = posterior["cn_posterior"]
        assert probs.shape == (data.n_bins, data.n_samples, 6)
        # Probabilities should sum to 1 across CN states
        np.testing.assert_allclose(
            probs.sum(axis=2), 1.0, atol=1e-5,
        )


# ---------------------------------------------------------------------------
# Call module tetraploid classification
# ---------------------------------------------------------------------------


class TestTetraploidClassification:
    """Tests for tetraploid detection in the call module."""

    def test_classify_sex_tetraploid_female(self):
        from gatk_sv_ploidy.call import _classify_sex

        assert _classify_sex(4, 0) == "TETRAPLOID_FEMALE"

    def test_classify_sex_tetraploid_male(self):
        from gatk_sv_ploidy.call import _classify_sex

        assert _classify_sex(2, 2) == "TETRAPLOID_MALE"

    def test_classify_sex_triple_x_y(self):
        from gatk_sv_ploidy.call import _classify_sex

        assert _classify_sex(3, 1) == "TRIPLE_X_Y"

    def test_classify_aneuploidy_tetraploid(self):
        from gatk_sv_ploidy.call import _classify_aneuploidy

        # All autosomes at CN=4
        aneu_map = {f"chr{i}": True for i in range(1, 23)}
        cn_map = {f"chr{i}": 4 for i in range(1, 23)}
        cn_map["chrX"] = 4
        cn_map["chrY"] = 0
        aneu_map["chrX"] = True
        aneu_map["chrY"] = False

        result = _classify_aneuploidy(aneu_map, cn_map, 4, 0)
        assert result == "TETRAPLOID"

    def test_classify_aneuploidy_tetrasomy_21(self):
        from gatk_sv_ploidy.call import _classify_aneuploidy

        aneu_map = {"chr21": True, "chr1": False, "chr2": False}
        cn_map = {"chr21": 4, "chr1": 2, "chr2": 2, "chrX": 2, "chrY": 0}
        result = _classify_aneuploidy(aneu_map, cn_map, 2, 0)
        assert result == "TETRASOMY_21"

    def test_classify_aneuploidy_normal_unchanged(self):
        from gatk_sv_ploidy.call import _classify_aneuploidy

        aneu_map = {"chr1": False, "chr21": False}
        cn_map = {"chr1": 2, "chr21": 2}
        result = _classify_aneuploidy(aneu_map, cn_map, 2, 0)
        assert result == "NORMAL"


# ---------------------------------------------------------------------------
# _util constants validation
# ---------------------------------------------------------------------------


class TestUtilConstants:
    """Validate the new constants in _util."""

    def test_max_genotype_states(self):
        from gatk_sv_ploidy._util import MAX_GENOTYPE_STATES
        assert MAX_GENOTYPE_STATES == 6

    def test_default_af_concentration(self):
        from gatk_sv_ploidy._util import DEFAULT_AF_CONCENTRATION
        assert DEFAULT_AF_CONCENTRATION > 0

    def test_default_af_weight(self):
        from gatk_sv_ploidy._util import DEFAULT_AF_WEIGHT
        assert DEFAULT_AF_WEIGHT > 0
