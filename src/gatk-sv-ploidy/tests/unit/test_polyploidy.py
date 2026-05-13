from __future__ import annotations

import logging

import numpy as np
import pandas as pd

from gatk_sv_ploidy import polyploidy
from gatk_sv_ploidy.models import _marginalized_af_log_lik_numpy


def _synthetic_ploidy_inputs() -> tuple[
    list[str],
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    sample_ids = [
        "haploid_sample",
        "diploid_sample",
        "triploid_sample",
        "tetraploid_sample",
    ]
    bin_chr = np.array([f"chr{i + 1}" for i in range(6)], dtype=object)
    site_total = np.full((6, 6, 4), 60, dtype=np.int32)
    site_alt = np.empty((6, 6, 4), dtype=np.int32)
    site_alt[:, :, 0] = np.array([0, 60, 0, 60, 0, 60], dtype=np.int32)
    site_alt[:, :, 1] = 30
    site_alt[:, :, 2] = np.array([20, 40, 20, 40, 20, 40], dtype=np.int32)
    site_alt[:, :, 3] = np.array([15, 45, 15, 45, 15, 45], dtype=np.int32)
    site_pop_af = np.full((6, 6), 0.5, dtype=np.float64)
    site_mask = np.ones((6, 6, 4), dtype=bool)
    return sample_ids, bin_chr, site_alt, site_total, site_pop_af, site_mask


def _classify_synthetic_ploidy() -> pd.DataFrame:
    sample_ids, bin_chr, site_alt, site_total, site_pop_af, site_mask = (
        _synthetic_ploidy_inputs()
    )
    return polyploidy.classify_polyploidy_from_site_data(
        sample_ids=sample_ids,
        bin_chr=bin_chr,
        site_alt=site_alt,
        site_total=site_total,
        site_pop_af=site_pop_af,
        site_mask=site_mask,
        haploidy_prior=0.2,
        triploidy_prior=0.2,
        tetraploidy_prior=0.2,
        min_informative_bins=3,
        min_informative_sites=18,
        pvalue_threshold=0.2,
        effect_size_threshold=0.0,
    )


def _reference_summed_af_log_lik_by_concentration(
    *,
    site_alt: np.ndarray,
    site_total: np.ndarray,
    site_pop_af: np.ndarray,
    site_mask: np.ndarray,
    cn_state: int,
    n_states: int,
    concentration_grid: np.ndarray,
) -> np.ndarray:
    return np.vstack([
        _marginalized_af_log_lik_numpy(
            site_alt,
            site_total,
            site_pop_af,
            site_mask,
            cn_state=cn_state,
            n_states=n_states,
            concentration=float(concentration),
        ).sum(axis=0)
        for concentration in concentration_grid
    ])


def test_fast_af_likelihood_matches_reference_for_masked_inputs() -> None:
    rng = np.random.default_rng(31415)
    site_total = rng.integers(20, 80, size=(4, 5, 3), dtype=np.int32)
    site_alt = rng.binomial(site_total, 0.45).astype(np.int32)
    site_mask = rng.random(size=(4, 5, 3)) > 0.25
    site_total[~site_mask] = 0
    site_alt[~site_mask] = 0
    site_pop_af_2d = rng.uniform(0.15, 0.85, size=(4, 5))
    sample_offsets = np.array([0.00, 0.03, -0.02], dtype=np.float64)
    site_pop_af_3d = np.clip(
        site_pop_af_2d[:, :, np.newaxis] + sample_offsets[np.newaxis, np.newaxis, :],
        0.05,
        0.95,
    )
    concentration_grid = np.array([4.0, 25.0, 80.0], dtype=np.float64)

    for site_pop_af in (site_pop_af_2d, site_pop_af_3d):
        for cn_state in (1, 2, 3, 4):
            expected = _reference_summed_af_log_lik_by_concentration(
                site_alt=site_alt,
                site_total=site_total,
                site_pop_af=site_pop_af,
                site_mask=site_mask,
                cn_state=cn_state,
                n_states=6,
                concentration_grid=concentration_grid,
            )
            observed = polyploidy._summed_af_log_lik_by_concentration(
                site_alt=site_alt,
                site_total=site_total,
                site_pop_af=site_pop_af,
                site_mask=site_mask,
                cn_state=cn_state,
                n_states=6,
                concentration_grid=concentration_grid,
            )
            assert np.allclose(observed, expected, rtol=1e-8, atol=1e-8)


def test_classify_polyploidy_from_site_data_classifies_cn1_to_cn4() -> None:
    results = _classify_synthetic_ploidy().set_index("sample")

    expected = {
        "haploid_sample": (1, "HAPLOID"),
        "diploid_sample": (2, "DIPLOID"),
        "triploid_sample": (3, "TRIPLOID"),
        "tetraploid_sample": (4, "TETRAPLOID"),
    }
    for sample_id, (expected_cn, expected_call) in expected.items():
        assert int(results.loc[sample_id, "autosomal_baseline_cn"]) == expected_cn
        assert results.loc[sample_id, "baseline_cn_call"] == expected_call
        assert float(results.loc[sample_id, "baseline_cn_posterior_probability"]) > 0.8

    assert float(results.loc["haploid_sample", "fraction_near_haploid_endpoints"]) == 1.0
    assert float(results.loc["triploid_sample", "fraction_near_triploid_thirds"]) == 1.0
    assert float(results.loc["tetraploid_sample", "fraction_near_tetraploid_quarters"]) == 1.0
    assert bool(results.loc["tetraploid_sample", "cn4_direct_peak_supported"]) is True
    assert float(results.loc["tetraploid_sample", "peak_evidence_weight"]) == 0.0
    cn4_quarter_fraction = float(
        results.loc["tetraploid_sample", "tetraploid_quarter_genotype_fraction"]
    )
    cn4_half_fraction = float(
        results.loc["tetraploid_sample", "tetraploid_half_genotype_fraction"]
    )
    assert cn4_quarter_fraction > 0.95
    assert cn4_half_fraction < 0.05
    assert (
        float(
            results.loc[
                "tetraploid_sample",
                "tetraploid_expected_quarter_genotype_fraction",
            ]
        ) > 0.45
    )
    assert (
        float(
            results.loc[
                "tetraploid_sample",
                "tetraploid_quarter_effective_site_fraction",
            ]
        ) > 0.5
    )
    posterior_sum = results[
        [
            "posterior_cn_1",
            "posterior_cn_2",
            "posterior_cn_3",
            "posterior_cn_4",
        ]
    ].sum(axis=1)
    assert np.allclose(posterior_sum, 1.0)


def test_classify_polyploidy_defaults_cn4_ambiguous_half_peaks_to_diploid() -> None:
    sample_ids = ["tetraploid_ambiguous"]
    bin_chr = np.array([f"chr{i + 1}" for i in range(8)], dtype=object)
    site_total = np.full((8, 8, 1), 60, dtype=np.int32)
    site_alt = np.full((8, 8, 1), 30, dtype=np.int32)
    site_pop_af = np.full((8, 8), 0.5, dtype=np.float64)
    site_mask = np.ones((8, 8, 1), dtype=bool)

    results = polyploidy.classify_polyploidy_from_site_data(
        sample_ids=sample_ids,
        bin_chr=bin_chr,
        site_alt=site_alt,
        site_total=site_total,
        site_pop_af=site_pop_af,
        site_mask=site_mask,
        haploidy_prior=0.001,
        triploidy_prior=0.001,
        tetraploidy_prior=0.49,
        min_informative_bins=4,
        min_informative_sites=32,
        pvalue_threshold=0.5,
        effect_size_threshold=0.0,
    ).set_index("sample")

    assert int(results.loc["tetraploid_ambiguous", "autosomal_baseline_cn"]) == 2
    assert results.loc["tetraploid_ambiguous", "baseline_cn_call"] == "DIPLOID"
    assert bool(results.loc["tetraploid_ambiguous", "cn4_direct_peak_supported"]) is False
    assert float(results.loc["tetraploid_ambiguous", "fraction_near_diploid_half"]) == 1.0
    assert float(results.loc["tetraploid_ambiguous", "fraction_near_tetraploid_quarters"]) == 0.0
    assert (
        float(
            results.loc[
                "tetraploid_ambiguous",
                "tetraploid_quarter_genotype_fraction",
            ]
        ) < 0.05
    )
    assert (
        float(
            results.loc[
                "tetraploid_ambiguous",
                "tetraploid_half_genotype_fraction",
            ]
        ) > 0.95
    )
    assert (
        float(
            results.loc[
                "tetraploid_ambiguous",
                "tetraploid_expected_quarter_genotype_fraction",
            ]
        ) > 0.45
    )
    assert (
        float(
            results.loc[
                "tetraploid_ambiguous",
                "tetraploid_quarter_observed_expected_ratio",
            ]
        ) < 0.1
    )


def test_tetraploid_direct_support_requires_raw_quarter_peaks() -> None:
    metrics = {
        "fraction_near_diploid_half": 1.0,
        "fraction_near_triploid_thirds": 0.0,
        "fraction_near_tetraploid_quarters": 0.0,
        "tetraploid_quarter_peak_fraction_advantage": -1.0,
        "tetraploid_quarter_genotype_fraction": 0.53,
        "tetraploid_quarter_genotype_advantage": 0.055,
        "tetraploid_quarter_effective_site_fraction": 0.29,
    }

    assert (
        polyploidy._has_direct_state_peak_support(
            metrics,
            4,
            min_haploid_endpoint_fraction=0.85,
            max_haploid_other_peak_fraction=0.05,
            min_triploid_peak_fraction_advantage=0.0,
            min_tetraploid_quarter_peak_fraction=0.10,
            min_tetraploid_quarter_peak_fraction_advantage=0.05,
        )
        is False
    )


def test_classify_polyploidy_defaults_low_fit_cn4_candidate_to_diploid() -> None:
    sample_ids = ["low_concentration_diploid_like"]
    bin_chr = np.array([f"chr{i + 1}" for i in range(8)], dtype=object)
    site_total = np.full((8, 8, 1), 100, dtype=np.int32)
    site_alt = np.full((8, 8, 1), 50, dtype=np.int32)
    site_pop_af = np.full((8, 8), 0.30, dtype=np.float64)
    site_mask = np.ones((8, 8, 1), dtype=bool)

    results = polyploidy.classify_polyploidy_from_site_data(
        sample_ids=sample_ids,
        bin_chr=bin_chr,
        site_alt=site_alt,
        site_total=site_total,
        site_pop_af=site_pop_af,
        site_mask=site_mask,
        af_concentration=2.489,
        af_concentration_grid=[2.489],
        haploidy_prior=0.001,
        triploidy_prior=0.001,
        tetraploidy_prior=0.49,
        min_diploid_het_prior=0.0,
        min_informative_bins=4,
        min_informative_sites=32,
        pvalue_threshold=0.9,
        effect_size_threshold=-100.0,
    ).set_index("sample")

    sample = "low_concentration_diploid_like"
    assert int(results.loc[sample, "autosomal_baseline_cn"]) == 2
    assert results.loc[sample, "baseline_cn_call"] == "DIPLOID"
    assert results.loc[sample, "baseline_cn_reason"] == (
        "cn4_empirical_peak_support_absent_diploid_default"
    )
    assert bool(results.loc[sample, "include_in_infer"]) is True
    assert bool(results.loc[sample, "cn4_posterior_supported"]) is True
    assert bool(results.loc[sample, "cn4_direct_peak_supported"]) is False
    assert bool(results.loc[sample, "cn4_af_fit_supported"]) is False
    assert float(results.loc[sample, "fraction_near_diploid_half"]) == 1.0
    assert float(results.loc[sample, "fraction_near_tetraploid_quarters"]) == 0.0
    assert float(results.loc[sample, "tetraploid_quarter_genotype_fraction"]) > 0.5
    assert float(results.loc[sample, "tetraploid_quarter_genotype_advantage"]) > 0.05


def test_classify_polyploidy_defaults_low_fit_cn3_candidate_to_diploid() -> None:
    sample_ids = ["low_fit_triploid_map"]
    bin_chr = np.array([f"chr{i + 1}" for i in range(10)], dtype=object)
    site_total = np.full((10, 10, 1), 120, dtype=np.int32)
    site_alt = np.empty((10, 10, 1), dtype=np.int32)
    site_alt[:, :4, 0] = 40
    site_alt[:, 4:, 0] = 60
    site_pop_af = np.full((10, 10), 0.30, dtype=np.float64)
    site_mask = np.ones((10, 10, 1), dtype=bool)

    results = polyploidy.classify_polyploidy_from_site_data(
        sample_ids=sample_ids,
        bin_chr=bin_chr,
        site_alt=site_alt,
        site_total=site_total,
        site_pop_af=site_pop_af,
        site_mask=site_mask,
        af_concentration=2.5,
        haploidy_prior=1e-12,
        triploidy_prior=0.49,
        tetraploidy_prior=1e-12,
        min_diploid_het_prior=0.0,
        min_informative_bins=4,
        min_informative_sites=50,
        pvalue_threshold=0.9,
        effect_size_threshold=-100.0,
    ).set_index("sample")

    sample = "low_fit_triploid_map"
    assert int(results.loc[sample, "candidate_baseline_cn"]) == 3
    assert int(results.loc[sample, "autosomal_baseline_cn"]) == 2
    assert results.loc[sample, "baseline_cn_call"] == "DIPLOID"
    assert results.loc[sample, "baseline_cn_reason"] == (
        "cn3_empirical_peak_support_absent_diploid_default"
    )
    assert bool(results.loc[sample, "include_in_infer"]) is True
    assert bool(results.loc[sample, "cn3_posterior_supported"]) is True
    assert bool(results.loc[sample, "cn3_direct_peak_supported"]) is False
    assert bool(results.loc[sample, "cn3_af_fit_supported"]) is False


def test_classify_polyploidy_marks_low_fit_cn3_candidate_undetermined() -> None:
    sample_ids = ["contamination_like_triploid_candidate"]
    bin_chr = np.array([f"chr{i + 1}" for i in range(10)], dtype=object)
    site_total = np.full((10, 10, 1), 120, dtype=np.int32)
    site_alt = np.empty((10, 10, 1), dtype=np.int32)
    site_alt[:, :8, 0] = 40
    site_alt[:, 8:, 0] = 60
    site_pop_af = np.full((10, 10), 0.5, dtype=np.float64)
    site_mask = np.ones((10, 10, 1), dtype=bool)

    results = polyploidy.classify_polyploidy_from_site_data(
        sample_ids=sample_ids,
        bin_chr=bin_chr,
        site_alt=site_alt,
        site_total=site_total,
        site_pop_af=site_pop_af,
        site_mask=site_mask,
        haploidy_prior=0.001,
        triploidy_prior=0.49,
        tetraploidy_prior=0.001,
        min_informative_bins=4,
        min_informative_sites=50,
        pvalue_threshold=0.9,
        effect_size_threshold=-100.0,
        min_triploid_peak_fraction=0.90,
    ).set_index("sample")

    sample = "contamination_like_triploid_candidate"
    assert int(results.loc[sample, "autosomal_baseline_cn"]) == 0
    assert results.loc[sample, "baseline_cn_call"] == "UNDETERMINED"
    assert results.loc[sample, "baseline_cn_reason"] == "cn3_af_shape_misfit"
    assert bool(results.loc[sample, "include_in_infer"]) is False
    assert bool(results.loc[sample, "cn3_posterior_supported"]) is True
    assert bool(results.loc[sample, "cn3_af_fit_supported"]) is False
    assert int(results.loc[sample, "candidate_baseline_cn"]) == 3


def test_classify_polyploidy_from_site_data_handles_overdispersed_diploid_af() -> None:
    sample_ids = ["overdispersed_diploid"]
    bin_chr = np.array([f"chr{i + 1}" for i in range(8)], dtype=object)
    site_total = np.full((8, 8, 1), 60, dtype=np.int32)
    site_alt = np.empty((8, 8, 1), dtype=np.int32)
    site_alt[:, :, 0] = np.array(
        [24, 36, 25, 35, 23, 37, 26, 34],
        dtype=np.int32,
    )
    site_pop_af = np.full((8, 8), 0.5, dtype=np.float64)
    site_mask = np.ones((8, 8, 1), dtype=bool)

    results = polyploidy.classify_polyploidy_from_site_data(
        sample_ids=sample_ids,
        bin_chr=bin_chr,
        site_alt=site_alt,
        site_total=site_total,
        site_pop_af=site_pop_af,
        site_mask=site_mask,
        min_informative_bins=4,
        min_informative_sites=32,
        effect_size_threshold=0.0,
    ).set_index("sample")

    assert int(results.loc["overdispersed_diploid", "autosomal_baseline_cn"]) == 2
    assert results.loc["overdispersed_diploid", "baseline_cn_call"] == "DIPLOID"
    assert float(results.loc["overdispersed_diploid", "posterior_cn_3"]) < 0.5
    assert bool(results.loc["overdispersed_diploid", "cn3_direct_peak_supported"]) is False


def test_privacy_safe_polyploidy_diagnostics_do_not_log_sample_ids(caplog) -> None:
    sample_ids, bin_chr, site_alt, site_total, site_pop_af, site_mask = (
        _synthetic_ploidy_inputs()
    )
    results = _classify_synthetic_ploidy()
    metrics = polyploidy.build_polyploidy_diagnostic_metrics(
        sample_ids=sample_ids,
        bin_chr=bin_chr,
        site_alt=site_alt,
        site_total=site_total,
        site_pop_af=site_pop_af,
        site_mask=site_mask,
        results_df=results,
        min_diploid_het_prior=0.1,
    )

    caplog.set_level(logging.INFO, logger=polyploidy.logger.name)
    polyploidy.log_privacy_safe_polyploidy_diagnostics(
        results,
        pvalue_threshold=0.2,
        effect_size_threshold=0.0,
        metrics_df=metrics,
    )

    log_text = "\n".join(record.getMessage() for record in caplog.records)
    assert "PRIVACY_SAFE_POLYPLOIDY" in log_text
    assert "autosomal_baseline_cn_counts" in log_text
    assert "posterior_threshold_counts state=cn4" in log_text
    assert "tetraploid_ambiguity_counts" in log_text
    assert "af_peak_counts" in log_text
    for sample_id in sample_ids:
        assert sample_id not in log_text


def test_polyploidy_progress_messages_do_not_log_sample_ids(monkeypatch, caplog) -> None:
    sample_ids, bin_chr, site_alt, site_total, site_pop_af, site_mask = (
        _synthetic_ploidy_inputs()
    )
    progress_calls = []

    def fake_tqdm(iterable, **kwargs):
        progress_calls.append(kwargs)
        return iterable

    monkeypatch.setattr(polyploidy, "tqdm", fake_tqdm)
    caplog.set_level(logging.INFO, logger=polyploidy.logger.name)

    polyploidy.classify_polyploidy_from_site_data(
        sample_ids=sample_ids,
        bin_chr=bin_chr,
        site_alt=site_alt,
        site_total=site_total,
        site_pop_af=site_pop_af,
        site_mask=site_mask,
        af_concentration_grid=[10.0],
        haploidy_prior=0.2,
        triploidy_prior=0.2,
        tetraploidy_prior=0.2,
        min_informative_bins=3,
        min_informative_sites=18,
        pvalue_threshold=0.2,
        effect_size_threshold=0.0,
        show_progress=True,
        log_progress=True,
    )

    descriptions = [str(call.get("desc", "")) for call in progress_calls]
    assert descriptions == [
        "CN1 AF likelihood",
        "CN2 AF likelihood",
        "CN3 AF likelihood",
        "CN4 AF likelihood",
        "Baseline CN peak evidence",
        "Classifying samples",
    ]

    log_text = "\n".join(record.getMessage() for record in caplog.records)
    assert "Polyploidy progress: prepared aggregate AF matrix" in log_text
    assert "baseline_cn_states=1,2,3,4" in log_text
    assert "Polyploidy progress: classification complete" in log_text
    for sample_id in sample_ids:
        assert sample_id not in log_text
        assert all(sample_id not in description for description in descriptions)


def test_polyploidy_main_writes_manifest_and_results(tmp_path, monkeypatch) -> None:
    sample_ids, bin_chr, site_alt, site_total, site_pop_af, site_mask = (
        _synthetic_ploidy_inputs()
    )

    depth_df = pd.DataFrame(
        {
            "Chr": bin_chr,
            "Start": np.arange(0, 6000, 1000),
            "End": np.arange(1000, 7000, 1000),
            sample_ids[0]: [1.0] * 6,
            sample_ids[1]: [1.0] * 6,
            sample_ids[2]: [1.0] * 6,
            sample_ids[3]: [1.0] * 6,
        },
        index=[f"chr{i + 1}:{i * 1000}-{(i + 1) * 1000}" for i in range(6)],
    )
    depth_path = tmp_path / "preprocessed_depth.tsv"
    depth_df.to_csv(depth_path, sep="\t")

    site_data_path = tmp_path / "site_data.npz"
    np.savez_compressed(
        site_data_path,
        site_alt=site_alt,
        site_total=site_total,
        site_pop_af=site_pop_af,
        site_mask=site_mask,
        sample_ids=np.asarray(sample_ids, dtype=object),
        bin_chr=bin_chr,
    )

    output_dir = tmp_path / "polyploidy_out"
    monkeypatch.setattr(
        "sys.argv",
        [
            "polyploidy",
            "--input",
            str(depth_path),
            "--site-data",
            str(site_data_path),
            "--output-dir",
            str(output_dir),
            "--min-informative-bins",
            "3",
            "--min-informative-sites",
            "18",
            "--pvalue-threshold",
            "0.2",
            "--effect-size-threshold",
            "0.0",
            "--haploidy-prior",
            "0.2",
            "--triploidy-prior",
            "0.2",
            "--tetraploidy-prior",
            "0.2",
            "--no-progress",
        ],
    )

    polyploidy.main()

    manifest_df = pd.read_csv(output_dir / "sample_autosomal_baseline_cn.tsv", sep="\t")
    results_df = pd.read_csv(output_dir / "polyploidy_test_results.tsv", sep="\t")

    manifest = dict(zip(manifest_df["sample"], manifest_df["autosomal_baseline_cn"]))
    calls = dict(zip(results_df["sample"], results_df["baseline_cn_call"]))

    assert manifest == {
        "haploid_sample": 1,
        "diploid_sample": 2,
        "triploid_sample": 3,
        "tetraploid_sample": 4,
    }
    assert calls == {
        "haploid_sample": "HAPLOID",
        "diploid_sample": "DIPLOID",
        "triploid_sample": "TRIPLOID",
        "tetraploid_sample": "TETRAPLOID",
    }
    assert "include_in_infer" in manifest_df.columns
    assert manifest_df["include_in_infer"].all()


def test_polyploidy_main_writes_diagnostics(tmp_path, monkeypatch) -> None:
    sample_ids, bin_chr, site_alt, site_total, site_pop_af, site_mask = (
        _synthetic_ploidy_inputs()
    )

    depth_df = pd.DataFrame(
        {
            "Chr": bin_chr,
            "Start": np.arange(0, 6000, 1000),
            "End": np.arange(1000, 7000, 1000),
            sample_ids[0]: [1.0] * 6,
            sample_ids[1]: [1.0] * 6,
            sample_ids[2]: [1.0] * 6,
            sample_ids[3]: [1.0] * 6,
        },
        index=[f"chr{i + 1}:{i * 1000}-{(i + 1) * 1000}" for i in range(6)],
    )
    depth_path = tmp_path / "preprocessed_depth.tsv"
    depth_df.to_csv(depth_path, sep="\t")

    site_data_path = tmp_path / "site_data.npz"
    np.savez_compressed(
        site_data_path,
        site_alt=site_alt,
        site_total=site_total,
        site_pop_af=site_pop_af,
        site_mask=site_mask,
        sample_ids=np.asarray(sample_ids, dtype=object),
        bin_chr=bin_chr,
    )

    output_dir = tmp_path / "polyploidy_out"
    monkeypatch.setattr(
        "sys.argv",
        [
            "polyploidy",
            "--input",
            str(depth_path),
            "--site-data",
            str(site_data_path),
            "--output-dir",
            str(output_dir),
            "--min-informative-bins",
            "3",
            "--min-informative-sites",
            "18",
            "--pvalue-threshold",
            "0.2",
            "--effect-size-threshold",
            "0.0",
            "--haploidy-prior",
            "0.2",
            "--triploidy-prior",
            "0.2",
            "--tetraploidy-prior",
            "0.2",
            "--diagnostics",
            "--diagnostic-sample-limit",
            "4",
            "--diagnostic-max-sites-per-sample",
            "10",
            "--no-progress",
        ],
    )

    polyploidy.main()

    diagnostics_dir = output_dir / "diagnostics"
    metrics_df = pd.read_csv(
        diagnostics_dir / "polyploidy_diagnostic_metrics.tsv",
        sep="\t",
    ).set_index("sample")
    raw_site_df = pd.read_csv(
        diagnostics_dir / "polyploidy_raw_af_sites.tsv.gz",
        sep="\t",
    )

    assert float(metrics_df.loc["haploid_sample", "fraction_near_haploid_endpoints"]) == 1.0
    assert float(metrics_df.loc["diploid_sample", "fraction_near_diploid_half"]) == 1.0
    assert float(metrics_df.loc["triploid_sample", "fraction_near_triploid_thirds"]) == 1.0
    assert float(metrics_df.loc["tetraploid_sample", "fraction_near_tetraploid_quarters"]) == 1.0
    diagnostic_cn4_quarter_fraction = float(
        metrics_df.loc["tetraploid_sample", "tetraploid_quarter_genotype_fraction"]
    )
    assert diagnostic_cn4_quarter_fraction > 0.95
    assert set(raw_site_df["sample"]) == set(sample_ids)
    assert raw_site_df.groupby("sample").size().max() <= 10
    assert (diagnostics_dir / "polyploidy_diagnostic_metrics.png").exists()
    assert (diagnostics_dir / "polyploidy_raw_af_monoploid_profiles.png").exists()
    assert (diagnostics_dir / "polyploidy_raw_af_diploid_profiles.png").exists()
    assert (diagnostics_dir / "polyploidy_raw_af_triploid_profiles.png").exists()
    assert (diagnostics_dir / "polyploidy_raw_af_tetraploid_profiles.png").exists()


def test_write_polyploidy_diagnostics_includes_undetermined_af_profiles(
    tmp_path,
    monkeypatch,
) -> None:
    sample_ids = ["diploid_sample", "undetermined_sample"]
    metrics_df = pd.DataFrame(
        {
            "sample": sample_ids,
            "autosomal_baseline_cn": [2, 0],
            "baseline_cn_call": ["DIPLOID", "UNDETERMINED"],
            "baseline_cn_posterior_probability": [0.98, np.nan],
            "candidate_baseline_cn_posterior_probability": [0.98, 0.91],
            "diagnostic_informative_sites": [120, 140],
            "posterior_cn_1": [0.01, 0.01],
            "posterior_cn_3": [0.01, 0.91],
            "posterior_cn_4": [0.0, 0.03],
            "cn1_direct_peak_supported": [False, False],
            "cn3_direct_peak_supported": [False, True],
            "cn4_direct_peak_supported": [False, False],
            "fraction_near_diploid_half": [0.95, 0.35],
            "fraction_near_triploid_thirds": [0.02, 0.28],
            "fraction_near_tetraploid_quarters": [0.01, 0.03],
            "fraction_near_diploid_compatible": [0.98, 0.62],
            "fraction_near_haploid_endpoints": [0.01, 0.02],
            "fraction_near_any_modeled_peak": [0.99, 0.74],
        }
    )

    def fake_build_polyploidy_diagnostic_metrics(**_kwargs) -> pd.DataFrame:
        return metrics_df.copy()

    def fake_collect_diagnostic_site_rows(*, selected_samples, **_kwargs) -> pd.DataFrame:
        rows = []
        for sample_id in selected_samples:
            rows.append(
                {
                    "sample": sample_id,
                    "chr": "chr1",
                    "autosomal_bin_index": 0,
                    "site_index": 0,
                    "site_alt": 20,
                    "site_total": 60,
                    "observed_af": 1.0 / 3.0,
                    "folded_observed_af": 1.0 / 3.0,
                }
            )
        return pd.DataFrame(rows)

    monkeypatch.setattr(
        polyploidy,
        "build_polyploidy_diagnostic_metrics",
        fake_build_polyploidy_diagnostic_metrics,
    )
    monkeypatch.setattr(
        polyploidy,
        "_collect_diagnostic_site_rows",
        fake_collect_diagnostic_site_rows,
    )

    polyploidy.write_polyploidy_diagnostics(
        output_dir=str(tmp_path),
        sample_ids=sample_ids,
        bin_chr=np.array(["chr1"], dtype=object),
        site_alt=np.zeros((1, 1, 2), dtype=np.int32),
        site_total=np.ones((1, 1, 2), dtype=np.int32),
        site_pop_af=np.full((1, 1), 0.5, dtype=np.float64),
        site_mask=np.ones((1, 1, 2), dtype=bool),
        results_df=metrics_df,
        min_diploid_het_prior=0.1,
        af_window=0.08,
        sample_limit=4,
        max_sites_per_sample=10,
        show_progress=False,
    )

    diagnostics_dir = tmp_path / "diagnostics"
    assert (diagnostics_dir / "polyploidy_raw_af_diploid_profiles.png").exists()
    assert (diagnostics_dir / "polyploidy_raw_af_undetermined_profiles.png").exists()
