from __future__ import annotations

import logging

import numpy as np
import pandas as pd

from gatk_sv_ploidy import triploidy


def _synthetic_triploidy_inputs() -> tuple[list[str], np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    sample_ids = ["diploid_sample", "triploid_sample"]
    bin_chr = np.array(["chr1", "chr2", "chr3", "chr4"], dtype=object)
    site_total = np.full((4, 4, 2), 60, dtype=np.int32)
    site_alt = np.empty((4, 4, 2), dtype=np.int32)
    site_alt[:, :, 0] = 30
    site_alt[:, :, 1] = np.array(
        [
            [20, 40, 20, 40],
            [20, 40, 20, 40],
            [20, 40, 20, 40],
            [20, 40, 20, 40],
        ],
        dtype=np.int32,
    )
    site_pop_af = np.full((4, 4), 0.5, dtype=np.float64)
    site_mask = np.ones((4, 4, 2), dtype=bool)
    return sample_ids, bin_chr, site_alt, site_total, site_pop_af, site_mask


def test_classify_triploidy_from_site_data_distinguishes_diploid_and_triploid() -> None:
    sample_ids, bin_chr, site_alt, site_total, site_pop_af, site_mask = _synthetic_triploidy_inputs()

    results = triploidy.classify_triploidy_from_site_data(
        sample_ids=sample_ids,
        bin_chr=bin_chr,
        site_alt=site_alt,
        site_total=site_total,
        site_pop_af=site_pop_af,
        site_mask=site_mask,
        triploidy_prior=0.5,
        min_informative_bins=2,
        min_informative_sites=8,
        pvalue_threshold=0.05,
        effect_size_threshold=0.0,
    ).set_index("sample")

    assert int(results.loc["diploid_sample", "autosomal_baseline_cn"]) == 2
    assert results.loc["diploid_sample", "triploidy_call"] == "DIPLOID"
    assert float(results.loc["diploid_sample", "triploidy_log_lik_ratio"]) < 0.0

    assert int(results.loc["triploid_sample", "autosomal_baseline_cn"]) == 3
    assert results.loc["triploid_sample", "triploidy_call"] == "TRIPLOID"
    assert float(results.loc["triploid_sample", "triploidy_log_lik_ratio"]) > 0.0
    assert float(results.loc["triploid_sample", "triploidy_posterior_probability"]) > 0.95


def test_classify_triploidy_from_site_data_handles_overdispersed_diploid_af() -> None:
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

    results = triploidy.classify_triploidy_from_site_data(
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
    assert results.loc["overdispersed_diploid", "triploidy_call"] == "DIPLOID"
    assert float(results.loc["overdispersed_diploid", "triploidy_posterior_probability"]) < 0.5
    assert float(results.loc["overdispersed_diploid", "triploidy_log_bayes_factor"]) < 0.0


def test_privacy_safe_triploidy_diagnostics_do_not_log_sample_ids(caplog) -> None:
    sample_ids, bin_chr, site_alt, site_total, site_pop_af, site_mask = _synthetic_triploidy_inputs()
    results = triploidy.classify_triploidy_from_site_data(
        sample_ids=sample_ids,
        bin_chr=bin_chr,
        site_alt=site_alt,
        site_total=site_total,
        site_pop_af=site_pop_af,
        site_mask=site_mask,
        triploidy_prior=0.5,
        min_informative_bins=2,
        min_informative_sites=8,
        pvalue_threshold=0.05,
        effect_size_threshold=0.0,
    )
    metrics = triploidy.build_triploidy_diagnostic_metrics(
        sample_ids=sample_ids,
        bin_chr=bin_chr,
        site_alt=site_alt,
        site_total=site_total,
        site_pop_af=site_pop_af,
        site_mask=site_mask,
        results_df=results,
        min_diploid_het_prior=0.1,
    )

    caplog.set_level(logging.INFO, logger=triploidy.logger.name)
    triploidy.log_privacy_safe_triploidy_diagnostics(
        results,
        pvalue_threshold=0.05,
        effect_size_threshold=0.0,
        metrics_df=metrics,
    )

    log_text = "\n".join(record.getMessage() for record in caplog.records)
    assert "PRIVACY_SAFE_TRIPLOIDY" in log_text
    assert "autosomal_baseline_cn_counts" in log_text
    assert "af_peak_counts" in log_text
    for sample_id in sample_ids:
        assert sample_id not in log_text


def test_triploidy_main_writes_manifest_and_results(tmp_path, monkeypatch) -> None:
    sample_ids, bin_chr, site_alt, site_total, site_pop_af, site_mask = _synthetic_triploidy_inputs()

    depth_df = pd.DataFrame(
        {
            "Chr": bin_chr,
            "Start": [0, 1000, 2000, 3000],
            "End": [1000, 2000, 3000, 4000],
            sample_ids[0]: [1.0, 1.0, 1.0, 1.0],
            sample_ids[1]: [1.0, 1.0, 1.0, 1.0],
        },
        index=[
            "chr1:0-1000",
            "chr2:1000-2000",
            "chr3:2000-3000",
            "chr4:3000-4000",
        ],
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

    output_dir = tmp_path / "triploidy_out"
    monkeypatch.setattr(
        "sys.argv",
        [
            "triploidy",
            "--input",
            str(depth_path),
            "--site-data",
            str(site_data_path),
            "--output-dir",
            str(output_dir),
            "--min-informative-bins",
            "2",
            "--min-informative-sites",
            "8",
            "--pvalue-threshold",
            "0.05",
            "--effect-size-threshold",
            "0.0",
            "--triploidy-prior",
            "0.5",
        ],
    )

    triploidy.main()

    manifest_df = pd.read_csv(output_dir / "sample_autosomal_baseline_cn.tsv", sep="\t")
    results_df = pd.read_csv(output_dir / "triploidy_test_results.tsv", sep="\t")

    manifest = dict(zip(manifest_df["sample"], manifest_df["autosomal_baseline_cn"]))
    calls = dict(zip(results_df["sample"], results_df["triploidy_call"]))

    assert manifest == {"diploid_sample": 2, "triploid_sample": 3}
    assert calls == {"diploid_sample": "DIPLOID", "triploid_sample": "TRIPLOID"}


def test_triploidy_main_writes_diagnostics(tmp_path, monkeypatch) -> None:
    sample_ids, bin_chr, site_alt, site_total, site_pop_af, site_mask = _synthetic_triploidy_inputs()

    depth_df = pd.DataFrame(
        {
            "Chr": bin_chr,
            "Start": [0, 1000, 2000, 3000],
            "End": [1000, 2000, 3000, 4000],
            sample_ids[0]: [1.0, 1.0, 1.0, 1.0],
            sample_ids[1]: [1.0, 1.0, 1.0, 1.0],
        },
        index=[
            "chr1:0-1000",
            "chr2:1000-2000",
            "chr3:2000-3000",
            "chr4:3000-4000",
        ],
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

    output_dir = tmp_path / "triploidy_out"
    monkeypatch.setattr(
        "sys.argv",
        [
            "triploidy",
            "--input",
            str(depth_path),
            "--site-data",
            str(site_data_path),
            "--output-dir",
            str(output_dir),
            "--min-informative-bins",
            "2",
            "--min-informative-sites",
            "8",
            "--pvalue-threshold",
            "0.05",
            "--effect-size-threshold",
            "0.0",
            "--triploidy-prior",
            "0.5",
            "--diagnostics",
            "--diagnostic-sample-limit",
            "2",
            "--diagnostic-max-sites-per-sample",
            "10",
        ],
    )

    triploidy.main()

    diagnostics_dir = output_dir / "diagnostics"
    metrics_df = pd.read_csv(
        diagnostics_dir / "triploidy_diagnostic_metrics.tsv",
        sep="\t",
    ).set_index("sample")
    raw_site_df = pd.read_csv(
        diagnostics_dir / "triploidy_raw_af_sites.tsv.gz",
        sep="\t",
    )

    assert float(metrics_df.loc["diploid_sample", "fraction_near_diploid_half"]) == 1.0
    assert float(metrics_df.loc["triploid_sample", "fraction_near_triploid_thirds"]) == 1.0
    assert set(raw_site_df["sample"]) == set(sample_ids)
    assert raw_site_df.groupby("sample").size().max() <= 10
    assert (diagnostics_dir / "triploidy_diagnostic_metrics.png").exists()
    assert (diagnostics_dir / "triploidy_raw_af_profiles.png").exists()
    assert (diagnostics_dir / "triploidy_raw_af_diploid_profiles.png").exists()
