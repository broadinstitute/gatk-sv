from __future__ import annotations

import json
import sys

import pandas as pd
import pytest

from gatk_sv_ploidy import call, eval as eval_module, infer

pytestmark = pytest.mark.integration


def test_medium_real_fixture_runs_infer_call_and_eval(
    medium_fixture_root,
    medium_expected_truth,
    medium_expected_sex,
    tmp_path,
    monkeypatch,
) -> None:
    infer_out = tmp_path / "infer"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "gatk-sv-ploidy infer",
            "--input",
            str(medium_fixture_root / "preprocessed_depth.tsv"),
            "--site-data",
            str(medium_fixture_root / "site_data.npz"),
            "--output-dir",
            str(infer_out),
            "--device",
            "cpu",
            "--max-iter",
            "12",
            "--log-freq",
            "6",
            "--no-early-stopping",
        ],
    )
    infer.main()

    chromosome_stats = pd.read_csv(infer_out / "chromosome_stats.tsv", sep="\t")
    assert chromosome_stats["sample"].nunique() == len(medium_expected_truth)
    assert sorted(chromosome_stats["chromosome"].unique().tolist()) == [
        "chr13",
        "chr18",
        "chr21",
        "chrX",
        "chrY",
    ]
    assert chromosome_stats["mean_cn_probability"].between(0.0, 1.0).all()

    call_out = tmp_path / "call"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "gatk-sv-ploidy call",
            "--chrom-stats",
            str(infer_out / "chromosome_stats.tsv"),
            "--truth-json",
            str(medium_fixture_root / "truth.json"),
            "--output-dir",
            str(call_out),
        ],
    )
    call.main()

    pred_df = pd.read_csv(call_out / "aneuploidy_type_predictions.tsv", sep="\t")
    assert pred_df["sample"].tolist() == sorted(medium_expected_truth)
    assert set(pred_df["predicted_aneuploidy_type"]) == {"NORMAL"}
    assert set(pred_df["true_aneuploidy_type"]) == {"NORMAL"}
    assert pred_df.set_index("sample")["sex"].to_dict() == medium_expected_sex

    eval_out = tmp_path / "eval"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "gatk-sv-ploidy eval",
            "--predictions",
            str(call_out / "aneuploidy_type_predictions.tsv"),
            "--truth-json",
            str(medium_fixture_root / "truth.json"),
            "--output-dir",
            str(eval_out),
        ],
    )
    eval_module.main()

    report = (eval_out / "metrics_report.txt").read_text()
    merged = pd.read_csv(eval_out / "predictions_with_truth.tsv", sep="\t")
    assert "Overall Accuracy" in report
    assert merged["sample"].tolist() == pred_df["sample"].tolist()
    assert merged.set_index("sample")["sex"].to_dict() == medium_expected_sex
    assert json.loads((medium_fixture_root / "manifest.json").read_text())["n_bins"] == 69
