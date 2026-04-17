from __future__ import annotations

import sys

import pandas as pd
import pytest

from gatk_sv_ploidy import call, infer, preprocess

pytestmark = [pytest.mark.integration, pytest.mark.large]


def test_large_real_fixture_runs_preprocess_and_infer(
    large_fixture_root,
    large_expected_truth,
    large_expected_sex,
    tmp_path,
    monkeypatch,
) -> None:
    preprocess_out = tmp_path / "preprocess"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "gatk-sv-ploidy preprocess",
            "--input",
            str(large_fixture_root / "raw_depth.tsv.gz"),
            "--output-dir",
            str(preprocess_out),
            "--viable-only",
            "--skip-bin-filter",
        ],
    )
    preprocess.main()

    preprocessed = pd.read_csv(preprocess_out / "preprocessed_depth.tsv", sep="\t", index_col=0)
    assert sorted(preprocessed["Chr"].unique().tolist()) == [
        "chr13",
        "chr18",
        "chr21",
        "chrX",
        "chrY",
    ]
    assert preprocessed.shape[0] >= 300
    assert len([c for c in preprocessed.columns if c not in {"Chr", "Start", "End", "source_file"}]) == len(large_expected_truth)

    infer_out = tmp_path / "infer"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "gatk-sv-ploidy infer",
            "--input",
            str(preprocess_out / "preprocessed_depth.tsv"),
            "--output-dir",
            str(infer_out),
            "--device",
            "cpu",
            "--max-iter",
            "15",
            "--log-freq",
            "5",
            "--no-early-stopping",
        ],
    )
    infer.main()

    chromosome_stats = pd.read_csv(infer_out / "chromosome_stats.tsv", sep="\t")
    assert chromosome_stats["sample"].nunique() == len(large_expected_truth)
    assert chromosome_stats["chromosome"].nunique() == 5
    assert chromosome_stats["mean_cn_probability"].between(0.0, 1.0).all()
    assert (infer_out / "inference_artifacts.npz").exists()

    call_out = tmp_path / "call"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "gatk-sv-ploidy call",
            "--chrom-stats",
            str(infer_out / "chromosome_stats.tsv"),
            "--truth-json",
            str(large_fixture_root / "truth.json"),
            "--output-dir",
            str(call_out),
        ],
    )
    call.main()

    pred_df = pd.read_csv(call_out / "aneuploidy_type_predictions.tsv", sep="\t")
    assert pred_df.set_index("sample")["predicted_aneuploidy_type"].to_dict() == large_expected_truth
    assert pred_df.set_index("sample")["sex"].to_dict() == large_expected_sex
