from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from types import ModuleType

import numpy as np
import pandas as pd
import pytest

from gatk_sv_ploidy import call, eval as eval_module, infer

pytestmark = pytest.mark.integration

_METADATA_COLUMNS = {"Bin", "Chr", "Start", "End", "BinLengthBp", "source_file"}
_SYNTHETIC_COUNT_SCALE = 250000.0
_SYNTHETIC_SEED = 17


def _load_synthetic_builder_module() -> ModuleType:
    builder_path = (
        Path(__file__).resolve().parents[4] /
        "scratch" /
        "synthetic_ploidy" /
        "build_synthetic_cohort.py"
    )
    spec = importlib.util.spec_from_file_location(
        "synthetic_ploidy_builder",
        builder_path,
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load synthetic builder from {builder_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _write_synthetic_cohort_from_medium_fixture(
    medium_fixture_root: Path,
    output_dir: Path,
) -> dict[str, str]:
    builder = _load_synthetic_builder_module()
    reference_df = pd.read_csv(
        medium_fixture_root / "preprocessed_depth.tsv",
        sep="\t",
        index_col=0,
    )
    sample_columns = [
        column for column in reference_df.columns if column not in _METADATA_COLUMNS
    ]
    reference_df["BinLengthBp"] = reference_df["End"] - reference_df["Start"]
    scaled_counts = np.rint(
        reference_df[sample_columns].to_numpy(dtype=np.float64) * _SYNTHETIC_COUNT_SCALE
    ).astype(np.int64)
    reference_df.loc[:, sample_columns] = scaled_counts

    template = builder.get_neutral_template(reference_df)
    sample_specs = builder.build_sample_specs()
    synthetic_depth = builder.synthesize_depth_matrix(
        template,
        sample_specs,
        np.random.default_rng(_SYNTHETIC_SEED),
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    synthetic_depth.to_csv(output_dir / "preprocessed_depth.tsv", sep="\t", index=True)
    (output_dir / "observation_type.txt").write_text("raw\n", encoding="ascii")

    truth, sex_truth = builder.build_truth_maps(sample_specs)
    (output_dir / "truth.json").write_text(
        json.dumps(truth, indent=2, sort_keys=True),
        encoding="ascii",
    )
    (output_dir / "sex_truth.json").write_text(
        json.dumps(sex_truth, indent=2, sort_keys=True),
        encoding="ascii",
    )
    return truth


def test_synthetic_fixture_absorbs_genomewide_depth_shifts_into_sample_depth(
    medium_fixture_root: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cohort_dir = tmp_path / "synthetic"
    truth = _write_synthetic_cohort_from_medium_fixture(medium_fixture_root, cohort_dir)

    infer_out = tmp_path / "infer"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "gatk-sv-ploidy infer",
            "--input",
            str(cohort_dir / "preprocessed_depth.tsv"),
            "--output-dir",
            str(infer_out),
            "--device",
            "cpu",
            "--depth-space",
            "raw",
            "--obs-likelihood",
            "negative_binomial",
            "--max-iter",
            "150",
            "--log-freq",
            "50",
            "--cn-inference-draws",
            "30",
            "--no-early-stopping",
        ],
    )
    infer.main()

    chromosome_stats = pd.read_csv(infer_out / "chromosome_stats.tsv", sep="\t")
    chromosome_indexed = chromosome_stats.set_index(["sample", "chromosome"])

    assert not (infer_out / "sample_baseline_autosome_ploidy.tsv").exists()

    for chromosome in ("chr13", "chr18", "chr21"):
        assert int(chromosome_indexed.loc[("synthetic_normal_01", chromosome), "copy_number"]) == 2
        assert int(chromosome_indexed.loc[("synthetic_triploid", chromosome), "copy_number"]) == 2
        assert int(chromosome_indexed.loc[("synthetic_tetraploid", chromosome), "copy_number"]) == 2

    normal_depth = float(
        chromosome_indexed.loc[("synthetic_normal_01", "chr13"), "sample_depth_map"]
    )
    triploid_depth = float(
        chromosome_indexed.loc[("synthetic_triploid", "chr13"), "sample_depth_map"]
    )
    tetraploid_depth = float(
        chromosome_indexed.loc[("synthetic_tetraploid", "chr13"), "sample_depth_map"]
    )
    assert triploid_depth > normal_depth * 1.3
    assert tetraploid_depth > normal_depth * 1.8

    call_out = tmp_path / "call"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "gatk-sv-ploidy call",
            "--chrom-stats",
            str(infer_out / "chromosome_stats.tsv"),
            "--truth-json",
            str(cohort_dir / "truth.json"),
            "--output-dir",
            str(call_out),
        ],
    )
    call.main()

    pred_df = pd.read_csv(call_out / "aneuploidy_type_predictions.tsv", sep="\t")
    pred_by_sample = pred_df.set_index("sample")
    assert pred_by_sample.loc["synthetic_normal_01", "predicted_aneuploidy_type"] == "NORMAL"
    assert pred_by_sample.loc["synthetic_trisomy21", "predicted_aneuploidy_type"] == "TRISOMY_21"
    assert pred_by_sample.loc["synthetic_xxy", "predicted_aneuploidy_type"] == "KLINEFELTER"
    assert pred_by_sample.loc["synthetic_triploid", "predicted_aneuploidy_type"] == "NORMAL"
    assert pred_by_sample.loc["synthetic_tetraploid", "predicted_aneuploidy_type"] == "NORMAL"
    assert pred_by_sample.loc["synthetic_triploid", "sample_depth_ratio"] > 1.3
    assert pred_by_sample.loc["synthetic_tetraploid", "sample_depth_ratio"] > 1.8
    assert pred_by_sample.loc["synthetic_tetraploid", "global_cn_scale_factor"] == pytest.approx(1.0)

    eval_out = tmp_path / "eval"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "gatk-sv-ploidy eval",
            "--predictions",
            str(call_out / "aneuploidy_type_predictions.tsv"),
            "--truth-json",
            str(cohort_dir / "truth.json"),
            "--output-dir",
            str(eval_out),
        ],
    )
    eval_module.main()

    report = (eval_out / "metrics_report.txt").read_text(encoding="ascii")
    merged = pd.read_csv(eval_out / "predictions_with_truth.tsv", sep="\t")
    assert "Overall Accuracy" in report
    assert merged["sample"].tolist() == pred_df["sample"].tolist()
    assert merged.set_index("sample").loc["synthetic_tetraploid", "predicted_aneuploidy_type"] == "NORMAL"
    assert pred_df.set_index("sample")["true_aneuploidy_type"].to_dict() == truth
