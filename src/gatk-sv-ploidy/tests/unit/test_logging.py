from __future__ import annotations

import argparse
import logging
import re

from gatk_sv_ploidy._logging import tool_logging_context


def test_tool_logging_writes_stderr_and_log_file(tmp_path, capsys) -> None:
    args = argparse.Namespace(
        input="/sensitive/project/cohort_depth.tsv",
        output_dir=str(tmp_path / "out"),
        highlight_sample="sensitive_sample_id",
        site_depth=["/sensitive/project/sample_a.sd.txt"],
        max_iter=3,
        device="cpu",
    )
    output_dir = tmp_path / "out"

    with tool_logging_context(
        tool_name="unit",
        output_dir=output_dir,
        args=args,
        random_seeds={"numpy": 42},
    ) as logger:
        logger.info("unit test payload")

    stderr_text = capsys.readouterr().err
    log_path = output_dir / "unit.log"
    log_text = log_path.read_text(encoding="utf-8")

    assert "unit test payload" in stderr_text
    assert "unit test payload" in log_text
    assert re.search(r"^\d{4}-\d{2}-\d{2}T.*Z INFO ", log_text, re.MULTILINE)
    assert "/sensitive/project" not in log_text
    assert "sensitive_sample_id" not in log_text
    assert "cohort_depth.tsv" in log_text
    assert "sample_a.sd.txt" in log_text
    assert "Random seeds" in log_text


def test_tool_logging_replaces_previous_handlers(tmp_path) -> None:
    args = argparse.Namespace(output_dir=str(tmp_path), device="cpu")

    with tool_logging_context(tool_name="first", output_dir=tmp_path, args=args):
        pass
    with tool_logging_context(tool_name="second", output_dir=tmp_path, args=args):
        logging.getLogger("gatk_sv_ploidy.second").info("second payload")

    first_log = (tmp_path / "first.log").read_text(encoding="utf-8")
    second_log = (tmp_path / "second.log").read_text(encoding="utf-8")

    assert "second payload" not in first_log
    assert "second payload" in second_log
