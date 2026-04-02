from __future__ import annotations

import os
from types import SimpleNamespace

import pandas as pd
import pysam
import pytest

from gatk_sv_compare.cli import AnalysisConfig, build_parser, main
from gatk_sv_compare.modules.base import AnalysisModule


def test_build_parser_exposes_phase1_commands() -> None:
    parser = build_parser()
    commands = parser._subparsers._group_actions[0].choices
    assert {"validate", "preprocess", "analyze", "run"}.issubset(commands)


def test_cli_validate_returns_zero_for_clean_vcf(make_vcf, capsys) -> None:
    vcf_path = make_vcf(
        file_name="ok.vcf",
        records=[
            "chr1\t100\tvar1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=25\tGT:GQ:ECN\t0/1:60:2\t0/0:50:2",
        ],
    )
    exit_code = main(["validate", "--vcf", str(vcf_path)])
    out = capsys.readouterr().out
    assert exit_code == 0
    assert "VCF:" in out


def test_cli_validate_fix_writes_output(make_vcf, tmp_path, capsys) -> None:
    vcf_path = make_vcf(
        file_name="fixable.vcf",
        records=[
            "chr1\t100\tvar1\tN\t<DEL>\t.\tPASS\tSVLEN=25\tGT:GQ:ECN\t0/1:60:2\t0/0:50:2",
        ],
    )
    requested_out_path = tmp_path / "fixed.vcf"
    out_path = tmp_path / "fixed.vcf.gz"

    exit_code = main(["validate", "--vcf", str(vcf_path), "--fix", "--out", str(requested_out_path)])

    out = capsys.readouterr().out
    assert exit_code == 0
    assert f"Wrote fixed VCF: {out_path}" in out
    assert out_path.exists()
    assert (tmp_path / "fixed.vcf.gz.tbi").exists()


def test_cli_validate_fix_writes_bgzip_output_and_index(make_vcf, tmp_path, capsys) -> None:
    vcf_path = make_vcf(
        file_name="fixable_bgzip.vcf",
        records=[
            "chr1\t100\tvar1\tN\t<INS:ME:ALU>\t.\tPASS\tSVLEN=25\tGT:GQ:ECN\t0/1:60:2\t0/0:50:2",
        ],
    )
    out_path = tmp_path / "fixed.vcf.gz"

    exit_code = main(["validate", "--vcf", str(vcf_path), "--fix", "--out", str(out_path)])

    out = capsys.readouterr().out
    assert exit_code == 0
    assert "Wrote fixed VCF:" in out
    assert out_path.exists()
    assert (tmp_path / "fixed.vcf.gz.tbi").exists()


def test_cli_validate_fix_defaults_output_path_for_vcf(make_vcf, capsys) -> None:
    vcf_path = make_vcf(
        file_name="autofix.vcf",
        records=[
            "chr1\t100\tvar1\tN\t<DEL>\t.\tPASS\tSVLEN=25\tGT:GQ:ECN\t0/1:60:2\t0/0:50:2",
        ],
    )
    expected_out = vcf_path.with_name("autofix.fixed.vcf.gz")

    exit_code = main(["validate", "--vcf", str(vcf_path), "--fix"])

    out = capsys.readouterr().out
    assert exit_code == 0
    assert f"Wrote fixed VCF: {expected_out}" in out
    assert expected_out.exists()
    assert expected_out.with_name(expected_out.name + ".tbi").exists()


def test_cli_validate_fix_defaults_output_path_for_vcfgz(make_vcf, capsys) -> None:
    vcf_path = make_vcf(
        file_name="autofix_source.vcf",
        records=[
            "chr1\t100\tvar1\tN\tN]chr2:200]\t.\tPASS\t.\tGT:GQ:ECN\t0/1:60:2\t0/0:50:2",
        ],
    )
    gz_path = vcf_path.with_suffix(".vcf.gz")
    pysam.tabix_compress(str(vcf_path), str(gz_path), force=True)
    vcf_path.unlink()
    expected_out = gz_path.with_name("autofix_source.fixed.vcf.gz")

    exit_code = main(["validate", "--vcf", str(gz_path), "--fix"])

    out = capsys.readouterr().out
    assert exit_code == 0
    assert f"Wrote fixed VCF: {expected_out}" in out
    assert expected_out.exists()
    assert expected_out.with_name(expected_out.name + ".tbi").exists()


def test_cli_validate_fix_reports_blocking_check(make_vcf, tmp_path, capsys) -> None:
    vcf_path = make_vcf(
        file_name="unfixable_cli.vcf",
        records=[
            "chr1\t100\tvar1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=25\tECN\t2\t2",
        ],
    )

    exit_code = main(["validate", "--vcf", str(vcf_path), "--fix"])

    out = capsys.readouterr().out
    assert exit_code == 1
    assert "Fix mode aborted:" in out
    assert "MISSING_GT" in out
    assert "Examples:" in out
    assert "[var1]" in out


def test_cli_validate_fix_suggests_ploidy_table_for_missing_ecn(make_vcf, tmp_path, capsys) -> None:
    vcf_path = make_vcf(
        file_name="missing_ecn_hint.vcf",
        records=[
            "chr1\t100\tvar1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=25\tGT:GQ\t0/1:60\t0/0:50",
        ],
    )

    exit_code = main(["validate", "--vcf", str(vcf_path), "--fix"])

    out = capsys.readouterr().out
    assert exit_code == 1
    assert "MISSING_ECN" in out
    assert "--ploidy-table" in out


def test_cli_validate_fix_repairs_missing_ecn_with_ploidy_table(make_vcf, tmp_path, capsys) -> None:
    vcf_path = make_vcf(
        file_name="missing_ecn_cli.vcf",
        records=[
            "chr1\t100\tvar1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=25\tGT:GQ\t0/1:60\t0/0:50",
        ],
    )
    ploidy_table = tmp_path / "ploidy.tsv"
    ploidy_table.write_text("SAMPLE\tchr1\nS1\t2\nS2\t2\n")
    out_path = tmp_path / "missing_ecn_fixed.vcf.gz"

    exit_code = main(["validate", "--vcf", str(vcf_path), "--fix", "--ploidy-table", str(ploidy_table), "--out", str(out_path)])

    out = capsys.readouterr().out
    assert exit_code == 0
    assert f"Wrote fixed VCF: {out_path}" in out
    with pysam.VariantFile(str(out_path)) as vcf:
        record = next(iter(vcf))
        assert record.samples["S1"]["ECN"] == 2
        assert record.samples["S2"]["ECN"] == 2


def test_cli_preprocess_builds_config(tmp_path, monkeypatch, capsys) -> None:
    captured = {}

    def fake_run_preprocess(config: AnalysisConfig):
        captured["config"] = config
        return tmp_path / "annotated_a.vcf.gz", tmp_path / "annotated_b.vcf.gz"

    contig_list = tmp_path / "contigs.list"
    contig_list.write_text("chr1\nchr2\n")
    reference_dict = tmp_path / "ref.dict"
    reference_dict.write_text("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:2000\n")

    monkeypatch.setattr("gatk_sv_compare.cli.run_preprocess", fake_run_preprocess)

    exit_code = main([
        "preprocess",
        "--vcf-a",
        str(tmp_path / "a.vcf.gz"),
        "--vcf-b",
        str(tmp_path / "b.vcf.gz"),
        "--reference-dict",
        str(reference_dict),
        "--contig-list",
        str(contig_list),
        "--output-dir",
        str(tmp_path / "out"),
        "--num-workers",
        "2",
    ])

    captured_io = capsys.readouterr()
    out = captured_io.out
    err = captured_io.err
    assert exit_code == 0
    assert "annotated_a=" in out
    assert "Loaded preprocess inputs: 2 contigs" in err
    assert "using 2 worker(s)" in err
    assert captured["config"].contigs == ["chr1", "chr2"]
    assert captured["config"].contig_lengths == {"chr1": 1000, "chr2": 2000}
    assert captured["config"].n_workers == 2


def test_cli_preprocess_restricts_to_single_contig(tmp_path, monkeypatch, capsys) -> None:
    captured = {}

    def fake_run_preprocess(config: AnalysisConfig):
        captured["config"] = config
        return tmp_path / "annotated_a.vcf.gz", tmp_path / "annotated_b.vcf.gz"

    contig_list = tmp_path / "contigs.list"
    contig_list.write_text("chr1\nchr22\n")
    reference_dict = tmp_path / "ref.dict"
    reference_dict.write_text("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr22\tLN:2000\n")

    monkeypatch.setattr("gatk_sv_compare.cli.run_preprocess", fake_run_preprocess)

    exit_code = main([
        "preprocess",
        "--vcf-a",
        str(tmp_path / "a.vcf.gz"),
        "--vcf-b",
        str(tmp_path / "b.vcf.gz"),
        "--reference-dict",
        str(reference_dict),
        "--contig-list",
        str(contig_list),
        "--contig",
        "chr22",
        "--output-dir",
        str(tmp_path / "out"),
    ])

    captured_io = capsys.readouterr()
    err = captured_io.err
    assert exit_code == 0
    assert "Loaded preprocess inputs: 1 contigs" in err
    assert captured["config"].contigs == ["chr22"]
    assert captured["config"].contig_lengths == {"chr22": 2000}


def test_cli_preprocess_rejects_unknown_requested_contig(tmp_path) -> None:
    contig_list = tmp_path / "contigs.list"
    contig_list.write_text("chr1\nchr2\n")
    reference_dict = tmp_path / "ref.dict"
    reference_dict.write_text("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:2000\n")

    with pytest.raises(ValueError, match="Requested contig chr22"):
        main([
            "preprocess",
            "--vcf-a",
            str(tmp_path / "a.vcf.gz"),
            "--vcf-b",
            str(tmp_path / "b.vcf.gz"),
            "--reference-dict",
            str(reference_dict),
            "--contig-list",
            str(contig_list),
            "--contig",
            "chr22",
            "--output-dir",
            str(tmp_path / "out"),
        ])


def test_cli_preprocess_defaults_to_auto_parallel_workers(tmp_path, monkeypatch, capsys) -> None:
    captured = {}

    def fake_run_preprocess(config: AnalysisConfig):
        captured["config"] = config
        return tmp_path / "annotated_a.vcf.gz", tmp_path / "annotated_b.vcf.gz"

    contig_list = tmp_path / "contigs.list"
    contig_list.write_text("chr1\nchr2\nchr3\nchr4\nchr5\n")
    reference_dict = tmp_path / "ref.dict"
    reference_dict.write_text(
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:2000\n@SQ\tSN:chr3\tLN:3000\n@SQ\tSN:chr4\tLN:4000\n@SQ\tSN:chr5\tLN:5000\n"
    )

    monkeypatch.setattr("gatk_sv_compare.cli.run_preprocess", fake_run_preprocess)

    exit_code = main([
        "preprocess",
        "--vcf-a",
        str(tmp_path / "a.vcf.gz"),
        "--vcf-b",
        str(tmp_path / "b.vcf.gz"),
        "--reference-dict",
        str(reference_dict),
        "--contig-list",
        str(contig_list),
        "--output-dir",
        str(tmp_path / "out"),
    ])

    captured_io = capsys.readouterr()
    err = captured_io.err
    expected_workers = min(5, os.cpu_count() or 1, 4)

    assert exit_code == 0
    assert f"using {expected_workers} worker(s)" in err
    assert captured["config"].n_workers == expected_workers


class DemoModule(AnalysisModule):
    ran = False

    @property
    def name(self) -> str:
        return "demo"

    def run(self, data, config: AnalysisConfig) -> None:
        del data, config
        DemoModule.ran = True


class SharedDemoModule(AnalysisModule):
    ran = False

    @property
    def name(self) -> str:
        return "shared_demo"

    @property
    def requires_shared_samples(self) -> bool:
        return True

    def run(self, data, config: AnalysisConfig) -> None:
        del data, config
        SharedDemoModule.ran = True


class PedDemoModule(AnalysisModule):
    ran = False

    @property
    def name(self) -> str:
        return "ped_demo"

    @property
    def requires_ped_file(self) -> bool:
        return True

    def run(self, data, config: AnalysisConfig) -> None:
        del data, config
        PedDemoModule.ran = True


def test_cli_analyze_runs_selected_modules_and_skips_unmet_requirements(make_vcf, tmp_path, monkeypatch, capsys) -> None:
    DemoModule.ran = False
    SharedDemoModule.ran = False
    PedDemoModule.ran = False
    captured = {}

    vcf_a = make_vcf(file_name="a.vcf", records=["chr1\t100\ta1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=50\tGT:GQ:ECN\t0/1:60:2\t0/0:50:2"])
    vcf_b = make_vcf(file_name="b.vcf", records=["chr1\t110\tb1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=55\tGT:GQ:ECN\t0/1:60:2\t0/0:50:2"])

    def fake_aggregate(config: AnalysisConfig):
        captured["config"] = config
        return SimpleNamespace(
            sites_a=pd.DataFrame([{"variant_id": "a1"}]),
            sites_b=pd.DataFrame([{"variant_id": "b1"}]),
            matched_pairs=pd.DataFrame(),
            shared_samples=[],
            label_a=config.vcf_a_label,
            label_b=config.vcf_b_label,
        )

    monkeypatch.setattr("gatk_sv_compare.cli.ALL_MODULES", [DemoModule, SharedDemoModule, PedDemoModule])
    monkeypatch.setattr("gatk_sv_compare.cli.aggregate", fake_aggregate)

    exit_code = main([
        "analyze",
        "--vcf-a",
        str(vcf_a),
        "--vcf-b",
        str(vcf_b),
        "--label-a",
        "CallsetA",
        "--label-b",
        "CallsetB",
        "--modules",
        "demo,shared_demo,ped_demo",
        "--output-dir",
        str(tmp_path / "out"),
        "--pass-only",
        "--context-overlap",
        "0.75",
        "--enable-site-match-table",
        "--per-sample-counts-table",
    ])

    captured_io = capsys.readouterr()
    out = captured_io.out
    err = captured_io.err
    assert exit_code == 0
    assert DemoModule.ran is True
    assert SharedDemoModule.ran is False
    assert PedDemoModule.ran is False
    assert "modules_ran=demo" in out
    assert "modules_skipped=shared_demo,ped_demo" in out
    assert "Skipping module shared_demo: no shared samples" in err
    assert "Skipping module ped_demo: no --ped provided" in err
    assert captured["config"].modules == ["demo", "shared_demo", "ped_demo"]
    assert captured["config"].contigs == ["chr1"]
    assert captured["config"].pass_only is True
    assert captured["config"].context_overlap == pytest.approx(0.75)
    assert captured["config"].enable_site_match_table is True
    assert captured["config"].per_sample_counts_table is True


def test_cli_analyze_rejects_unknown_modules(make_vcf, tmp_path) -> None:
    vcf_a = make_vcf(file_name="a.vcf", records=["chr1\t100\ta1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=50\tGT:GQ:ECN\t0/1:60:2\t0/0:50:2"])
    vcf_b = make_vcf(file_name="b.vcf", records=["chr1\t110\tb1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=55\tGT:GQ:ECN\t0/1:60:2\t0/0:50:2"])

    with pytest.raises(ValueError, match="Unknown module"):
        main([
            "analyze",
            "--vcf-a",
            str(vcf_a),
            "--vcf-b",
            str(vcf_b),
            "--modules",
            "does_not_exist",
            "--output-dir",
            str(tmp_path / "out"),
        ])


def test_cli_run_executes_preprocess_then_analyze(tmp_path, monkeypatch, capsys) -> None:
    DemoModule.ran = False
    captured = {}

    def fake_run_preprocess(config: AnalysisConfig):
        captured["preprocess_config"] = config
        return tmp_path / "annotated_a.vcf.gz", tmp_path / "annotated_b.vcf.gz"

    def fake_aggregate(config: AnalysisConfig):
        captured["analyze_config"] = config
        return SimpleNamespace(
            sites_a=pd.DataFrame([{"variant_id": "a1"}]),
            sites_b=pd.DataFrame([{"variant_id": "b1"}]),
            matched_pairs=pd.DataFrame(),
            shared_samples=["S1"],
            label_a=config.vcf_a_label,
            label_b=config.vcf_b_label,
        )

    contig_list = tmp_path / "contigs.list"
    contig_list.write_text("chr1\nchr2\n")
    reference_dict = tmp_path / "ref.dict"
    reference_dict.write_text("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:2000\n")

    monkeypatch.setattr("gatk_sv_compare.cli.ALL_MODULES", [DemoModule])
    monkeypatch.setattr("gatk_sv_compare.cli.run_preprocess", fake_run_preprocess)
    monkeypatch.setattr("gatk_sv_compare.cli.aggregate", fake_aggregate)

    exit_code = main([
        "run",
        "--vcf-a",
        str(tmp_path / "input_a.vcf.gz"),
        "--vcf-b",
        str(tmp_path / "input_b.vcf.gz"),
        "--label-a",
        "CallsetA",
        "--label-b",
        "CallsetB",
        "--reference-dict",
        str(reference_dict),
        "--contig-list",
        str(contig_list),
        "--output-dir",
        str(tmp_path / "out"),
        "--context-overlap",
        "0.25",
        "--modules",
        "demo",
    ])

    captured_io = capsys.readouterr()
    out = captured_io.out
    err = captured_io.err
    assert exit_code == 0
    assert DemoModule.ran is True
    assert "annotated_a=" in out
    assert "modules_ran=demo" in out
    assert "Starting end-to-end run with 2 contig(s)" in err
    assert captured["preprocess_config"].contigs == ["chr1", "chr2"]
    assert captured["analyze_config"].vcf_a_path == tmp_path / "annotated_a.vcf.gz"
    assert captured["analyze_config"].vcf_b_path == tmp_path / "annotated_b.vcf.gz"
    assert captured["analyze_config"].contig_lengths == {"chr1": 1000, "chr2": 2000}
    assert captured["analyze_config"].context_overlap == pytest.approx(0.25)

