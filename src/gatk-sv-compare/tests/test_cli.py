from __future__ import annotations

from gatk_sv_compare.cli import AnalysisConfig, build_parser, main


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

    out = capsys.readouterr().out
    assert exit_code == 0
    assert "annotated_a=" in out
    assert captured["config"].contigs == ["chr1", "chr2"]
    assert captured["config"].contig_lengths == {"chr1": 1000, "chr2": 2000}
    assert captured["config"].n_workers == 2

