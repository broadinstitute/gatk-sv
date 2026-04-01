from __future__ import annotations

import pytest

from gatk_sv_compare.cli import build_parser, main


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


def test_cli_preprocess_is_stubbed() -> None:
    with pytest.raises(NotImplementedError):
        main([
            "preprocess",
            "--vcf-a",
            "a.vcf.gz",
            "--vcf-b",
            "b.vcf.gz",
            "--reference-dict",
            "ref.dict",
            "--contig-list",
            "contigs.list",
            "--output-dir",
            "out",
        ])
