from __future__ import annotations

import pysam

from gatk_sv_compare.config import ValidateConfig
from gatk_sv_compare.validate import render_summary, validate_and_fix, validate_vcf
from gatk_sv_compare.vcf_format import PipelineStage


def test_validate_vcf_returns_summary(make_vcf) -> None:
    vcf_path = make_vcf(
        file_name="cleaned.vcf",
        records=[
            "chr1\t100\tvar1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=100\tGT:GQ:ECN\t0/1:60:2\t0/0:50:2",
        ],
    )
    summary = validate_vcf(ValidateConfig(vcf_path=vcf_path))
    assert summary.stage is PipelineStage.CLEANED
    assert summary.record_count == 1
    assert summary.has_errors is False
    assert "Stage: cleaned" in render_summary(summary)


def test_validate_vcf_reports_errors(make_vcf) -> None:
    vcf_path = make_vcf(
        file_name="bad.vcf",
        records=[
            "chr1\t100\tvar1\tN\tA]chr1:200]\t.\t.\tSVLEN=-10\tGQ:ECN\t60:2\t50:2",
        ],
    )
    summary = validate_vcf(ValidateConfig(vcf_path=vcf_path))
    assert summary.has_errors is True
    check_ids = {issue.check_id for issue in summary.issues}
    assert "MISSING_SVTYPE" in check_ids
    assert "BREAKEND_NOTATION" not in check_ids


def test_render_summary_deduplicates_detail_lines_and_reports_counts(make_vcf) -> None:
    vcf_path = make_vcf(
        file_name="repeated_issues.vcf",
        records=[
            "chr1\t100\tvar1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=250000\tGT:GQ:ECN\t0/1:60:2\t0/0:50:2",
            "chr1\t200\tvar2\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=250001\tGT:GQ:ECN\t0/1:60:2\t0/0:50:2",
        ],
    )

    summary = validate_vcf(ValidateConfig(vcf_path=vcf_path))
    rendered = render_summary(summary)

    assert rendered.count("IMPLAUSIBLE_SVLEN") == 2
    assert rendered.count("WARN IMPLAUSIBLE_SVLEN") == 2
    assert "- WARN IMPLAUSIBLE_SVLEN [var1]: SVLEN exceeds 20% of chromosome length and may be artifactual" in rendered
    assert "[var2]: SVLEN exceeds 20% of chromosome length and may be artifactual" not in rendered
    assert "Issue counts:" in rendered
    assert "- WARN IMPLAUSIBLE_SVLEN: 2" in rendered


def test_validate_and_fix_rewrites_fixable_issues(make_vcf, tmp_path) -> None:
    vcf_path = make_vcf(
        file_name="fixable.vcf",
        extra_header_lines=["##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Partner chromosome\">"],
        records=[
            "chr1\t100\tvar1\tN\t<INS>\t.\t.\tSVTYPE=INS;SVLEN=-25;END=150;CHR2=chr2\tGT:GQ:ECN\t0/1:120:2\t0/0:50:2",
        ],
    )
    requested_out_path = tmp_path / "fixed.vcf"
    out_path = tmp_path / "fixed.vcf.gz"

    result = validate_and_fix(vcf_path, requested_out_path)

    assert result.wrote_output is True
    assert result.fixed_summary is not None
    assert result.out_path == out_path
    fixed_check_ids = {issue.check_id for issue in result.fixed_summary.issues}
    assert "SVLEN_SIGN" not in fixed_check_ids
    assert "INS_END_MISMATCH" not in fixed_check_ids
    assert "EMPTY_FILTER" not in fixed_check_ids
    assert "GQ_RANGE" not in fixed_check_ids
    assert "CHR2_ON_NON_BND" not in fixed_check_ids
    assert out_path.exists()
    assert (tmp_path / "fixed.vcf.gz.tbi").exists()
    with pysam.VariantFile(str(out_path)) as vcf:
        record = next(iter(vcf))
        assert set(record.filter.keys()) == {"PASS"}
        assert int(record.info["SVLEN"]) == 25
        assert record.stop == 100
        assert "CHR2" not in record.info
        assert record.samples[0]["GQ"] == 12


def test_validate_and_fix_blocks_unfixable_errors(make_vcf, tmp_path) -> None:
    vcf_path = make_vcf(
        file_name="unfixable.vcf",
        records=[
            "chr1\t100\tvar1\tN\t<DEL>\t.\tPASS\tSVLEN=25\tGT:GQ:ECN\t0/1:60:2\t0/0:50:2",
        ],
    )
    out_path = tmp_path / "blocked.vcf"

    result = validate_and_fix(vcf_path, out_path)

    assert result.wrote_output is False
    assert result.fixed_summary is None
    assert out_path.exists() is False
    assert (tmp_path / "blocked.vcf.gz").exists() is False


def test_validate_and_fix_writes_bgzip_and_tabix_for_gz_output(make_vcf, tmp_path) -> None:
    vcf_path = make_vcf(
        file_name="fixable_for_bgzip.vcf",
        records=[
            "chr1\t100\tvar1\tN\t<INS>\t.\t.\tSVTYPE=INS;SVLEN=-25;END=150\tGT:GQ:ECN\t0/1:120:2\t0/0:50:2",
        ],
    )
    out_path = tmp_path / "fixed.vcf.gz"

    result = validate_and_fix(vcf_path, out_path)

    assert result.wrote_output is True
    assert out_path.exists()
    assert (tmp_path / "fixed.vcf.gz.tbi").exists()
    with pysam.VariantFile(str(out_path)) as vcf:
        record = next(iter(vcf))
        assert int(record.info["SVLEN"]) == 25
        assert record.samples[0]["GQ"] == 12
