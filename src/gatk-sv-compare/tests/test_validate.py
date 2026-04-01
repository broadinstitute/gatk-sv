from __future__ import annotations

from gatk_sv_compare.config import ValidateConfig
from gatk_sv_compare.validate import render_summary, validate_vcf
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
