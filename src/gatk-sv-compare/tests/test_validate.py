from __future__ import annotations

import pysam

from gatk_sv_compare.config import ValidateConfig
from gatk_sv_compare.validate import render_fix_result, render_summary, validate_and_fix, validate_vcf
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


def test_validate_vcf_reports_missing_ecn_as_error(make_vcf) -> None:
    vcf_path = make_vcf(
        file_name="missing_ecn_error.vcf",
        records=[
            "chr1\t100\tvar1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=25\tGT:GQ\t0/1:60\t0/0:50",
        ],
    )

    summary = validate_vcf(ValidateConfig(vcf_path=vcf_path))

    check_ids = {issue.check_id for issue in summary.issues if issue.severity == "ERROR"}
    assert "MISSING_ECN" in check_ids


def test_render_summary_deduplicates_detail_lines_and_reports_counts(make_vcf) -> None:
    vcf_path = make_vcf(
        file_name="repeated_issues.vcf",
        records=[
            "chr1\t100\tvar1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-25\tGT:GQ:ECN\t0/1:60:2\t0/0:50:2",
            "chr1\t200\tvar2\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-30\tGT:GQ:ECN\t0/1:60:2\t0/0:50:2",
        ],
    )

    summary = validate_vcf(ValidateConfig(vcf_path=vcf_path))
    rendered = render_summary(summary)

    assert rendered.count("SVLEN_SIGN") == 2
    assert rendered.count("WARN SVLEN_SIGN") == 2
    assert "- WARN SVLEN_SIGN [var1]: SVLEN is negative" in rendered
    assert "[var2]: SVLEN is negative" not in rendered
    assert "Issue counts:" in rendered
    assert "- WARN SVLEN_SIGN: 2" in rendered


def test_validate_and_fix_rewrites_fixable_issues(make_vcf, tmp_path) -> None:
    vcf_path = make_vcf(
        file_name="fixable.vcf",
        records=[
            "chr1\t100\tvar1\tN\t<DEL>\t.\tPASS\tSVLEN=25\tGT:GQ:ECN\t0/1:60:2\t0/0:50:2",
        ],
    )
    requested_out_path = tmp_path / "fixed.vcf"
    out_path = tmp_path / "fixed.vcf.gz"

    result = validate_and_fix(vcf_path, requested_out_path)

    assert result.wrote_output is True
    assert result.fixed_summary is not None
    assert result.out_path == out_path
    fixed_check_ids = {issue.check_id for issue in result.fixed_summary.issues}
    assert "MISSING_SVTYPE" not in fixed_check_ids
    assert out_path.exists()
    assert (tmp_path / "fixed.vcf.gz.tbi").exists()
    with pysam.VariantFile(str(out_path)) as vcf:
        record = next(iter(vcf))
        assert str(record.info["SVTYPE"]) == "DEL"
        assert set(record.filter.keys()) == {"PASS"}
        assert int(record.info["SVLEN"]) == 25


def test_validate_and_fix_canonicalizes_breakend_non_lossily(make_vcf, tmp_path) -> None:
    vcf_path = make_vcf(
        file_name="breakend_fixable.vcf",
        records=[
            "chr1\t100\tvar1\tN\tN]chr2:200]\t.\tPASS\t.\tGT:GQ:ECN\t0/1:60:2\t0/0:50:2",
        ],
    )
    out_path = tmp_path / "breakend_fixed.vcf.gz"

    result = validate_and_fix(vcf_path, out_path)

    assert result.wrote_output is True
    assert result.fixed_summary is not None
    assert result.fixed_summary.has_errors is False
    with pysam.VariantFile(str(out_path)) as vcf:
        record = next(iter(vcf))
        assert record.alts == ("<CTX>",)
        assert str(record.info["SVTYPE"]) == "CTX"
        assert str(record.info["CHR2"]) == "chr2"
        assert int(record.info["END2"]) == 200
        assert int(record.stop) == 100
        assert str(record.info["ORIGINAL_ALT"]) == "N]chr2:200]"


def test_validate_and_fix_blocks_unfixable_errors(make_vcf, tmp_path) -> None:
    vcf_path = make_vcf(
        file_name="unfixable.vcf",
        records=[
            "chr1\t100\tvar1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=25\tECN\t2\t2",
        ],
    )
    out_path = tmp_path / "blocked.vcf"

    result = validate_and_fix(vcf_path, out_path)

    assert result.wrote_output is False
    assert result.fixed_summary is None
    assert out_path.exists() is False
    assert (tmp_path / "blocked.vcf.gz").exists() is False
    assert "MISSING_GT" in render_fix_result(result)


def test_validate_and_fix_repairs_missing_ecn_with_ploidy_table(make_vcf, tmp_path) -> None:
    vcf_path = make_vcf(
        file_name="missing_ecn_fixable.vcf",
        records=[
            "chr1\t100\tvar1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=25\tGT:GQ\t0/1:60\t0/0:50",
        ],
    )
    ploidy_table = tmp_path / "ploidy.tsv"
    ploidy_table.write_text("SAMPLE\tchr1\nS1\t2\nS2\t1\n")
    out_path = tmp_path / "fixed_ecn.vcf.gz"

    result = validate_and_fix(vcf_path, out_path, ploidy_table_path=ploidy_table)

    assert result.wrote_output is True
    assert result.fixed_summary is not None
    assert result.fixed_summary.has_errors is False
    with pysam.VariantFile(str(out_path)) as vcf:
        record = next(iter(vcf))
        assert "ECN" in record.format
        assert record.samples["S1"]["ECN"] == 2
        assert record.samples["S2"]["ECN"] == 1


def test_validate_and_fix_reports_missing_ecn_without_ploidy_table(make_vcf, tmp_path) -> None:
    vcf_path = make_vcf(
        file_name="missing_ecn_unfixable.vcf",
        records=[
            "chr1\t100\tvar1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=25\tGT:GQ\t0/1:60\t0/0:50",
        ],
    )
    out_path = tmp_path / "blocked_ecn.vcf"

    result = validate_and_fix(vcf_path, out_path)

    assert result.wrote_output is False
    assert result.fixed_summary is None
    assert "MISSING_ECN" in render_fix_result(result)


def test_validate_and_fix_converts_cn_multiallelic_record_to_cnv(make_vcf, tmp_path) -> None:
    vcf_path = make_vcf(
        file_name="cn_multiallelic_fixable.vcf",
        records=[
            "chr1\t100\tvar1\tN\t<CN0>,<CN1>,<CN2>,<CN3>\t.\tPASS\tSVTYPE=DUP;SVLEN=10\tGT:GQ:ECN\t0/1:40:2\t0/0:35:2",
        ],
    )
    out_path = tmp_path / "cnv_fixed.vcf.gz"

    result = validate_and_fix(vcf_path, out_path)

    assert result.wrote_output is True
    assert result.fixed_summary is not None
    assert result.fixed_summary.has_errors is False
    with pysam.VariantFile(str(out_path)) as vcf:
        record = next(iter(vcf))
        assert record.alts == ("<CNV>",)
        assert str(record.info["SVTYPE"]) == "CNV"
        assert record.samples["S1"]["GT"] == (None, None)
        assert record.samples["S2"]["GT"] == (None, None)


def test_validate_and_fix_converts_sparse_cn_multiallelic_record_to_cnv(make_vcf, tmp_path) -> None:
    vcf_path = make_vcf(
        file_name="sparse_cn_multiallelic_fixable.vcf",
        records=[
            "chr1\t100\tvar1\tN\t<CN0>,<CN2>,<CN3>\t.\tPASS\tSVTYPE=DUP;SVLEN=10\tGT:GQ:ECN\t0/1:40:2\t0/0:35:2",
        ],
    )
    out_path = tmp_path / "sparse_cnv_fixed.vcf.gz"

    result = validate_and_fix(vcf_path, out_path)

    assert result.wrote_output is True
    assert result.fixed_summary is not None
    assert result.fixed_summary.has_errors is False
    with pysam.VariantFile(str(out_path)) as vcf:
        record = next(iter(vcf))
        assert record.alts == ("<CNV>",)
        assert str(record.info["SVTYPE"]) == "CNV"
        assert record.samples["S1"]["GT"] == (None, None)
        assert record.samples["S2"]["GT"] == (None, None)


def test_validate_and_fix_missing_ecn_does_not_surface_cn_multiallelic_as_blocking(make_vcf, tmp_path) -> None:
    """MULTI_ALLELIC_NON_CNV (fixable) must not appear in the blocking-error list
    when MISSING_ECN is the actual reason the fix was aborted."""
    vcf_path = make_vcf(
        file_name="cn_multiallelic_missing_ecn.vcf",
        records=[
            # <CN0>,<CN2>,<CN3> alts are fixable; GT:GQ without ECN is the real blocker
            "chr1\t100\tvar1\tN\t<CN0>,<CN2>,<CN3>\t.\tPASS\tSVTYPE=DUP;SVLEN=10\tGT:GQ\t0/1:40\t0/0:35",
        ],
    )
    out_path = tmp_path / "blocked_cn_ecn.vcf.gz"

    result = validate_and_fix(vcf_path, out_path)

    assert result.wrote_output is False
    assert result.fixed_summary is None
    rendered = render_fix_result(result)
    # MISSING_ECN is the actual blocker and must appear in the "Fix mode aborted" line
    assert "Fix mode aborted: unresolved critical checks remain: MISSING_ECN" in rendered
    # MULTI_ALLELIC_NON_CNV is fixable; it must not appear in the Examples / blocking section
    blocking_section = rendered.split("Fix mode aborted:")[-1]
    assert "MULTI_ALLELIC_NON_CNV" not in blocking_section
    assert result.unresolved_error_counts().get("MULTI_ALLELIC_NON_CNV", 0) == 0


def test_validate_and_fix_still_blocks_non_cn_multiallelic_record(make_vcf, tmp_path) -> None:
    vcf_path = make_vcf(
        file_name="non_cn_multiallelic_unfixable.vcf",
        records=[
            "chr1\t100\tvar1\tN\t<DEL>,<DUP>\t.\tPASS\tSVTYPE=DUP;SVLEN=10\tGT:GQ:ECN\t0/1:40:2\t0/0:35:2",
        ],
    )
    out_path = tmp_path / "non_cn_multiallelic_fixed.vcf.gz"

    result = validate_and_fix(vcf_path, out_path)

    assert result.wrote_output is False
    rendered = render_fix_result(result)
    assert "MULTI_ALLELIC_NON_CNV" in rendered
    assert "Examples:" in rendered
    assert "[var1]" in rendered


def test_validate_and_fix_writes_bgzip_and_tabix_for_gz_output(make_vcf, tmp_path) -> None:
    vcf_path = make_vcf(
        file_name="fixable_for_bgzip.vcf",
        records=[
            "chr1\t100\tvar1\tN\t<INS:ME:ALU>\t.\tPASS\tSVLEN=25\tGT:GQ:ECN\t0/1:60:2\t0/0:50:2",
        ],
    )
    out_path = tmp_path / "fixed.vcf.gz"

    result = validate_and_fix(vcf_path, out_path)

    assert result.wrote_output is True
    assert out_path.exists()
    assert (tmp_path / "fixed.vcf.gz.tbi").exists()
    with pysam.VariantFile(str(out_path)) as vcf:
        record = next(iter(vcf))
        assert str(record.info["SVTYPE"]) == "INS"
