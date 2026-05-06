from __future__ import annotations

import pysam

from gatk_sv_compare.vcf_format import PipelineStage, check_record, detect_pipeline_stage, filter_values


def test_detect_pipeline_stage_final_annotated(make_vcf) -> None:
    vcf_path = make_vcf(
        file_name="final.vcf",
        final_annotated=True,
        records=[
            "chr1\t100\tvar1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=50;N_BI_GENOS=2;N_HOMREF=1;N_HET=1;N_HOMALT=0;gnomad_v4.1_sv_AF=0.1\tGT:GQ:ECN:OGQ:SL\t0/1:60:2:55:1.0\t0/0:50:2:49:1.2",
        ],
    )
    with pysam.VariantFile(str(vcf_path)) as vcf:
        assert detect_pipeline_stage(vcf) is PipelineStage.FINAL_ANNOTATED


def test_check_record_reports_ins_end_mismatch(make_vcf) -> None:
    vcf_path = make_vcf(
        file_name="ins_mismatch.vcf",
        records=[
            "chr1\t100\tvar1\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=10;END=150\tGT:GQ:ECN\t0/1:40:2\t0/0:35:2",
        ],
    )
    with pysam.VariantFile(str(vcf_path)) as vcf:
        record = next(iter(vcf))
    issues = check_record(record, contig_length=1_000_000)
    assert any(issue.check_id == "INS_END_MISMATCH" and issue.severity == "INFO" for issue in issues)


def test_check_record_flags_multiallelic_non_cnv(make_vcf) -> None:
    vcf_path = make_vcf(
        file_name="multi.vcf",
        records=[
            "chr1\t100\tvar1\tN\t<CN0>,<CN1>,<CN2>,<CN3>\t.\tPASS\tSVTYPE=DUP;SVLEN=10\tGT:GQ:ECN\t0/1:40:2\t0/0:35:2",
        ],
    )
    with pysam.VariantFile(str(vcf_path)) as vcf:
        record = next(iter(vcf))
    issues = check_record(record)
    assert any(issue.check_id == "MULTI_ALLELIC_NON_CNV" for issue in issues)


def test_check_record_does_not_emit_cnv_no_gt_info(make_vcf) -> None:
    vcf_path = make_vcf(
        file_name="cnv.vcf",
        records=[
            "chr1\t100\tcnv1\tN\t<CNV>\t.\tMULTIALLELIC\tSVTYPE=CNV;SVLEN=1000\tGT:GQ:ECN\t./.:40:2\t./.:35:2",
        ],
    )
    with pysam.VariantFile(str(vcf_path)) as vcf:
        record = next(iter(vcf))
    issues = check_record(record)
    assert all(issue.check_id != "CNV_NO_GT" for issue in issues)


def test_filter_values_treats_empty_filter_as_empty_set(make_vcf) -> None:
    vcf_path = make_vcf(
        file_name="empty_filter.vcf",
        records=[
            "chr1\t100\tvar1\tN\t<DEL>\t.\t.\tSVTYPE=DEL;SVLEN=10\tGT:GQ:ECN\t0/1:40:2\t0/0:35:2",
        ],
    )
    with pysam.VariantFile(str(vcf_path)) as vcf:
        record = next(iter(vcf))

    assert filter_values(record) == set()


def test_check_record_handles_missing_cpx_type_header(make_vcf) -> None:
    vcf_path = make_vcf(
        file_name="cpx_missing_header.vcf",
        records=[
            "chr1\t100\tvar1\tN\t<CPX>\t.\tPASS\tSVTYPE=CPX;SVLEN=10;END=150;CPX_TYPE=dDUP\tGT:GQ:ECN\t0/1:40:2\t0/0:35:2",
        ],
    )
    text = vcf_path.read_text()
    vcf_path.write_text(text.replace('##INFO=<ID=CPX_TYPE,Number=1,Type=String,Description="Complex subtype">\n', ''))

    with pysam.VariantFile(str(vcf_path)) as vcf:
        record = next(iter(vcf))

    issues = check_record(record)

    assert isinstance(issues, list)


def test_check_record_allows_999_scale_gq(make_vcf) -> None:
    vcf_path = make_vcf(
        file_name="gq_999.vcf",
        records=[
            "chr1\t100\tvar1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=50\tGT:GQ:ECN\t0/1:600:2\t0/0:300:2",
        ],
    )
    with pysam.VariantFile(str(vcf_path)) as vcf:
        record = next(iter(vcf))

    issues = check_record(record)

    assert any(issue.check_id == "GQ_RANGE" and issue.severity == "WARN" for issue in issues)


def test_check_record_allows_str_without_ecn(make_vcf) -> None:
    vcf_path = make_vcf(
        file_name="str.vcf",
        extra_header_lines=[
            "##INFO=<ID=RU,Number=1,Type=String,Description=\"Repeat unit\">",
            "##INFO=<ID=PERIOD,Number=1,Type=Integer,Description=\"Repeat period\">",
            "##INFO=<ID=LOCUS,Number=1,Type=String,Description=\"STR locus\">",
            "##FORMAT=<ID=REPCN,Number=.,Type=Integer,Description=\"Repeat copy number\">",
        ],
        records=[
            "chr1\t100\tstr1\tN\t<STR>\t.\tPASS\tSVTYPE=STR;END=150;RU=AC;PERIOD=2;LOCUS=L1\tGT:REPCN\t0/1:10,11\t0/0:10,10",
        ],
    )
    with pysam.VariantFile(str(vcf_path)) as vcf:
        record = next(iter(vcf))

    issues = check_record(record)

    assert all(issue.check_id not in {"UNKNOWN_SVTYPE", "MISSING_ECN"} for issue in issues)
