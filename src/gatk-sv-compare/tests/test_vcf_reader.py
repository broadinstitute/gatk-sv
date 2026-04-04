from __future__ import annotations

import numpy as np
import pytest

from gatk_sv_compare.vcf_reader import iter_contig


def test_iter_contig_uses_precomputed_counts_and_concordance(make_vcf) -> None:
    vcf_path = make_vcf(
        file_name="annotated.vcf",
        final_annotated=True,
        extra_header_lines=[
            "##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description=\"Calling algorithms\">",
            "##INFO=<ID=EVIDENCE,Number=.,Type=String,Description=\"Evidence types\">",
            "##INFO=<ID=STATUS,Number=1,Type=String,Description=\"Match status\">",
            "##INFO=<ID=TRUTH_VID,Number=1,Type=String,Description=\"Truth variant id\">",
            "##INFO=<ID=OVERLAP_FRAC_SEGDUP,Number=1,Type=Float,Description=\"Segmental duplication overlap fraction\">",
            "##INFO=<ID=NUM_END_OVERLAPS_SEGDUP,Number=1,Type=Integer,Description=\"Segmental duplication end overlap count\">",
            "##INFO=<ID=VAR_PPV,Number=1,Type=Float,Description=\"Variant PPV\">",
        ],
        records=[
            "chr1\t100\tvar1\tN\t<INS:ME:ALU>\t.\tPASS\tSVTYPE=INS;SVLEN=310;ALGORITHMS=melt,wham;EVIDENCE=SR,BAF;N_BI_GENOS=2;N_HOMREF=1;N_HET=1;N_HOMALT=0;STATUS=MATCHED;TRUTH_VID=truth1;OVERLAP_FRAC_SEGDUP=0.1;NUM_END_OVERLAPS_SEGDUP=1;VAR_PPV=0.95;gnomad_v4.1_sv_AF=0.1\tGT:GQ:ECN:OGQ:SL\t0/1:60:2:55:1.0\t0/0:50:2:49:1.2",
        ],
    )

    records = list(iter_contig(vcf_path, "chr1", extract_concordance=True))
    assert len(records) == 1
    record = records[0]
    assert record.svtype == "INS:MEI"
    assert record.n_hom_ref == 1
    assert record.n_het == 1
    assert record.status == "MATCHED"
    assert record.truth_vid == "truth1"
    assert record.genomic_context == "segdup"
    assert record.algorithms == ("melt", "wham")
    assert record.evidence_bucket == "SR"
    assert record.concordance_metrics is not None
    assert record.concordance_metrics["VAR_PPV"] == pytest.approx(0.95)


def test_iter_contig_falls_back_to_gt_counts_and_subsets_gq(make_vcf) -> None:
    vcf_path = make_vcf(
        file_name="cleaned.vcf",
        sample_names=["S1", "S2", "S3"],
        records=[
            "chr1\t100\tvar1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=100\tGT:GQ:ECN\t0/0:10:2\t0/1:20:2\t1/1:30:2",
        ],
    )

    records = list(iter_contig(vcf_path, "chr1", extract_gq=True, sample_indices=np.asarray([0, 2], dtype=int)))
    assert len(records) == 1
    record = records[0]
    assert record.n_bi_genos == 3
    assert record.n_hom_ref == 1
    assert record.n_het == 1
    assert record.n_hom_alt == 1
    assert record.gq_array is not None
    assert record.gq_array.tolist() == [10.0, 30.0]


def test_iter_contig_keeps_raw_999_scale_gq_values(make_vcf) -> None:
    vcf_path = make_vcf(
        file_name="scaled_gq.vcf",
        sample_names=["S1", "S2", "S3"],
        records=[
            "chr1\t100\tvar1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=100\tGT:GQ:ECN\t0/0:100:2\t0/1:200:2\t1/1:300:2",
        ],
    )

    records = list(iter_contig(vcf_path, "chr1", extract_gq=True, sample_indices=np.asarray([0, 2], dtype=int)))

    assert len(records) == 1
    record = records[0]
    assert record.gq_array is not None
    assert record.gq_array.tolist() == [100.0, 300.0]


def test_iter_contig_handles_cnv_frequency(make_vcf) -> None:
    vcf_path = make_vcf(
        file_name="cnv.vcf",
        extra_header_lines=[
            "##INFO=<ID=CN_NONREF_FREQ,Number=1,Type=Float,Description=\"CNV nonref frequency\">",
            "##INFO=<ID=CN_NUMBER,Number=1,Type=Integer,Description=\"CNV count\">",
        ],
        records=[
            "chr1\t100\tcnv1\tN\t<CNV>\t.\tMULTIALLELIC\tSVTYPE=CNV;SVLEN=1000;CN_NONREF_FREQ=0.25;CN_NUMBER=2\tGT:GQ:ECN\t./.:40:2\t./.:35:2",
        ],
    )

    records = list(iter_contig(vcf_path, "chr1"))
    assert len(records) == 1
    record = records[0]
    assert record.svtype == "CNV"
    assert record.af == 0.25
    assert record.cn_nonref_freq == 0.25
    assert record.n_bi_genos == 2


def test_iter_contig_selects_highest_qualifying_context_overlap(make_vcf) -> None:
    vcf_path = make_vcf(
        file_name="contexts.vcf",
        extra_header_lines=[
            "##INFO=<ID=OVERLAP_FRAC_SEGDUP,Number=1,Type=Float,Description=\"Segmental duplication overlap fraction\">",
            "##INFO=<ID=OVERLAP_FRAC_SIMPLE_REPEAT,Number=1,Type=Float,Description=\"Simple repeat overlap fraction\">",
            "##INFO=<ID=NUM_END_OVERLAPS_SEGDUP,Number=1,Type=Integer,Description=\"Segmental duplication end overlap count\">",
            "##INFO=<ID=NUM_END_OVERLAPS_SIMPLE_REPEAT,Number=1,Type=Integer,Description=\"Simple repeat end overlap count\">",
        ],
        records=[
            "chr1\t100\tctx1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=100;OVERLAP_FRAC_SEGDUP=0.8;OVERLAP_FRAC_SIMPLE_REPEAT=0.6;NUM_END_OVERLAPS_SEGDUP=2;NUM_END_OVERLAPS_SIMPLE_REPEAT=2\tGT:GQ:ECN\t0/1:40:2\t0/0:35:2",
        ],
    )

    records = list(iter_contig(vcf_path, "chr1"))

    assert len(records) == 1
    assert records[0].genomic_context == "segdup"


def test_iter_contig_applies_context_overlap_threshold_for_span_variants(make_vcf) -> None:
    vcf_path = make_vcf(
        file_name="thresholds.vcf",
        extra_header_lines=[
            "##INFO=<ID=OVERLAP_FRAC_SEGDUP,Number=1,Type=Float,Description=\"Segmental duplication overlap fraction\">",
            "##INFO=<ID=NUM_END_OVERLAPS_SEGDUP,Number=1,Type=Integer,Description=\"Segmental duplication end overlap count\">",
        ],
        records=[
            "chr1\t100\tctx2\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=100;OVERLAP_FRAC_SEGDUP=0.4;NUM_END_OVERLAPS_SEGDUP=2\tGT:GQ:ECN\t0/1:40:2\t0/0:35:2",
        ],
    )

    default_records = list(iter_contig(vcf_path, "chr1"))
    relaxed_records = list(iter_contig(vcf_path, "chr1", context_overlap=0.4))

    assert default_records[0].genomic_context == "none"
    assert relaxed_records[0].genomic_context == "segdup"
