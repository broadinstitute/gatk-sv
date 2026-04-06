from __future__ import annotations

import pandas as pd

from gatk_sv_compare.aggregate import aggregate, build_matched_pairs, _apply_min_size_filter
from gatk_sv_compare.config import AnalysisConfig


def test_build_matched_pairs_uses_truth_vid_symmetry() -> None:
    sites_a = pd.DataFrame(
        [
            {
                "variant_id": "a1",
                "contig": "chr1",
                "svtype": "DEL",
                "af": 0.1,
                "status": "MATCHED",
                "truth_vid": "b1",
            }
        ]
    )
    sites_b = pd.DataFrame(
        [
            {
                "variant_id": "b1",
                "contig": "chr1",
                "svtype": "DEL",
                "af": 0.2,
                "status": "MATCHED",
                "truth_vid": "a1",
            }
        ]
    )

    matched = build_matched_pairs(sites_a, sites_b)
    assert matched.shape[0] == 1
    assert matched.loc[0, "variant_id_a"] == "a1"
    assert matched.loc[0, "variant_id_b"] == "b1"


def test_aggregate_builds_site_tables_and_matched_pairs(tmp_path, make_vcf) -> None:
    extra_headers = [
        "##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description=\"Calling algorithms\">",
        "##INFO=<ID=EVIDENCE,Number=.,Type=String,Description=\"Evidence types\">",
        "##INFO=<ID=STATUS,Number=1,Type=String,Description=\"Match status\">",
        "##INFO=<ID=TRUTH_VID,Number=1,Type=String,Description=\"Truth variant id\">",
        "##INFO=<ID=OVERLAP_FRAC_SEGDUP,Number=1,Type=Float,Description=\"Segmental duplication overlap fraction\">",
    ]
    vcf_a = make_vcf(
        file_name="a.vcf",
        sample_names=["S1", "S2", "S3"],
        extra_header_lines=extra_headers,
        records=[
            "chr1\t100\ta1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=100;ALGORITHMS=manta,wham;EVIDENCE=RD,PE;STATUS=MATCHED;TRUTH_VID=b1;OVERLAP_FRAC_SEGDUP=0.8\tGT:GQ:ECN\t0/0:10:2\t0/1:20:2\t1/1:30:2",
            "chr1\t200\ta2\tN\t<CNV>\t.\tMULTIALLELIC\tSVTYPE=CNV;SVLEN=1000;ALGORITHMS=cnmops;EVIDENCE=RD,PE,SR,BAF\tGT:GQ:ECN\t./.:40:2\t./.:35:2\t./.:25:2",
        ],
    )
    vcf_b = make_vcf(
        file_name="b.vcf",
        sample_names=["S2", "S3", "S4"],
        extra_header_lines=extra_headers,
        records=[
            "chr1\t110\tb1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=110;ALGORITHMS=manta;EVIDENCE=RD,PE;STATUS=MATCHED;TRUTH_VID=a1\tGT:GQ:ECN\t0/1:18:2\t0/1:22:2\t0/0:40:2",
        ],
    )
    config = AnalysisConfig(
        vcf_a_path=vcf_a,
        vcf_b_path=vcf_b,
        output_dir=tmp_path / "out",
        contigs=["chr1"],
        n_workers=1,
        vcf_a_label="A",
        vcf_b_label="B",
    )

    data = aggregate(config)

    assert data.sites_a.shape[0] == 2
    assert data.sites_b.shape[0] == 1
    assert data.shared_samples == ["S2", "S3"]
    assert data.sample_indices_a.tolist() == [1, 2]
    assert data.sample_indices_b.tolist() == [0, 1]
    assert data.matched_pairs.shape[0] == 1
    assert bool(data.sites_a.loc[data.sites_a["variant_id"] == "a2", "in_filtered_pass_view"].iloc[0]) is False
    assert data.sites_a.loc[data.sites_a["variant_id"] == "a1", "genomic_context"].iloc[0] == "segdup"
    assert data.sites_a.loc[data.sites_a["variant_id"] == "a1", "algorithms"].iloc[0] == "manta,wham"
    assert data.sites_a.loc[data.sites_a["variant_id"] == "a1", "evidence_bucket"].iloc[0] == "RD,PE"
    assert data.sites_a.loc[data.sites_a["variant_id"] == "a2", "evidence_bucket"].iloc[0] == "RD,PE,SR"
    assert (tmp_path / "out" / "aggregate" / "sites_a.chr1.parquet").exists()
    assert (tmp_path / "out" / "aggregate" / "matched_pairs.parquet").exists()


def test_apply_min_size_filter_removes_small_keeps_none() -> None:
    sites = pd.DataFrame({
        "variant_id": ["small_del", "big_del", "no_size_bnd"],
        "svlen": [500, 6000, None],
    })
    result = _apply_min_size_filter(sites, min_size=5000)
    assert list(result["variant_id"]) == ["big_del", "no_size_bnd"]


def test_apply_min_size_filter_empty_dataframe() -> None:
    sites = pd.DataFrame(columns=["variant_id", "svlen"])
    result = _apply_min_size_filter(sites, min_size=5000)
    assert result.empty


def test_aggregate_with_min_size_filters_before_matching(tmp_path, make_vcf) -> None:
    extra_headers = [
        "##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description=\"Calling algorithms\">",
        "##INFO=<ID=EVIDENCE,Number=.,Type=String,Description=\"Evidence types\">",
        "##INFO=<ID=STATUS,Number=1,Type=String,Description=\"Match status\">",
        "##INFO=<ID=TRUTH_VID,Number=1,Type=String,Description=\"Truth variant id\">",
    ]
    # a1 (SVLEN=100) should be filtered; a2 (SVLEN=6000) should remain
    vcf_a = make_vcf(
        file_name="a.vcf",
        sample_names=["S1", "S2"],
        extra_header_lines=extra_headers,
        records=[
            "chr1\t100\ta1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=100;ALGORITHMS=manta;EVIDENCE=RD;STATUS=MATCHED;TRUTH_VID=b1\tGT:GQ:ECN\t0/1:10:2\t0/0:20:2",
            "chr1\t500\ta2\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=6000;ALGORITHMS=depth;EVIDENCE=RD;STATUS=UNMATCHED\tGT:GQ:ECN\t0/1:30:2\t0/0:40:2",
        ],
    )
    vcf_b = make_vcf(
        file_name="b.vcf",
        sample_names=["S1", "S2"],
        extra_header_lines=extra_headers,
        records=[
            "chr1\t110\tb1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=110;ALGORITHMS=manta;EVIDENCE=RD;STATUS=MATCHED;TRUTH_VID=a1\tGT:GQ:ECN\t0/1:15:2\t0/0:25:2",
        ],
    )
    config = AnalysisConfig(
        vcf_a_path=vcf_a,
        vcf_b_path=vcf_b,
        output_dir=tmp_path / "out",
        contigs=["chr1"],
        n_workers=1,
        min_size=5000,
    )
    data = aggregate(config)
    # a1 (svlen=100) should have been filtered out
    assert data.sites_a.shape[0] == 1
    assert data.sites_a["variant_id"].iloc[0] == "a2"
    # b1 (svlen=110) should also be filtered
    assert data.sites_b.shape[0] == 0
    # No matches possible after filtering
    assert data.matched_pairs.shape[0] == 0
