from __future__ import annotations

import pandas as pd
import pytest

from gatk_sv_compare.aggregate import aggregate
from gatk_sv_compare.config import AnalysisConfig
from gatk_sv_compare.modules.genotype_exact_match import (
    GenotypeExactMatchModule,
    build_exact_match_table,
)


def test_build_exact_match_table_computes_expected_rates(module_test_context) -> None:
    table = build_exact_match_table(module_test_context.data, module_test_context.config)

    assert set(table["variant_id_a"]) == {"a_del", "a_ins"}
    del_row = table.loc[table["variant_id_a"] == "a_del"].iloc[0]
    ins_row = table.loc[table["variant_id_a"] == "a_ins"].iloc[0]
    assert int(del_row["n_compared"]) == 3
    assert del_row["exact_match_rate"] == pytest.approx(1.0 / 3.0)
    assert del_row["homref_to_het_rate"] == pytest.approx(1.0 / 3.0)
    assert del_row["homalt_to_homref_rate"] == pytest.approx(1.0 / 3.0)
    assert ins_row["exact_match_rate"] == pytest.approx(1.0)


def test_genotype_exact_match_module_writes_outputs(module_test_context) -> None:
    module = GenotypeExactMatchModule()
    module_test_context.config.enable_site_match_table = True

    module.run(module_test_context.data, module_test_context.config)

    output_dir = module_test_context.config.output_dir / "genotype_exact_match"
    per_site = pd.read_csv(output_dir / "tables" / "genotype_match_rates.per_site.tsv.gz", sep="\t")
    summary = pd.read_csv(output_dir / "tables" / "genotype_match_rates.tsv.gz", sep="\t")
    assert not per_site.empty
    assert not summary.empty
    assert "pair_id" not in per_site.columns
    assert {"svtype_CallsetA", "size_bucket_CallsetA", "af_bucket_CallsetA", "evidence_bucket_CallsetA", "svtype_CallsetB", "size_bucket_CallsetB", "af_bucket_CallsetB", "evidence_bucket_CallsetB"}.issubset(per_site.columns)
    assert not (output_dir / "exact_match.by_type.png").exists()
    assert not (output_dir / "exact_match.by_context.png").exists()


def test_genotype_exact_match_site_table_disabled_by_default(module_test_context) -> None:
    module = GenotypeExactMatchModule()

    module.run(module_test_context.data, module_test_context.config)

    output_dir = module_test_context.config.output_dir / "genotype_exact_match"
    assert not (output_dir / "tables" / "genotype_match_rates.per_site.tsv.gz").exists()
    assert (output_dir / "tables" / "genotype_match_rates.tsv.gz").exists()


def test_build_exact_match_table_compares_str_repcn(make_vcf, tmp_path) -> None:
    extra_headers = [
        "##INFO=<ID=STATUS,Number=1,Type=String,Description=\"Match status\">",
        "##INFO=<ID=TRUTH_VID,Number=1,Type=String,Description=\"Truth variant id\">",
        "##INFO=<ID=RU,Number=1,Type=String,Description=\"Repeat unit\">",
        "##INFO=<ID=PERIOD,Number=1,Type=Integer,Description=\"Repeat period\">",
        "##INFO=<ID=LOCUS,Number=1,Type=String,Description=\"STR locus\">",
        "##FORMAT=<ID=REPCN,Number=.,Type=Integer,Description=\"Repeat copy number\">",
    ]
    vcf_a = make_vcf(
        file_name="str_a.vcf",
        sample_names=["S1", "S2", "S3"],
        extra_header_lines=extra_headers,
        records=[
            "chr1\t100\ta_str\tN\t<STR>\t.\tPASS\tSVTYPE=STR;END=150;RU=AC;PERIOD=2;LOCUS=L1;STATUS=MATCHED;TRUTH_VID=b_str\tGT:REPCN\t0/1:10,10\t0/1:8,9\t0/1:7,7",
        ],
    )
    vcf_b = make_vcf(
        file_name="str_b.vcf",
        sample_names=["S1", "S2", "S3"],
        extra_header_lines=extra_headers,
        records=[
            "chr1\t100\tb_str\tN\t<STR>\t.\tPASS\tSVTYPE=STR;END=150;RU=AC;PERIOD=2;LOCUS=L1;STATUS=MATCHED;TRUTH_VID=a_str\tGT:REPCN\t0/1:10,10\t0/1:9,8\t0/1:7,8",
        ],
    )
    config = AnalysisConfig(
        vcf_a_path=vcf_a,
        vcf_b_path=vcf_b,
        output_dir=tmp_path / "out",
        contigs=["chr1"],
        n_workers=1,
    )
    data = aggregate(config)

    table = build_exact_match_table(data, config)

    assert table.shape[0] == 1
    row = table.iloc[0]
    assert row["svtype_a"] == "STR"
    assert int(row["n_compared"]) == 3
    assert row["exact_match_rate"] == pytest.approx(2.0 / 3.0)
