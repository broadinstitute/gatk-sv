from __future__ import annotations

from pathlib import Path

import pandas as pd

from gatk_sv_compare.modules.counts_per_genome import CountsPerGenomeModule, collect_per_sample_counts


def test_collect_per_sample_counts_counts_alt_sites_and_alleles(module_test_context) -> None:
    counts = collect_per_sample_counts(module_test_context.config.vcf_a_path, pass_only=True)

    del_s3 = counts.loc[(counts["sample"] == "S3") & (counts["svtype"] == "DEL")].iloc[0]
    assert int(del_s3["sites"]) == 1
    assert int(del_s3["alleles"]) == 2
    assert counts.loc[counts["svtype"] == "CNV"].empty


def test_counts_per_genome_module_writes_outputs(module_test_context) -> None:
    module = CountsPerGenomeModule()
    module_test_context.config.per_sample_counts_table = True

    module.run(module_test_context.data, module_test_context.config)

    output_dir = module_test_context.config.output_dir / "counts_per_genome"
    table_a = pd.read_csv(output_dir / "tables" / "per_sample_counts.CallsetA.tsv.gz", sep="\t")
    table_b = pd.read_csv(output_dir / "tables" / "per_sample_counts.CallsetB.tsv.gz", sep="\t")
    assert not table_a.empty
    assert not table_b.empty
    assert (output_dir / "sites_per_genome.by_type.png").exists()
    assert (output_dir / "alleles_per_genome.by_type.png").exists()


def test_counts_per_genome_table_disabled_by_default(module_test_context) -> None:
    module = CountsPerGenomeModule()

    module.run(module_test_context.data, module_test_context.config)

    output_dir = module_test_context.config.output_dir / "counts_per_genome"
    assert not (output_dir / "tables" / "per_sample_counts.CallsetA.tsv.gz").exists()


def test_collect_per_sample_counts_treats_empty_filter_as_pass(make_vcf, tmp_path: Path) -> None:
    vcf_path = make_vcf(
        file_name="empty_filter_pass.vcf",
        final_annotated=True,
        sample_names=["S1", "S2"],
        records=[
            "chr1\t100\tv1\tN\t<DEL>\t.\t.\tSVTYPE=DEL;SVLEN=100;AC=1;AF=0.25;AN=4\tGT:GQ:ECN\t0/1:50:2\t0/0:45:2",
            "chr1\t200\tv2\tN\t<DUP>\t.\tMULTIALLELIC\tSVTYPE=DUP;SVLEN=200;AC=1;AF=0.25;AN=4\tGT:GQ:ECN\t0/1:55:2\t0/0:40:2",
        ],
    )

    counts = collect_per_sample_counts(vcf_path, pass_only=True)

    assert set(counts["svtype"]) == {"DEL"}
    del_s1 = counts.loc[(counts["sample"] == "S1") & (counts["svtype"] == "DEL")].iloc[0]
    assert int(del_s1["sites"]) == 1
    assert int(del_s1["alleles"]) == 1
