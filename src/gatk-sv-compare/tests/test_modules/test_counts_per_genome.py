from __future__ import annotations

import pandas as pd

from gatk_sv_compare.modules.counts_per_genome import CountsPerGenomeModule, collect_per_sample_counts


def test_collect_per_sample_counts_counts_alt_sites_and_alleles(module_test_context) -> None:
    counts = collect_per_sample_counts(module_test_context.config.vcf_a_path, pass_only=True)

    del_s3 = counts.loc[(counts["sample"] == "S3") & (counts["svtype"] == "DEL")].iloc[0]
    assert int(del_s3["sites"]) == 1
    assert int(del_s3["alleles"]) == 2


def test_counts_per_genome_module_writes_outputs(module_test_context) -> None:
    module = CountsPerGenomeModule()

    module.run(module_test_context.data, module_test_context.config)

    output_dir = module_test_context.config.output_dir / "counts_per_genome"
    table_a = pd.read_csv(output_dir / "tables" / "per_sample_counts.CallsetA.tsv", sep="\t")
    table_b = pd.read_csv(output_dir / "tables" / "per_sample_counts.CallsetB.tsv", sep="\t")
    assert not table_a.empty
    assert not table_b.empty
    assert (output_dir / "sites_per_genome.by_type.CallsetA.png").exists()
    assert (output_dir / "alleles_per_genome.by_type.CallsetB.png").exists()
