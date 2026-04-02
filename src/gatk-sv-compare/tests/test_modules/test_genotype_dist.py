from __future__ import annotations

import pandas as pd

from gatk_sv_compare.modules.genotype_dist import GenotypeDistModule, build_hwe_table


def test_build_hwe_table_returns_expected_rows(module_test_context) -> None:
    table = build_hwe_table(module_test_context.data.sites_a, pass_only=True)

    assert set(table["variant_id"]) == {"a_del", "a_ins"}
    assert set(table["hwe_class"]).issubset({"pass", "nominal", "bonferroni"})
    assert (table[["aa", "ab", "bb"]].sum(axis=1).round(6) == 1.0).all()


def test_genotype_dist_module_writes_outputs(module_test_context) -> None:
    module = GenotypeDistModule()

    module.run(module_test_context.data, module_test_context.config)

    output_dir = module_test_context.config.output_dir / "genotype_dist"
    table_a = pd.read_csv(output_dir / "tables" / "hwe_stats.CallsetA.tsv", sep="\t")
    table_b = pd.read_csv(output_dir / "tables" / "hwe_stats.CallsetB.tsv", sep="\t")
    assert not table_a.empty
    assert not table_b.empty
    assert (output_dir / "ternary.all.CallsetA.png").exists()
    assert (output_dir / "carrier_freq_vs_af.CallsetB.png").exists()
