from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd

from gatk_sv_compare.modules.genotype_dist import GenotypeDistModule, build_hwe_table, summarize_hwe_by_bucket
from gatk_sv_compare.plot_utils import plot_ternary


def test_build_hwe_table_returns_expected_rows(module_test_context) -> None:
    table = build_hwe_table(module_test_context.data.sites_a, pass_only=True)

    assert set(table["variant_id"]) == {"a_del", "a_ins"}
    assert set(table["hwe_class"]).issubset({"pass", "nominal", "bonferroni"})
    assert (table[["aa", "ab", "bb"]].sum(axis=1).round(6) == 1.0).all()


def test_genotype_dist_module_writes_outputs(module_test_context) -> None:
    module = GenotypeDistModule()

    module.run(module_test_context.data, module_test_context.config)

    output_dir = module_test_context.config.output_dir / "genotype_dist"
    table = pd.read_csv(output_dir / "tables" / "hwe_stats.tsv.gz", sep="\t")
    assert not table.empty
    assert {"n_variants_CallsetA", "frac_pass_CallsetA", "frac_nominal_CallsetA", "frac_bonferroni_CallsetA", "n_variants_CallsetB", "frac_pass_CallsetB", "frac_nominal_CallsetB", "frac_bonferroni_CallsetB"}.issubset(table.columns)
    assert (output_dir / "ternary.all.CallsetA.png").exists()
    assert (output_dir / "carrier_freq_vs_af.CallsetB.png").exists()


def test_summarize_hwe_by_bucket_writes_label_specific_columns(module_test_context) -> None:
    table = build_hwe_table(module_test_context.data.sites_a, pass_only=True)
    summary = summarize_hwe_by_bucket(table, "CallsetA")

    assert not summary.empty
    assert {"n_variants_CallsetA", "frac_pass_CallsetA", "mean_af_CallsetA"}.issubset(summary.columns)


def test_plot_ternary_labels_corners() -> None:
    fig, ax = plt.subplots()

    plot_ternary(ax, [0.5], [0.25], [0.25], ["#000000"])

    labels = {text.get_text() for text in ax.texts}
    assert {"REF", "HET", "HOMVAR"}.issubset(labels)
    plt.close(fig)
