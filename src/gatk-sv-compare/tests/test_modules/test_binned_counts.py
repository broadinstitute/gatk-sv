from __future__ import annotations

import pandas as pd

from gatk_sv_compare.modules.binned_counts import BinnedCountsModule, build_combined_binned_counts, summarize_binned_counts


def test_summarize_binned_counts_respects_pass_view(module_test_context) -> None:
    sites = module_test_context.data.sites_a

    summary = summarize_binned_counts(sites, pass_only=True)

    assert set(summary["svtype"]) == {"DEL", "INS:MEI", "CNV"}
    assert int(summary["n_variants"].sum()) == 3
    assert int(summary["n_matched"].sum()) == 2
    assert int(summary["n_unmatched"].sum()) == 1


def test_binned_counts_module_writes_outputs(module_test_context) -> None:
    module = BinnedCountsModule()

    module.run(module_test_context.data, module_test_context.config)

    output_dir = module_test_context.config.output_dir / "binned_counts"
    counts = pd.read_csv(output_dir / "counts.tsv.gz", sep="\t")
    assert (output_dir / "counts.parquet").exists()
    assert not counts.empty
    assert {"n_variants_CallsetA", "n_matched_CallsetA", "n_unmatched_CallsetA", "n_variants_CallsetB", "n_matched_CallsetB", "n_unmatched_CallsetB"}.issubset(counts.columns)


def test_build_combined_binned_counts_merges_both_callsets(module_test_context) -> None:
    counts = build_combined_binned_counts(module_test_context.data.sites_a, module_test_context.data.sites_b, "CallsetA", "CallsetB")

    assert not counts.empty
    assert int(counts["n_variants_CallsetA"].sum()) == len(module_test_context.data.sites_a)
    none_mask = counts["svtype"] == "CNV"
    none_mask &= counts["size_bucket"] == "2.5-10kb"
    none_mask &= counts["af_bucket"] == "10-50%"
    none_mask &= counts["genomic_context"] == "none"
    assert none_mask.any()
