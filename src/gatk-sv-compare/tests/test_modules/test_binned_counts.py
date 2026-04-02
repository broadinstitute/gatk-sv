from __future__ import annotations

import pandas as pd

from gatk_sv_compare.modules.binned_counts import BinnedCountsModule, summarize_binned_counts


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
    counts_a = pd.read_csv(output_dir / "counts_a.tsv", sep="\t")
    counts_b = pd.read_csv(output_dir / "counts_b.tsv", sep="\t")
    assert (output_dir / "counts_a.parquet").exists()
    assert (output_dir / "counts_b.parquet").exists()
    assert not counts_a.empty
    assert not counts_b.empty
