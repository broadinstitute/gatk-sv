from __future__ import annotations

import pandas as pd

from gatk_sv_compare.modules.genotype_concordance import (
    GenotypeConcordanceModule,
    summarize_concordance_metrics,
    _extract_concordance_rows,
)


def test_extract_concordance_rows_reads_available_info_metrics(module_test_context) -> None:
    rows = _extract_concordance_rows(
        module_test_context.config.vcf_a_path,
        module_test_context.data.sites_a,
        module_test_context.data.label_a,
    )

    assert not rows.empty
    assert "var_ppv" in rows.columns
    assert set(rows["label"]) == {"CallsetA"}

    summary = summarize_concordance_metrics(rows)
    assert not summary.empty
    assert "var_ppv" in set(summary["metric"])


def test_genotype_concordance_module_writes_outputs(module_test_context) -> None:
    module = GenotypeConcordanceModule()

    module.run(module_test_context.data, module_test_context.config)

    output_dir = module_test_context.config.output_dir / "genotype_concordance"
    table = pd.read_csv(output_dir / "tables" / "concordance_metrics.tsv", sep="\t")
    assert not table.empty
    assert (output_dir / "var_ppv.by_type.CallsetA.png").exists()
    assert (output_dir / "genotype_concordance.by_type.CallsetB.png").exists()
