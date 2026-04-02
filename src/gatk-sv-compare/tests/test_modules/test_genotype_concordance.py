from __future__ import annotations

import pandas as pd
import pytest

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
        "a",
    )

    assert not rows.empty
    assert "var_ppv" in rows.columns
    assert set(rows["label"]) == {"CallsetA"}
    assert set(rows["source"]) == {"a"}

    summary = summarize_concordance_metrics(rows)
    assert not summary.empty
    assert "var_ppv" in set(summary["metric"])
    assert "mean_value_a" in summary.columns


def test_genotype_concordance_module_writes_outputs(module_test_context) -> None:
    module = GenotypeConcordanceModule()

    module.run(module_test_context.data, module_test_context.config)

    output_dir = module_test_context.config.output_dir / "genotype_concordance"
    table = pd.read_csv(output_dir / "tables" / "concordance_metrics.tsv.gz", sep="\t")
    assert not table.empty
    assert (output_dir / "concordance_metrics.by_type.png").exists()
    assert (output_dir / "exact_match.by_type.png").exists()
    assert (output_dir / "exact_match.by_context.png").exists()
    assert {"n_sites_CallsetA", "mean_value_CallsetA", "median_value_CallsetA", "n_sites_CallsetB", "mean_value_CallsetB", "median_value_CallsetB"}.issubset(table.columns)
    assert "none" in set(table["genomic_context"])


def test_summarize_concordance_metrics_uses_cnv_concordance_for_cnv() -> None:
    metrics = pd.DataFrame(
        [
            {
                "source": "a",
                "label": "CallsetA",
                "variant_id": "cnv1",
                "svtype": "CNV",
                "size_bucket": "2.5-10kb",
                "af_bucket": "1-10%",
                "genomic_context": "none",
                "genotype_concordance": 1.0,
                "non_ref_genotype_concordance": 1.0,
                "var_ppv": 1.0,
                "var_sensitivity": 1.0,
                "cnv_concordance": 0.42,
            }
        ]
    )

    summary = summarize_concordance_metrics(metrics)
    row_mask = summary["metric"] == "genotype_concordance"
    row_mask &= summary["svtype"] == "CNV"
    row_mask &= summary["genomic_context"] == "none"
    row = summary.loc[row_mask].iloc[0]

    assert row["mean_value_a"] == pytest.approx(0.42)
    assert not ((summary["svtype"] == "CNV") & (summary["metric"] == "var_ppv")).any()
