from __future__ import annotations

import pandas as pd
import pytest

from gatk_sv_compare.modules.genotype_quality import GenotypeQualityModule, summarize_gq


def test_summarize_gq_collects_alt_genotype_statistics(module_test_context) -> None:
    summary = summarize_gq(module_test_context.config.vcf_a_path, pass_only=True)

    overall = summary.loc[summary["group"] == "overall"].iloc[0]
    assert int(overall["n"]) == 3
    assert overall["median_gq"] == pytest.approx(65.0)


def test_genotype_quality_module_writes_outputs(module_test_context) -> None:
    module = GenotypeQualityModule()

    module.run(module_test_context.data, module_test_context.config)

    output_dir = module_test_context.config.output_dir / "genotype_quality"
    summary_a = pd.read_csv(output_dir / "tables" / "gq_summary.CallsetA.tsv", sep="\t")
    summary_b = pd.read_csv(output_dir / "tables" / "gq_summary.CallsetB.tsv", sep="\t")
    assert not summary_a.empty
    assert not summary_b.empty
    assert (output_dir / "gq_histogram.overall.CallsetA.png").exists()
    assert (output_dir / "gq_histogram.overall.CallsetB.png").exists()
