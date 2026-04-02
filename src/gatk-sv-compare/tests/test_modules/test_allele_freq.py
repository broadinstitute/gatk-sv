from __future__ import annotations

import pandas as pd
import pytest

from gatk_sv_compare.modules.allele_freq import AlleleFreqModule, build_af_correlation_table


def test_build_af_correlation_table_uses_matched_pairs(module_test_context) -> None:
    table = build_af_correlation_table(module_test_context.data)

    overall = table.loc[table["group"] == "overall"].iloc[0]
    assert int(overall["n_matched"]) == 2
    assert overall["mean_abs_diff"] == pytest.approx(0.0833, abs=1e-3)


def test_allele_freq_module_writes_outputs(module_test_context) -> None:
    module = AlleleFreqModule()

    module.run(module_test_context.data, module_test_context.config)

    output_dir = module_test_context.config.output_dir / "allele_freq"
    stats = pd.read_csv(output_dir / "tables" / "af_correlation_stats.tsv", sep="\t")
    assert set(stats["group"]) >= {"overall", "DEL", "INS:MEI"}
    assert (output_dir / "af_correlation.overall.png").exists()
    assert (output_dir / "af_correlation.by_type.png").exists()
