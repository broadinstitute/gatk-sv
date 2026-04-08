from __future__ import annotations

from dataclasses import replace

import matplotlib.pyplot as plt
import pandas as pd
import pytest
from matplotlib.figure import Figure

from gatk_sv_compare.modules.allele_freq import AlleleFreqModule, _plot_grouped_af_correlation, build_af_correlation_table
from gatk_sv_compare.plot_utils import plot_scatter_af


def test_build_af_correlation_table_uses_matched_pairs(module_test_context) -> None:
    table = build_af_correlation_table(module_test_context.data)

    overall = table.loc[table["group"] == "overall"].iloc[0]
    assert int(overall["n_matched"]) == 2
    assert overall["mean_abs_diff"] == pytest.approx(0.0833, abs=1e-3)
    assert any(group.startswith("algorithm:") for group in table["group"])
    assert any(group.startswith("evidence:") for group in table["group"])


def test_build_af_correlation_table_applies_pass_only_to_both_sides(module_test_context) -> None:
    sites_a = module_test_context.data.sites_a.copy()
    sites_b = module_test_context.data.sites_b.copy()
    matched_pairs = module_test_context.data.matched_pairs.copy()

    sites_a.loc[sites_a["variant_id"] == "a_del", ["svtype", "in_filtered_pass_view"]] = ["BND", False]
    sites_b.loc[sites_b["variant_id"] == "b_del", ["svtype", "in_filtered_pass_view"]] = ["BND", False]
    matched_pairs.loc[matched_pairs["variant_id_a"] == "a_del", "svtype_a"] = "BND"
    matched_pairs.loc[matched_pairs["variant_id_b"] == "b_del", "svtype_b"] = "BND"

    filtered_data = replace(module_test_context.data, sites_a=sites_a, sites_b=sites_b, matched_pairs=matched_pairs)

    table = build_af_correlation_table(filtered_data, pass_only=True)

    assert "BND" not in set(table["group"])
    overall = table.loc[table["group"] == "overall"].iloc[0]
    assert int(overall["n_matched"]) == 1


def test_allele_freq_module_writes_outputs(module_test_context) -> None:
    module = AlleleFreqModule()

    module.run(module_test_context.data, module_test_context.config)

    output_dir = module_test_context.config.output_dir / "allele_freq"
    stats = pd.read_csv(output_dir / "tables" / "af_correlation_stats.tsv.gz", sep="\t")
    assert set(stats["group"]) >= {"overall", "DEL", "INS:MEI"}
    assert {"pearson_r", "pearson_p", "spearman_rho", "spearman_p"}.issubset(stats.columns)
    assert (output_dir / "af_correlation.overall.png").exists()
    assert (output_dir / "af_correlation.by_type.png").exists()
    assert (output_dir / "af_correlation.by_size.png").exists()
    assert (output_dir / "af_correlation.by_context.png").exists()
    assert (output_dir / "af_correlation.by_evidence.png").exists()
    assert (output_dir / "af_correlation.by_algorithm.png").exists()


def test_plot_scatter_af_trend_stays_on_diagonal_for_perfect_data() -> None:
    x_values = pd.Series([index / 49.0 for index in range(50)], dtype=float)
    y_values = x_values.copy()
    fig, ax = plt.subplots()

    plot_scatter_af(ax, x_values, y_values, "A", "B")

    trend_line = ax.lines[-1]
    assert trend_line.get_xdata() == pytest.approx(trend_line.get_ydata(), abs=1e-6)
    plt.close(fig)


def test_plot_grouped_af_correlation_includes_ctx_svtype(tmp_path, monkeypatch) -> None:
    matched_pairs = pd.DataFrame(
        [
            {
                "svtype_a": "CTX",
                "af_a": 0.10,
                "af_b": 0.15,
                "size_bucket_a": "N/A",
                "genomic_context_a": "none",
                "evidence_bucket_a": "PE,SR",
                "algorithms_a": "manta",
            },
            {
                "svtype_a": "DEL",
                "af_a": 0.20,
                "af_b": 0.18,
                "size_bucket_a": "100-500bp",
                "genomic_context_a": "none",
                "evidence_bucket_a": "RD,PE",
                "algorithms_a": "manta",
            },
        ]
    )
    captured = {}
    original_savefig = Figure.savefig

    def spy_savefig(self, *args, **kwargs):
        captured["titles"] = [axis.get_title() for axis in self.axes]
        return original_savefig(self, *args, **kwargs)

    monkeypatch.setattr(Figure, "savefig", spy_savefig)

    _plot_grouped_af_correlation(matched_pairs, "svtype_a", tmp_path / "af_by_type.png", "CallsetA", "CallsetB")

    assert "CTX" in captured["titles"]
