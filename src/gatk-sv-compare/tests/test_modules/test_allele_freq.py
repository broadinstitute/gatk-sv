from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from gatk_sv_compare.modules.allele_freq import AlleleFreqModule, _af_bucket_axis_limits, _af_bucket_axis_ticks, build_af_correlation_table
from gatk_sv_compare.plot_utils import plot_scatter_af


def test_build_af_correlation_table_uses_matched_pairs(module_test_context) -> None:
    table = build_af_correlation_table(module_test_context.data)

    overall = table.loc[table["group"] == "overall"].iloc[0]
    assert int(overall["n_matched"]) == 2
    assert overall["mean_abs_diff"] == pytest.approx(0.0833, abs=1e-3)
    assert any(group.startswith("algorithm:") for group in table["group"])


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
    assert (output_dir / "af_correlation.by_af.png").exists()
    assert (output_dir / "af_correlation.by_context.png").exists()
    assert (output_dir / "af_correlation.by_algorithm.png").exists()


def test_plot_scatter_af_trend_stays_on_diagonal_for_perfect_data() -> None:
    x_values = np.linspace(0.0, 1.0, num=50)
    y_values = x_values.copy()
    fig, ax = plt.subplots()

    plot_scatter_af(ax, x_values, y_values, "A", "B")

    trend_line = ax.lines[-1]
    assert np.allclose(trend_line.get_xdata(), trend_line.get_ydata(), atol=1e-6)
    plt.close(fig)


def test_af_bucket_axis_limits_use_bucket_bounds_when_available() -> None:
    frame = pd.DataFrame({"af_a": [0.02, 0.08], "af_b": [0.03, 0.09]})

    limits = _af_bucket_axis_limits("1-10%", frame)

    assert limits == pytest.approx((0.01, 0.10))


def test_af_bucket_axis_limits_fall_back_to_data_window_for_ac1() -> None:
    frame = pd.DataFrame({"af_a": [0.11, 0.17], "af_b": [0.12, 0.18]})

    lower, upper = _af_bucket_axis_limits("AC=1", frame)

    assert lower == pytest.approx(0.11)
    assert upper == pytest.approx(0.18)


def test_af_bucket_axis_ticks_span_panel_limits() -> None:
    ticks = _af_bucket_axis_ticks(0.01, 0.10)

    np.testing.assert_allclose(ticks, np.asarray([0.01, 0.055, 0.10]))
