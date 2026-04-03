from __future__ import annotations

import numpy as np
import pandas as pd

from gatk_sv_compare.modules.overall_counts import OverallCountsModule, _af_bin_edges, _normalized_histogram


def test_af_bin_edges_use_discrete_low_ac_bins() -> None:
    sites = pd.DataFrame(
        {
            "af": [0.01, 0.02, 0.03, 0.2, 0.5],
            "ac": [1, 2, 3, 21, 50],
            "an": [100, 100, 100, 100, 100],
            "svtype": ["DEL", "DEL", "DUP", "DEL", "DUP"],
        }
    )

    edges = _af_bin_edges(sites, low_ac_max=20, high_bin_count=4)

    expected_prefix = np.asarray([0.005, 0.015, 0.025, 0.035, 0.045], dtype=float)
    np.testing.assert_allclose(edges[:5], expected_prefix)
    np.testing.assert_allclose(edges[19:22], np.asarray([0.195, 0.205, 0.40375], dtype=float))
    assert np.all(np.diff(edges) > 0)
    assert edges[-1] == 1.0


def test_overall_counts_module_writes_expected_plots(module_test_context) -> None:
    module = OverallCountsModule()

    module.run(module_test_context.data, module_test_context.config)

    output_dir = module_test_context.config.output_dir / "overall_counts"
    expected = [
        "sv_count_by_type.CallsetA.png",
        "sv_count_by_type.CallsetB.png",
        "sv_count_by_algorithm.CallsetA.png",
        "sv_count_by_algorithm.CallsetB.png",
        "size_distribution.CallsetA.png",
        "size_distribution.CallsetB.png",
        "size_distribution.by_algorithm.CallsetA.png",
        "size_distribution.by_algorithm.CallsetB.png",
        "af_distribution.CallsetA.png",
        "af_distribution.CallsetB.png",
        "af_distribution.by_algorithm.CallsetA.png",
        "af_distribution.by_algorithm.CallsetB.png",
        "genomic_context_by_type.CallsetA.png",
        "genomic_context_by_type.CallsetB.png",
        "genomic_context_by_algorithm.CallsetA.png",
        "genomic_context_by_algorithm.CallsetB.png",
    ]
    for name in expected:
        assert (output_dir / name).exists()


def test_normalized_histogram_returns_proportions() -> None:
    counts, edges = _normalized_histogram(np.asarray([0.1, 0.1, 0.3, 0.7], dtype=float), np.asarray([0.0, 0.2, 0.5, 1.0], dtype=float))

    np.testing.assert_allclose(counts, np.asarray([0.5, 0.25, 0.25]))
    np.testing.assert_allclose(edges, np.asarray([0.0, 0.2, 0.5, 1.0]))
