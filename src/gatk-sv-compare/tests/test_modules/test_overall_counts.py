from __future__ import annotations

import numpy as np
import pandas as pd

from gatk_sv_compare.modules.overall_counts import OverallCountsModule, _af_bin_edges


def test_af_bin_edges_use_discrete_low_ac_bins() -> None:
    sites = pd.DataFrame(
        {
            "af": [0.01, 0.02, 0.03, 0.2, 0.5],
            "ac": [1, 2, 3, 21, 50],
            "svtype": ["DEL", "DEL", "DUP", "DEL", "DUP"],
        }
    )

    edges = _af_bin_edges(sites, low_ac_max=20, high_bin_count=4)

    expected_prefix = np.array(
        [
            0.01 / np.sqrt(2.0),
            np.sqrt(0.01 * 0.02),
            np.sqrt(0.02 * 0.03),
            np.sqrt(0.03 * 0.2),
        ]
    )
    np.testing.assert_allclose(edges[:4], expected_prefix)
    assert np.all(np.diff(edges) > 0)
    assert edges[-1] == 1.0


def test_overall_counts_module_writes_expected_plots(module_test_context) -> None:
    module = OverallCountsModule()

    module.run(module_test_context.data, module_test_context.config)

    output_dir = module_test_context.config.output_dir / "overall_counts"
    expected = [
        "sv_count_by_type.CallsetA.png",
        "sv_count_by_type.CallsetB.png",
        "size_distribution.CallsetA.png",
        "size_distribution.CallsetB.png",
        "af_distribution.CallsetA.png",
        "af_distribution.CallsetB.png",
        "genomic_context_by_type.CallsetA.png",
        "genomic_context_by_type.CallsetB.png",
    ]
    for name in expected:
        assert (output_dir / name).exists()
