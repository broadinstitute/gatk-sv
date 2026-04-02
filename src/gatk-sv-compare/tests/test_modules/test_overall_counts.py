from __future__ import annotations

from gatk_sv_compare.modules.overall_counts import OverallCountsModule


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
