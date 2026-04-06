from __future__ import annotations

import numpy as np

from gatk_sv_compare.modules.size_signatures import (
    SizeSignaturesModule,
    flag_implausible_variants,
    quantify_retrotransposon_peaks,
    summarize_mei_subtypes,
)


def test_size_signature_helpers_return_expected_summaries(module_test_context) -> None:
    peaks = quantify_retrotransposon_peaks(np.array([280.0, 300.0, 320.0, 6000.0, 6100.0, 6200.0]))
    implausible = flag_implausible_variants(module_test_context.data.sites_a, {"chr1": 10000}, max_chrom_fraction=0.2)
    mei = summarize_mei_subtypes(module_test_context.data.sites_a)

    assert int(peaks["n_insertions"]) == 6
    assert not implausible.empty
    assert set(mei["subtype"]) == {"<INS:ME:ALU>"}


def test_size_signatures_module_writes_outputs(module_test_context) -> None:
    module = SizeSignaturesModule()

    module.run(module_test_context.data, module_test_context.config)

    output_dir = module_test_context.config.output_dir / "size_signatures"
    assert (output_dir / "implausible_variants.CallsetA.tsv.gz").exists()
    assert (output_dir / "mei_subtype_summary.CallsetB.tsv.gz").exists()
    assert not (output_dir / "tables").exists()
    assert not (output_dir / "retrotransposon_peaks.tsv.gz").exists()
    assert not (output_dir / "size_distribution_comparison.tsv.gz").exists()
    assert not (output_dir / "retrotransposon_peaks.CallsetA.png").exists()
    assert not (output_dir / "ins_size_overlay.png").exists()
