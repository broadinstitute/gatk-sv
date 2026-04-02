from __future__ import annotations

import numpy as np
import pandas as pd

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
    peaks = pd.read_csv(output_dir / "tables" / "retrotransposon_peaks.tsv", sep="\t")
    ks_table = pd.read_csv(output_dir / "tables" / "size_distribution_comparison.tsv", sep="\t")
    assert set(peaks["label"]) == {"CallsetA", "CallsetB"}
    assert ks_table.loc[0, "svtype"] == "INS+INS:MEI"
    assert (output_dir / "tables" / "implausible_variants.CallsetA.tsv").exists()
    assert (output_dir / "tables" / "mei_subtype_summary.CallsetB.tsv").exists()
    assert (output_dir / "retrotransposon_peaks.CallsetA.png").exists()
    assert (output_dir / "ins_size_overlay.png").exists()
