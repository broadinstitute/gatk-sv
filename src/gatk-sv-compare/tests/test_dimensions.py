from __future__ import annotations

from gatk_sv_compare.dimensions import categorize_variant, is_filtered_pass, normalize_svtype


def test_normalize_svtype_maps_mei_insertions() -> None:
    assert normalize_svtype("INS", "<INS:ME:ALU>") == "INS:MEI"
    assert normalize_svtype("INS", "<INS>") == "INS"


def test_is_filtered_pass_includes_multiallelic() -> None:
    assert is_filtered_pass({"PASS"}) is True
    assert is_filtered_pass({"MULTIALLELIC"}) is True
    assert is_filtered_pass({"HIGH_NCR"}) is False


def test_categorize_variant_uses_cnv_frequency_and_context() -> None:
    category = categorize_variant(
        {
            "svtype": "CNV",
            "alt_allele": "<CNV>",
            "svlen": 5000,
            "af": None,
            "ac": 0,
            "cn_nonref_freq": 0.02,
            "genomic_context": "segdup",
        }
    )
    assert category.svtype == "CNV"
    assert category.af_bucket == "1-10%"
    assert category.size_bucket == "2.5-10kb"
    assert category.genomic_context == "segdup"
