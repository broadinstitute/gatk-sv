from __future__ import annotations

import pandas as pd

from gatk_sv_compare.dimensions import categorize_variant, explode_algorithm_buckets, is_filtered_pass, normalize_algorithms, normalize_evidence_bucket, normalize_svtype, ordered_algorithms, ordered_contexts, ordered_evidence_buckets, ordered_plot_af_buckets, ordered_plot_evidence_buckets, ordered_plot_size_buckets


def test_normalize_svtype_maps_mei_insertions() -> None:
    assert normalize_svtype("INS", "<INS:ME:ALU>") == "INS:MEI"
    assert normalize_svtype("INS", "<INS>") == "INS"


def test_is_filtered_pass_accepts_pass_or_empty_filter() -> None:
    assert is_filtered_pass({"PASS"}) is True
    assert is_filtered_pass(set()) is True
    assert is_filtered_pass({"MULTIALLELIC"}) is False
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


def test_ordered_contexts_always_includes_none_bucket() -> None:
    assert ordered_contexts(["segdup"]) == ["simple_repeat", "segdup", "repeatmasker", "none"]


def test_ordered_plot_af_buckets_excludes_unknown() -> None:
    assert ordered_plot_af_buckets(["1-10%", "unknown", "AC=1"]) == ["AC=1", "1-10%"]


def test_ordered_plot_size_buckets_excludes_unknown_and_na() -> None:
    assert ordered_plot_size_buckets(["2.5-10kb", "unknown", "N/A", "<100bp"]) == ["<100bp", "2.5-10kb"]


def test_normalize_algorithms_and_ordered_algorithms_handle_multi_value_inputs() -> None:
    assert normalize_algorithms(("manta", "wham", "manta")) == ("manta", "wham")
    assert ordered_algorithms(["wham", "unknown", "manta"]) == ["manta", "wham", "unknown"]


def test_explode_algorithm_buckets_duplicates_multi_algorithm_rows() -> None:
    frame = pd.DataFrame([{"variant_id": "v1", "algorithms": "manta,wham"}, {"variant_id": "v2", "algorithms": "melt"}])

    expanded = explode_algorithm_buckets(frame)

    assert expanded["variant_id"].tolist() == ["v1", "v1", "v2"]
    assert expanded["algorithm"].tolist() == ["manta", "wham", "melt"]


def test_normalize_and_order_evidence_buckets() -> None:
    assert normalize_evidence_bucket(("SR", "RD", "BAF", "PE")) == "RD,BAF,PE,SR"
    assert normalize_evidence_bucket("BAF,SR") == "BAF,SR"
    assert normalize_evidence_bucket("BAF") == "BAF"
    assert ordered_evidence_buckets(["unknown", "PE,SR", "RD", "RD,PE"]) == ["RD", "RD,PE", "PE,SR", "unknown"]
    assert ordered_plot_evidence_buckets(["unknown", "PE,SR", "RD"]) == ["RD", "PE,SR"]
