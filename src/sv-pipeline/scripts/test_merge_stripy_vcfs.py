import importlib.util
from pathlib import Path

import pysam
import pytest


MODULE_PATH = Path(__file__).with_name("merge_stripy_vcfs.py")
MODULE_SPEC = importlib.util.spec_from_file_location("merge_stripy_vcfs", MODULE_PATH)
assert MODULE_SPEC is not None
assert MODULE_SPEC.loader is not None
merge_stripy_vcfs = importlib.util.module_from_spec(MODULE_SPEC)
MODULE_SPEC.loader.exec_module(merge_stripy_vcfs)


def _write_main_vcf(path: Path, samples):
    header = pysam.VariantHeader()
    header.contigs.add("chr1", length=1000)
    header.add_line('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Variant type">')
    header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    for sample in samples:
        header.add_sample(sample)

    with pysam.VariantFile(path, "w", header=header) as vcf:
        record = vcf.new_record(contig="chr1", start=9, stop=10, alleles=("A", "C"), id="main1")
        record.samples[samples[0]]["GT"] = (0, 1)
        vcf.write(record)


def _write_stripy_vcf(
        path: Path,
        samples,
        site_id: str,
        start: int,
        dp_by_sample,
        filter_ids=None,
        str_filter_by_sample=None,
):
    filter_ids = list(filter_ids or ["PASS"])
    header = pysam.VariantHeader()
    header.contigs.add("chr1", length=1000)
    header.add_line('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Variant type">')
    merge_stripy_vcfs.update_header(header)
    for filter_id in filter_ids:
        if filter_id != "PASS":
            header.add_meta("FILTER", items=[("ID", filter_id), ("Description", f"Test filter {filter_id}")])
    for sample in samples:
        header.add_sample(sample)

    with pysam.VariantFile(path, "w", header=header) as vcf:
        record = vcf.new_record(contig="chr1", start=start, stop=start + 10, alleles=("N", "<STR>"), id=site_id)
        record.filter.clear()
        for filter_id in filter_ids:
            record.filter.add(filter_id)
        record.info["RU"] = "CAG"
        record.info["PERIOD"] = 3
        record.info["LOCUS"] = "ATXN1"
        record.info["SVTYPE"] = "STR"
        record.info["DISEASES"] = ("SCA1",)
        for sample, depth in dp_by_sample.items():
            record.samples[sample]["GT"] = (None, None)
            record.samples[sample]["DP"] = depth
            record.samples[sample]["REPCN"] = (10.0, 11.0)
            if str_filter_by_sample is not None and sample in str_filter_by_sample:
                record.samples[sample]["STR_FILTER"] = str_filter_by_sample[sample]
        vcf.write(record)


def test_process_ignores_stripy_samples_missing_from_main_header(tmp_path):
    main_vcf_path = tmp_path / "main.vcf"
    stripy_mixed_path = tmp_path / "stripy_mixed.vcf"
    stripy_excluded_path = tmp_path / "stripy_excluded.vcf"
    output_path = tmp_path / "merged.vcf"

    _write_main_vcf(main_vcf_path, ["sample_a"])
    _write_stripy_vcf(
        stripy_mixed_path,
        ["sample_a", "sample_b"],
        site_id="str_keep",
        start=99,
        dp_by_sample={"sample_a": 12, "sample_b": 40},
    )
    _write_stripy_vcf(
        stripy_excluded_path,
        ["sample_b"],
        site_id="str_drop",
        start=199,
        dp_by_sample={"sample_b": 25},
    )

    with pysam.VariantFile(main_vcf_path) as main_vcf, \
            pysam.VariantFile(stripy_mixed_path) as stripy_mixed, \
            pysam.VariantFile(stripy_excluded_path) as stripy_excluded:
        stripy_vcfs = [stripy_mixed, stripy_excluded]
        selected_samples_by_vcf = merge_stripy_vcfs.validate_sample_ids(
            main_header=main_vcf.header,
            stripy_vcfs=stripy_vcfs,
        )

        assert selected_samples_by_vcf == [["sample_a"], []]

        output_header = merge_stripy_vcfs.build_output_header(
            main_header=main_vcf.header,
            stripy_vcfs=stripy_vcfs,
            selected_samples_by_vcf=selected_samples_by_vcf,
        )
        with pysam.VariantFile(output_path, "w", header=output_header) as out_vcf:
            merge_stripy_vcfs._process(
                main_vcf=main_vcf,
                stripy_inputs=stripy_vcfs,
                selected_samples_by_vcf=selected_samples_by_vcf,
                out_vcf=out_vcf,
            )

    with pysam.VariantFile(output_path) as merged_vcf:
        assert "STR" in merged_vcf.header.alts
        assert list(merged_vcf.header.samples) == ["sample_a"]
        records = list(merged_vcf)

    assert [record.id for record in records] == ["main1", "str_keep"]
    assert records[1].samples["sample_a"]["DP"] == 12


def test_validate_sample_ids_still_rejects_duplicate_retained_samples(tmp_path):
    main_vcf_path = tmp_path / "main.vcf"
    stripy_one_path = tmp_path / "stripy_one.vcf"
    stripy_two_path = tmp_path / "stripy_two.vcf"

    _write_main_vcf(main_vcf_path, ["sample_a"])
    _write_stripy_vcf(
        stripy_one_path,
        ["sample_a", "sample_b"],
        site_id="str_one",
        start=99,
        dp_by_sample={"sample_a": 12, "sample_b": 40},
    )
    _write_stripy_vcf(
        stripy_two_path,
        ["sample_a"],
        site_id="str_two",
        start=199,
        dp_by_sample={"sample_a": 25},
    )

    with pysam.VariantFile(main_vcf_path) as main_vcf, \
            pysam.VariantFile(stripy_one_path) as stripy_one, \
            pysam.VariantFile(stripy_two_path) as stripy_two:
        with pytest.raises(ValueError, match="appears in multiple STRipy VCFs"):
            merge_stripy_vcfs.validate_sample_ids(
                main_header=main_vcf.header,
                stripy_vcfs=[stripy_one, stripy_two],
            )


def test_process_merges_site_filters_and_populates_str_filter(tmp_path):
    stripy_pass_path = tmp_path / "stripy_pass.vcf"
    stripy_filtered_path = tmp_path / "stripy_filtered.vcf"
    output_path = tmp_path / "merged.vcf"

    _write_stripy_vcf(
        stripy_pass_path,
        ["sample_a"],
        site_id="str_shared",
        start=99,
        dp_by_sample={"sample_a": 12},
    )
    _write_stripy_vcf(
        stripy_filtered_path,
        ["sample_b"],
        site_id="str_shared",
        start=99,
        dp_by_sample={"sample_b": 25},
        filter_ids=["LOW_DP"],
        str_filter_by_sample={"sample_b": ["Low depth (<10 reads)"]},
    )

    with pysam.VariantFile(stripy_pass_path) as stripy_pass, \
            pysam.VariantFile(stripy_filtered_path) as stripy_filtered:
        stripy_vcfs = [stripy_pass, stripy_filtered]
        selected_samples_by_vcf = merge_stripy_vcfs.validate_sample_ids(
            main_header=None,
            stripy_vcfs=stripy_vcfs,
        )

        output_header = merge_stripy_vcfs.build_output_header(
            main_header=None,
            stripy_vcfs=stripy_vcfs,
            selected_samples_by_vcf=selected_samples_by_vcf,
        )
        with pysam.VariantFile(output_path, "w", header=output_header) as out_vcf:
            merge_stripy_vcfs._process(
                main_vcf=None,
                stripy_inputs=stripy_vcfs,
                selected_samples_by_vcf=selected_samples_by_vcf,
                out_vcf=out_vcf,
            )

    with pysam.VariantFile(output_path) as merged_vcf:
        records = list(merged_vcf)

    assert len(records) == 1
    assert list(records[0].filter.keys()) == ["LOW_DP"]
    assert list(records[0].samples["sample_a"]["STR_FILTER"]) == ["PASS"]
    assert list(records[0].samples["sample_b"]["STR_FILTER"]) == ["Low depth (<10 reads)"]
