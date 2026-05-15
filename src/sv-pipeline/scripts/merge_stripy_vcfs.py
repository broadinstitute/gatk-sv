#!/bin/env python

import argparse
from contextlib import ExitStack
import pysam
import sys
from typing import Dict, Iterable, List, Optional, Sequence, Text, Tuple


INFO_FIELDS = ["RU", "PERIOD", "LOCUS", "SVTYPE", "DISEASES"]
FORMAT_FIELDS = ["REPCN", "REPCI1", "REPCI2", "OUTLIER", "ZSCORE", "DP", "STR_FILTER"]
StripyRecordKey = Tuple[Text, int, int, Optional[Text], Tuple[Text, ...]]


def _add_header_line_if_absent(header: pysam.VariantHeader,
                               header_collection: Iterable[Text],
                               key: Text,
                               line: Text) -> None:
    if key not in header_collection:
        header.add_line(line)


def update_header(header: pysam.VariantHeader) -> None:
    """
    Adds needed STRipy header lines. These should correspond to the FORMAT_FIELDS and INFO_FIELDS variables.
    """
    _add_header_line_if_absent(header, header.formats, "GT", '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    _add_header_line_if_absent(header, header.info, "END", '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">')
    _add_header_line_if_absent(header, header.info, "RU", '##INFO=<ID=RU,Number=1,Type=String,Description="Repeat unit in the reference orientation">')
    _add_header_line_if_absent(header, header.info, "PERIOD", '##INFO=<ID=PERIOD,Number=1,Type=Integer,Description="Length of the repeat unit">')
    _add_header_line_if_absent(header, header.info, "DISEASES", '##INFO=<ID=DISEASES,Number=.,Type=String,Description="Associated disease symbols for this STR locus (| separated)">')
    _add_header_line_if_absent(header, header.info, "LOCUS", '##INFO=<ID=LOCUS,Number=1,Type=String,Description="Gene/locus identifier from STRipy">')
    _add_header_line_if_absent(header, header.formats, "REPCN", '##FORMAT=<ID=REPCN,Number=2,Type=Float,Description="Number of repeat units spanned by each allele">')
    _add_header_line_if_absent(header, header.formats, "REPCI1", '##FORMAT=<ID=REPCI1,Number=2,Type=Integer,Description="95% CI min,max on repeat counts of first allele">')
    _add_header_line_if_absent(header, header.formats, "REPCI2", '##FORMAT=<ID=REPCI2,Number=2,Type=Integer,Description="95% CI min,max on repeat counts of second allele">')
    _add_header_line_if_absent(header, header.formats, "OUTLIER", '##FORMAT=<ID=OUTLIER,Number=2,Type=Integer,Description="Allelic population outlier flags (0/1) assigned by STRipy">')
    _add_header_line_if_absent(header, header.formats, "ZSCORE", '##FORMAT=<ID=ZSCORE,Number=2,Type=Float,Description="Allelic population Z-scores assigned by STRipy">')
    _add_header_line_if_absent(header, header.formats, "DP", '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total depth at STR site">')
    _add_header_line_if_absent(header, header.formats, "STR_FILTER", '##FORMAT=<ID=STR_FILTER,Number=.,Type=String,Description="Filter status assigned by STRipy">')


def _copy_stripy_header_metadata(output_header: pysam.VariantHeader,
                                 stripy_header: pysam.VariantHeader) -> None:
    for contig in stripy_header.contigs:
        if contig not in output_header.contigs:
            contig_length = stripy_header.contigs[contig].length
            if contig_length is None:
                output_header.contigs.add(contig)
            else:
                output_header.contigs.add(contig, length=contig_length)

    for filter_name in stripy_header.filters:
        if filter_name != "PASS" and filter_name not in output_header.filters:
            filter_description = stripy_header.filters[filter_name].description or "STRipy filter"
            output_header.add_meta("FILTER", items=[("ID", filter_name), ("Description", filter_description)])


def build_output_header(main_header: Optional[pysam.VariantHeader],
                        stripy_vcfs: Sequence[pysam.VariantFile],
                        selected_samples_by_vcf: Sequence[Sequence[Text]]) -> pysam.VariantHeader:
    if main_header is None:
        output_header = stripy_vcfs[0].header.copy()
    else:
        output_header = main_header.copy()

    for stripy_vcf, selected_samples in zip(stripy_vcfs, selected_samples_by_vcf):
        _copy_stripy_header_metadata(output_header=output_header, stripy_header=stripy_vcf.header)
        for sample in selected_samples:
            if sample not in output_header.samples:
                output_header.add_sample(sample)

    update_header(header=output_header)
    return output_header


def validate_sample_ids(main_header: Optional[pysam.VariantHeader],
                        stripy_vcfs: Sequence[pysam.VariantFile]) -> List[List[Text]]:
    selected_samples_by_vcf: List[List[Text]] = []
    seen_samples = set()
    for stripy_vcf in stripy_vcfs:
        if len(stripy_vcf.header.samples) == 0:
            raise ValueError("Expected at least 1 sample in STRipy header but got 0")
        selected_samples: List[Text] = []
        for sample in stripy_vcf.header.samples:
            if main_header is not None and sample not in main_header.samples:
                continue
            if sample in seen_samples:
                raise ValueError(f"Sample {sample} appears in multiple STRipy VCFs")
            seen_samples.add(sample)
            selected_samples.append(sample)
        selected_samples_by_vcf.append(selected_samples)
    return selected_samples_by_vcf


def _build_stripy_record_key(stripy_record: pysam.VariantRecord) -> StripyRecordKey:
    return (
        stripy_record.contig,
        stripy_record.start,
        stripy_record.stop,
        stripy_record.id,
        tuple(stripy_record.alleles or ()),
    )


def _new_stripy_output_record(stripy_record: pysam.VariantRecord,
                              out_vcf: pysam.VariantFile) -> pysam.VariantRecord:
    record = out_vcf.new_record(contig=stripy_record.contig,
                                start=stripy_record.start,
                                stop=stripy_record.stop,
                                alleles=stripy_record.alleles,
                                id=stripy_record.id,
                                qual=stripy_record.qual,
                                filter=list(stripy_record.filter.keys()))
    for key in INFO_FIELDS:
        if key in stripy_record.info:
            record.info[key] = stripy_record.info[key]
    return record


def _get_stripy_site_filter_values(stripy_record: pysam.VariantRecord) -> List[Text]:
    filter_values = list(stripy_record.filter.keys())
    return filter_values or ["PASS"]


def _normalize_str_filter_values(filter_values: Optional[Sequence[Text]]) -> List[Text]:
    if filter_values is None:
        return []
    if isinstance(filter_values, str):
        return [filter_values]
    return [str(value) for value in filter_values if value is not None]


def _merge_stripy_site_filters(existing_record: pysam.VariantRecord,
                               stripy_record: pysam.VariantRecord) -> None:
    merged_filters = [filter_name for filter_name in existing_record.filter.keys() if filter_name != "PASS"]
    for filter_name in _get_stripy_site_filter_values(stripy_record):
        if filter_name != "PASS" and filter_name not in merged_filters:
            merged_filters.append(filter_name)

    existing_record.filter.clear()
    if merged_filters:
        for filter_name in merged_filters:
            existing_record.filter.add(filter_name)
    else:
        existing_record.filter.add("PASS")


def _get_sample_str_filter_values(stripy_record: pysam.VariantRecord,
                                  sample: Text) -> List[Text]:
    stripy_sample = stripy_record.samples[sample]
    if "STR_FILTER" in stripy_sample:
        str_filter_values = _normalize_str_filter_values(stripy_sample["STR_FILTER"])
        if str_filter_values:
            return str_filter_values
    return _get_stripy_site_filter_values(stripy_record)


def _validate_stripy_site_metadata(existing_record: pysam.VariantRecord,
                                   stripy_record: pysam.VariantRecord) -> None:
    for key in INFO_FIELDS:
        existing_value = existing_record.info[key] if key in existing_record.info else None
        incoming_value = stripy_record.info[key] if key in stripy_record.info else None
        if existing_value != incoming_value:
            raise ValueError(
                f"Conflicting STRipy site metadata for {stripy_record.contig}:{stripy_record.pos}-{stripy_record.stop} "
                f"({stripy_record.id}): {key} differs across inputs"
            )


def _copy_stripy_sample_fields(record: pysam.VariantRecord,
                               stripy_record: pysam.VariantRecord,
                               samples: Sequence[Text]) -> None:
    for sample in samples:
        record.samples[sample]['GT'] = (None, None)
        stripy_sample = stripy_record.samples[sample]
        for key in FORMAT_FIELDS:
            if key == "STR_FILTER":
                continue
            if key in stripy_sample:
                record.samples[sample][key] = stripy_sample[key]
        record.samples[sample]["STR_FILTER"] = _get_sample_str_filter_values(stripy_record=stripy_record, sample=sample)


def _write_stripy_records(stripy_inputs: Sequence[pysam.VariantFile],
                          selected_samples_by_vcf: Sequence[Sequence[Text]],
                          out_vcf: pysam.VariantFile) -> None:
    merged_records: Dict[StripyRecordKey, pysam.VariantRecord] = {}
    for stripy_vcf, samples in zip(stripy_inputs, selected_samples_by_vcf):
        if not samples:
            continue
        for stripy_record in stripy_vcf:
            key = _build_stripy_record_key(stripy_record)
            if key not in merged_records:
                merged_records[key] = _new_stripy_output_record(stripy_record=stripy_record, out_vcf=out_vcf)
            else:
                _validate_stripy_site_metadata(existing_record=merged_records[key], stripy_record=stripy_record)
                _merge_stripy_site_filters(existing_record=merged_records[key], stripy_record=stripy_record)
            _copy_stripy_sample_fields(record=merged_records[key], stripy_record=stripy_record, samples=samples)

    for record in merged_records.values():
        out_vcf.write(record)


def _process(main_vcf: Optional[pysam.VariantFile],
             stripy_inputs: Sequence[pysam.VariantFile],
             selected_samples_by_vcf: Sequence[Sequence[Text]],
             out_vcf: pysam.VariantFile) -> None:
    """
    Copies the main VCF records when provided and appends STRipy records merged by site across inputs.
    """
    if main_vcf is not None:
        for record in main_vcf:
            translated_record = record.copy()
            translated_record.translate(out_vcf.header)
            out_vcf.write(translated_record)
    _write_stripy_records(
        stripy_inputs=stripy_inputs,
        selected_samples_by_vcf=selected_samples_by_vcf,
        out_vcf=out_vcf,
    )


def _read_stripy_vcf_paths(file_list_path: Text) -> List[Text]:
    with open(file_list_path, "r", encoding="utf-8") as file_list:
        stripy_vcf_paths = [line.strip() for line in file_list if line.strip()]
    if not stripy_vcf_paths:
        raise ValueError(f"No STRipy VCFs found in list file: {file_list_path}")
    return stripy_vcf_paths


def _parse_arguments(argv: List[Text]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merges one or more single- or multi-sample STRipy VCFs into a GATK-SV call set, or into a STRipy-only VCF",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--main-vcf", type=str,
                        help="Optional GATK-SV VCF; when provided, STRipy samples absent from this header are ignored")
    parser.add_argument("--stripy-vcfs-list", type=str, required=True,
                        help="Text file listing one STRipy VCF path per line")
    parser.add_argument("--out", type=str, required=True,
                        help="Output VCF (unsorted)")
    if len(argv) <= 1:
        parser.parse_args(["--help"])
        sys.exit(0)
    parsed_arguments = parser.parse_args(argv[1:])
    return parsed_arguments


def main(argv: Optional[List[Text]] = None):
    if argv is None:
        argv = sys.argv
    arguments = _parse_arguments(argv)
    stripy_vcf_paths = _read_stripy_vcf_paths(arguments.stripy_vcfs_list)

    with ExitStack() as exit_stack:
        main_vcf = exit_stack.enter_context(pysam.VariantFile(arguments.main_vcf)) if arguments.main_vcf else None
        stripy_vcfs = [exit_stack.enter_context(pysam.VariantFile(path)) for path in stripy_vcf_paths]
        main_header = main_vcf.header if main_vcf is not None else None
        selected_samples_by_vcf = validate_sample_ids(main_header=main_header, stripy_vcfs=stripy_vcfs)
        output_header = build_output_header(
            main_header=main_header,
            stripy_vcfs=stripy_vcfs,
            selected_samples_by_vcf=selected_samples_by_vcf,
        )
        with pysam.VariantFile(arguments.out, mode='w', header=output_header) as out_vcf:
            _process(
                main_vcf=main_vcf,
                stripy_inputs=stripy_vcfs,
                selected_samples_by_vcf=selected_samples_by_vcf,
                out_vcf=out_vcf,
            )


if __name__ == "__main__":
    main()
