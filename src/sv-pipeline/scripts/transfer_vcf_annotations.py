#!/bin/env python

import argparse
import pysam
import sys
from typing import Optional, List, Text, Set, Dict, Tuple


def create_header(header_in: pysam.VariantHeader,
                  header_ann: pysam.VariantHeader,
                  transfer_infos: Set[Text],
                  transfer_formats: Set[Text]) -> pysam.VariantHeader:
    """
    Adds annotation metadata lines to the input header. Also checks that the fields exist in the annotation vcf.

    Parameters
    ----------
    header_in: pysam.VariantHeader
        input header
    header_ann: pysam.VariantHeader
        annotation vcf header
    transfer_infos: Set[Text]
        set of info fields to annotate
    transfer_formats: Set[Text]
        set of format fields to annotate

    Returns
    -------
    header: pysam.VariantHeader
        output header
    """
    header = pysam.VariantHeader()
    for sample in header_in.samples:
        header.add_sample(sample)
    for line in header_in.records:
        header.add_line(str(line))

    # Check annotation fields exist
    for k in transfer_infos:
        if k not in header_ann.info:
            raise ValueError(f"INFO field {k} not found in annotation vcf header")
    for k in transfer_formats:
        if k not in header_ann.formats:
            raise ValueError(f"FORMAT field {k} not found in annotation vcf header")

    # Add annotation headerlines
    for line in header_ann.records:
        if len(line.attrs) > 0 and 'ID' in line.keys() and \
                (line['ID'] in transfer_infos or line['ID'] in transfer_formats):
            header.add_line(str(line))
    return header


def get_annotations(ann_vcf: pysam.VariantFile,
                    header_out: pysam.VariantHeader,
                    transfer_infos: Set[Text],
                    transfer_formats: Set[Text]) -> Tuple[Dict, Dict]:
    """
    Gathers info/format metadata from a source vcf

    Parameters
    ----------
    ann_vcf: pysam.VariantFile
        Source vcf with annotations
    header_out: pysam.VariantHeader
        Output vcf header
    transfer_infos: Set[Text]
        set of info fields to annotate
    transfer_formats: Set[Text]
        set of format fields to annotate

    Returns
    -------
    (dict, dict)
        pair of nested dicts with INFO and FORMAT fields: {vid: {info: value}}, {vid: {sample: {format: value}}}
    """
    sample_set = set(header_out.samples).intersection(set(ann_vcf.header.samples))
    if len(sample_set) == 0 and len(transfer_formats) > 0:
        raise ValueError(f"There are FORMAT fields to transfer but no samples in common between the vcfs")
    info_dict = dict()
    format_dict = dict()
    for record in ann_vcf:
        if record.id in info_dict:
            raise ValueError(f"Encountered duplicate variant id: {record.id}")
        info_dict[record.id] = {k: record.info.get(k, None) for k in transfer_infos}
        format_dict[record.id] = \
            {s: {k: record.samples[s].get(k, None) for k in transfer_formats} for s in record.samples if s in sample_set}
    return info_dict, format_dict


def annotate(record: pysam.VariantRecord,
             vcf_out: pysam.VariantFile,
             info_dict: Tuple[Dict, Dict],
             format_dict: Tuple[Dict, Dict]) -> pysam.VariantRecord:
    """
    Annotates record with INFO and FORMAT fields

    Parameters
    ----------
    record: pysam.VariantRecord
        input record
    vcf_out: pysam.VariantFile
        output vcf object
    info_dict: Set[Text]
        info field dict from get_annotations: {vid: {info: value}}
    format_dict: Set[Text]
        format field dict from get_annotations: {vid: {sample: {format: value}}}

    Returns
    -------
    header: pysam.VariantRecord
        annotated record
    """
    new_record = vcf_out.new_record(alleles=record.alleles, contig=record.contig, filter=record.filter, id=record.id,
                                    info=record.info, start=record.start, qual=record.qual, stop=record.stop)
    for sample, dest in new_record.samples.items():
        source = record.samples[sample]
        for key, val in source.items():
            dest[key] = val
    if record.id in info_dict:
        for key, val in info_dict[record.id].items():
            if val is not None:
                new_record.info[key] = val
    if record.id in format_dict:
        for sample_id, formats in format_dict[record.id].items():
            sample_data = new_record.samples[sample_id]
            for key, val in formats.items():
                if val is not None:
                    sample_data[key] = val
    return new_record


def __parse_arg_list(arg: Text) -> List[Text]:
    if arg is None:
        return list()
    else:
        return arg.split(',')


def __parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Annotate the given vcf with fields from another, matching by variant id",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("vcf", type=str,
                        help="VCF to annotate")
    parser.add_argument("--out", type=str, required=True,
                        help="Output VCF")
    parser.add_argument("--ann-vcf", type=str, required=True,
                        help="Vcf with source annotations")
    parser.add_argument("--formats", type=str,
                        help="Comma-delimited list of FORMAT fields to transfer (e.g. SR_GT,SR_GQ)")
    parser.add_argument("--infos", type=str,
                        help="Comma-delimited list of INFO fields to transfer (e.g. AC,AN,AF)")
    if len(argv) <= 1:
        parser.parse_args(["--help"])
        sys.exit(0)
    parsed_arguments = parser.parse_args(argv[1:])
    return parsed_arguments


def main(argv: Optional[List[Text]] = None):
    if argv is None:
        argv = sys.argv
    arguments = __parse_arguments(argv)
    transfer_formats = set(__parse_arg_list(arguments.formats))
    transfer_infos = set(__parse_arg_list(arguments.infos))
    if len(transfer_infos) + len(transfer_formats) == 0:
        raise ValueError('No INFO or FORMAT fields were specified')

    # convert vcf header and records
    with pysam.VariantFile(arguments.vcf) as vcf_in, pysam.VariantFile(arguments.ann_vcf) as ann_in:
        header = create_header(
            header_in=vcf_in.header,
            header_ann=ann_in.header,
            transfer_formats=transfer_formats,
            transfer_infos=transfer_infos
        )
        print("Retrieving annotations...")
        info_dict, format_dict = get_annotations(ann_vcf=ann_in,
                                                 header_out=header,
                                                 transfer_infos=transfer_infos,
                                                 transfer_formats=transfer_formats)
        with pysam.VariantFile(arguments.out, mode='w', header=header) as vcf_out:
            print("Applying annotations...")
            for record in vcf_in:
                record_out = annotate(record=record, vcf_out=vcf_out, info_dict=info_dict, format_dict=format_dict)
                vcf_out.write(record_out)


if __name__ == "__main__":
    main()
