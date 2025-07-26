#!/bin/env python

import argparse
import pysam
import sys
import gzip
from math import floor
from typing import Any, List, Text, Dict, Optional, Set

GQ_FIELDS = ["GQ", "PE_GQ", "SR_GQ", "RD_GQ"]


def _parse_bnd_ends(vcf_path: Text) -> Dict[Text, int]:
    """
    Since pysam automatically changes invalid END fields (i.e. when less than the start position), they must
    be parsed manually.

    Parameters
    ----------
    vcf_path: Text
        input vcf path

    Returns
    -------
    header: Dict[Text, int]
        map from variant ID to END position
    """
    bnd_end_dict = dict()
    with gzip.open(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            columns = line.split('\t', 8)
            vid = columns[2]
            info = columns[7]
            if 'SVTYPE=BND' not in info and 'SVTYPE=CTX' not in info:
                continue
            info_tokens = info.split(';')
            end_field_list = [x for x in info_tokens if x.startswith("END=")]
            if len(end_field_list) > 0:
                end = int(end_field_list[0].replace("END=", ""))
            else:
                # Special case where END and POS happen to be equal
                end = int(columns[1])
            bnd_end_dict[vid] = end
    return bnd_end_dict


def _parse_ploidy_table(path: Text) -> Dict[Text, Dict[Text, int]]:
    """
    Parses tsv of sample ploidy values.

    Parameters
    ----------
    path: Text
        table path

    Returns
    -------
    header: Dict[Text, Dict[Text, int]]
        map of sample to contig to ploidy, i.e. Dict[sample][contig] = ploidy
    """
    ploidy_dict = dict()
    with open(path, 'r') as f:
        header = f.readline().strip().split('\t')
        for line in f:
            tokens = line.strip().split('\t')
            sample = tokens[0]
            ploidy_dict[sample] = {header[i]: int(tokens[i]) for i in range(1, len(header))}
    return ploidy_dict


def update_header(header: pysam.VariantHeader,
                  add_sr_pos: bool) -> None:
    """
    Ingests the given header, removes specified fields, and adds necessary fields.

    Parameters
    ----------
    header: pysam.VariantHeader
        input header
    add_sr_pos: bool
        add SR1POS and SR2POS INFO fields
    """
    header.add_line('##FORMAT=<ID=ECN,Number=1,Type=Integer,Description="Expected copy number for ref genotype">')
    # Add these just in case (no effect if they exist)
    header.add_line('##INFO=<ID=END2,Number=1,Type=Integer,Description="Second position">')
    header.add_line('##INFO=<ID=CHR2,Number=1,Type=String,Description="Second contig">')
    header.add_line('##INFO=<ID=BOTHSIDES_SUPPORT,Number=0,Type=Flag,Description="Variant has read-level support for both sides of breakpoint">')
    header.add_line('##INFO=<ID=HIGH_SR_BACKGROUND,Number=0,Type=Flag,Description="High number of SR splits in background samples indicating messy region">')
    if add_sr_pos:
        header.add_line('##INFO=<ID=SR1POS,Number=1,Type=Integer,Description="Split read position at start">')
        header.add_line('##INFO=<ID=SR2POS,Number=1,Type=Integer,Description="Split read position at end">')


def rescale_gq(record):
    for sample in record.samples:
        for gq_field in GQ_FIELDS:
            if gq_field in record.samples[sample] and record.samples[sample][gq_field] is not None:
                record.samples[sample][gq_field] = floor(record.samples[sample][gq_field] / 10)


def convert(record: pysam.VariantRecord,
            bnd_end_dict: Optional[Dict[Text, int]],
            ploidy_dict: Dict[Text, Dict[Text, int]],
            scale_down_gq: bool,
            bothside_pass_vid_set: Set[Text],
            background_fail_vid_set: Set[Text],
            add_sr_pos: bool) -> pysam.VariantRecord:
    """
    Converts a record from svtk to gatk style. This includes updating END/END2 and adding
    necessary fields such as ECN.

    Parameters
    ----------
    record: pysam.VariantRecord
        svtk-style record
    bnd_end_dict: Optional[Dict[Text, int]]
        map from BND variant ID to END coordinate
    ploidy_dict: Dict[Text, Dict[Text, int]]
        map from sample to contig to ploidy
    scale_down_gq: bool
        scale GQs to 0-99 range
    bothside_pass_vid_set: Set[Text]
        set of variant IDs with bothside SR pass flag
    background_fail_vid_set: Set[Text]
        set of variant Ids with background SR fail flag
    add_sr_pos: bool
        add SR1POS and SR2POS INFO fields

    Returns
    -------
    header: pysam.VariantRecord
        gatk-style record
    """

    def is_null(val):
        return val is None or val == "."

    svtype = record.info['SVTYPE']
    contig = record.contig
    # Version of htsjdk currently in gatk only supports these base alleles
    if record.ref in ['A', 'C', 'G', 'T', 'a', 'c', 'g', 't', 'N', 'n']:
        record.ref = record.ref
    else:
        record.ref = 'N'
    if svtype == 'BND' or svtype == 'CTX':
        record.info['END2'] = bnd_end_dict[record.id] if bnd_end_dict is not None \
            else record.info.get('END2', record.stop)
    # Add SR1POS/SR2POS, CPX not supported
    if add_sr_pos and svtype != 'CPX':
        record.info['SR1POS'] = record.pos
        if svtype == 'BND' or svtype == 'CTX':
            record.info['SR2POS'] = record.info['END2']
        else:
            record.info['SR2POS'] = record.stop
    # Fix this weird edge case (may be from CPX review workflow)
    if svtype == 'INV' and '<CPX>' in record.alleles[1]:
        svtype = 'CPX'
        record.info['SVTYPE'] = svtype
    is_ddup = svtype == 'CPX' and 'dDUP' in record.info.get('CPX_TYPE', '')
    if svtype == 'BND' or svtype == 'INS' or svtype == 'CTX' or is_ddup:
        record.stop = record.pos + 1
    if svtype == 'CPX' or svtype == 'INS':
        if 'CHR2' in record.info:
            record.info.pop('CHR2')
        if 'END2' in record.info:
            record.info.pop('END2')
    # Delete empty INFO fields (GATK does not like "." for non-String types)
    keys = record.info.keys()
    for k in keys:
        val = record.info[k]
        if is_null(val) or (isinstance(val, tuple) and len(val) == 1 and is_null(val[0])):
            del record.info[k]
    # Add SR flag to INFO field
    if record.id in bothside_pass_vid_set:
        record.info['BOTHSIDES_SUPPORT'] = True
    if record.id in background_fail_vid_set:
        record.info['HIGH_SR_BACKGROUND'] = True
    # copy FORMAT fields
    for sample, genotype in record.samples.items():
        genotype['ECN'] = ploidy_dict[sample][contig]
    if scale_down_gq:
        rescale_gq(record)
    return record


def parse_last_column(path):
    if path is None:
        return set()
    with open(path) as f:
        return set([line.strip().split('\t')[-1] for line in f])


def _process(vcf_in: pysam.VariantFile,
             vcf_out: pysam.VariantFile,
             arguments: Dict[Text, Any]) -> None:
    """"
    Master function for processing the given input vcf and writing output

    Parameters
    ----------
    vcf_in: pysam.VariantFile
        input vcf
    vcf_out: pysam.VariantFile
        output vcf
    arguments: Dict[Text, Any]
        commandline arguments

    Returns
    -------
    header: pysam.VariantRecord
        record with ECN fields added"""
    if arguments.fix_end:
        bnd_end_dict = _parse_bnd_ends(arguments.vcf)
    else:
        bnd_end_dict = None
    ploidy_dict = _parse_ploidy_table(arguments.ploidy_table)

    bothside_pass_vid_set = parse_last_column(arguments.bothside_pass_list)
    background_fail_vid_set = parse_last_column(arguments.background_fail_list)

    for record in vcf_in:
        out = convert(record=record, bnd_end_dict=bnd_end_dict,
                      ploidy_dict=ploidy_dict, scale_down_gq=arguments.scale_down_gq,
                      bothside_pass_vid_set=bothside_pass_vid_set,
                      background_fail_vid_set=background_fail_vid_set,
                      add_sr_pos=arguments.add_sr_pos)
        vcf_out.write(out)


def _parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Convert a SVTK-style SV VCF to GATK-style",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--vcf", type=str, required=True,
                        help="GATK VCF")
    parser.add_argument("--out", type=str, required=True,
                        help="Output VCF")
    parser.add_argument("--ploidy-table", type=str, required=True,
                        help="Tab-delimited table of sample ploidies. The table should have a header row where the "
                             "first column is SAMPLE, and the remaining columns are contig names. For each row "
                             "thereafter, the first column is the sample name, and remaining columns are the contig "
                             "ploidy values for that sample.")
    parser.add_argument("--fix-end", action='store_true',
                        help="Fix END tags and assign END2 to END")
    parser.add_argument("--scale-down-gq", action='store_true',
                        help="Scales all GQs down from [0-999] to [0-99]")
    parser.add_argument("--bothside-pass-list", type=str,
                        help="Path to bothside SR pass flag variant list")
    parser.add_argument("--background-fail-list", type=str,
                        help="Path to background SR fail flag variant list")
    parser.add_argument("--add-sr-pos", action='store_true',
                        help="Add SR1POS and SR2POS INFO fields")
    if len(argv) <= 1:
        parser.parse_args(["--help"])
        sys.exit(0)
    parsed_arguments = parser.parse_args(argv[1:])
    return parsed_arguments


def main(argv: Optional[List[Text]] = None):
    if argv is None:
        argv = sys.argv
    arguments = _parse_arguments(argv)

    # convert vcf header and records
    with pysam.VariantFile(arguments.vcf) as vcf_in:
        update_header(
            header=vcf_in.header,
            add_sr_pos=arguments.add_sr_pos
        )
        with pysam.VariantFile(arguments.out, mode='w', header=vcf_in.header) as vcf_out:
            _process(vcf_in, vcf_out, arguments)


if __name__ == "__main__":
    main()
