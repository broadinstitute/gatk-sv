#!/bin/env python

import argparse
import pysam
import sys
import gzip
from math import floor
from typing import Any, List, Text, Dict, Optional

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


def update_header(header: pysam.VariantHeader) -> None:
    """
    Ingests the given header, removes specified fields, and adds necessary fields.

    Parameters
    ----------
    header: pysam.VariantHeader
        input header
    """
    header.add_line('##FORMAT=<ID=ECN,Number=1,Type=Integer,Description="Expected copy number for ref genotype">')
    # Add these just in case (no effect if they exist)
    header.add_line('##INFO=<ID=END2,Number=1,Type=Integer,Description="Second position">')
    header.add_line('##INFO=<ID=CHR2,Number=1,Type=String,Description="Second contig">')


def rescale_gq(record):
    for sample in record.samples:
        for gq_field in GQ_FIELDS:
            if gq_field in record.samples[sample] and record.samples[sample][gq_field] is not None:
                record.samples[sample][gq_field] = floor(record.samples[sample][gq_field] / 10)


def convert(record: pysam.VariantRecord,
            bnd_end_dict: Optional[Dict[Text, int]],
            ploidy_dict: Dict[Text, Dict[Text, int]],
            scale_down_gq: bool) -> pysam.VariantRecord:
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
    # Fix this weird edge case (may be from CPX review workflow)
    if svtype == 'INV' and '<CPX>' in record.alleles[1]:
        svtype = 'CPX'
        record.info['SVTYPE'] = svtype
    is_ddup = svtype == 'CPX' and 'dDUP' in record.info.get('CPX_TYPE', '')
    if svtype == 'BND' or svtype == 'INS' or svtype == 'CTX' or is_ddup:
        record.stop = record.start + 1
        if is_ddup:
            # e.g. SOURCE=DUP_chrX:49151588-49151850
            source = record.info.get('SOURCE', None)
            if source is not None:
                tokens = source.split(':')
                chr2 = tokens[0].split('_')[-1]
                end2 = int(tokens[-1].split('-')[0])
                record.info['CHR2'] = chr2
                record.info['END2'] = end2
            else:
                # Sometimes SOURCE is not set (may be from CPX review workflow)
                record.info['CHR2'] = record.chrom
                record.info['END2'] = record.stop
    # Delete empty INFO fields (GATK does not like "." for non-String types)
    keys = record.info.keys()
    for k in keys:
        val = record.info[k]
        if is_null(val) or (isinstance(val, tuple) and len(val) == 1 and is_null(val[0])):
            del record.info[k]
    # copy FORMAT fields
    for sample, genotype in record.samples.items():
        genotype['ECN'] = ploidy_dict[sample][contig]
    if scale_down_gq:
        rescale_gq(record)
    return record


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

    for record in vcf_in:
        out = convert(record=record, bnd_end_dict=bnd_end_dict,
                      ploidy_dict=ploidy_dict, scale_down_gq=arguments.scale_down_gq)
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
            header=vcf_in.header
        )
        with pysam.VariantFile(arguments.out, mode='w', header=vcf_in.header) as vcf_out:
            _process(vcf_in, vcf_out, arguments)


if __name__ == "__main__":
    main()
