#!/bin/env python

import argparse
import pysam
import sys
from typing import Any, List, Text, Set, Dict, Optional


def create_header(header: pysam.VariantHeader,) -> pysam.VariantHeader:
    """
    Ingests the given header and removes specified fields, and subsets sequence dictionary to given contigs.

    Parameters
    ----------
    header_in: pysam.VariantHeader
        input header

    Returns
    -------
    header: pysam.VariantHeader
        svtk-style header
    """
    header.add_line('##FORMAT=<ID=ECN,Number=1,Type=Integer,Description="Expected copy number for ref genotype">')
    return header


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


def convert(record: pysam.VariantRecord,
            ploidy_dict: Dict[Text, Dict[Text, int]]) -> pysam.VariantRecord:
    """
    Converts a record from svtk to gatk style. This includes updating all GT fields with proper ploidy, and adding
    necessary fields such as ECN and CN.

    Parameters
    ----------
    record: pysam.VariantRecord
        input record
    vcf_out: pysam.VariantFile
        new vcf, to which the converted record will be written
    ploidy_dict: Dict[Text, Dict[Text, int]]
        map from sample to contig to ploidy

    Returns
    -------
    record: pysam.VariantRecord
        annotated record
    """
    # fix SVLEN, STRANDS, CHR2, and END2 where needed
    #if svtype == 'INS' or svtype == 'CPX':
    #    record.info['SVLEN'] = record.info['SVLEN']
    #if svtype == 'BND':
    #    record.info['STRANDS'] = record.info['STRANDS']
    #if svtype == 'BND' or svtype == 'CTX' or svtype == 'CPX':
    #    record.info['END2'] = bnd_end_dict[record.id] if bnd_end_dict is not None \
    #        else record.info.get('END2', record.stop)
    #if svtype == 'BND' or svtype == 'CTX':
    #    record.stop = record.start + 1
    contig = record.contig
    # copy FORMAT fields
    for sample, genotype in record.samples.items():
        genotype['ECN'] = ploidy_dict[sample][contig]
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
        record with CN and ECN fields added"""
    ploidy_dict = _parse_ploidy_table(arguments.ploidy_table)
    for record in vcf_in:
        out = convert(record=record, ploidy_dict=ploidy_dict)
        vcf_out.write(out)


def _parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Convert a GATK-style SV VCF to SVTK-style",
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
        header = create_header(vcf_in.header)
        with pysam.VariantFile(arguments.out, mode='w', header=header) as vcf_out:
            _process(vcf_in, vcf_out, arguments)


if __name__ == "__main__":
    main()
