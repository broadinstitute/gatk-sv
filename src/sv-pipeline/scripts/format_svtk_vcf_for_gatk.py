#!/bin/env python

import argparse
import pysam
import sys
import gzip
from typing import List, Text, Set, Dict, Optional


def __parse_bnd_ends(vcf_path: Text) -> Dict[Text, int]:
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
            if 'SVTYPE=BND' not in info:
                continue
            info_tokens = info.split(';')
            end = [x for x in info_tokens if x.startswith("END=")][0]
            bnd_end_dict[vid] = int(end.replace("END=", ""))
    return bnd_end_dict


def __parse_ploidy_table(path: Text) -> Dict[Text, Dict[Text, int]]:
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


def create_header(header_in: pysam.VariantHeader,
                  remove_infos: Set[Text],
                  remove_formats: Set[Text]) -> pysam.VariantHeader:
    """
    Ingests the given header, removes specified fields, and adds necessary fields.

    Parameters
    ----------
    header_in: pysam.VariantHeader
        input header
    remove_infos: Set[Text]
        set of info fields to remove
    remove_formats: Set[Text]
        set of format fields to remove

    Returns
    -------
    header: pysam.VariantHeader
        gatk-style header
    """
    header = pysam.VariantHeader()
    for sample in header_in.samples:
        header.add_sample(sample)
    for line in header_in.records:
        # remove fields
        if len(line.attrs) > 0 and 'ID' in line.keys() and (line['ID'] in remove_infos or line['ID'] in remove_formats):
            continue
        line_str = str(line)
        # remove source line
        if line_str.startswith('##source='):
            continue
        header.add_line(line_str)
    # new fields
    header.add_line('##INFO=<ID=END2,Number=1,Type=Integer,Description="Second breakend position">')
    header.add_line('##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number">')
    header.add_line('##FORMAT=<ID=ECN,Number=1,Type=Integer,Description="Expected copy number for ref genotype">')
    return header


def convert(record: pysam.VariantRecord,
            vcf_out: pysam.VariantFile,
            remove_infos: Set[Text],
            remove_formats: Set[Text],
            bnd_end_dict: Dict[Text, int],
            ploidy_dict: Dict[Text, Dict[Text, int]]) -> pysam.VariantRecord:
    """
    Converts a record from svtk to gatk style. This includes updating all GT fields with proper ploidy, and adding
    necessary fields such as ECN and CN.

    Parameters
    ----------
    record: pysam.VariantRecord
        svtk-style record
    vcf_out: pysam.VariantFile
        new vcf, to which the converted record will be written
    remove_infos: Set[Text]
        info fields to remove
    remove_formats: Set[Text]
        format fields to remove
    bnd_end_dict: Dict[Text, int]
        map from BND variant ID to END coordinate
    ploidy_dict: Dict[Text, Dict[Text, int]]
        map from sample to contig to ploidy

    Returns
    -------
    header: pysam.VariantRecord
        gatk-style record
    """
    # info fields we drop by default (unless needed for certain SV types)
    default_remove_infos = set(["SVLEN", "STRANDS", "CHR2", "END2"])
    svtype = record.info['SVTYPE']
    # Force symbolic BND alleles
    if svtype == 'BND':
        alleles = (record.alleles[0], '<BND>')
    else:
        alleles = record.alleles
    contig = record.contig
    new_record = vcf_out.new_record(contig=contig, start=record.start, stop=record.stop, alleles=alleles)
    new_record.id = record.id
    # copy INFO fields
    for key in record.info:
        if key not in default_remove_infos and key not in remove_infos:
            new_record.info[key] = record.info[key]
    # fix SVLEN, STRANDS, CHR2, and END2 where needed
    if svtype == 'INS':
        new_record.info['SVLEN'] = record.info['SVLEN']
    elif svtype == 'BND':
        new_record.info['STRANDS'] = record.info['STRANDS']
        new_record.info['CHR2'] = record.info['CHR2']
        new_record.info['END2'] = bnd_end_dict[record.id]
        new_record.stop = record.start + 1
    elif svtype == 'INV':
        new_record.info['STRANDS'] = record.info['STRANDS']
    # copy FORMAT fields
    for sample in record.samples:
        genotype = record.samples[sample]
        new_genotype = new_record.samples[sample]
        for key in genotype.keys():
            if key not in remove_formats:
                new_genotype[key] = genotype[key]
        new_genotype['ECN'] = ploidy_dict[sample][contig]
        if new_genotype['ECN'] == 0:
            new_genotype['GT'] = ()
            new_genotype['CN'] = 0
        elif new_genotype['ECN'] == 1:
            if svtype == 'DUP':
                new_genotype['CN'] = 1 + sum(genotype['GT'])
                new_genotype['GT'] = (None,)
            elif sum(genotype['GT']) == 0:
                new_genotype['GT'] = (0,)
                if svtype == 'DEL':
                    new_genotype['CN'] = 1
            else:
                new_genotype['GT'] = (1,)
                if svtype == 'DEL':
                    new_genotype['CN'] = 0
        else:
            if svtype == 'DUP':
                new_genotype['CN'] = 2 + sum(genotype['GT'])
                new_genotype['GT'] = (None, None)
            elif sum(genotype['GT']) == 0:
                new_genotype['GT'] = (0, 0)
                if svtype == 'DEL':
                    new_genotype['CN'] = 2
            elif sum(genotype['GT']) == 1:
                new_genotype['GT'] = (0, 1)
                if svtype == 'DEL':
                    new_genotype['CN'] = 1
            else:
                new_genotype['GT'] = (1, 1)
                if svtype == 'DEL':
                    new_genotype['CN'] = 0
        if svtype == 'CNV':
            new_genotype['GT'] = (None,) * new_genotype['ECN']
    return new_record


def __parse_arg_list(arg: Text) -> List[Text]:
    if arg is None:
        return set()
    else:
        return arg.split(',')


def __parse_arguments(argv: List[Text]) -> argparse.Namespace:
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
    parser.add_argument("--remove-formats", type=str,
                        help="Comma-delimited list of FORMAT fields to remove")
    parser.add_argument("--remove-infos", type=str,
                        help="Comma-delimited list of INFO fields to remove")
    if len(argv) <= 1:
        parser.parse_args(["--help"])
        sys.exit(0)
    parsed_arguments = parser.parse_args(argv[1:])
    return parsed_arguments


def main(argv: Optional[List[Text]] = None):
    if argv is None:
        argv = sys.argv
    arguments = __parse_arguments(argv)
    remove_formats = set(__parse_arg_list(arguments.remove_formats))
    remove_infos = set(__parse_arg_list(arguments.remove_infos))
    bnd_end_dict = __parse_bnd_ends(arguments.vcf)
    ploidy_dict = __parse_ploidy_table(arguments.ploidy_table)

    # convert vcf header and records
    vcf_in = pysam.VariantFile(arguments.vcf)
    header = create_header(
        header_in=vcf_in.header,
        remove_infos=remove_infos,
        remove_formats=remove_formats
    )
    vcf_out = pysam.VariantFile(arguments.out, mode='w', header=header)
    for record in vcf_in:
        vcf_out.write(convert(
            record=record,
            vcf_out=vcf_out,
            remove_infos=remove_infos,
            remove_formats=remove_formats,
            bnd_end_dict=bnd_end_dict,
            ploidy_dict=ploidy_dict
        ))
    vcf_in.close()
    vcf_out.close()


if __name__ == "__main__":
    main()
