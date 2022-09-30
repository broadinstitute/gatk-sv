#!/bin/env python

import argparse
import pysam
import sys
from typing import Any, List, Text, Set, Dict, Optional

_gt_sum_map = dict()


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


def create_header(header_in: pysam.VariantHeader) -> pysam.VariantHeader:
    """
    Ingests the given header, removes specified fields, and adds necessary fields.

    Parameters
    ----------
    header_in: pysam.VariantHeader
        input header

    Returns
    -------
    header: pysam.VariantHeader
        gatk-style header
    """
    header = pysam.VariantHeader()
    for sample in header_in.samples:
        header.add_sample(sample)
    # new fields
    header.add_line('##fileformat=VCFv4.2')
    header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header.add_line('##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number">')
    header.add_line('##FORMAT=<ID=ECN,Number=1,Type=Integer,Description="Expected copy number for ref genotype">')
    header.add_line('##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description="Source algorithms">')
    header.add_line('##INFO=<ID=CHR2,Number=1,Type=String,Description="Second contig">')
    header.add_line('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">')
    header.add_line('##INFO=<ID=END2,Number=1,Type=Integer,Description="Second position">')
    header.add_line('##INFO=<ID=STRANDS,Number=1,Type=String,Description="First and second strands">')
    header.add_line('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of affected segment on the reference">')
    header.add_line('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">')
    header.add_line('##contig=<ID=chr1,length=248956422,assembly=38>')
    header.add_line('##contig=<ID=chr2,length=242193529,assembly=38>')
    header.add_line('##contig=<ID=chr3,length=198295559,assembly=38>')
    header.add_line('##contig=<ID=chr4,length=190214555,assembly=38>')
    header.add_line('##contig=<ID=chr5,length=181538259,assembly=38>')
    header.add_line('##contig=<ID=chr6,length=170805979,assembly=38>')
    header.add_line('##contig=<ID=chr7,length=159345973,assembly=38>')
    header.add_line('##contig=<ID=chr8,length=145138636,assembly=38>')
    header.add_line('##contig=<ID=chr9,length=138394717,assembly=38>')
    header.add_line('##contig=<ID=chr10,length=133797422,assembly=38>')
    header.add_line('##contig=<ID=chr11,length=135086622,assembly=38>')
    header.add_line('##contig=<ID=chr12,length=133275309,assembly=38>')
    header.add_line('##contig=<ID=chr13,length=114364328,assembly=38>')
    header.add_line('##contig=<ID=chr14,length=107043718,assembly=38>')
    header.add_line('##contig=<ID=chr15,length=101991189,assembly=38>')
    header.add_line('##contig=<ID=chr16,length=90338345,assembly=38>')
    header.add_line('##contig=<ID=chr17,length=83257441,assembly=38>')
    header.add_line('##contig=<ID=chr18,length=80373285,assembly=38>')
    header.add_line('##contig=<ID=chr19,length=58617616,assembly=38>')
    header.add_line('##contig=<ID=chr20,length=64444167,assembly=38>')
    header.add_line('##contig=<ID=chr21,length=46709983,assembly=38>')
    header.add_line('##contig=<ID=chr22,length=50818468,assembly=38>')
    header.add_line('##contig=<ID=chrX,length=156040895,assembly=38>')
    header.add_line('##contig=<ID=chrY,length=57227415,assembly=38>')
    return header


def convert(record: pysam.VariantRecord,
            vcf_out: pysam.VariantFile,
            algorithm: Text,
            min_size: int,
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
    ploidy_dict: Dict[Text, Dict[Text, int]]
        map from sample to contig to ploidy

    Returns
    -------
    pysam.VariantRecord
        gatk-style record
    """
    svtype = record.info['SVTYPE']
    if not supported_type(record):
        return list()
    contig = record.contig
    if contig not in vcf_out.header.contigs:
        return list()
    # Force DUPs to insertions for consistency
    isDup = svtype == 'DUP'
    if isDup:
        svtype = 'INS'
    if 'SVLEN' in record.info and record.info['SVLEN'] is not None:
        svlen = record.info['SVLEN']
        if isinstance(svlen, tuple):
            svlen = svlen[0]
        svlen = abs(int(svlen))
        if svtype == 'INS':
            end = record.start + 1
        else:
            end = record.start + svlen
    else:
        svlen = record.stop - record.pos
        end = record.stop
    if svlen < min_size:
        return list()
    if end <= record.start:
        end = record.start + 1
    # Force symbolic allele
    alleles = ["N", f"<{svtype}>"]
    new_record = vcf_out.new_record(id=record.id, contig=contig, start=record.start, stop=end, alleles=alleles)
    new_record.info['ALGORITHMS'] = [algorithm]
    new_record.info['SVTYPE'] = svtype
    # fix SVLEN, STRANDS, CHR2, and END2 where needed
    if svtype == 'INS':
        new_record.info['SVLEN'] = svlen
    elif svtype == 'INV':
        new_record.info['STRANDS'] = '++'
    # copy FORMAT fields
    for sample in record.samples:
        genotype = record.samples[sample]
        new_genotype = new_record.samples[sample]
        ecn = ploidy_dict[sample][contig]
        new_genotype['ECN'] = ecn
        gt_sum = _cache_gt_sum(genotype['GT'])
        if svtype == 'DEL':
            new_genotype['CN'] = max(ecn - gt_sum, 0)
        if new_genotype['ECN'] == 0:
            new_genotype['GT'] = ()
        elif ecn == 1:
            if gt_sum == 0:
                new_genotype['GT'] = (0,)
            else:
                new_genotype['GT'] = (1,)
        else:
            if gt_sum == 0:
                new_genotype['GT'] = (0, 0)
            elif gt_sum == 1:
                new_genotype['GT'] = (0, 1)
            else:
                new_genotype['GT'] = (1, 1)
    out = [new_record]
    return out


def _cache_gt_sum(gt):
    s = _gt_sum_map.get(gt, None)
    if s is None:
        s = sum([1 for a in gt if a is not None and a > 0])
        _gt_sum_map[gt] = s
    return s


def supported_type(record: pysam.VariantRecord) -> bool:
    return record.info['SVTYPE'] in ['DEL', 'DUP', 'INS', 'INV']


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
    """
    ploidy_dict = _parse_ploidy_table(arguments.ploidy_table)

    for record in vcf_in:
        out = convert(
            record=record,
            vcf_out=vcf_out,
            algorithm=arguments.algorithm,
            min_size=arguments.min_size,
            ploidy_dict=ploidy_dict
        )
        for record in out:
            vcf_out.write(record)


def _parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Convert a PacBio-derived SV VCF to SVTK-style",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--vcf", type=str, required=True,
                        help="Input VCF. Supported callers include: pbsv, pav, sniffles")
    parser.add_argument("--out", type=str, required=True,
                        help="Output VCF")
    parser.add_argument("--algorithm", type=str, required=True,
                        help="Algorithm name")
    parser.add_argument("--min-size", type=int, default=25,
                        help="Min SV size")
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
