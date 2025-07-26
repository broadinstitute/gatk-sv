#!/bin/env python

import argparse
import pysam
import sys
from typing import Optional, List, Text, Set

_gt_map = dict()
_cnv_gt_map = dict()


def _cache_gt(gt):
    if gt is None:
        return 0, 0
    s = _gt_map.get(gt, None)
    if s is None:
        x = sum([1 for a in gt if a is not None and a > 0])
        if x == 0:
            s = (0, 0)
        elif x == 1:
            s = (0, 1)
        else:
            s = (1, 1)
        _gt_map[gt] = s
    return s


def _cache_cnv_gt(cn):
    if cn is None:
        return 0, 0
    s = _cnv_gt_map.get(cn, None)
    if s is None:
        # Split copies evenly between alleles, giving one more copy to the second allele if odd
        x = max(0, cn - 2)
        alt1 = min(x // 2, 4)
        alt2 = min(x - alt1, 4)
        _cnv_gt_map[cn] = (alt1, alt2)
    return s


def create_header(header_in: pysam.VariantHeader,
                  source: Text,
                  remove_infos: Set[Text],
                  remove_formats: Set[Text],
                  contigs: Set[Text]) -> pysam.VariantHeader:
    """
    Ingests the given header and removes specified fields, and subsets sequence dictionary to given contigs.

    Parameters
    ----------
    header_in: pysam.VariantHeader
        input header
    source: Text
        source application, for which a '##source=' line is added
    remove_infos: Set[Text]
        set of info fields to remove
    remove_formats: Set[Text]
        set of format fields to remove
    contigs: Set[Text]
        set of allowed contigs

    Returns
    -------
    header: pysam.VariantHeader
        svtk-style header
    """
    header = pysam.VariantHeader()
    for sample in header_in.samples:
        header.add_sample(sample)
    for line in header_in.records:
        # remove INFO/FORMAT fields
        if len(line.attrs) > 0 and 'ID' in line.keys() and \
                (line['ID'] in remove_infos or line['ID'] in remove_formats):
            continue
        # GATK may produce SVLEN with Number=.
        if 'ID' in line.keys() and line['ID'] == 'SVLEN':
            header.add_line("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of "
                            "affected segment on the reference\">")
            continue
        line_str = str(line)
        # subset contigs
        if line_str.startswith('##contig=') and line['ID'] not in contigs:
            continue
        header.add_line(line_str)
    # source line
    if source is not None:
        header.add_line("##source=" + source)
    return header


def convert(record: pysam.VariantRecord,
            vcf_out: pysam.VariantFile,
            remove_infos: Set[Text],
            remove_formats: Set[Text],
            set_pass: bool) -> pysam.VariantRecord:
    """
    Converts a record from gatk to svtk style. This includes updating all GT fields to diploid and reverting END
    tag values.

    Parameters
    ----------
    record: pysam.VariantRecord
        gatk-style record
    vcf_out: pysam.VariantFile
        new vcf, to which the converted record will be written
    remove_infos: Set[Text]
        info fields to remove
    remove_formats: Set[Text]
        format fields to remove
    set_pass: bool
        set empty FILTER statuses to PASS

    Returns
    -------
    header: pysam.VariantRecord
        svtk-style record
    """
    svtype = record.info['SVTYPE']
    alleles = record.alleles
    # Reset alleles
    if len(alleles) != 2:
        raise ValueError(f"Encountered non-biallelic site {record.id} with alleles {alleles}")
    if svtype == 'BND':
        # Ensure we aren't using breakend notation here, since it isn't supported in some modules
        alleles = ('N', '<BND>')
    elif svtype == 'CNV':
        # Prior to CleanVcf, all CNVs have <CNx> alleles
        alleles = ('N', '<CN0>', '<CN1>', '<CN2>', '<CN3>')
    else:
        alleles = ('N', alleles[1])
    contig = record.contig
    # Change filter to PASS if requested
    if set_pass and len(record.filter) == 0:
        filter = ('PASS',)
    else:
        filter = record.filter
    new_record = vcf_out.new_record(contig=contig, start=record.start, stop=record.stop, alleles=alleles, filter=filter)
    new_record.id = record.id
    # copy INFO fields
    for key in record.info:
        if key not in remove_infos:
            new_record.info[key] = record.info[key]
    if svtype == 'CNV':
        # Prior to CleanVcf, all mCNVs are DUP type
        new_record.info['SVTYPE'] = 'DUP'
    # svtk generally expects all records to have CHR2 assigned
    chr2 = record.info.get('CHR2', None)
    if chr2 is None:
        new_record.info['CHR2'] = contig
    # fix END, SVLEN, STRANDS
    if svtype == 'INS':
        new_record.info['SVLEN'] = record.info.get('SVLEN', -1)
        new_record.info['STRANDS'] = '+-'
        # END information is lost when setting POS=END, so we need to set it to SR2POS if it's available
        # Note that the END position is only important following SR breakpoint refinement in FilterBatch
        if 'SR2POS' in record.info and record.info['SR2POS'] is not None:
            new_record.stop = record.info['SR2POS']
    elif svtype == 'BND' or svtype == 'CTX':
        new_record.info['SVLEN'] = -1
    elif svtype == 'CPX':
        new_record.info['SVLEN'] = record.info.get('SVLEN', -1)
    elif svtype == 'DEL':
        new_record.info['STRANDS'] = '+-'
    elif svtype == 'DUP' or svtype == 'CNV':
        new_record.info['STRANDS'] = '-+'
    elif svtype == 'INV':
        new_record.info['STRANDS'] = record.info.get('STRANDS', None)

    for sample in record.samples:
        new_genotype = new_record.samples[sample]
        genotype = record.samples[sample]
        # copy FORMAT fields
        for key in genotype.keys():
            if key not in remove_formats:
                new_genotype[key] = genotype[key]
        # fix GT, always assuming diploid
        if svtype == 'CNV':
            new_genotype['GT'] = _cache_cnv_gt(genotype.get('CN', genotype.get('RD_CN')))
        else:
            new_genotype['GT'] = _cache_gt(genotype.get('GT', None))
    return new_record


def __read_contigs(path: Text) -> List[Text]:
    with open(path, 'r') as f:
        return [line.strip().split('\t')[0] for line in f]


def __parse_arg_list(arg: Text) -> List[Text]:
    if arg is None:
        return set()
    else:
        return arg.split(',')


def __parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Convert a GATK-style SV VCF for consumption by svtk. Not to be used after CleanVcf.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--vcf", type=str, required=True,
                        help="GATK VCF")
    parser.add_argument("--contigs", type=str, required=True,
                        help="List of contigs")
    parser.add_argument("--out", type=str, required=True,
                        help="Output VCF")
    parser.add_argument("--source", type=str,
                        help="Source application (adds source header line)")
    parser.add_argument("--remove-formats", type=str,
                        help="Comma-delimited list of FORMAT fields to remove")
    parser.add_argument("--remove-infos", type=str,
                        help="Comma-delimited list of INFO fields to remove")
    parser.add_argument("--set-pass", default=False, action='store_true',
                        help="Set empty FILTER fields (\".\") to PASS")
    if len(argv) <= 1:
        parser.parse_args(["--help"])
        sys.exit(0)
    parsed_arguments = parser.parse_args(argv[1:])
    return parsed_arguments


def main(argv: Optional[List[Text]] = None):
    if argv is None:
        argv = sys.argv
    arguments = __parse_arguments(argv)
    contigs = set(__read_contigs(arguments.contigs))
    remove_formats = set(__parse_arg_list(arguments.remove_formats))
    remove_infos = set(__parse_arg_list(arguments.remove_infos))

    # fields removed by default
    remove_formats.add('ECN')

    # convert vcf header and records
    with pysam.VariantFile(arguments.vcf) as vcf_in:
        header = create_header(
            header_in=vcf_in.header,
            source=arguments.source,
            remove_infos=remove_infos,
            remove_formats=remove_formats,
            contigs=contigs
        )
        with pysam.VariantFile(arguments.out, mode='w', header=header) as vcf_out:
            for record in vcf_in:
                vcf_out.write(convert(
                    record=record,
                    vcf_out=vcf_out,
                    remove_infos=remove_infos,
                    remove_formats=remove_formats,
                    set_pass=arguments.set_pass
                ))


if __name__ == "__main__":
    main()
