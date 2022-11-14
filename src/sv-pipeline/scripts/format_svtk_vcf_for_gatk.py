#!/bin/env python

import argparse
import pysam
import sys
import gzip
from typing import Any, List, Text, Set, Dict, Optional

_gt_sum_map = dict()


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
            if 'SVTYPE=BND' not in info and 'SVTYPE=CTX' not in info and 'SVTYPE=CPX' not in info:
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


def create_header(header_in: pysam.VariantHeader,
                  replace_ev_format: bool,
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
    header.add_line('##INFO=<ID=OSVTYPE,Number=1,Type=String,Description="Original SVTYPE">')
    header.add_line('##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number">')
    header.add_line('##FORMAT=<ID=ECN,Number=1,Type=Integer,Description="Expected copy number for ref genotype">')
    if replace_ev_format:
        header.add_line('##FORMAT=<ID=EV,Number=.,Type=String,Description="Classes of evidence supporting final '
                        'genotype">')
    return header


def convert(record: pysam.VariantRecord,
            vcf_out: pysam.VariantFile,
            remove_infos: Set[Text],
            remove_formats: Set[Text],
            bnd_end_dict: Optional[Dict[Text, int]],
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
    bnd_end_dict: Optional[Dict[Text, int]]
        map from BND variant ID to END coordinate
    ploidy_dict: Dict[Text, Dict[Text, int]]
        map from sample to contig to ploidy

    Returns
    -------
    header: pysam.VariantRecord
        gatk-style record
    """
    svtype = record.info['SVTYPE']
    # Force symbolic BND alleles
    if svtype == 'BND':
        alleles = (record.alleles[0], '<BND>')
    else:
        alleles = record.alleles
    contig = record.contig
    new_record = vcf_out.new_record(id=record.id, contig=contig, start=record.start, stop=record.stop, alleles=alleles,
                                    info={key: value for key, value in record.info.items() if key not in remove_infos})
    # fix SVLEN, STRANDS, CHR2, and END2 where needed
    if svtype == 'INS':
        new_record.info['SVLEN'] = record.info['SVLEN']
    elif svtype == 'BND' or svtype == 'CTX':
        if svtype == 'CTX':
            svtype = 'BND'
            new_record.info['OSVTYPE'] = record.info['SVTYPE']
            new_record.info['SVTYPE'] = svtype
            new_record.info['STRANDS'] = record.info.get('STRANDS', '++')
        else:
            new_record.info['STRANDS'] = record.info['STRANDS']
        new_record.info['CHR2'] = record.info['CHR2']
        new_record.info['END2'] = bnd_end_dict[record.id] if bnd_end_dict is not None else record.info.get('END2', record.pos + record.info['SVLEN'])
        new_record.stop = record.start + 1
    elif svtype == 'INV':
        new_record.info['STRANDS'] = record.info.get('STRANDS', '++')
    elif svtype == 'CTX':
        new_record.info['STRANDS'] = record.info.get('STRANDS', '++')
        new_record.info['OSVTYPE'] = record.info['SVTYPE']
        new_record.info['SVTYPE'] = 'BND'
    elif svtype == 'CPX':
        svtype = 'INV'
        new_record.info['STRANDS'] = record.info.get('STRANDS', '++')
        new_record.info['OSVTYPE'] = record.info['SVTYPE']
        new_record.info['SVTYPE'] = svtype
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
                new_genotype['CN'] = 1 + _cache_gt_sum(genotype['GT'])
                new_genotype['GT'] = (None,)
            elif _cache_gt_sum(genotype['GT']) == 0:
                new_genotype['GT'] = (0,)
                if svtype == 'DEL':
                    new_genotype['CN'] = 1
            else:
                new_genotype['GT'] = (1,)
                if svtype == 'DEL':
                    new_genotype['CN'] = 0
        else:
            gt_sum = _cache_gt_sum(genotype['GT'])
            if svtype == 'DUP':
                new_genotype['CN'] = 2 + gt_sum
                new_genotype['GT'] = (None, None)
            elif gt_sum == 0:
                new_genotype['GT'] = (0, 0)
                if svtype == 'DEL':
                    new_genotype['CN'] = 2
            elif gt_sum == 1:
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


def _cache_gt_sum(gt):
    s = _gt_sum_map.get(gt, None)
    if s is None:
        s = sum([1 for a in gt if a is not None and a > 0])
        _gt_sum_map[gt] = s
    return s


def add_cn_ecn(record: pysam.VariantRecord,
               vcf_out: pysam.VariantFile,
               ploidy_dict: Dict[Text, Dict[Text, int]]) -> pysam.VariantRecord:
    """"
    Only modifies records by adding CN and ECN INFO fields, e.g. for 'fixed' VCFs that just need
    this metadata for certain GATK tools such as SVCluster and SVConcordance

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
    header: pysam.VariantRecord
        record with CN and ECN fields added"""
    svtype = record.info['SVTYPE']
    contig = record.contig
    new_record = vcf_out.new_record(id=record.id, contig=contig, start=record.start, stop=record.stop,
                                    alleles=record.alleles, info=record.info)

    # copy FORMAT fields
    for sample in record.samples:
        genotype = record.samples[sample]
        new_genotype = new_record.samples[sample]
        for key in genotype.keys():
            new_genotype[key] = genotype[key]
        ecn = ploidy_dict[sample][contig]
        new_genotype['ECN'] = ecn
        if svtype == 'DEL':
            new_genotype['CN'] = max(0, ecn - _cache_gt_sum(genotype['GT']))
        elif svtype == 'DUP':
            new_genotype['CN'] = ecn + _cache_gt_sum(genotype['GT'])
        elif svtype == 'CNV':
            # Disambiguates non-existent and empty (i.e. ".") CN
            cn = genotype.get('CN', None)
            if cn is None:
                cn = ecn
            new_genotype['CN'] = cn
    return new_record


def filter_unsupported_type(record: pysam.VariantRecord) -> bool:
    svtype = record.info['SVTYPE']
    return svtype == 'CPX' or svtype == 'CTX'


def _parse_arg_list(arg: Text) -> List[Text]:
    if arg is None:
        return set()
    else:
        return arg.split(',')


def _process(vcf_in: pysam.VariantFile,
             vcf_out: pysam.VariantFile,
             arguments: Dict[Text, Any],
             vcf_filter: Optional[pysam.VariantFile] = None) -> None:
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
    vcf_filter: Optional[pysam.VariantFile]
        if provided, write filtered records to this vcf

    Returns
    -------
    header: pysam.VariantRecord
        record with CN and ECN fields added"""
    remove_formats = set(_parse_arg_list(arguments.remove_formats))
    remove_infos = set(_parse_arg_list(arguments.remove_infos))
    if not arguments.only_add_cn_fields and not arguments.use_end2:
        bnd_end_dict = _parse_bnd_ends(arguments.vcf)
    else:
        bnd_end_dict = None
    ploidy_dict = _parse_ploidy_table(arguments.ploidy_table)

    # info fields we drop by default (unless needed for certain SV types)
    default_remove_infos = set(["SVLEN", "STRANDS", "CHR2"])
    if bnd_end_dict is not None:
        default_remove_infos.add("END2")
    remove_infos = remove_infos.union(default_remove_infos)

    for record in vcf_in:
        if arguments.filter_unsupported_types and filter_unsupported_type(record):
            if vcf_filter is not None:
                vcf_filter.write(record)
        else:
            if arguments.only_add_cn_fields:
                out = add_cn_ecn(record=record, vcf_out=vcf_out, ploidy_dict=ploidy_dict)
            else:
                out = convert(
                    record=record,
                    vcf_out=vcf_out,
                    remove_infos=remove_infos,
                    remove_formats=remove_formats,
                    bnd_end_dict=bnd_end_dict,
                    ploidy_dict=ploidy_dict
                )
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
    parser.add_argument("--only-add-cn-fields", action='store_true',
                        help="Only add CN and ECN info fields. All other corrections are skipped.")
    parser.add_argument("--use-end2", action='store_true',
                        help="Use existing END2 fields rather than getting them from END")
    parser.add_argument("--filter-unsupported-types", action='store_true',
                        help="Filter CPX and CTX types, which are not currently supported by GATK")
    parser.add_argument("--filter-out", type=str,
                        help="Write any filtered variants to the specified VCF")
    parser.add_argument("--replace-ev-format", action='store_true',
                        help="Adds EV FORMAT field with unbounded Number to header")
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
    arguments = _parse_arguments(argv)
    remove_formats = set(_parse_arg_list(arguments.remove_formats))
    remove_infos = set(_parse_arg_list(arguments.remove_infos))

    # convert vcf header and records
    with pysam.VariantFile(arguments.vcf) as vcf_in:
        header = create_header(
            header_in=vcf_in.header,
            replace_ev_format=arguments.replace_ev_format,
            remove_infos=remove_infos,
            remove_formats=remove_formats
        )
        with pysam.VariantFile(arguments.out, mode='w', header=header) as vcf_out:
            vcf_filter = pysam.VariantFile(arguments.filter_out, mode='w', header=vcf_in.header) if \
                arguments.filter_out is not None else None
            _process(vcf_in, vcf_out, arguments, vcf_filter=vcf_filter)
            if vcf_filter is not None:
                vcf_filter.close()


if __name__ == "__main__":
    main()
