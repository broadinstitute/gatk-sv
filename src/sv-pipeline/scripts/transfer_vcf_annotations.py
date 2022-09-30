#!/bin/env python

import argparse
import re
import pysam
import sys
from typing import Optional, List, Text, Set, Dict, Tuple, Iterable

from more_itertools import peekable


class PeekIterator:

    def __init__(self, iterable):
        self.iterator = iter(iterable)
        self.peeked = deque()

    def __iter__(self):
        return self

    def __next__(self):
        if self.peeked:
            return self.peeked.popleft()
        return next(self.iterator)

    def peek(self, ahead=0):
        while len(self.peeked) <= ahead:
            self.peeked.append(next(self.iterator))
        return self.peeked[ahead]


def create_header(header_in: pysam.VariantHeader,
                  header_ann: pysam.VariantHeader,
                  info_keys: Set[Text],
                  format_keys: Set[Text]) -> pysam.VariantHeader:
    """
    Adds annotation metadata lines to the input header. Also checks that the fields exist in the annotation vcf.

    Parameters
    ----------
    header_in: pysam.VariantHeader
        input header
    header_ann: pysam.VariantHeader
        annotation vcf header
    info_keys: Set[Text]
        set of info fields to annotate
    format_keys: Set[Text]
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
    for k in info_keys:
        if k not in header_ann.info:
            raise ValueError(f"INFO field {k} not found in annotation vcf header")
    for k in format_keys:
        if k not in header_ann.formats:
            raise ValueError(f"FORMAT field {k} not found in annotation vcf header")

    # Add annotation headerlines
    for line in header_ann.records:
        if len(line.attrs) > 0 and 'ID' in line.keys() and \
                (line['ID'] in info_keys or line['ID'] in format_keys):
            header.add_line(str(line))
    return header


class ContigComparator(str):
    """
    Comparator for contig name sorting, equivalent to 'sort -V'
    """
    def __lt__(self, other):
        num_self = re.sub("[^0-9]", "", self)
        num_other = re.sub("[^0-9]", "", other)
        if num_self == "" and num_other == "":
            return str(self) < str(other)  # e.g. chrX, chrY
        elif num_self == "":
            return False
        elif num_other == "":
            return True
        else:
            return int(num_self) < int(num_other)


def pairs(vcf_a: pysam.VariantFile,
          vcf_b: pysam.VariantFile,
          header: pysam.VariantHeader) -> Iterable[Tuple]:
    """
    Zips records from vcf_a and vcf_b, matching by variant ID. Input vcfs must be sorted

    Parameters
    ----------
    vcf_a: pysam.VariantFile
        first vcf
    vcf_b: pysam.VariantFile
        second vcf
    header: pysam.VariantHeader
        header defining reference contigs

    Returns
    -------
    Iterable[Tuple]
        paired records (record_a, record_b); one record may be None if no match was found
    """

    def _cmp_coord(a_chrom: Text,
                   a_pos: int,
                   b_chrom: Text,
                   b_pos: int,
                   contig_indices: Dict[Text, int]) -> int:
        """
        Comparing function for two genomic coordinates
            a_chrom: Text
                First contig
            a_pos: int
                First position
            b_chrom: Text
                Second contig
            b_pos: int
                Second position
            contig_indices: Dict[Text, int]
                Specifies contig order
        """
        x_chrom_i = contig_indices[a_chrom]
        y_chrom_i = contig_indices[b_chrom]
        if x_chrom_i < y_chrom_i:
            return -1
        elif x_chrom_i > y_chrom_i:
            return 1
        elif a_pos < b_pos:
            return -1
        elif a_pos > b_pos:
            return 1
        else:
            return 0

    vcf_a = peekable(vcf_a)
    vcf_b = peekable(vcf_b)

    # Sometimes contigs aren't sorted properly in the header
    sorted_contig_list = [c for c in header.contigs]
    sorted_contig_list.sort(key=ContigComparator)
    contig_indices = {c: sorted_contig_list.index(c) for c in sorted_contig_list}

    chrom = None
    pos = None
    a_buff = dict()  # buffer records at the current coordinate (chrom, pos)
    b_buff = dict()
    while vcf_a and vcf_b:
        x = vcf_a.peek()
        y = vcf_b.peek()
        if chrom is None:
            cxy = _cmp_coord(x.chrom, x.pos, y.chrom, y.pos, contig_indices)
            if cxy <= 0:
                chrom = x.chrom
                pos = x.pos
            else:
                chrom = y.chrom
                pos = y.pos
        else:
            cxy = _cmp_coord(x.chrom, x.pos, y.chrom, y.pos, contig_indices)
            cx = _cmp_coord(x.chrom, x.pos, chrom, pos, contig_indices)
            cy = _cmp_coord(y.chrom, y.pos, chrom, pos, contig_indices)
            if cx < 0:
                raise ValueError(f"Position of {x.id} precedes current position {chrom}:{pos}. "
                                 f"Check that input vcf is sorted.")
            elif cy < 0:
                raise ValueError(f"Position of {y.id} precedes current position {chrom}:{pos}. "
                                 f"Check that annotations vcf is sorted.")
            if cx == 0:
                a_buff[x.id] = next(vcf_a)
            if cy == 0:
                b_buff[y.id] = next(vcf_b)
            if cx > 0 and cy > 0:
                for key, val in a_buff.items():
                    yield val, b_buff.get(key, None)
                a_buff = dict()
                b_buff = dict()
                if cxy <= 0:
                    chrom = x.chrom
                    pos = x.pos
                else:
                    chrom = y.chrom
                    pos = y.pos
    for key, val in a_buff.items():
        yield val, b_buff.get(key, None)
    for x in vcf_a:
        yield x, None
    for y in vcf_b:
        yield None, y


def annotate_record(in_record: pysam.VariantRecord,
                    ann_record: pysam.VariantRecord,
                    vcf_out: pysam.VariantFile,
                    info_keys: Set[Text],
                    format_keys: Set[Text]) -> pysam.VariantRecord:
    """
    Annotates record with INFO and FORMAT fields

    Parameters
    ----------
    in_record: pysam.VariantRecord
        base record
    ann_record: pysam.VariantRecord
        annotation source record
    vcf_out: pysam.VariantFile
        output vcf object
    info_keys: List[Text]
        list of format keys as ordered in format_dict
    format_keys: List[Text]
        list of format keys as ordered in format_dict

    Returns
    -------
    pysam.VariantRecord
        annotated record
    """
    new_record = vcf_out.new_record(alleles=in_record.alleles, contig=in_record.contig, filter=in_record.filter,
                                    id=in_record.id, info=in_record.info, start=in_record.start,
                                    qual=in_record.qual, stop=in_record.stop)
    sample_list = sorted(list(set(vcf_out.header.samples).intersection(set(ann_record.samples))))
    for sample, dest in new_record.samples.items():
        source = in_record.samples[sample]
        for key, val in source.items():
            dest[key] = val
    for key in info_keys:
        val = ann_record.info.get(key, None)
        if in_record.info['SVTYPE'] == 'INS' and key == 'TRUTH_AF':
            print(val)
        if val is not None and not (isinstance(val, Iterable) and val[0] is None):
            new_record.info[key] = val
    for sample_id in sample_list:
        ann_format = ann_record.samples[sample_id]
        new_sample_format = new_record.samples[sample_id]
        for key in format_keys:
            val = ann_format.get(key, None)
            if val is not None and not (isinstance(val, Iterable) and val[0] is None):
                new_sample_format[key] = val
    return new_record


def annotate_vcf(vcf_in: pysam.VariantFile,
                 vcf_ann: pysam.VariantFile,
                 vcf_out: pysam.VariantFile,
                 info_keys: Set[Text],
                 format_keys: Set[Text]) -> None:
    """
    Transfers annotations from one vcf to another, matching on variant ID. Result is written to the given output vcf.

    Parameters
    ----------
    vcf_in: pysam.VariantFile
        Base vcf to add annotations to
    vcf_ann: pysam.VariantFile
        Vcf containing source annotations
    vcf_out: pysam.VariantFile
        Output vcf
    info_keys: Set[Text]
        Annotation INFO keys
    format_keys: Set[Text]
        Annotation FORMAT keys
    """
    for in_rec, ann_rec in pairs(vcf_in, vcf_ann, vcf_out.header):
        if in_rec is None:
            continue
        if ann_rec is None:
            vcf_out.write(in_rec)
        else:
            if in_rec.id != ann_rec.id:
                raise ValueError(f"Record ids do not match: {in_rec.id} {ann_rec.id}")
            if in_rec.chrom != ann_rec.chrom or in_rec.pos != ann_rec.pos:
                raise ValueError(f"Record positions do not match for {in_rec.id}: "
                                 f"{in_rec.chrom}:{in_rec.pos} {ann_rec.chrom}:{ann_rec.pos}")
            vcf_out.write(annotate_record(in_record=in_rec, ann_record=ann_rec, vcf_out=vcf_out,
                                          info_keys=info_keys, format_keys=format_keys))


def _parse_arg_list(arg: Text) -> List[Text]:
    if arg is None:
        return list()
    else:
        return arg.split(',')


def _parse_arguments(argv: List[Text]) -> argparse.Namespace:
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
    arguments = _parse_arguments(argv)
    info_keys = set(_parse_arg_list(arguments.infos))
    format_keys = set(_parse_arg_list(arguments.formats))
    if len(info_keys) + len(format_keys) == 0:
        raise ValueError('No INFO or FORMAT fields were specified')

    with pysam.VariantFile(arguments.vcf) as vcf_in, pysam.VariantFile(arguments.ann_vcf) as vcf_ann:
        header = create_header(
            header_in=vcf_in.header,
            header_ann=vcf_ann.header,
            info_keys=info_keys,
            format_keys=format_keys
        )
        with pysam.VariantFile(arguments.out, mode='w', header=header) as vcf_out:
            annotate_vcf(vcf_in=vcf_in, vcf_ann=vcf_ann, vcf_out=vcf_out,
                         info_keys=info_keys, format_keys=format_keys)


if __name__ == "__main__":
    main()
