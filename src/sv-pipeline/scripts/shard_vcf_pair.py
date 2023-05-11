#!/bin/env python

import argparse
import os
import pysam
import re
import sys
from typing import Optional, List, Text, Dict, Tuple, Iterable

from more_itertools import peekable


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


def shard_vcfs(vcf_a: pysam.VariantFile,
               vcf_b: pysam.VariantFile,
               out_dir: Text,
               prefix_a: Text,
               prefix_b: Text,
               shard_size: int,
               drop_a: bool,
               drop_b: bool) -> None:
    """
    Transfers annotations from one vcf to another, matching on variant ID. Result is written to the given output vcf.

    Parameters
    ----------
    vcf_in: pysam.VariantFile
        Base vcf to add annotations to
    vcf_ann: pysam.VariantFile
        Vcf containing source annotations
    prefix_a: Text
        First vcf file prefix
    prefix_b: Text
        Second vcf file prefix
    format_keys: Set[Text]
        Annotation FORMAT keys
    """

    def _get_shard_files(dir: Text,
                         header_a: pysam.VariantHeader,
                         header_b: pysam.VariantHeader,
                         prefix_a: Text,
                         prefix_b: Text,
                         index: int) -> Text:
        path_a = os.path.join(dir, '{}_{:08d}.vcf.gz'.format(prefix_a, index))
        path_b = os.path.join(dir, '{}_{:08d}.vcf.gz'.format(prefix_b, index))
        return pysam.VariantFile(path_a, mode='w', header=header_a), \
            pysam.VariantFile(path_b, mode='w', header=header_b)

    n_a = 0  # record counts for the current shard
    n_b = 0
    shard_index = 0
    out_a, out_b = _get_shard_files(dir=out_dir, header_a=vcf_a.header, header_b=vcf_b.header,
                                    prefix_a=prefix_a, prefix_b=prefix_b, index=shard_index)
    for rec_a, rec_b in pairs(vcf_a, vcf_b, vcf_a.header):
        if rec_a is not None and not (drop_a and rec_b is None):
            out_a.write(rec_a)
            n_a += 1
        if rec_b is not None and not (drop_b and rec_a is None):
            out_b.write(rec_b)
            n_b += 1
        if max(n_a, n_b) >= shard_size:
            shard_index += 1
            n_a = 0
            n_b = 0
            out_a.close()
            out_b.close()
            out_a, out_b = _get_shard_files(dir=out_dir, header_a=vcf_a.header, header_b=vcf_b.header,
                                            prefix_a=prefix_a, prefix_b=prefix_b, index=shard_index)
    out_a.close()
    out_b.close()


def _parse_arg_list(arg: Text) -> List[Text]:
    if arg is None:
        return list()
    else:
        return arg.split(',')


def _parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Splits records from two vcfs with overlapping variant IDs, "
                    "ensuring any two records with the same ID end up in the same shard",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("vcf_a", type=str, help="First vcf")
    parser.add_argument("vcf_b", type=str, help="Second vcf")
    parser.add_argument("--out-dir", type=str, default="./", help="Output directory")
    parser.add_argument("--prefix-a", type=str, required=True, help="First vcf's filename prefix")
    parser.add_argument("--prefix-b", type=str, required=True, help="Second vcf's filename prefix")
    parser.add_argument("--shard-size", type=int, default=30000, help="Target shard size")
    parser.add_argument("--drop-a", action='store_true',
                        help="Drop records from vcf_a without matching IDs in vcf_b")
    parser.add_argument("--drop-b", action='store_true',
                        help="Drop records from vcf_b without matching IDs in vcf_a")
    if len(argv) <= 1:
        parser.parse_args(["--help"])
        sys.exit(0)
    parsed_arguments = parser.parse_args(argv[1:])
    return parsed_arguments


def main(argv: Optional[List[Text]] = None):
    if argv is None:
        argv = sys.argv
    arguments = _parse_arguments(argv)
    if arguments.shard_size <= 0:
        raise ValueError('Shard size must be positive')
    if arguments.prefix_a == arguments.prefix_b:
        raise ValueError('Prefixes cannot be equal')

    with pysam.VariantFile(arguments.vcf_a) as vcf_a, pysam.VariantFile(arguments.vcf_b) as vcf_b:
        shard_vcfs(
            vcf_a=vcf_a,
            vcf_b=vcf_b,
            out_dir=arguments.out_dir,
            prefix_a=arguments.prefix_a,
            prefix_b=arguments.prefix_b,
            shard_size=arguments.shard_size,
            drop_a=arguments.drop_a,
            drop_b=arguments.drop_b
        )


if __name__ == "__main__":
    main()
