#!/bin/env python

"""
This script creates a VCF from a Scramble output table
"""

import argparse
from collections import defaultdict, deque
import gzip
import logging
import os
import sys
from typing import Optional, List, Text

from intervaltree import IntervalTree
import pandas as pd
import pysam

"""
Creates a VCF from the Scramble output table. Applies critical filters to the raw calls, clusters redundant calls,
 and estimates variant length based on MEI type and alignment information.
"""

CIGAR_TUPLE_INS_INDEX = 1
CIGAR_TUPLE_DEL_INDEX = 2

FLOAT_COLUMNS = ["Alignment_Score", "Alignment_Percent_Length", "Alignment_Percent_Identity"]
INT_COLUMNS = ["pos", "end", "Clipped_Reads_In_Cluster", "Start_In_MEI", "Stop_In_MEI", "polyA_Position",
               "polyA_SupportingReads", "TSD_length"]


def make_header(reference_path, sample_id):
    """
    Creates the header
    Args:
        reference_path: str
            Local path to the reference FASTA file
        sample_id: str
            Sample name
    Returns:
        header: VariantHeader
            New VCF header
    """
    reference_filename = os.path.basename(reference_path)
    header_lines = [
        "##fileformat=VCFv4.3",
        f"##reference={reference_filename}",
        "##source=scramble",
        "##ALT=<ID=INS:ME:LINE1,Description=\"LINE1 element insertion\">",
        "##ALT=<ID=INS:ME:SVA,Description=\"SVA element insertion\">",
        "##ALT=<ID=INS:ME:ALU,Description=\"ALU element insertion\">",
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this "
        "record\">",
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
        "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT "
        "alleles\">",
        "##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description=\"Source algorithms\">",
        "##INFO=<ID=STRANDS,Number=1,Type=String,Description=\"Breakpoint strandedness [++,+-,-+,--]\">",
        "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate\">",
        "##INFO=<ID=MEI_START,Number=1,Type=Integer,Description=\"Start of alignment to canonical MEI "
        "sequence\">",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    ]
    reference_index_path = reference_path + ".fai"
    with open(reference_index_path) as f:
        for line in f:
            tokens = line.strip().split('\t')
            seq_name = tokens[0]
            seq_len = int(tokens[1])
            header_lines.append(f"##contig=<ID={seq_name},length={seq_len}>")
    header = pysam.VariantHeader()
    for line in header_lines:
        header.add_line(line)
    header.add_sample(sample_id)
    return header


def read_table(file_lines, ref_path, cluster_distance, alu_size, sva_size, l1_size):
    """
    Parses input table of Scramble calls. Performs call deduplication and SVLEN estimation.
    Args:
        file_lines: list of str
            Table file lines
        ref_path: str
            Local path to reference FASTA file
        cluster_distance: int
            Distance in bp for clustering redundant calls
        alu_size: int
            Size in bp of ALU elements
        sva_size: int
            Size in bp of SVA elements
        l1_size: int
            Size in bp of LINE1 elements
    Returns:
        pandas DataFrame
            Unfiltered calls
    """
    family_lengths = {
        "l1": l1_size,
        "sva": sva_size,
        "alu": alu_size
    }

    # See discussion for description: https://github.com/broadinstitute/gatk-sv/discussions/734
    def _calculate_svlen_one_sided(r):
        if r['Clipped_Side'] == 'right':
            if r['Insertion_Direction'] == 'Minus':
                return r['Stop_In_MEI']
            else:
                return family_lengths[r['MEI_Family']] - r['Start_In_MEI']
        else:
            if r['Insertion_Direction'] == 'Minus':
                return family_lengths[r['MEI_Family']] - r['Start_In_MEI']
            else:
                return r['Stop_In_MEI']

    # See discussion for description: https://github.com/broadinstitute/gatk-sv/discussions/734
    def _calculate_svlen_two_sided(r1, r2):
        if r1['Clipped_Side'] == 'right':
            r_right = r1
            r_left = r2
        else:
            r_right = r2
            r_left = r1
        if r_right['Insertion_Direction'] == 'Minus':
            return r_right['Stop_In_MEI'] - r_left['Start_In_MEI']
        else:
            return r_left['Stop_In_MEI'] - r_right['Start_In_MEI']

    def _cast_numeric(row, columns, type):
        for col in columns:
            if row[col] in ['None Found', 'NA']:
                row[col] = None
            else:
                row[col] = type(row[col])

    def _record_sort_key(record, ref_sequences_index_dict):
        if record['chrom'] not in ref_sequences_index_dict:
            raise ValueError(f"Record contig \"{record['chrom']}\" not found in reference dictionary")
        return [ref_sequences_index_dict[record['chrom']], record['pos'], record['end']]

    with pysam.FastaFile(ref_path) as f:
        ref_sequences_index_dict = {seq: i for i, seq in enumerate(f.references)}

    # Buffer of active items within cluster_distance
    buffer = deque()
    data = deque()
    columns = None
    for line in file_lines:
        if not line:
            continue
        record = line.strip().split('\t')
        if columns is None:
            columns = record
            # Remove leading '#'
            columns[0] = columns[0][1:]
            continue
        row = {columns[i]: record[i] for i in range(len(record))}

        row['sides'] = 1
        _cast_numeric(row, FLOAT_COLUMNS, float)
        _cast_numeric(row, INT_COLUMNS, int)
        new_buffer = deque()
        found_match = False
        for item in buffer:
            if found_match:
                new_buffer.append(item)
            elif item['chrom'] != row['chrom'] \
                    or abs(int(item['pos']) - int(row['pos'])) > cluster_distance:
                item['svlen'] = _calculate_svlen_one_sided(item)
                data.appendleft(item)
            elif item['Clipped_Side'] != row['Clipped_Side'] \
                    and item['Insertion_Direction'] == row['Insertion_Direction']:
                row['pos'] = item['pos']
                row['svlen'] = _calculate_svlen_two_sided(item, row)
                data.append(row)
                found_match = True
            else:
                new_buffer.appendleft(item)
        if not found_match:
            new_buffer.appendleft(row)
        buffer = new_buffer
    for item in buffer:
        item['svlen'] = _calculate_svlen_one_sided(item)
        data.append(item)
    # Sort records in memory
    data = sorted(list(data), key=lambda x: _record_sort_key(record=x, ref_sequences_index_dict=ref_sequences_index_dict))
    return pd.DataFrame(data=data)


def fails_indel_filter(record, samfile, mei_trees, indel_window, min_reads, min_ins_size, max_ins_size,
                       min_del_size, max_del_size):
    """
    Checks if a given record lies within a reference MEI and is sufficiently close to a nearby indel,
    based on local read alignments, which is a common error mode for Scramble. The check for indels is to count
    the number of insertions and deletions occurring anywhere in the local vicinity of the breakpoint and determine
    if either count exceeds the specified threshold. Note that indel read counts are tallied across the entire window,
    as identifying individual indels within simple repeats is very difficult (CIGAR indels tend to be scattered). Also,
    the indels are tallied over all reads overlapping the window, so the full search window is indel_window + read_length.
    Args:
        record: VariantRecord
            Variant in question
        samefile: AlignmentFile
            SAM/BAM/CRAM file opened with pysam
        mei_trees: dict of IntervalTree
            Reference MEI intervals
        indel_window: int
            Window in bp to look for indels around the insertion site
        min_reads: int
            Required number of reads containing an indel to filter
        min_ins_size: int
            Minimum insertion length to look for in the read CIGAR
        max_ins_size: int
            Maximum insertion length to look for in the read CIGAR
        min_del_size: int
            Minimum deletion length to look for in the read CIGAR
        max_del_size: int
            Maximum deletion length to look for in the read CIGAR
    Returns:
        True if the call should be filtered
    """
    def _get_num_indels(reads, min_reads, min_ins_size, max_ins_size, min_del_size, max_del_size):
        ins_count = 0
        del_count = 0
        for read in reads:
            if read.cigartuples is None:
                continue
            for tup in read.cigartuples:
                if tup[0] == CIGAR_TUPLE_INS_INDEX and min_ins_size <= tup[1] <= max_ins_size:
                    ins_count += 1
                    if ins_count >= min_reads:
                        return ins_count
                elif tup[0] == CIGAR_TUPLE_DEL_INDEX and min_del_size <= tup[1] <= max_del_size:
                    del_count += 1
                    if del_count >= min_reads:
                        return del_count
        return max(ins_count, del_count)

    if min_reads < 1:
        raise ValueError("Minimum number of indel reads must be positive")
    if min_ins_size < 1:
        raise ValueError("Minimum insertion size must be positive")
    if max_ins_size < 1:
        raise ValueError("Maximum insertion size must be positive")
    if min_del_size < 1:
        raise ValueError("Minimum del size must be positive")
    if max_del_size < 1:
        raise ValueError("Maximum del size must be positive")
    if record.chrom not in mei_trees or len(mei_trees[record.chrom].overlap(record.pos, record.pos + 1)) == 0:
        return False
    start = max(record.pos - indel_window, 1)
    end = record.pos + indel_window
    num_indels = _get_num_indels(reads=samfile.fetch(record.chrom, start, end), min_reads=min_reads,
                                 min_ins_size=min_ins_size, max_ins_size=max_ins_size,
                                 min_del_size=min_del_size, max_del_size=max_del_size)
    return num_indels >= min_reads


def filter_and_write_vcf(vcf, samfile, df, sample_id, enable_indel_filter, mei_trees, del_filter_trees,
                         indel_window, indel_reads, min_ins_size, max_ins_size, min_del_size, max_del_size):
    """
    Iterates through calls in the given DataFrame, creating VCF variant records. Applies indel and SV deletion filters.
    Writes resulting records to the VCF. All variants are genotyped as diploid/heterozygous.
    Args:
        vcf: VariantFile
            Output VCF
        samfile: AlignmentFile
            SAM/BAM/CRAM file opened with pysam
        df: DataFrame
            Scramble calls
        sample_id: str
            Sample name
        enable_indel_filter: bool
            Indel filter flag
        mei_trees: dict of IntervalTree
            Reference MEI intervals
        del_filter_trees: dict of IntervalTree
            Deletion calls, used for filtering
        indel_window: int
            Window in bp to look for indels around the insertion site
        indel_reads: int
            Required number of reads containing an indel to filter
        min_ins_size: int
            Minimum insertion length to look for in the read CIGAR
        max_ins_size: int
            Maximum insertion length to look for in the read CIGAR
        min_del_size: int
            Minimum deletion length to look for in the read CIGAR
        max_del_size: int
            Maximum deletion length to look for in the read CIGAR
    Returns:
        None
    """
    for index, row in df.iterrows():
        if row['MEI_Family'] == 'alu':
            allele = '<INS:ME:ALU>'
        elif row['MEI_Family'] == 'sva':
            allele = '<INS:ME:SVA>'
        elif row['MEI_Family'] == 'l1':
            allele = '<INS:ME:LINE1>'
        else:
            raise ValueError(f"Unrecognized MEI_Family: {row['MEI_Family']}")
        record = vcf.new_record(
            contig=row['chrom'],
            start=row['pos'],
            stop=row['pos'] + 1,
            alleles=['N', allele],
            id=f"scramble_{sample_id}_{row['chrom']}_{index}"
        )
        if enable_indel_filter and fails_indel_filter(record=record, samfile=samfile, mei_trees=mei_trees,
                                                      indel_window=indel_window, min_reads=indel_reads,
                                                      min_ins_size=min_ins_size, max_ins_size=max_ins_size,
                                                      min_del_size=min_del_size, max_del_size=max_del_size):
            print(f"Filtering site at {row['chrom']}:{row['pos']} due to nearby INDEL")
            continue
        if record.chrom in del_filter_trees and \
                len(del_filter_trees[record.chrom].overlap(record.pos, record.pos + 1)) > 0:
            print(f"Filtering site at {row['chrom']}:{row['pos']} due to nearby DEL")
            continue
        record.info['SVTYPE'] = 'INS'
        record.info['SVLEN'] = row['svlen']
        record.info['ALGORITHMS'] = ['scramble']
        record.samples[sample_id]['GT'] = (0, 1)
        vcf.write(record)


def create_trees_from_bed_records(path, padding):
    """
    Creates a set of interval trees (one per contig) from a bed file at the given path. Adds padding to each interval.
    Args:
        path: str
            Path to reference MEI bed file
        padding: int
            Padding to apply to both ends of every interval, in bp
    Returns:
        Dict[contig str] -> IntervalTree
    """
    open_fn = get_open_read_function(path)
    with open_fn(path) as f:
        trees = defaultdict(IntervalTree)
        for line in f:
            record = line.strip().split('\t')
            contig = record[0]
            start = int(record[1])
            end = int(record[2])
            trees[contig].addi(start - padding, end + padding)
        return trees


def add_del_ends_to_trees(vcf, trees, padding):
    """
    Iterates through a VCF of SV calls and adds all deletions to a set of interval trees.
    Args:
        vcf: VariantFile
            SV VCF
        trees:Dict[contig] -> IntervalTree
            Interval trees
        padding: int
            Padding to apply to both ends of every deletion interval, in bp
    Returns:
        Original trees with added intervals (not a copy)
    """
    for record in vcf:
        if record.info.get("SVTYPE", "") == "DEL":
            if record.chrom not in trees:
                trees[record.chrom] = IntervalTree()
            trees[record.chrom].addi(record.pos - padding, record.pos + padding)
            trees[record.chrom].addi(record.stop - padding, record.stop + padding)
    return trees


def get_open_read_function(path):
    """
    Helper for opening a path that may be gzipped, ending in .gz
    Args:
        path: str
            File path
    Returns:
        Method to open the file
    """
    if path.endswith('.gz'):
        return gzip_open_read_text
    else:
        return open


def gzip_open_read_text(path):
    """
    Helper for opening a gzipped file
    Args:
        path: str
            File path
    Returns:
        File handle
    """
    return gzip.open(path, "rt")


def preprocess_scramble_table(path):
    """
    Parses the Scramble calls table.
    Args:
        path: str
            Scramble tsv table path
    Returns: list[str]
        List of file lines, including header
    """
    open_fn = get_open_read_function(path)
    with open_fn(path) as f_in:
        header = "#chrom\tpos\tend\t" + f_in.readline()
        output = [header]
        for line in f_in:
            tokens = line.strip().split('\t')
            chrom, pos = tokens[0].split(':')
            pos = int(pos)
            output.append(f"{chrom}\t{pos}\t{pos + 1}\t" + line)
        return output


def __parse_arguments(argv: List[Text]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert Scramble output into a VCF",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--table", type=str, required=True,
                        help="Scramble variant table (may be gzipped)")
    parser.add_argument("--alignments-file", type=str, required=True,
                        help="BAM/CRAM file, must be indexed")
    parser.add_argument("--mei-bed", type=str, required=True,
                        help="Bed file containing MEI intervals from the reference")
    parser.add_argument("--input-vcf", type=str, required=True,
                        help="Input VCF (e.g.from DRAGEN-SV or Manta)")
    parser.add_argument("--sample", type=str, required=True,
                        help="Sample ID")
    parser.add_argument("--reference", type=str, required=True,
                        help="Reference fasta, must be indexed")
    parser.add_argument("--out", type=str, required=True,
                        help="Output VCF")
    parser.add_argument("--alu-size", type=int, default=282,
                        help="Alu size")
    parser.add_argument("--sva-size", type=int, default=1362,
                        help="SVA size")
    parser.add_argument("--l1-size", type=int, default=6023,
                        help="LINE1 size")
    parser.add_argument("--cluster-distance", type=int, default=300,
                        help="Maximum distance to collapse pairs of break-end calls")
    parser.add_argument("--del-filter-window", type=int, default=50,
                        help="Window for SV deletion breakend filtering")
    parser.add_argument("--disable-indel-filter", action='store_true', default=False,
                        help="Enable filtering of calls near indels")
    parser.add_argument("--indel-window", type=int, default=200,
                        help="Window for indel filtering")
    parser.add_argument("--indel-reads", type=int, default=3,
                        help="Min number of indel reads required for indel filter")
    parser.add_argument("--small-ins-min-size", type=int, default=5,
                        help="Min small insertion size required for indel filter")
    parser.add_argument("--small-ins-max-size", type=int, default=50,
                        help="Max small insertion size allowed for indel filter")
    parser.add_argument("--small-del-min-size", type=int, default=5,
                        help="Min small deletion size required for indel filter")
    parser.add_argument("--small-del-max-size", type=int, default=1000,
                        help="Max small deletion size allowed for indel filter")
    parser.add_argument("--mei-padding", type=int, default=0,
                        help="Reference MEI padding in bp")
    parser.add_argument("--log-level", required=False, default="INFO", help="Specify level of logging information")
    if len(argv) <= 1:
        parser.parse_args(["--help"])
        sys.exit(0)
    parsed_arguments = parser.parse_args(argv[1:])
    return parsed_arguments


def main(argv: Optional[List[Text]] = None):
    if argv is None:
        argv = sys.argv
    arguments = __parse_arguments(argv)

    # Set logging level from --log-level input
    log_level = arguments.log_level
    numeric_level = getattr(logging, log_level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % log_level)
    logging.basicConfig(level=numeric_level, format='%(asctime)s - %(levelname)s - %(message)s')

    header = make_header(reference_path=arguments.reference, sample_id=arguments.sample)
    logging.info("Loading Scramble table...")
    bed_lines = preprocess_scramble_table(path=arguments.table)
    df = read_table(file_lines=bed_lines,
                    ref_path=arguments.reference,
                    cluster_distance=arguments.cluster_distance,
                    alu_size=arguments.alu_size,
                    sva_size=arguments.sva_size,
                    l1_size=arguments.l1_size)
    logging.info("Loading MEI bed...")
    mei_trees = create_trees_from_bed_records(arguments.mei_bed, padding=arguments.mei_padding)
    logging.info("Loading deletions...")
    with pysam.VariantFile(arguments.input_vcf) as f_vcf:
        del_filter_trees = dict()
        add_del_ends_to_trees(vcf=f_vcf, trees=del_filter_trees, padding=arguments.del_filter_window)
    logging.info("Writing vcf...")
    with pysam.VariantFile(arguments.out, "w", header=header) as vcf, \
            pysam.AlignmentFile(arguments.alignments_file, reference_filename=arguments.reference) as samfile:
        filter_and_write_vcf(vcf=vcf, samfile=samfile, df=df, sample_id=arguments.sample,
                             enable_indel_filter=not arguments.disable_indel_filter,
                             mei_trees=mei_trees,
                             del_filter_trees=del_filter_trees,
                             indel_window=arguments.indel_window, indel_reads=arguments.indel_reads,
                             min_ins_size=arguments.small_ins_min_size,
                             max_ins_size=arguments.small_ins_max_size,
                             min_del_size=arguments.small_del_min_size,
                             max_del_size=arguments.small_del_max_size)
    pysam.tabix_index(arguments.out, preset="vcf", force=True)


if __name__ == "__main__":
    main()
