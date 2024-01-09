#!/bin/env python

"""
This script creates a VCF from a Scramble output table
"""

import argparse
import os
import sys
import gzip
from typing import Optional, List, Text

import pandas as pd
import pysam

FLOAT_COLUMNS = ["Alignment_Score", "Alignment_Percent_Length", "Alignment_Percent_Identity"]
INT_COLUMNS = ["pos", "end", "Clipped_Reads_In_Cluster", "Start_In_MEI", "Stop_In_MEI", "polyA_Position",
               "polyA_SupportingReads", "TSD_length"]


def make_header(reference_path, sample_id):

    reference_filename = os.path.basename(reference_path)
    header_lines = [
        f"##fileformat=VCFv4.3",
        f"##reference={reference_filename}",
        f"##source=scramble",
        f"##ALT=<ID=INS:ME:LINE1,Description=\"LINE1 element insertion\">",
        f"##ALT=<ID=INS:ME:SVA,Description=\"SVA element insertion\">",
        f"##ALT=<ID=INS:ME:ALU,Description=\"ALU element insertion\">",
        f"##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this "
        f"record\">",
        f"##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
        f"##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT "
        f"alleles\">",
        f"##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description=\"Source algorithms\">",
        f"##INFO=<ID=STRANDS,Number=1,Type=String,Description=\"Breakpoint strandedness [++,+-,-+,--]\">",
        f"##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate\">",
        f"##INFO=<ID=MEI_START,Number=1,Type=Integer,Description=\"Start of alignment to canonical MEI "
        f"sequence\">",
        f"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
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


def read_table(path, cluster_distance, alu_size, sva_size, l1_size):

    family_lengths = {
        "l1": l1_size,
        "sva": sva_size,
        "alu": alu_size
    }

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

    if path.endswith('.gz'):
        open_fn = gzip.open
        decoder = lambda x: str(x, 'utf-8')
    else:
        open_fn = open
        decoder = lambda x: x
    with open_fn(path) as f:
        # Buffer of active items within cluster_distance
        buffer = list()
        data = list()
        columns = None
        for line in f:
            tokens = decoder(line).strip().split('\t')
            if columns is None:
                columns = tokens
                continue
            row = {columns[i]: tokens[i] for i in range(len(tokens))}
            row['sides'] = 1
            _cast_numeric(row, FLOAT_COLUMNS, float)
            _cast_numeric(row, INT_COLUMNS, int)
            new_buffer = list()
            found_match = False
            for item in buffer:
                if found_match:
                    new_buffer.append(item)
                elif item['chrom'] != row['chrom'] \
                        or abs(int(item['pos']) - int(row['pos'])) > cluster_distance:
                    item['svlen'] = _calculate_svlen_one_sided(item)
                    data.insert(0, item)
                elif item['Clipped_Side'] != row['Clipped_Side'] \
                        and item['Insertion_Direction'] == row['Insertion_Direction']:
                    row['pos'] = item['pos']
                    row['svlen'] = _calculate_svlen_two_sided(item, row)
                    data.append(row)
                    found_match = True
                else:
                    new_buffer.insert(0, item)
            if not found_match:
                new_buffer.insert(0, row)
            buffer = new_buffer
        for item in buffer:
            item['svlen'] = _calculate_svlen_one_sided(item)
            data.append(item)
        return pd.DataFrame(data=data)


def write_vcf(vcf, df, sample_id):
    for index, row in df.iterrows():
        if row['MEI_Family'] == 'alu':
            allele = '<INS:ME:ALU>'
        elif row['MEI_Family'] == 'sva':
            allele = '<INS:ME:SVA>'
        elif row['MEI_Family'] == 'l1':
            allele = '<INS:ME:L1>'
        else:
            raise ValueError(f"Unrecognized MEI_Family: {row['MEI_Family']}")
        record = vcf.new_record(
            contig=row['chrom'],
            start=row['pos'],
            stop=row['pos'],
            alleles=['N', allele],
            id=f"scramble_{sample_id}_{row['chrom']}_{index}"
        )
        record.info['SVTYPE'] = 'INS'
        record.info['SVLEN'] = row['svlen']
        record.info['ALGORITHMS'] = ['scramble']
        record.samples[sample_id]['GT'] = (0, 1)
        vcf.write(record)

def __parse_arguments(argv: List[Text]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert Scramble output into a VCF",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--table", type=str, required=True,
                        help="Scramble variant table (may be gzipped)")
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
    if len(argv) <= 1:
        parser.parse_args(["--help"])
        sys.exit(0)
    parsed_arguments = parser.parse_args(argv[1:])
    return parsed_arguments


def main(argv: Optional[List[Text]] = None):
    if argv is None:
        argv = sys.argv
    arguments = __parse_arguments(argv)
    header = make_header(reference_path=arguments.reference, sample_id=arguments.sample)
    df = read_table(path=arguments.table,
                    cluster_distance=arguments.cluster_distance,
                    alu_size=arguments.alu_size,
                    sva_size=arguments.sva_size,
                    l1_size=arguments.l1_size)
    with pysam.VariantFile(arguments.out, 'w', header=header) as f:
        write_vcf(vcf=f, df=df, sample_id=arguments.sample)


if __name__ == "__main__":
    main()
