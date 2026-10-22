#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Fixes alleles for multi-allelic CNVs. Optionally drops records with invalid coordinates.
"""

import argparse
import logging
import pysam
import sys


def get_new_header(header):
    header_list = str(header).split('\n')
    new_header_lines = list()
    for line in header_list:
        if line.startswith('##INFO=<ID=MULTIALLELIC,'):
            # Remove MULTIALLELIC field (legacy)
            continue
        elif line.startswith('#CHROM'):
            # Exclude samples line
            continue
        elif not line:
            # Skip empty lines
            continue
        else:
            new_header_lines.append(line)
    new_header = pysam.VariantHeader()
    for s in header.samples:
        new_header.add_sample(s)
    for line in new_header_lines:
        new_header.add_line(line)
    return new_header


def valid_coord(chrom, coord, ref_contig_length_dict):
    if chrom not in ref_contig_length_dict:
        raise ValueError(f"Contig not found in reference: {chrom}")
    result = 1 <= coord <= ref_contig_length_dict[chrom]
    if not result:
        logging.warning(f"Invalid coordinate {chrom}:{coord}")
    return result


def get_startswith_match_in_list(substr, str_list):
    matches = [t.replace(substr, '') for t in str_list if t.startswith(substr)]
    if len(matches) == 0:
        return None
    elif len(matches) == 1:
        if matches[0] == ".":
            return None
        else:
            return matches[0]
    else:
        raise ValueError(f"More than one match for {substr} in {str_list}")


def valid_record_coordinates(record_tokens, ref_contig_length_dict):
    chrom = record_tokens[0]
    pos = int(record_tokens[1])
    info_tokens = record_tokens[7].split(';')
    end = get_startswith_match_in_list("END=", info_tokens)
    if end is not None:
        end = int(end)
    chr2 = get_startswith_match_in_list("CHR2=", info_tokens)
    end2 = get_startswith_match_in_list("END2=", info_tokens)
    if end2 is not None:
        end2 = int(end2)
    if not valid_coord(chrom=chrom,
                       coord=pos,
                       ref_contig_length_dict=ref_contig_length_dict):
        return False
    elif not valid_coord(chrom=chrom,
                         coord=end,
                         ref_contig_length_dict=ref_contig_length_dict):
        return False
    elif end2 is not None:
        # For END2 check: use CHR2 if available, otherwise default to CHROM
        if chr2 is None:
            chr2 = chrom
        if not valid_coord(chrom=chr2,
                           coord=end2,
                           ref_contig_length_dict=ref_contig_length_dict):
            return False
    return True


def process_record(record, drop_invalid_coords, ref_contig_length_dict):
    # Fix multi-allelic CNV alts (legacy)
    if record.alts[0].startswith('<CN'):
        record.alts = ('<CNV>',)
        record.info['SVTYPE'] = 'CNV'
    # Remove MULTIALLELIC field (legacy)
    record.info.pop('MULTIALLELIC')
    # Since pysam messes with some of the formatting (i.e. END limitation) we parse the string and replace
    record_tokens = str(record).strip().split('\t')
    if drop_invalid_coords:
        if not valid_record_coordinates(record_tokens=record_tokens, ref_contig_length_dict=ref_contig_length_dict):
            logging.warning(f"Dropping record {record_tokens[2]}")
            return None
    format_keys = record_tokens[8].split(':')
    gt_index = format_keys.index('GT')
    new_record_tokens = record_tokens[:9]
    for format in record_tokens[9:]:
        format_tokens = format.split(':')
        # Reset multiallelic CNV genotypes to ./.
        if record.info['SVTYPE'] == 'CNV':
            gt_tokens = format_tokens[gt_index].split('/')
            format_tokens[gt_index] = "/".join(["." for _ in gt_tokens])
        new_record_tokens.append(':'.join(format_tokens))
    return '\t'.join(new_record_tokens)


def parse_reference_fai(path):
    if not path:
        raise ValueError("Reference index required if using --drop-invalid-coords")
    with open(path) as f:
        result = dict()
        for line in f:
            tokens = line.strip().split('\t')
            chrom = tokens[0]
            length = int(tokens[1])
            result[chrom] = length
        return result


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--vcf", type=str, required=True,
                        help="Path to input VCF from GATK-SV GenotypeBatch module")
    parser.add_argument("--drop-invalid-coords", action='store_true',
                        help="Drop records with invalid coordinates, i.e. POS/END/END2 greater than chromosome length")
    parser.add_argument("--reference-fai", type=str, required=False,
                        help="Reference fasta index (.fai), required only if using --drop-invalid-coords")
    args = parser.parse_args()
    logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

    if args.drop_invalid_coords:
        ref_contig_length_dict = parse_reference_fai(args.reference_fai)
    else:
        ref_contig_length_dict = None

    with pysam.VariantFile(args.vcf) as vcf:
        new_header = get_new_header(vcf.header)
        sys.stdout.write(str(new_header))
        for record in vcf:
            new_record = process_record(record=record,
                                        drop_invalid_coords=args.drop_invalid_coords,
                                        ref_contig_length_dict=ref_contig_length_dict)
            if new_record is not None:
                sys.stdout.write(new_record + "\n")


if __name__ == '__main__':
    main()
