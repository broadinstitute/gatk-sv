#!/usr/bin/env python
# -*- coding: utf-8 -*-


import argparse
import gzip
import pysam


def load_introns_for_contig(intron_file, target_contig):
    introns = []
    is_gzipped = intron_file.endswith('.gz')
    opener = gzip.open if is_gzipped else open

    with opener(intron_file, 'rt') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 4:
                continue
            chrom, start, end = fields[0], int(float(fields[2])), int(float(fields[3]))

            if chrom != target_contig:
                continue

            if start > end:
                start, end = end, start
            introns.append((start, end))

    introns.sort()
    return introns


def find_matching_introns(var_start, var_end, introns, max_distance=8):
    for intron_start, intron_end in introns:
        if abs(var_start - intron_start) + abs(var_end - intron_end) <= max_distance:
            return True

        if intron_start > var_end + max_distance:
            break

    return False


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('vcf', help='Input VCF file')
    parser.add_argument('intron_reference', help='Intron reference file (can be gzipped)')
    parser.add_argument('contig', help='Contig/chromosome to process')
    parser.add_argument('output', help='Output VCF file (will be bgzipped)')
    parser.add_argument(
        '--max-distance',
        type=int,
        default=8,
        help='Maximum combined breakpoint distance to consider a match (default: 8)'
    )

    args = parser.parse_args()

    introns = load_introns_for_contig(args.intron_reference, args.contig)

    with pysam.VariantFile(args.vcf, 'r') as fin:
        header = fin.header
        header.add_line("##FILTER=<ID=RETRO_DEL,Description=\"Deletion is close to an intron (breakpoint distance <= 8bp), so is likely a retrotransposon deletion\">")
        with pysam.VariantFile(args.output, 'w', header=header) as fo:
            for record in fin:
                if record.info.get('SVTYPE') == 'DEL':
                    var_start = record.pos
                    var_end = record.stop
                    if find_matching_introns(var_start, var_end, introns, args.max_distance):
                        record.filter.add('RETRO_DEL')
                fo.write(record)


if __name__ == '__main__':
    main()
