#!/usr/bin/env python
# -*- coding: utf-8 -*-


import argparse
import gzip
import pysam


def iter_introns(intron_file):
    """Yield (contig, start, end) for every reference row, with the coordinates ordered."""
    is_gzipped = intron_file.endswith('.gz')
    opener = gzip.open if is_gzipped else open

    with opener(intron_file, 'rt') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 4:
                continue
            contig, start, end = fields[0], int(float(fields[2])), int(float(fields[3]))
            yield (contig, end, start) if start > end else (contig, start, end)


def load_introns_for_contig(intron_file, contig):
    return sorted((start, end) for c, start, end in iter_introns(intron_file) if c == contig)


def load_introns_by_contig(intron_file):
    """Introns grouped by contig, each sorted - for a whole-genome VCF in one pass."""
    by_contig = {}
    for contig, start, end in iter_introns(intron_file):
        by_contig.setdefault(contig, []).append((start, end))

    for introns in by_contig.values():
        introns.sort()
    return by_contig


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
    parser.add_argument('output', help='Output VCF file (will be bgzipped)')
    parser.add_argument('--contig', default=None,
                        help='Contig to restrict intron matching to (required unless --all-contigs)')
    parser.add_argument('--all-contigs', action='store_true',
                        help='Match every record against the introns of its own contig, for callers '
                             'that hand over a whole-genome VCF instead of scattering per contig. '
                             'Equivalent to running once per contig, in one pass.')
    parser.add_argument(
        '--max-distance',
        type=int,
        default=8,
        help='Maximum combined breakpoint distance to consider a match (default: 8)'
    )

    args = parser.parse_args()
    if bool(args.contig) == args.all_contigs:
        parser.error('exactly one of --contig / --all-contigs is required')

    if args.all_contigs:
        introns_by_contig = load_introns_by_contig(args.intron_reference)
    else:
        # --contig keeps #953's behavior verbatim: the named contig's introns are matched against
        # every record in the file. Callers that hand over a multi-contig VCF want --all-contigs.
        shared_introns = load_introns_for_contig(args.intron_reference, args.contig)

    with pysam.VariantFile(args.vcf, 'r') as fin:
        header = fin.header
        header.add_line("##FILTER=<ID=RETRO_DEL,Description=\"Deletion is close to an intron (breakpoint distance <= 8bp), so is likely a retrotransposon deletion\">")
        with pysam.VariantFile(args.output, 'w', header=header) as fo:
            for record in fin:
                if record.info.get('SVTYPE') == 'DEL':
                    introns = introns_by_contig.get(record.chrom, []) if args.all_contigs \
                        else shared_introns
                    var_start = record.pos
                    var_end = record.stop
                    if find_matching_introns(var_start, var_end, introns, args.max_distance):
                        record.filter.add('RETRO_DEL')
                fo.write(record)


if __name__ == '__main__':
    main()
