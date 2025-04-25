#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

"""
Standardize a VCF of SV calls.

Each record corresponds to a single SV breakpoint and will have the following
INFO fields, with specified constraints:
  SVTYPE:  SV type [DEL,DUP,INV,BND]
  CHR2:    Secondary chromosome [Must be lexicographically greater than CHROM]
  END:     SV end position (or position on CHR2 in translocations)
  STRANDS: Breakpoint strandedness [++,+-,-+,--]
  SVLEN:   SV length (-1 if translocation)
  ALGORITHMS:  Source algorithm
"""

import argparse
import sys
import pkg_resources
from svtk.standardize import VCFStandardizer
from pysam import VariantFile


def main(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        prog='svtk standardize',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf', help='Raw VCF.')
    parser.add_argument('fout', help='Standardized VCF.')
    parser.add_argument('source', help='Source algorithm. '
                        '[delly,dragen,lumpy,manta,wham,melt,scramble]')
    parser.add_argument('-p', '--prefix', help='If provided, variant names '
                        'will be overwritten with this prefix.')
    parser.add_argument('--include-reference-sites', action='store_true',
                        default=False, help='Include records where all '
                        'samples are called 0/0 or ./.')
    parser.add_argument('--standardizer', help='Path to python file with '
                        'custom standardizer definition. (Not yet supported.)')
    parser.add_argument('--contigs', type=argparse.FileType('r'),
                        help='Reference fasta index (.fai). If provided, '
                        'contigs in index will be used in VCF header. '
                        'Otherwise all GRCh37 contigs will be used in header. '
                        'Variants on contigs not in provided list will be '
                        'removed.')
    parser.add_argument('--min-size', type=int, default=50,
                        help='Minimum SV size to report [50].')
    parser.add_argument('--call-null-sites', action='store_true',
                        default=False,
                        help='Call sites with null genotypes (./.). Generally '
                        'useful when an algorithm has been run on a single '
                        'sample and has only reported variant sites.')
    parser.add_argument('--sample-names', type=str, default=None,
                        help='Comma-delimited list of sample names to use in '
                             'header [use existing].')

    # Print help if no arguments specified
    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    vcf = VariantFile(args.vcf)
    # Parse new sample names if provided
    if args.sample_names:
        sample_names_list = args.sample_names.split(',')
    else:
        sample_names_list = vcf.header.samples

    # Add contigs to header if provided
    if args.contigs:
        template = pkg_resources.resource_filename(
            'svtk', 'data/no_contigs_template.vcf')
        # pysam can no longer open up a VCF header with FORMAT but no samples, so copy template to temporary file
        # and add samples, then open and return header
        header = VCFStandardizer.get_header_from_template(template, sample_names_list)
        contig_line = '##contig=<ID={contig},length={length}>'
        for line in args.contigs:
            contig, length = line.split()[:2]
            header.add_line(contig_line.format(contig=contig, length=length))
    # Use GRCh37 by default
    else:
        template = pkg_resources.resource_filename(
            'svtk', 'data/GRCh37_template.vcf')
        # pysam can no longer open up a VCF header with FORMAT but no samples, so copy template to temporary file
        # and add samples, then open and return header
        header = VCFStandardizer.get_header_from_template(template, sample_names_list)

    # Tag source in header
    meta = '##FORMAT=<ID={0},Number=1,Type=Integer,Description="Called by {1}">'
    meta = meta.format(args.source, args.source.capitalize())
    header.add_line(meta)
    header.add_line('##source={0}'.format(args.source))

    fout = VariantFile(args.fout, mode='w', header=header)

    standardizer = VCFStandardizer.create(
        args.source, vcf, fout, sample_names_list,
        args.prefix, args.min_size, args.include_reference_sites,
        args.call_null_sites)

    for record in standardizer.standardize_vcf():
        fout.write(record)

    fout.close()
    vcf.close()


if __name__ == '__main__':
    main(sys.argv[1:])
