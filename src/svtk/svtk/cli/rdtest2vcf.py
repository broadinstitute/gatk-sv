#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

"""
Convert an RdTest-formatted bed to the standard VCF format.

Each record corresponds to a single CNV interval and will have the following
INFO fields, with specified constraints:
  SVTYPE:  SV type [DEL,DUP]
  CHR2:    Secondary chromosome; set to same as CHROM
  END:     SV end position
  STRANDS: Breakpoint strandedness [DEL:+-,DUP:-+]
  SVLEN:   SV length
  ALGORITHMS:  Tagged with "depth"
"""

import argparse
import sys
import os
import pkg_resources
from collections import namedtuple
import pysam
from svtk.standardize import VCFStandardizer


def RdtestParser(bed):
    CNV = namedtuple('CNV', 'chrom start end name samples svtype'.split())

    for line in bed:
        if line.startswith('#'):
            continue
        data = line.strip().split()

        chrom = data[0]
        start = int(data[1])
        end = int(data[2])
        name = data[3]
        samples = data[4].split(',')
        svtype = data[5].upper()

        yield CNV(chrom, start, end, name, samples, svtype)


def rdtest2vcf(bed, vcf):
    for cnv in RdtestParser(bed):
        record = vcf.new_record()
        record.chrom = cnv.chrom

        if cnv.start < 1:
            record.pos = 1
        else:
            record.pos = cnv.start
        record.id = cnv.name

        record.ref = 'N'
        record.alts = ('<{0}>'.format(cnv.svtype), )

        record.filter.add('PASS')

        # Add required INFO fields
        record.info['SVTYPE'] = cnv.svtype
        record.info['CHR2'] = cnv.chrom
        record.stop = cnv.end
        record.info['SVLEN'] = cnv.end - cnv.start
        record.info['ALGORITHMS'] = ['depth']
        if cnv.svtype == 'DEL':
            record.info['STRANDS'] = '+-'
        elif cnv.svtype == 'DUP':
            record.info['STRANDS'] = '-+'

        # Seed with ref genotypes
        for sample in vcf.header.samples:
            record.samples[sample]['GT'] = (0, 0)
            record.samples[sample]['depth'] = 0

        # Call any samples with variant as heterozygous
        called = 0
        for sample in cnv.samples:
            if sample in vcf.header.samples:
                called += 1
                record.samples[sample]['GT'] = (0, 1)
                record.samples[sample]['depth'] = 1

        if called > 0:
            vcf.write(record)


def main(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        prog='svtk rdtest2vcf',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('bed', type=argparse.FileType('r'),
                        help='RdTest-formatted bed file. '
                        '(chrom, start, end, name, samples, svtype)')
    parser.add_argument('samples', help='List of all samples present in '
                        'variant callset.')
    parser.add_argument('fout', help='Standardized VCF. Will be compressed '
                        'with bgzip and tabix indexed if filename ends with '
                        '.gz')
    parser.add_argument('--contigs', type=argparse.FileType('r'),
                        help='Reference fasta index (.fai). If provided, '
                        'contigs in index will be used in VCF header. '
                        'Otherwise all GRCh37 contigs will be used in header.')

    # Print help if no arguments specified
    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    # Get list of samples
    with open(args.samples) as slist:
        samples = sorted([s.strip() for s in slist.readlines()])

    # Add contigs to header if provided
    if args.contigs:
        template = pkg_resources.resource_filename(
            'svtk', 'data/no_contigs_template.vcf')
        # pysam can no longer open up a VCF header with FORMAT but no samples, so copy template to temporary file
        # and add samples, then open and return header
        header = VCFStandardizer.get_header_from_template(template, samples)
        contig_line = '##contig=<ID={contig},length={length}>'
        for line in args.contigs:
            contig, length = line.split()[:2]
            header.add_line(contig_line.format(**locals()))
    # Use GRCh37 by default
    else:
        template = pkg_resources.resource_filename(
            'svtk', 'data/GRCh37_template.vcf')
        # pysam can no longer open up a VCF header with FORMAT but no samples, so copy template to temporary file
        # and add samples, then open and return header
        header = VCFStandardizer.get_header_from_template(template, samples)
        header = template.header

    # Tag source in header
    meta = ('##FORMAT=<ID=depth,Number=1,Type=Integer,'
            'Description="Called by read-depth algorithms">')
    header.add_line(meta)
    header.add_line('##source=depth')

    if args.fout.endswith('.vcf.gz'):
        fname = os.path.splitext(args.fout)[0]
    elif args.fout.endswith('.vcf'):
        fname = args.fout
    else:
        msg = 'Invalid VCF filename; must end with .vcf or .vcf.gz: {0}'
        msg = msg.format(args.fout)
        raise ValueError(msg)

    fout = pysam.VariantFile(fname, mode='w', header=header)

    rdtest2vcf(args.bed, fout)
    fout.close()

    # TODO: do this with subprocess so we don't have to write to disk twice
    if args.fout.endswith('.gz'):
        pysam.tabix_compress(fname, args.fout)
        pysam.tabix_index(args.fout, preset='vcf')
        os.remove(fname)


if __name__ == '__main__':
    main(sys.argv[1:])
