#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""

"""

import argparse
import sys
import pandas as pd
import pysam
import svtk.utils as svu


def final_filter(vcf, bed, chrom=None):
    bed = bed.set_index('name')
    new_variant = ((bed.operation == 'add-mosaic-variant') |
                   (bed.operation == 'add-variant'))

    filter_names = bed.loc[~new_variant].index.values

    for record in vcf:
        if record.id in filter_names:
            op = bed.loc[record.id, 'operation']
            samples = bed.loc[record.id, 'samples'].split(',')

            if op == 'add-samples':
                for sample in samples:
                    record.samples[sample]['GT'] = (0, 1)
            elif op == 'flag-mosaic':
                for sample in samples:
                    record.samples[sample]['GT'] = (0, 1)
                record.info['MOSAIC'] = samples
            elif op == 'remove-variant':
                continue

        yield record

    for new_name in bed.loc[new_variant].index.values:
        # Stick to variants on same chromosome
        data = bed.loc[new_name]
        if chrom is not None and data.chrom != chrom:
            continue

        new_record = record.copy()
        for sample in new_record.samples.keys():
            svu.set_null(new_record, sample)

        samples = data['samples'].split(',')
        for sample in samples:
            new_record.samples[sample]['GT'] = (0, 1)

        new_record.id = new_name
        new_record.chrom = data.chrom
        new_record.pos = data.start
        new_record.alts = ('<{0}>'.format(data.svtype), )
        new_record.info['CHR2'] = data.chrom
        new_record.stop = data.end

        new_record.info['SVTYPE'] = data.svtype
        if data.svtype == 'DEL':
            new_record.info['STRANDS'] = '+-'
        else:
            new_record.info['STRANDS'] = '-+'
        new_record.info['SVLEN'] = int(data.end - data.start)
        new_record.info['SOURCES'] = data.sources.split(',')

        new_record.info['MEMBERS'] = (new_name, )

        op = bed.loc[new_name, 'operation']
        if op == 'add-mosaic-variant':
            if len(samples) == 1:
                new_record.info['MOSAIC'] = samples
            else:
                parents = [s for s in samples if s.endswith(
                    'fa') or s.endswith('mo')]
                new_record.info['MOSAIC'] = parents

        yield new_record


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf')
    parser.add_argument('filter_bed')
    parser.add_argument('fout')
    parser.add_argument('--chrom')
    args = parser.parse_args()

    bed = pd.read_table(args.filter_bed)

    vcf = pysam.VariantFile(args.vcf)

    header = vcf.header
    header.add_line(
        '##INFO=<ID=MOSAIC,Number=.,Type=String,Description="Samples predicted to harbor somatic or germline mosaicism">')

    if args.fout in '- stdout'.split():
        fout = pysam.VariantFile(sys.stdout, 'w', header=header)
    else:
        fout = pysam.VariantFile(args.fout, 'w', header=header)

    for record in final_filter(vcf, bed, args.chrom):
        # Apply size filter
        if record.info['SVLEN'] >= 50 or record.info['SVLEN'] == -1:
            fout.write(record)


if __name__ == '__main__':
    main()
