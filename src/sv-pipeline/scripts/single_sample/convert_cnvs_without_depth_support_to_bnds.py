#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
After genotyping, identifies CNV variants over a size threshold without depth support but with PESR evidence and changes
them to BNDs. This task is ordinarily accomplished by Module 3 but in some cases (ie the single sample pipeline) that is skipped.

Record-level decisions (drop / convert-to-BND) consider all "proband" samples given on the command line. In the
single-sample pipeline that is just the case; in trio de novo mode it is case + provided parents, so that a variant
supported by depth or PE/SR evidence in ANY proband is retained.
"""

import argparse
import csv
import pysam
import sys
from svtk.famfile import parse_famfile


def has_depth_support_autosome(record, sample):
    return record.samples[sample].get('RD_CN') is not None and \
        ((record.info['SVTYPE'] == 'DUP' and record.samples[sample].get('RD_CN') > 2) or
            (record.info['SVTYPE'] == 'DEL' and record.samples[sample].get('RD_CN') < 2))


def has_sr_or_pe_support(record, sample):
    return (record.samples[sample].get('PE_GT') is not None and record.samples[sample].get('PE_GT') > 0) \
        or (record.samples[sample].get('SR_GT') is not None and record.samples[sample].get('SR_GT') > 0)


def has_depth_support_allosome(record, sample, samples_with_same_sex):
    if record.samples[sample].get('RD_CN') is None:
        return False
    cns = [record.samples[s].get('RD_CN') for s in samples_with_same_sex]
    cns = [cn for cn in cns if cn is not None]
    if len(cns) == 0:
        return False
    cns.sort()
    median_cn = cns[int((len(cns) + 1) / 2)]
    return record.samples[sample].get('RD_CN') is not None and \
        ((record.info['SVTYPE'] == 'DUP' and record.samples[sample].get('RD_CN') > median_cn) or
         (record.info['SVTYPE'] == 'DEL' and record.samples[sample].get('RD_CN') < median_cn))


def read_contigs_list(contigs_list):
    contigs = []
    with open(contigs_list) as fin:
        reader = csv.reader(fin, delimiter='\t')
        for row in reader:
            contig = row[-1]
            if contig not in contigs:
                contigs.append(contig)

    return contigs


def convert_record_to_bnd(record):
    if record.info['SVTYPE'] == 'DEL':
        record.info['STRANDS'] = '+-'
    elif record.info['SVTYPE'] == 'DUP':
        record.info['STRANDS'] = '-+'
    record.info['SVTYPE'] = 'BND'
    record.info['CHR2'] = record.chrom
    record.info['END2'] = record.stop
    record.info['SVLEN'] = record.stop - record.start
    record.stop = record.pos
    record.alts = ['<BND>']


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf', help='Module 04 genotyped PESR vcf')
    parser.add_argument('allosome_contigs_file')
    parser.add_argument('famfile', type=argparse.FileType('r'))
    parser.add_argument(
        'proband_samples', nargs='+',
        help='Case sample and any provided parent sample names (space-separated).')
    parser.add_argument(
        'min_size', help='minumum size at which to apply conversions', type=int)
    parser.add_argument('-o', '--outfile',
                        help='Output file [default: stdout]')

    args = parser.parse_args()

    if args.vcf in '- stdin'.split():
        vcf = pysam.VariantFile(sys.stdin)
    else:
        vcf = pysam.VariantFile(args.vcf)
    header = vcf.header

    proband_samples = [s for s in args.proband_samples if s in vcf.header.samples]
    if not proband_samples:
        sys.exit("Error: none of the proband samples {} are present in the VCF "
                 "samples {}".format(args.proband_samples, vcf.header.samples))
    min_size = args.min_size

    if args.outfile is None:
        fout = pysam.VariantFile(sys.stdout, 'w', header=header)
    else:
        out = args.outfile
        fout = pysam.VariantFile(out, 'w', header=header)

    allosome_contigs = read_contigs_list(args.allosome_contigs_file)

    fam = parse_famfile(args.famfile)
    samples_by_sex = {'1': [s for s in fam.samples if fam.samples[s].sex == '1'],
                      '2': [s for s in fam.samples if fam.samples[s].sex == '2']}

    def has_depth_support(record, sample):
        """Depth support for one sample on one record (autosome or allosome)."""
        if record.contig not in allosome_contigs:
            return has_depth_support_autosome(record, sample)
        return has_depth_support_allosome(record, sample, samples_by_sex[fam.samples[sample].sex])

    for record in vcf:
        svtype = record.info['SVTYPE']
        if (svtype == 'DEL' or svtype == 'DUP') and record.info['SVLEN'] >= min_size:
            any_pesr_support = any(has_sr_or_pe_support(record, s) for s in proband_samples)
            all_depth_missing = all(
                record.samples[s].get('RD_CN') is None for s in proband_samples)
            if all_depth_missing and not any_pesr_support:
                sys.stderr.write("Record {} has a missing depth genotype and no PE/SR support in any proband; dropping\n".format(record.id))
                continue
            no_depth_support = not any(has_depth_support(record, s) for s in proband_samples)
            if no_depth_support and any_pesr_support:
                convert_record_to_bnd(record)
        fout.write(record)


if __name__ == '__main__':
    main()
