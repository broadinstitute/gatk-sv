#!/usr/bin/env python

"""
Annotates variants with no-call rate
"""

import argparse
import sys
import pysam


_gt_status_map = dict()
_gt_sum_map = dict()


def _cache_gt_sum(gt):
    s = _gt_sum_map.get(gt, None)
    if s is None:
        s = sum([1 for a in gt if a is not None and a > 0])
        _gt_sum_map[gt] = s
    return s


def _is_no_call(gt):
    if gt not in _gt_status_map:
        _gt_status_map[gt] = gt is None or all(e is None for e in gt)
    return _gt_status_map[gt]


def annotate_ncr(vcf, fout):
    n_samples = float(len(fout.header.samples))
    for record in vcf:
        if record.info['SVTYPE'] == 'CNV':
            continue
        gt_list = record.samples.values()
        n_no_call = sum([_is_no_call(gt['GT']) for gt in gt_list])
        an = sum([len(gt.get('GT', list())) for gt in gt_list])
        ac = sum([_cache_gt_sum(gt['GT']) for gt in gt_list])
        af = ac / float(an) if an > 0 else None
        record.info['NCR'] = n_no_call / n_samples
        record.info['AN'] = an
        record.info['AC'] = ac
        record.info['AF'] = af
        fout.write(record)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--vcf', help='Input vcf (defaults to stdin).')
    parser.add_argument('--out', help='Output file (defaults to stdout)')

    args = parser.parse_args()

    if args.vcf is None:
        vcf = pysam.VariantFile(sys.stdin)
    else:
        vcf = pysam.VariantFile(args.vcf)

    header = vcf.header
    header.add_line('##INFO=<ID=NCR,Number=1,Type=Float,Description="Fraction of no-call genotypes">')
    header.add_line('##INFO=<ID=AC,Number=1,Type=Integer,Description="Allele count">')
    header.add_line('##INFO=<ID=AN,Number=1,Type=Integer,Description="Allele number">')
    header.add_line('##INFO=<ID=AF,Number=1,Type=Float,Description="Allele frequency">')
    if args.out is None:
        fout = pysam.VariantFile(sys.stdout, 'w', header=header)
    else:
        fout = pysam.VariantFile(args.out, 'w', header=header)

    annotate_ncr(vcf, fout)
    fout.close()


if __name__ == '__main__':
    main()
