#!/usr/bin/env python

"""
Annotates variants with no-call rate
"""

import argparse
import sys
import pysam


_gt_status_map = dict()
_gt_non_ref_map = dict()


def _is_non_ref(gt):
    s = _gt_non_ref_map.get(gt, None)
    if s is None:
        s = any([a is not None and a > 0 for a in gt])
        _gt_non_ref_map[gt] = s
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
        n_non_ref = sum([_is_non_ref(gt['GT']) for gt in gt_list])
        record.info['N_NO_CALL'] = n_no_call
        record.info['N_NON_REF'] = n_non_ref
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
    header.add_line('##INFO=<ID=N_NO_CALL,Number=1,Type=Integer,Description="Number of no-call genotypes">')
    header.add_line('##INFO=<ID=N_NON_REF,Number=1,Type=Integer,Description="Number of non-ref genotypes">')
    if args.out is None:
        fout = pysam.VariantFile(sys.stdout, 'w', header=header)
    else:
        fout = pysam.VariantFile(args.out, 'w', header=header)

    annotate_ncr(vcf, fout)
    fout.close()


if __name__ == '__main__':
    main()
