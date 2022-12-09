#!/usr/bin/env python

"""
Annotates variants with no-call rate and SL statistics
"""

import argparse
import sys
import pysam
from typing import Any, List, Text, Set, Dict, Optional


_gt_status_map = dict()
_gt_non_ref_or_no_call_map = dict()
_gt_non_ref_map = dict()


def _is_non_ref_or_no_call(gt):
    s = _gt_non_ref_or_no_call_map.get(gt, None)
    if s is None:
        s = any([a is not None and a > 0 for a in gt]) or all([a is None for a in gt])
        _gt_non_ref_or_no_call_map[gt] = s
    return s


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


def annotate_ncr(vcf, fout, ploidy_dict):
    n_samples = float(len(fout.header.samples))
    if n_samples == 0:
        raise ValueError("Cannot annotate sites-only vcf")
    for record in vcf:
        if record.info['SVTYPE'] == 'CNV':
            continue
        gt_items = record.samples.items()
        n_no_call = sum([_is_no_call(gt['GT']) for s, gt in gt_items if ploidy_dict[s][record.chrom] > 0])
        sl = [float(gt['SL']) for s, gt in gt_items
              if _is_non_ref_or_no_call(gt['GT']) and ploidy_dict[s][record.chrom] > 0]
        record.info['NCN'] = n_no_call
        record.info['NCR'] = n_no_call / n_samples
        record.info['SL_MEAN'] = sum(sl) / len(sl) if len(sl) > 0 else None
        record.info['SL_MAX'] = max(sl) if len(sl) > 0 else None
        record.info['SL_MIN'] = min(sl) if len(sl) > 0 else None

        sl = [float(gt['SL']) for _, gt in gt_items if _is_non_ref(gt['GT'])]
        record.info['SL_MEAN_UNFILTERED'] = sum(sl) / len(sl) if len(sl) > 0 else None
        record.info['SL_MAX_UNFILTERED'] = max(sl) if len(sl) > 0 else None
        record.info['SL_MIN_UNFILTERED'] = min(sl) if len(sl) > 0 else None
        fout.write(record)


def _parse_ploidy_table(path: Text) -> Dict[Text, Dict[Text, int]]:
    """
    Parses tsv of sample ploidy values.

    Parameters
    ----------
    path: Text
        table path

    Returns
    -------
    header: Dict[Text, Dict[Text, int]]
        map of sample to contig to ploidy, i.e. Dict[sample][contig] = ploidy
    """
    ploidy_dict = dict()
    with open(path, 'r') as f:
        header = f.readline().strip().split('\t')
        for line in f:
            tokens = line.strip().split('\t')
            sample = tokens[0]
            ploidy_dict[sample] = {header[i]: int(tokens[i]) for i in range(1, len(header))}
    return ploidy_dict


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--vcf', type=str, help='Input vcf (defaults to stdin).')
    parser.add_argument('--out', type=str, help='Output file (defaults to stdout)')
    parser.add_argument('--ploidy-table', type=str, required=True, help='Ploidy table tsv')

    args = parser.parse_args()

    if args.vcf is None:
        vcf = pysam.VariantFile(sys.stdin)
    else:
        vcf = pysam.VariantFile(args.vcf)

    header = vcf.header
    header.add_line('##INFO=<ID=NCN,Number=1,Type=Integer,Description="Number of no-call genotypes">')
    header.add_line('##INFO=<ID=NCR,Number=1,Type=Float,Description="Rate of no-call genotypes">')
    header.add_line('##INFO=<ID=SL_MEAN,Number=1,Type=Float,Description="Mean SL of no-call and non-ref genotypes">')
    header.add_line('##INFO=<ID=SL_MIN,Number=1,Type=Float,Description="Min SL of no-call and non-ref genotypes">')
    header.add_line('##INFO=<ID=SL_MAX,Number=1,Type=Float,Description="Max SL of no-call and non-ref genotypes">')
    header.add_line('##INFO=<ID=SL_MEAN_UNFILTERED,Number=1,Type=Float,Description="Mean SL of non-ref genotypes">')
    header.add_line('##INFO=<ID=SL_MIN_UNFILTERED,Number=1,Type=Float,Description="Min SL of non-ref genotypes">')
    header.add_line('##INFO=<ID=SL_MAX_UNFILTERED,Number=1,Type=Float,Description="Max SL of non-ref genotypes">')
    if args.out is None:
        fout = pysam.VariantFile(sys.stdout, 'w', header=header)
    else:
        fout = pysam.VariantFile(args.out, 'w', header=header)

    ploidy_dict = _parse_ploidy_table(args.ploidy_table)
    annotate_ncr(vcf, fout, ploidy_dict)
    fout.close()


if __name__ == '__main__':
    main()
