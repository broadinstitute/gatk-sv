#!/bin/env python

import argparse
import sys
import pysam
from typing import Any, List, Text, Set, Dict, Optional


_gt_status_map = dict()
_gt_non_ref_or_no_call_map = dict()
_gt_non_ref_map = dict()
_cnv_types = ('DEL', 'DUP')


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


def annotate_ncr(record, n_samples, ploidy_dict):
    gt_items = record.samples.items()
    n_no_call = sum([_is_no_call(gt['GT']) for s, gt in gt_items if ploidy_dict[s][record.chrom] > 0])
    sl = [float(gt['SL']) for s, gt in gt_items
          if _is_non_ref_or_no_call(gt['GT']) and ploidy_dict[s][record.chrom] > 0]
    record.info['NCN'] = n_no_call
    record.info['NCR'] = n_no_call / n_samples
    record.info['SL_MEAN'] = sum(sl) / len(sl) if len(sl) > 0 else None
    record.info['SL_MAX'] = max(sl) if len(sl) > 0 else None
    record.info['SL_MIN'] = min(sl) if len(sl) > 0 else None


def apply_filter(record, sl_threshold):
    record.info['MINSL'] = sl_threshold
    if sl_threshold is None:
        return
    for gt in record.samples.values():
        if _is_non_ref(gt['GT']) and gt['SL'] >= sl_threshold:
            gt['GT'] = tuple(None for _ in gt['GT'])


def get_threshold(record, sl_thresholds, med_size, large_size):
    svtype = record.info['SVTYPE']
    if svtype in _cnv_types:
        svlen = record.info['SVLEN']
        if svlen < med_size:
            size_index = 0
        elif svlen < large_size:
            size_index = 1
        else:
            size_index = 2
        return sl_thresholds[svtype][size_index]
    else:
        return sl_thresholds[svtype][0]


def process(vcf, fout, ploidy_dict, thresholds, args):
    n_samples = float(len(fout.header.samples))
    if n_samples == 0:
        raise ValueError("This is a sites-only vcf")
    k = 0
    for record in vcf:
        sl_threshold = get_threshold(record, thresholds, args.medium_size, args.large_size)
        apply_filter(record, sl_threshold)
        annotate_ncr(record, n_samples, ploidy_dict)
        fout.write(record)


def _create_threshold_dict(args):
    return {
        'DEL': [args.small_del_threshold, args.medium_del_threshold, args.large_del_threshold],
        'DUP': [args.small_dup_threshold, args.medium_dup_threshold, args.large_dup_threshold],
        'INS': [args.ins_threshold],
        'INV': [args.inv_threshold],
        'BND': [args.bnd_threshold],
        'CPX': [args.cpx_threshold],
        'CTX': [args.ctx_threshold],
        'CNV': [None]
    }


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


def _parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Apply class-specific genotype filtering using SL scores and annotate NCR",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--vcf', type=str, help='Input vcf (defaults to stdin)')
    parser.add_argument('--out', type=str, help='Output file (defaults to stdout)')
    parser.add_argument('--ploidy-table', type=str, required=True, help='Ploidy table tsv')
    parser.add_argument("--medium-size", type=float, default=500,
                        help="Min size for medium DEL/DUP")
    parser.add_argument("--large-size", type=float, default=1000,
                        help="Min size for large DEL/DUP")
    parser.add_argument("--small-del-threshold", type=float,
                        help="Threshold SL for small DELs")
    parser.add_argument("--medium-del-threshold", type=float,
                        help="Threshold SL for medium DELs")
    parser.add_argument("--large-del-threshold", type=float,
                        help="Threshold SL for large DELs")
    parser.add_argument("--small-dup-threshold", type=float,
                        help="Threshold SL for small DUPs")
    parser.add_argument("--medium-dup-threshold", type=float,
                        help="Threshold SL for medium DUPs")
    parser.add_argument("--large-dup-threshold", type=float,
                        help="Threshold SL for large DUPs")
    parser.add_argument("--ins-threshold", type=float, default=1,
                        help="Threshold SL for INS")
    parser.add_argument("--inv-threshold", type=float, default=1,
                        help="Threshold SL for INV")
    parser.add_argument("--bnd-threshold", type=float,
                        help="Threshold SL for BND")
    parser.add_argument("--cpx-threshold", type=float,
                        help="Threshold SL for CPX")
    parser.add_argument("--ctx-threshold", type=float,
                        help="Threshold SL for CTX")
    if len(argv) <= 1:
        parser.parse_args(["--help"])
        sys.exit(0)
    parsed_arguments = parser.parse_args(argv[1:])
    return parsed_arguments


def main(argv: Optional[List[Text]] = None):
    if argv is None:
        argv = sys.argv
    args = _parse_arguments(argv)

    if args.vcf is None:
        vcf = pysam.VariantFile(sys.stdin)
    else:
        vcf = pysam.VariantFile(args.vcf)

    header = vcf.header
    header.add_line('##INFO=<ID=NCN,Number=1,Type=Integer,Description="Number of no-call genotypes">')
    header.add_line('##INFO=<ID=NCR,Number=1,Type=Float,Description="Rate of no-call genotypes">')
    header.add_line('##INFO=<ID=SL_MEAN,Number=1,Type=Float,Description="Mean SL of filtered and non-ref genotypes">')
    header.add_line('##INFO=<ID=SL_MIN,Number=1,Type=Float,Description="Min SL of filtered and non-ref genotypes">')
    header.add_line('##INFO=<ID=SL_MAX,Number=1,Type=Float,Description="Max SL of filtered and non-ref genotypes">')
    if args.out is None:
        fout = pysam.VariantFile(sys.stdout, 'w', header=header)
    else:
        fout = pysam.VariantFile(args.out, 'w', header=header)

    ploidy_dict = _parse_ploidy_table(args.ploidy_table)
    thresholds = _create_threshold_dict(args)
    process(vcf, fout, ploidy_dict, thresholds, args)
    fout.close()


if __name__ == "__main__":
    main()
