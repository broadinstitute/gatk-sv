#!/bin/env python

import argparse
import sys
import pysam
import math
from numpy import median
from typing import List, Text, Dict, Optional

_gt_no_call_map = dict()
_gt_hom_var_map = dict()
_gt_ref_map = dict()
_gt_to_filter_status_map = dict()
_gt_to_filter_gt_map = dict()
_cnv_types = ('DEL', 'DUP')
_HIGH_NCR_FILTER = "HIGH_NCR"


def _is_no_call(gt):
    s = _gt_no_call_map.get(gt, None)
    if s is None:
        s = all([a is None for a in gt])
        _gt_no_call_map[gt] = s
    return s


def _is_hom_ref(gt):
    s = _gt_ref_map.get(gt, None)
    if s is None:
        s = all([a is not None and a == 0 for a in gt])
        _gt_ref_map[gt] = s
    return s


def _is_hom_var(gt):
    s = _gt_hom_var_map.get(gt, None)
    if s is None:
        s = all([a is not None and a > 0 for a in gt])
        _gt_hom_var_map[gt] = s
    return s


def _gt_to_filter_status(gt, fails_filter):
    if not fails_filter:
        return 'pass'
    s = _gt_to_filter_status_map.get(gt, None)
    if s is None:
        if _is_hom_ref(gt):
            s = 'homRefFail'
        elif _is_hom_var(gt):
            s = 'homVarFail'
        elif _is_no_call(gt):
            s = 'noCallFail'
        else:
            s = 'hetFail'
        _gt_to_filter_status_map[gt] = s
    return s


def _filter_gt(gt, allele):
    s = _gt_to_filter_gt_map.get(gt, None)
    if s is None:
        s = tuple(allele for _ in gt)
        _gt_to_filter_gt_map[gt] = s
    return s


def recalculate_gq(sl, scale_factor, is_hom_ref, upper, lower, shift):
    sign = -1 if is_hom_ref else 1
    sl = (sign * min(max(sl, lower), upper)) + shift
    return int(math.floor(-10 * math.log10(1 / (math.pow(scale_factor, sl) + 1))))


def _fails_filter(sl, gt_is_ref, cutoff):
    if cutoff is None:
        return False
    if gt_is_ref:
        sl = -sl
    return sl < cutoff


def recal_qual_score(record):
    """
    Recalibrate quality score for a single variant
    """
    quals = []
    for s in [s for s in record.samples]:
        gt = record.samples[s]['GT']
        if _is_hom_ref(gt) or _is_no_call(gt):
            continue
        elif _is_hom_var(gt):
            quals.append(99)
        else:
            quals.append(record.samples[s]['GQ'])

    if len(quals) > 0:
        return int(median(quals))


def _apply_filter(record, sl_threshold, ploidy_dict, apply_hom_ref, ncr_threshold,
                  keep_gq, gq_scale_factor, upper_sl_cap, lower_sl_cap, sl_shift, max_gq):
    record.info['MINSL'] = sl_threshold
    allele = 0 if apply_hom_ref else None
    sl_list = []
    n_samples = 0
    n_no_call = 0
    for s, gt in record.samples.items():
        # Always set this to avoid weird values from pysam
        gt['GT_FILTER'] = '.'
        # Skip over ploidy 0 or if SL is missing
        if ploidy_dict[s][record.chrom] == 0:
            continue
        # Count every sample where ploidy > 0 for NCR
        n_samples += 1
        sl = gt.get('SL', None)
        if sl is None:
            continue
        gt_is_ref = _is_hom_ref(gt['GT'])
        # SL metrics are over non-ref / no-call samples
        if not gt_is_ref:
            sl_list.append(sl)
        if not keep_gq:
            # Recalculate GQ values based on SL
            gt['GQ'] = min(recalculate_gq(sl=sl, scale_factor=gq_scale_factor, is_hom_ref=gt_is_ref,
                                          upper=upper_sl_cap, lower=lower_sl_cap, shift=sl_shift), max_gq)
        # Check and apply filter
        fails_filter = _fails_filter(sl, gt_is_ref, sl_threshold)
        gt['GT_FILTER'] = _gt_to_filter_status(gt['GT'], fails_filter)
        if fails_filter:
            gt['GT'] = _filter_gt(gt['GT'], allele)
        if record.info['SVTYPE'] != 'CNV' and _is_no_call(gt['GT']):
            n_no_call += 1
    # Annotate metrics
    record.info['NCN'] = n_no_call
    record.info['NCR'] = n_no_call / n_samples if n_samples > 0 else None
    if ncr_threshold is not None and record.info['NCR'] is not None and record.info['NCR'] >= ncr_threshold:
        record.filter.add(_HIGH_NCR_FILTER)
    if len(record.filter) == 0:
        record.filter.add('PASS')
    n_sl = len(sl_list)
    record.info['SL_MEAN'] = sum(sl_list) / n_sl if n_sl > 0 else None
    record.info['SL_MAX'] = max(sl_list) if n_sl > 0 else None
    record.qual = recal_qual_score(record)
    # Clean out AF metrics since they're no longer valid
    for field in ['AC', 'AF']:
        if field in record.info:
            del record.info[field]


def get_threshold(record, cutoff_rules):
    svtype = record.info['SVTYPE']
    svlen = record.info.get('SVLEN', None)
    
    for rule in cutoff_rules:
        if rule['sv_type'] != svtype:
            continue
        
        if rule['min_size'] is not None and (svlen is None or svlen < rule['min_size']):
            continue
        
        if rule['max_size'] is not None and (svlen is None or svlen >= rule['max_size']):
            continue
        
        return rule['sl_cutoff']
    
    return None


def process(vcf, fout, ploidy_dict, cutoff_rules, args):
    n_samples = float(len(fout.header.samples))
    if n_samples == 0:
        raise ValueError("This is a sites-only vcf")
    for record in vcf:
        sl_threshold = get_threshold(record=record, cutoff_rules=cutoff_rules)
        _apply_filter(record=record, sl_threshold=sl_threshold,
                      ploidy_dict=ploidy_dict, apply_hom_ref=args.apply_hom_ref,
                      ncr_threshold=args.ncr_threshold,
                      keep_gq=args.keep_gq, gq_scale_factor=args.gq_scale_factor,
                      upper_sl_cap=args.upper_sl_cap, lower_sl_cap=args.lower_sl_cap, sl_shift=args.sl_shift,
                      max_gq=args.max_gq)
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


def _parse_sl_cutoff_table(path: Text) -> List[Dict[Text, float]]:
    """
    Parses tsv of SL cutoff values.

    Parameters
    ----------
    path: Text
        table path with columns: sv_type, min_size, max_size, sl_cutoff

    Returns
    -------
    List[Dict[Text, float]]
        list of cutoff rules
    """
    cutoff_rules = []
    with open(path, 'r') as f:
        header = f.readline().strip().split('\t')
        for line in f:
            tokens = line.strip().split('\t')
            rule = {
                'sv_type': tokens[0],
                'min_size': None if tokens[1] == '-1' else float(tokens[1]),
                'max_size': None if tokens[2] == '-1' else float(tokens[2]),
                'sl_cutoff': float(tokens[3])
            }
            cutoff_rules.append(rule)
    return cutoff_rules


def _parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Apply class-specific genotype filtering using SL scores and annotates NCR",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--vcf', type=str, help='Input vcf (defaults to stdin)')
    parser.add_argument('--out', type=str, help='Output file (defaults to stdout)')
    parser.add_argument('--ploidy-table', type=str, required=True, help='Ploidy table tsv')
    parser.add_argument('--sl-cutoff-table', type=str, help='TSV table with columns: sv_type, min_size, max_size, sl_cutoff')

    parser.add_argument('--apply-hom-ref', action='store_true',
                        help='Set filtered genotypes to hom-ref (0/0) instead of no-call (./.)')
    parser.add_argument("--ncr-threshold", type=float,
                        help=f"If provided, adds {_HIGH_NCR_FILTER} filter to records with no-call rates equal to or "
                             f"exceeding this threshold")

    parser.add_argument('--keep-gq', action='store_true',
                        help='Do not replace GQ. Default: Recalculate GQ = '
                             'math.floor(-10 * math.log10(1 / (math.pow(gq_scale_factor, SL) + 1))),'
                             'where SL=-SL is applied when hom-ref')
    parser.add_argument("--gq-scale-factor", type=float, default=(0.52 / 0.48),
                        help="Scale factor if using --replace-gq")
    parser.add_argument("--upper-sl-cap", type=float, default=200,
                        help="SL=min(SL, upper_cap) is applied for GQ calculations")
    parser.add_argument("--lower-sl-cap", type=float, default=-200,
                        help="SL=max(SL, lower_cap) is applied for GQ calculations")
    parser.add_argument("--sl-shift", type=float, default=100,
                        help="SL=SL+shift is applied for GQ calculations")
    parser.add_argument("--max-gq", type=int, default=99, help="GQ cap")


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
    header.add_line('##INFO=<ID=SL_MAX,Number=1,Type=Float,Description="Max SL of filtered and non-ref genotypes">')
    header.add_line('##FORMAT=<ID=GT_FILTER,Number=1,Type=String,Description="Genotype filter status">')
    header.add_line(f"##FILTER=<ID={_HIGH_NCR_FILTER},Description=\"Unacceptably high rate of no-call GTs\">")
    if args.out is None:
        fout = pysam.VariantFile(sys.stdout, 'w', header=header)
    else:
        fout = pysam.VariantFile(args.out, 'w', header=header)

    ploidy_dict = _parse_ploidy_table(args.ploidy_table)
    cutoff_rules = []
    if args.sl_cutoff_table:
        cutoff_rules = _parse_sl_cutoff_table(args.sl_cutoff_table)
    process(vcf, fout, ploidy_dict, cutoff_rules, args)
    fout.close()


if __name__ == "__main__":
    main()
