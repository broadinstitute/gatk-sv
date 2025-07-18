#!/bin/env python

import argparse
import sys
import pysam
from typing import List, Text, Dict, Optional

_gt_no_call_map = dict()
_gt_hom_var_map = dict()
_gt_ref_map = dict()
_HIGH_NCR_FILTER = "HIGH_NCR"
_REFERENCE_ARTIFACT_FILTER = "REFERENCE_ARTIFACT"
_REFERENCE_ARTIFACT_THRESHOLD = 0.99
SVTYPES = ["DEL", "DUP", "INS", "CNV", "INV", "CPX", "CTX", "BND"]
GT_UPDATE = {(None,): (None, None), (0,): (0, 0), (1,): (0, 1)}


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


def _apply_filters(record, ploidy_dict, ncr_threshold, filter_reference_artifacts, remove_zero_carrier_sites):
    # remove previous NCR annotations and filters prior to recalculating
    if _HIGH_NCR_FILTER in record.filter:
        del record.filter[_HIGH_NCR_FILTER]
    if _REFERENCE_ARTIFACT_FILTER in record.filter:
        del record.filter[_REFERENCE_ARTIFACT_FILTER]
    for field in ['NCR', 'NCN']:
        if field in record.info:
            record.info.pop(field)

    if record.info['SVTYPE'] == 'CNV':
        return True  # skip checks and annotations for mCNVs due to lack of GT

    n_samples = 0
    n_no_call = 0
    n_hom_var = 0
    has_carrier = False
    for s, gt in record.samples.items():
        # set all GTs for biallelic SVs back to diploid
        if len(gt['GT']) == 1:
            orig_gt = gt['GT']
            gt['GT'] = GT_UPDATE[orig_gt]
        # Count every sample where ploidy > 0
        if ploidy_dict[s][record.chrom] == 0:
            continue
        n_samples += 1
        # Count no-calls and hom vars
        if _is_no_call(gt['GT']):
            n_no_call += 1
        elif _is_hom_var(gt['GT']):
            n_hom_var += 1
            has_carrier = True
        elif not _is_hom_ref(gt['GT']):
            has_carrier = True

    # hard filter sites with no carriers
    if remove_zero_carrier_sites and not has_carrier:
        return False

    # Annotate metrics and apply filters
    record.info['NCN'] = n_no_call
    record.info['NCR'] = n_no_call / n_samples if n_samples > 0 else None
    if ncr_threshold is not None and record.info['NCR'] is not None and record.info['NCR'] >= ncr_threshold:
        record.filter.add(_HIGH_NCR_FILTER)
        if 'PASS' in record.filter:
            del record.filter['PASS']

    if filter_reference_artifacts and n_hom_var / n_samples > _REFERENCE_ARTIFACT_THRESHOLD:
        record.filter.add(_REFERENCE_ARTIFACT_FILTER)
        if 'PASS' in record.filter:
            del record.filter['PASS']

    # if filter is . set to PASS
    if len(record.filter) == 0:
        record.filter.add('PASS')

    # Clean out AF metrics since they're no longer valid
    for field in ['AC', 'AF']:
        if field in record.info:
            del record.info[field]
    return True


def _increment_counter(record, counter, prev_chrom):
    if prev_chrom is not None and prev_chrom != record.chrom:
        counter = {svtype: 0 for svtype in SVTYPES}
    counter[record.info['SVTYPE']] += 1
    return counter, record.chrom


def _new_variant_id(record, counter, cohort_id, shard_str):
    svtype = record.info['SVTYPE']
    return f"{cohort_id}.{svtype}_{record.chrom}_{shard_str}{str(counter[svtype])}"


def process(vcf, fout, ploidy_dict, args):
    if args.cohort_id is not None:
        id_map = open(f"{args.cohort_id}.vid_map.tsv", 'w')
        counter = {svtype: 0 for svtype in SVTYPES}
        prev_chrom = None
        shard_str = ""
        if args.shard_index is not None:
            shard_str = f"shard{args.shard_index}_"
    n_samples = float(len(fout.header.samples))
    if n_samples == 0:
        raise ValueError("This is a sites-only vcf")
    for record in vcf:
        keep = _apply_filters(record=record,
                              ploidy_dict=ploidy_dict,
                              ncr_threshold=args.ncr_threshold,
                              filter_reference_artifacts=args.filter_reference_artifacts,
                              remove_zero_carrier_sites=args.remove_zero_carrier_sites)
        if keep:
            # re-set all variant IDs if cohort ID is provided to be part of variant IDs
            if args.cohort_id is not None:
                counter, prev_chrom = _increment_counter(record, counter, prev_chrom)
                new_id = _new_variant_id(record, counter, args.cohort_id, shard_str)
                id_map.write(f"{record.id}\t{new_id}\n")
                record.id = new_id
            fout.write(record)
    if args.cohort_id is not None:
        id_map.close()


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
        description="Apply class-specific genotype filtering using SL scores and annotates NCR",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--vcf', type=str, help='Input vcf (defaults to stdin)')
    parser.add_argument('--out', type=str, help='Output file (defaults to stdout)')
    parser.add_argument('--ploidy-table', type=str, required=True, help='Ploidy table tsv')
    parser.add_argument("--ncr-threshold", type=float,
                        help=f"If provided, adds {_HIGH_NCR_FILTER} filter to records with no-call rates equal to or "
                             f"exceeding this threshold")
    parser.add_argument("--filter-reference-artifacts", action='store_true', default=False,
                        help=f"If provided, adds {_REFERENCE_ARTIFACT_FILTER} filter to records that are hom alt in "
                             f">{_REFERENCE_ARTIFACT_THRESHOLD*100}% of samples>")
    parser.add_argument("--remove-zero-carrier-sites", action='store_true', default=False,
                        help="If provided, hard filters sites with zero carriers>")
    parser.add_argument("--cohort-id", type=str, required=False,
                        help="If provided, rename all variant IDs, using cohort-id as base")
    parser.add_argument("--shard-index", type=str, required=False,
                        help="Provide if sharding within chromosomes for unique variant IDs")
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
    if args.ncr_threshold is not None:
        header.add_line(f"##FILTER=<ID={_HIGH_NCR_FILTER},Description=\"Unacceptably high rate of no-call GTs\">")
    if args.filter_reference_artifacts:
        header.add_line(f"##FILTER=<ID={_REFERENCE_ARTIFACT_FILTER},""Description=\"Likely reference artifact sites"
                        " that are homozygous alternative in over 99% of the samples\">")
    if args.out is None:
        fout = pysam.VariantFile(sys.stdout, 'w', header=header)
    else:
        fout = pysam.VariantFile(args.out, 'w', header=header)

    ploidy_dict = _parse_ploidy_table(args.ploidy_table)
    process(vcf, fout, ploidy_dict, args)
    fout.close()


if __name__ == "__main__":
    main()
