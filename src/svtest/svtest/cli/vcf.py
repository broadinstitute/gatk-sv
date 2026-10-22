#!/usr/bin/env python

"""
Collect vcf metrics. Writes metrics to stdout.

Metrics:
  <prefix>_vcf_<type>_count   : variants of a type in the test set
  <prefix>_vcf_<type>_tp      : variants with at least one matching variant in the baseline set (if baseline-vcf provided)
  <prefix>_vcf_<type>_fp      : variants with no matching variants in the baseline set (if baseline-vcf provided)
  <prefix>_vcf_<type>_fn      : variants in baseline set that had no matching variants in the test set (if baseline-vcf provided)
  <prefix>_vcf_<type>_size_X_Y     : variants with size >= X and < Y
  <prefix>_vcf_<type>_size_gte_X   : variants with size >= X
  <prefix>_vcf_<type>_vargq_X_Y    : variants with varGQ >= X and < Y (if varGQ present)
  <prefix>_vcf_<type>_vargq_gte_X  : variants with varGQ >= X (if varGQ present)
  <prefix>_vcf_<type>_af_X_Y     : variants with allele frequency >= X and < Y (if EV present)
  <prefix>_vcf_<type>_af_gte_X   : variants with allele frequency >= X (if EV present)
  <prefix>_vcf_<type>_ac_1       : singleton variants (if EV present)
  <prefix>_vcf_<type>_evidence_<evidence> : variants supported by evidence type (if EVIDENCE present)

"""

import argparse
import sys
import svtest.utils.IntervalUtils as iu
import svtest.utils.TestUtils as tu
import svtest.utils.VCFUtils as vu
import svtest.utils.IOUtils as iou
from pysam import VariantFile
import pandas as pd

VCF_METRIC_STR = "_vcf_"

# Size bins
SIZES = [500, 5000, 100000]

# varGQ bins
VARGQ_BINS = [2, 200, 400, 600, 800, 999]

# AF bins
AF_BINS = [0.01, 0.1, 0.5]

# Valid evidence types
EVIDENCE_TYPES = ["RD", "BAF", "PE", "SR"]

# Accepted "passing" filters
PASSING_FILTERS = ["PASS", "BOTHSIDES_SUPPORT",
                   "MULTIALLELIC"]

INVALID_CHR2_STR = "invalid_chr2"
INVALID_END_STR = "invalid_end"

BED_FILE_HEADER_CHAR = "#"
BED_FILE_CHROM_COL = "chrom"
BED_FILE_START_COL = "start"
BED_FILE_END_COL = "end"
BED_FILE_SVTYPE_COL = "svtype"


def main(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        prog='svtest vcf',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('test_vcf', type=str)
    parser.add_argument('contig_list', type=str)
    parser.add_argument('sample_list', type=str)
    parser.add_argument('types', type=str,
                        help='Comma-delimited list of variant types (case-sensitive)')
    parser.add_argument('metric_prefix', type=str)
    parser.add_argument('--baseline-vcf', type=str,
                        help='Baseline vcf to provide evaluation metrics against')
    parser.add_argument('--baseline-bed', type=str,
                        help='Baseline bed file to provide evaluation metrics against. Must have header beginning with "' +
                             BED_FILE_HEADER_CHAR + '" and the following columns: "' +
                        '", "'.join([BED_FILE_CHROM_COL, BED_FILE_START_COL, BED_FILE_END_COL, BED_FILE_SVTYPE_COL]) + '"')
    parser.add_argument('--min-reciprocal-overlap', type=float, default=0.5,
                        help='Minimum reciprocal overlap for validation metrics [0.5]')
    parser.add_argument('--padding', type=int, default=50,
                        help='Interval padding for validation metrics [50]')
    parser.add_argument('--max-warnings', type=int, default=50,
                        help='Maximum number of records to print warnings for [50]')
    parser.add_argument('--fp-file', type=str, default=None,
                        help='Write false positives to file')
    parser.add_argument('--fn-file', type=str, default=None,
                        help='Write false negatives to file')
    parser.add_argument('--fp-pass-file', type=str, default=None,
                        help='Write PASS false positives to file')
    parser.add_argument('--fn-pass-file', type=str, default=None,
                        help='Write PASS false negatives to file')

    # Print help if no arguments specified
    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)
    if (args.baseline_vcf is None and args.baseline_bed is None) and (args.fp_file is not None or args.fn_file is not None):
        raise ValueError(
            "FP and FN files cannot be generated if --baseline-vcf and --baseline-bed aren't specified")
    if args.baseline_vcf is not None and args.baseline_bed is not None:
        raise ValueError(
            "Cannot specify both --baseline-vcf and --baseline-bed")
    types_list = args.types.split(',')

    contigs = iou.read_contig_list(args.contig_list)
    samples = iou.read_samples_list(args.sample_list)
    metrics, fp_intervals, fn_intervals, fp_intervals_pass, fn_intervals_pass = get_metrics(args.test_vcf,
                                                                                            args.baseline_vcf,
                                                                                            args.baseline_bed,
                                                                                            contigs,
                                                                                            types_list,
                                                                                            args.min_reciprocal_overlap,
                                                                                            args.padding,
                                                                                            samples,
                                                                                            args.metric_prefix,
                                                                                            args.max_warnings)

    # Write metrics
    write_metrics(metrics)
    if args.fp_file is not None and fp_intervals is not None:
        write_intervals(args.fp_file, fp_intervals)
    if args.fn_file is not None and fn_intervals is not None:
        write_intervals(args.fn_file, fn_intervals)
    if args.fp_pass_file is not None and fp_intervals_pass is not None:
        write_intervals(args.fp_pass_file, fp_intervals_pass)
    if args.fn_pass_file is not None and fn_intervals_pass is not None:
        write_intervals(args.fn_pass_file, fn_intervals_pass)


def write_metrics(metrics):
    for key in metrics:
        sys.stdout.write("%s\t%s\n" % (key, str(metrics[key])))


def write_intervals(path, intervals):
    with open(path, 'w') as f:
        for type in intervals:
            for contig in intervals[type]:
                for interval in intervals[type][contig]:
                    line = "\t".join(
                        [contig, str(interval[0]), str(interval[1]), type]) + "\n"
                    f.write(line)


def get_metrics(ftest, fbase_vcf, fbase_bed, contigs, variant_types, min_ro, padding, samples, metric_prefix, max_warnings):

    test_vcf = VariantFile(ftest)

    check_header(test_vcf, samples)

    genotyped = check_if_genotyped(test_vcf)
    has_vargq = check_if_vargq(test_vcf)
    collect_evidence = check_if_evidence(test_vcf)
    test_records = list(test_vcf.fetch())

    unfiltered_variant_type_counts = get_count_by_type(
        test_records, variant_types)

    pass_filter_set = set(PASSING_FILTERS)
    pass_records = [r for r in test_records if (
        "PASS" in r.filter or len(set(r.filter) - pass_filter_set) == 0)]

    error_counts = count_errors(test_records, contigs, max_warnings)

    variant_type_counts = get_count_by_type(pass_records, variant_types)
    size_counts = get_distributions_by_type(
        pass_records, variant_types, "SVLEN", SIZES, exclude_types=['BND'])

    metrics = add_error_count_metrics({}, error_counts, metric_prefix)

    if fbase_vcf is not None:
        base_vcf = VariantFile(fbase_vcf)
        if genotyped != check_if_genotyped(base_vcf):
            raise ValueError(
                "One of the vcfs seems to be genotyped but the other does not")
        if has_vargq != check_if_vargq(base_vcf):
            raise ValueError(
                "One of the vcfs has the varGQ field but the other does not")
        if collect_evidence != check_if_evidence(base_vcf):
            raise ValueError(
                "One of the vcfs has the EVIDENCE field but the other does not")
        base_records = list(base_vcf.fetch())
        test_tree = iu.create_trees_from_records(
            test_records, variant_types, contigs)
        test_pass_tree = iu.create_trees_from_records(
            pass_records, variant_types, contigs)
        base_tree = iu.create_trees_from_records(
            base_records, variant_types, contigs)
        base_pass_records = [r for r in base_records if (
            "PASS" in r.filter or len(set(r.filter) - pass_filter_set) == 0)]
        base_pass_tree = iu.create_trees_from_records(
            base_pass_records, variant_types, contigs)
    elif fbase_bed is not None:
        base_records = parse_bed_file(fbase_bed)
        test_tree = iu.create_trees_from_records(
            test_records, variant_types, contigs)
        test_pass_tree = iu.create_trees_from_records(
            pass_records, variant_types, contigs)
        base_tree = iu.create_trees_from_bed_records(
            base_records, variant_types, contigs)
        base_pass_tree = None
    else:
        test_tree = None
        test_pass_tree = None
        base_tree = None
        base_pass_tree = None

    if base_tree is not None:
        metrics, fp_intervals, fn_intervals = add_evaluation_metrics(
            metrics, test_tree, base_tree, variant_types, min_ro, padding, metric_prefix)
    else:
        fp_intervals = None
        fn_intervals = None

    if base_pass_tree is not None:
        metrics, fp_intervals_pass, fn_intervals_pass = add_evaluation_metrics(
            metrics, test_pass_tree, base_pass_tree, variant_types, min_ro, padding,
            metric_prefix, metric_suffix="_pass")
    else:
        fp_intervals_pass = None
        fn_intervals_pass = None

    if genotyped:
        allele_frequencies, num_singletons = get_allele_frequency_counts(
            pass_records, test_vcf.header, variant_types)
    if has_vargq:
        vargq_counts = get_distributions_by_type(
            pass_records, variant_types, "varGQ", VARGQ_BINS)
    if collect_evidence:
        evidence_counts = collect_evidence_fields(pass_records, variant_types)

    for type in variant_types:
        metrics[metric_prefix + VCF_METRIC_STR + type +
                "_count"] = unfiltered_variant_type_counts[type]
        metrics[metric_prefix + VCF_METRIC_STR + type +
                "_pass_count"] = variant_type_counts[type]
        if type != 'BND':
            metrics = add_binned_metrics(
                size_counts, SIZES, type, metrics, metric_prefix, "pass_size")
        if genotyped:
            metrics = add_binned_metrics(
                allele_frequencies, AF_BINS, type, metrics, metric_prefix, "pass_af")
            if type in num_singletons:
                metrics[metric_prefix + VCF_METRIC_STR +
                        type + "_pass_ac_1"] = num_singletons[type]
        if has_vargq:
            metrics = add_binned_metrics(
                vargq_counts, VARGQ_BINS, type, metrics, metric_prefix, "pass_vargq")
        if collect_evidence:
            metrics = add_metrics_from_dict(
                evidence_counts, type, metrics, metric_prefix, "pass_evidence")

    return metrics, fp_intervals, fn_intervals, fp_intervals_pass, fn_intervals_pass


def add_error_count_metrics(metrics, error_counts, metric_prefix):
    for key in error_counts:
        metrics[metric_prefix + VCF_METRIC_STR + key] = error_counts[key]
    return metrics


def add_evaluation_metrics(metrics, test_tree, base_tree, variant_types, min_ro, padding,
                           metric_prefix, metric_suffix=""):
    tp_test = {}
    fp_test = {}
    fp_intervals_test = {}
    tp_base = {}
    fp_base = {}
    fp_intervals_base = {}
    for type in variant_types:
        tp_test_, fp_test_, fp_intervals_test_ = iu.evaluate_tree(
            test_tree[type], base_tree[type], min_ro, padding=padding)
        tp_test[type] = tp_test_
        fp_test[type] = fp_test_
        fp_intervals_test[type] = fp_intervals_test_
        tp_base_, fp_base_, fp_intervals_base_ = iu.evaluate_tree(
            base_tree[type], test_tree[type], min_ro, padding=padding)
        tp_base[type] = tp_base_
        fp_base[type] = fp_base_
        fp_intervals_base[type] = fp_intervals_base_

    tp_test_by_type = sum_counts_over_contigs(tp_test)
    fp_test_by_type = sum_counts_over_contigs(fp_test)
    fp_base_by_type = sum_counts_over_contigs(fp_base)

    for type in variant_types:
        metrics[metric_prefix + VCF_METRIC_STR + type +
                "_tp" + metric_suffix] = tp_test_by_type[type]
        metrics[metric_prefix + VCF_METRIC_STR + type +
                "_fp" + metric_suffix] = fp_test_by_type[type]
        metrics[metric_prefix + VCF_METRIC_STR + type +
                "_fn" + metric_suffix] = fp_base_by_type[type]
    return metrics, fp_intervals_test, fp_intervals_base


def get_allele_frequency_counts(records, header, variant_types):
    num_samples = float(len(header.samples))
    allele_freq = {}
    num_singletons = {}
    types_set = set(variant_types)
    # Don't calculate MCNV AF since non-ref alleles cannot be determined without chromosome ploidy
    af_types = types_set - set(["CNV"])
    for type in af_types:
        allele_freq[type] = []
        num_singletons[type] = 0
    for record in records:
        type = vu.get_sv_type(record, types_set)
        if type not in af_types:
            continue
        af = 0
        for sample in record.samples.values():
            for val in sample["GT"]:
                if val is not None and val > 0:
                    af += 1
        if af == 1:
            num_singletons[type] += 1
        allele_freq[type].append(af / num_samples)
    allele_freq_counts = {}
    num_bins = len(AF_BINS)
    for type in af_types:
        allele_freq_counts[type] = [0] * (num_bins + 1)
        for val in allele_freq[type]:
            idx = get_distribution_index(val, AF_BINS, num_bins)
            allele_freq_counts[type][idx] += 1
    return allele_freq_counts, num_singletons


def count_errors(records, contigs, max_warnings):
    contigs_set = set(contigs)
    error_counts = {INVALID_CHR2_STR: 0, INVALID_END_STR: 0}
    print_warnings = True
    for record in records:
        warnings_maxed = sum(error_counts.values()) > max_warnings
        if warnings_maxed and print_warnings:
            sys.stderr.write(
                "Max of %d warnings have been given, the rest will be suppressed.\n" % max_warnings)
            print_warnings = False
        error_counts = check_record(
            error_counts, record, contigs_set, print_warnings)
    return error_counts


def check_record(error_counts, record, contigs_set, warn):
    if "CHR2" in record.info and record.info["CHR2"] not in contigs_set:
        if warn:
            sys.stderr.write("Invalid CHR2 value: %s\n" % record.info["CHR2"])
        error_counts[INVALID_CHR2_STR] += 1
    if not valid_end_field(record):
        if warn:
            sys.stderr.write(
                "Position was not less than END in record %s\n" % record.id)
        error_counts[INVALID_END_STR] += 1
    return error_counts


def valid_end_field(record):
    return record.pos <= record.stop


def check_if_genotyped(vcf):
    return check_header_format_field(vcf, "EV")


def check_if_vargq(vcf):
    return check_header_info_field(vcf, "varGQ")


def check_if_evidence(vcf):
    return check_header_info_field(vcf, "EVIDENCE")


def check_header_info_field(vcf, name):
    if name in vcf.header.info:
        return True
    return False


def check_header_format_field(vcf, name):
    if name in vcf.header.formats:
        return True
    return False


def check_header(vcf, expected_samples):
    vcf_samples = vcf.header.samples
    tu.test_sets_equal(vcf_samples, expected_samples, item_str="sample",
                       name_a="VCF header", name_b="samples list")


def collect_evidence_fields(records, variant_types):
    evidence_counts = {}
    for variant_type in variant_types:
        evidence_counts[variant_type] = {}
        for evidence_type in EVIDENCE_TYPES:
            evidence_counts[variant_type][evidence_type] = 0
    variant_types_set = set(variant_types)
    evidence_types_set = set(EVIDENCE_TYPES)
    for record in records:
        variant_type = vu.get_sv_type(record, variant_types_set)
        evidence_types = vu.get_evidence_types(record, evidence_types_set)
        for evidence_type in evidence_types:
            evidence_counts[variant_type][evidence_type] += 1
    return evidence_counts


def add_binned_metrics(metric_counts, bins, type, metrics, prefix, name):
    if type not in metric_counts:
        return metrics
    for i in range(len(bins)):
        if i == 0:
            lb = 0
        else:
            lb = bins[i - 1]
        key = prefix + VCF_METRIC_STR + type + "_" + \
            name + "_" + str(lb) + "_" + str(bins[i])
        metrics[key] = metric_counts[type][i]
    key = prefix + VCF_METRIC_STR + type + "_" + name + "_gte_" + str(bins[-1])
    metrics[key] = metric_counts[type][-1]
    return metrics


def add_metrics_from_dict(dict, type, metrics, prefix, name):
    for key in dict[type]:
        metric_name = prefix + VCF_METRIC_STR + type + "_" + name + "_" + key
        metrics[metric_name] = dict[type][key]
    return metrics


def get_count_by_type(records, variant_types):
    counts = {}
    types_set = set(variant_types)
    for type in variant_types:
        counts[type] = 0
    for record in records:
        type = vu.get_sv_type(record, types_set)
        counts[type] += 1
    return counts


def get_distributions_by_type(records, variant_types, field, bins, exclude_types=[]):
    num_bins = len(bins)
    counts = {}
    types_set = set(variant_types)
    for type in variant_types:
        counts[type] = [0] * (num_bins + 1)
    for record in records:
        type = vu.get_sv_type(record, types_set)
        if type not in exclude_types:
            print(record.id, file=sys.stderr)
            print(field, file=sys.stderr)
            val = vu.get_info_field(record, field, singularize=True)
            idx = get_distribution_index(val, bins, num_bins)
            counts[type][idx] += 1
    return counts


def get_distribution_index(val, bins, num_bins):
    for i in range(num_bins):
        if val < bins[i]:
            return i
    return len(bins)


def sum_counts_over_contigs(x):
    result = {}
    for type in x:
        result[type] = 0
        for contig in x[type]:
            result[type] += x[type][contig]
    return result


def parse_bed_file(path):
    df = pd.read_csv(path, delimiter='\t')
    return df[[BED_FILE_HEADER_CHAR + BED_FILE_CHROM_COL, BED_FILE_START_COL, BED_FILE_END_COL, BED_FILE_SVTYPE_COL]].values


if __name__ == '__main__':
    main()
