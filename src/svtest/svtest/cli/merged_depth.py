#!/usr/bin/env python

"""
Collect merged depth file metrics. Writes SVTYPE counts to stdout.

Note --test-hits and --baseline-hits files can be generated with bedtools intersect:
    bedtools intersect -wa -u -f 0.5 -r -a ${test_bed} -b ${baseline_bed} | cut -f4 > overlap.test.list
    bedtools intersect -wa -u -f 0.5 -r -b ${test_bed} -a ${baseline_bed} | cut -f4 > overlap.base.list

Metrics:
 merged_depth_<TYPE>_count : total variant count
 merged_depth_<TYPE>_tp    : count of variants in test set that had at least one matching variant in the baseline set
                             (if baseline-bed specified)
 merged_depth_<TYPE>_fp    : count of variants in test set that had no matching variants in the baseline set
                             (if baseline-bed specified)
 merged_depth_<TYPE>_tp    : count of variants in baseline set that had no matching variants in the test set
                             (if baseline-bed specified)
 merged_depth_<TYPE>_X_Y     : variants with size >= X and < Y
 merged_depth_<TYPE>_gte_X   : variants with size >= X

"""

import gzip
import argparse
import sys
import svtest.utils.TestUtils as tu
import svtest.utils.IOUtils as iou

# For size distribution
SIZES = [1000, 2000, 3000, 4000, 5000, 10000, 100000]

KEY_PREFIX = "merged_depth_"


def main(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        prog='svtest merged-depth',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('test_bed', type=str)
    parser.add_argument('contig_list', type=str)
    parser.add_argument('type', type=str)
    parser.add_argument('--baseline-bed', type=str, default=None,
                        help="Baseline bed file to evaluate against")
    parser.add_argument('--test-hits', type=str,
                        help="List of test record ids that overlap baseline set (required if using --baseline-bed)")
    parser.add_argument('--baseline-hits', type=str,
                        help="List of baseline record ids that overlap test set (required if using --baseline-bed)")

    # Print help if no arguments specified
    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    if (bool(args.baseline_bed) ^ bool(args.test_hits)) or \
       (bool(args.baseline_bed) ^ bool(args.baseline_hits)) or \
       (bool(args.test_hits) ^ bool(args.baseline_hits)):
        raise ValueError(
            "Inconsistent arguments specified: --baseline-bed, --test-hits, and --baseline-hits must be specified together.")

    contigs = iou.read_contig_list(args.contig_list)

    # Read file
    with gzip.open(args.test_bed, mode='rb') as ftest:
        if args.baseline_bed is None:
            metrics = get_metrics(ftest, None, contigs,
                                  args.type, args.test_hits, args.baseline_hits)
        else:
            with gzip.open(args.baseline_bed, mode='rb') as fbase:
                metrics = get_metrics(
                    ftest, fbase, contigs, args.type, args.test_hits, args.baseline_hits)

    # Write metrics
    write_metrics(metrics)


def write_metrics(metrics):
    for key in metrics:
        sys.stdout.write("%s\t%s\n" % (key, str(metrics[key])))


def get_metrics(ftest, fbase, contigs, type, test_hits_path, base_hits_path):
    test_header, test_ids, size_counts, num_test_records = parse_bed(
        ftest, type, contigs)
    metrics = {
        KEY_PREFIX + type + "_count": num_test_records
    }
    metrics = get_size_metrics(metrics, size_counts, type)
    if fbase is not None:
        metrics = get_baseline_metrics(
            metrics, fbase, test_hits_path, base_hits_path, test_header, test_ids, num_test_records, type, contigs)
    return metrics


def get_size_metrics(metrics, size_counts, type):
    for i in range(len(SIZES)):
        if i == 0:
            lb = 0
        else:
            lb = SIZES[i - 1]
        key = KEY_PREFIX + type + "_" + str(lb) + "_" + str(SIZES[i])
        metrics[key] = size_counts[i]
    key = KEY_PREFIX + type + "_gte_" + str(SIZES[-1])
    metrics[key] = size_counts[-1]
    return metrics


def get_baseline_metrics(metrics, fbase, test_hits_path, base_hits_path, test_header, test_ids, num_test_records, type, contigs):
    base_header, base_ids, _, num_baseline_records = parse_bed(
        fbase, type, contigs)
    tu.test_sets_equal(test_header, base_header, item_str="header column",
                       name_a="test file header", name_b="baseline file header")
    if len(base_header) != len(test_header):
        raise ValueError('Files have different column header sizes')

    with open(test_hits_path, mode='r') as f:
        test_hits = f.read().splitlines()
    with open(base_hits_path, mode='r') as f:
        base_hits = f.read().splitlines()
    check_hit_ids(test_hits, test_ids, "test")
    check_hit_ids(base_hits, base_ids, "baseline")

    tp_test = len(test_hits)
    fp_test = num_test_records - tp_test
    fp_base = num_baseline_records - len(base_hits)
    metrics[KEY_PREFIX + type + "_tp"] = tp_test
    metrics[KEY_PREFIX + type + "_fp"] = fp_test
    metrics[KEY_PREFIX + type + "_fn"] = fp_base
    return metrics


def parse_bed(file, type, contigs):
    header = file.readline().decode().strip().split('\t')
    n_cols = len(header)
    num_records = 0
    num_sizes = len(SIZES)
    size_counts = [0] * (num_sizes + 1)
    contigs_set = set(contigs)
    ids = []
    for line in file:
        num_records += 1
        tokens = line.decode().strip().split('\t')
        check_record(tokens, n_cols, type, contigs_set)
        s = get_size_distribution_index(tokens, num_sizes)
        size_counts[s] += 1
        ids.append(tokens[3])
    return header, ids, size_counts, num_records


def check_hit_ids(hits, ids, name):
    unknown_ids = set(hits) - set(ids)
    if len(unknown_ids) > 0:
        raise ValueError("Unknown %s record ids: %s" % (name, str(ids)))


def get_size_distribution_index(tokens, num_sizes):
    start = int(tokens[1])
    end = int(tokens[2])
    interval_size = end - start
    for i in range(num_sizes):
        if interval_size < SIZES[i]:
            return i
    return len(SIZES)


def check_record(columns, n_cols, type, contigs):
    tu.test_iterable_size(columns, n_cols)
    tu.test_column_in_iterable(columns, 0, contigs)
    tu.test_is_int(columns, 1)
    tu.test_is_int(columns, 2)
    tu.test_column_equals(columns, 5, type)


def sum_counts_over_contigs(tree_counts):
    return sum([tree_counts[contig] for contig in tree_counts])


if __name__ == '__main__':
    main()
