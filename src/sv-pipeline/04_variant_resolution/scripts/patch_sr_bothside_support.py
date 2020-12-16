#!/bin/python

import sys
from collections import defaultdict


def count_vids(list_path):
    counts = defaultdict(lambda: 0)
    with open(list_path, 'r') as f_list:
        for path in f_list:
            with open(path.strip(), 'r') as f:
                for vid in f:
                    counts[vid.strip()] += 1
    return counts


def count_sr_pass(path, n):
    counts = defaultdict(lambda: 0)
    with open(path, 'r') as f:
        for line in f:
            tokens = line.strip().split('\t')
            n_support = round(float(tokens[0]) * n)
            vid = tokens[-1]
            counts[vid] = n_support
    return counts


NON_REF_VIDS_LIST = sys.argv[1]
BOTHSIDE_PASS_FILE = sys.argv[2]
NUM_BATCHES = int(sys.argv[3])

non_ref_counts = count_vids(NON_REF_VIDS_LIST)
bothside_pass_counts = count_sr_pass(BOTHSIDE_PASS_FILE, NUM_BATCHES)

with open(BOTHSIDE_PASS_FILE, 'r') as f:
    for line in f:
        tokens = line.strip().split('\t')
        vid = tokens[-1]
        bothside_pass_count = bothside_pass_counts[vid]
        if bothside_pass_count == 0:
            continue
        non_ref_count = non_ref_counts[vid]
        if non_ref_count == 0:
            continue
        fraction_support = min(1., bothside_pass_count / float(non_ref_count))
        sys.stdout.write("{}\t{}\n".format(fraction_support, "\t".join(tokens[1:])))
