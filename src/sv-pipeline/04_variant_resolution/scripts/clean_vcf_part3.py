#!/bin/python

import argparse
from collections import defaultdict
from os import mkdir, path


def count_variants(infile):
    variant_counts = defaultdict(int)
    with open(infile, 'r') as IN:
        for line in IN:
            var_id = line.strip().split('\t')[0]
            variant_counts[var_id] += 1
    return dict(sorted(variant_counts.items(), key=lambda item: item[1], reverse=True))


def assign_shards(variant_counts, max_samples):
    shard_assignments = {}
    shard_number = 0
    sample_counter = 0
    first = True
    for variant in variant_counts.keys():
        if not first and (sample_counter + variant_counts[variant] > max_samples):
            shard_number += 1
            sample_counter = 0
        shard_assignments[variant] = shard_number
        sample_counter += variant_counts[variant]
        first = False
    return shard_number, shard_assignments


def create_shards(infile, shard_assignments, num_shards):
    if not path.isdir("./shards"):
        mkdir("./shards")
    with open(infile, 'r') as IN:
        for line in IN:
            var_id = line.strip().split('\t')[0]
            shard = shard_assignments[var_id]
            shard_file = f"shards/out.{shard}_{num_shards}.txt"
            with open(shard_file, 'a') as OUT:
                OUT.write(line)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("combined_file", help="rd_cn_revise file with variant ID, sample ID, and CN columns")
    parser.add_argument("-s", "--max-samples",
                        help="Maximum number of variant x sample entries in a shard (default = 7,000)",
                        default=7000, type=int)
    args = parser.parse_args()

    variant_counts = count_variants(args.combined_file)
    num_shards, shard_assignments = assign_shards(variant_counts, args.max_samples)
    create_shards(args.combined_file, shard_assignments, num_shards)


if __name__ == "__main__":
    main()
