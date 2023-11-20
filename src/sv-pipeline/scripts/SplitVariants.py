#!/bin/python

import pandas as pd
import csv
import os
import argparse


def process_bed_file(input_bed, N, bca=True):
    condition_prefixes = {
        'gt5kb': {
            'condition': lambda line: (line[4] == 'DEL' or line[4] == 'DUP') and (int(line[2]) - int(line[1]) >= 5000)},
        'lt5kb': {
            'condition': lambda line: (line[4] == 'DEL' or line[4] == 'DUP') and (int(line[2]) - int(line[1]) < 5000)},
        'bca': {'condition': lambda line: bca and (line[4] != 'DEL' and line[4] != 'DUP' and line[4] != 'INS')},
        'ins': {'condition': lambda line: bca and line[4] == 'INS'}
    }

    current_lines = {prefix: [] for prefix in condition_prefixes.keys()}
    current_counts = {prefix: 0 for prefix in condition_prefixes.keys()}
    current_suffixes = {prefix: 'a' for prefix in condition_prefixes.keys()}

    with open(input_bed, 'r') as infile:
        for line in infile:
            line = line.strip().split('\t')

            for prefix, conditions in condition_prefixes.items():
                if conditions['condition'](line):
                    current_lines[prefix].append('\t'.join(line))
                    current_counts[prefix] += 1

                    if current_counts[prefix] == N:
                        output_suffix = current_suffixes[prefix].rjust(6, 'a')
                        output_file = f"{prefix}.{output_suffix}.bed"
                        with open(output_file, 'w') as outfile:
                            outfile.write('\n'.join(current_lines[prefix]))

                        print(f"File {output_file} written.")
                        current_lines[prefix] = []
                        current_counts[prefix] = 0
                        current_suffixes[prefix] = increment_suffix(current_suffixes[prefix])

    # Handle remaining lines after the loop
    for prefix, lines in current_lines.items():
        if lines:
            output_suffix = current_suffixes[prefix].rjust(6, 'a')
            output_file = f"{prefix}.{output_suffix}.bed"
            with open(output_file, 'w') as outfile:
                outfile.write('\n'.join(lines))

            print(f"File {output_file} written.")


def increment_suffix(suffix):
    alphabet = 'abcdefghijklmnopqrstuvwxyz'
    if suffix == 'z' * 6:
        return 'a' * 6
    else:
        index = alphabet.index(suffix[0])
        next_char = alphabet[(index + 1) % 26]
        return next_char + suffix[1:]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bed", help="Path to input bed file")
    parser.add_argument("--n", help="number of variants per file")
    parser.add_argument("--bca", default="FALSE", help="")
    args = parser.parse_args()
    process_bed_file(args.bed, args.n, args.bca)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()
