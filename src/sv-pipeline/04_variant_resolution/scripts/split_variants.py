#!/bin/python
import argparse
import logging


def process_bed_file(input_bed, n_per_split, bca=True):
    SVTYPE_FIELD=4
    END_POS=2
    START_POS=1

    condition_prefixes = {
        'gt5kb': {
            'condition': lambda line: (line[SVTYPE_FIELD] == 'DEL' or line[SVTYPE_FIELD] == 'DUP') and (int(line[END_POS]) - int(line[START_POS]) >= 5000)},
        'lt5kb': {
            'condition': lambda line: (line[SVTYPE_FIELD] == 'DEL' or line[SVTYPE_FIELD] == 'DUP') and (int(line[END_POS]) - int(line[START_POS]) < 5000)},
        'bca': {'condition': lambda line: bca and (line[SVTYPE_FIELD] != 'DEL' and line[SVTYPE_FIELD] != 'DUP' and line[SVTYPE_FIELD] != 'INS')},
        'ins': {'condition': lambda line: bca and line[SVTYPE_FIELD] == 'INS'}
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

                    if current_counts[prefix] == n_per_split:
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

            logging.info(f"File '{output_file}' written.")


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
        "--bed", help="Path to input bed file", required=True)
    parser.add_argument("--n", help="number of variants per file",required=True)
    parser.add_argument("--bca", default=False, help="If there are ", action='store_true')
    args = parser.parse_args()
    process_bed_file(args.bed, args.n, args.bca)

if __name__ == '__main__':
    main()
