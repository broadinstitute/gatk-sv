#!/bin/python
import argparse
import logging


def process_bed_file(input_bed, n_per_split, bca=True):
    SVTYPE_FIELD = 5
    END_FIELD = 2
    START_FIELD = 1

    # Check the conditions to generate prefixes for the output files
    condition_prefixes = {
        'gt5kb': {'condition': lambda line: (line[SVTYPE_FIELD] == 'DEL' or line[SVTYPE_FIELD] == 'DUP') and (int(line[END_FIELD]) - int(line[START_FIELD]) >= 5000)},
        'lt5kb': {'condition': lambda line: (line[SVTYPE_FIELD] == 'DEL' or line[SVTYPE_FIELD] == 'DUP') and (int(line[END_FIELD]) - int(line[START_FIELD]) < 5000)},
        'bca': {'condition': lambda line: bca and line[SVTYPE_FIELD] not in ['DEL', 'DUP', 'INS']},
        'ins': {'condition': lambda line: bca and line[SVTYPE_FIELD] == 'INS'}
    }

    # Create trackers for the current file information
    current_lines = {prefix: [] for prefix in condition_prefixes.keys()}
    current_counts = {prefix: 0 for prefix in condition_prefixes.keys()}
    current_suffixes = {prefix: 'a' for prefix in condition_prefixes.keys()}

    with open(input_bed, 'r') as infile:
        for line in infile:
            line = line.strip('\n').split('\t')
            # This line swaps the last two columns so the sample names are in the fifth column and SV type in the last
            line[4], line[5] = line[5], line[4]
            for prefix, conditions in condition_prefixes.items():
                # If a line matches a condition add it to the appropriate file
                if conditions['condition'](line):
                    current_lines[prefix].append('\t'.join(line))
                    current_counts[prefix] += 1
                    # If a file has met the number of records per file create a new file with the next suffix and write
                    # the current line to that new file
                    if current_counts[prefix] == n_per_split:
                        output_suffix = current_suffixes[prefix].rjust(6, 'a')
                        output_file = f"{prefix}.{output_suffix}.bed"
                        with open(output_file, 'w') as outfile:
                            outfile.write('\n'.join(current_lines[prefix]))
                        # Keep track of which files have been written after reaching the max number of files
                        logging.info(f"File '{output_file}' written.")
                        # Update the tracking information
                        current_lines[prefix] = []
                        current_counts[prefix] = 0
                        current_suffixes[prefix] = increment_suffix(current_suffixes[prefix])
    # Handle the samples after files with the given number of lines per file have been written
    for prefix, lines in current_lines.items():
        if lines:
            output_suffix = current_suffixes[prefix].rjust(6, 'a')
            output_file = f"{prefix}.{output_suffix}.bed"
            with open(output_file, 'w') as outfile:
                outfile.write('\n'.join(lines))
            logging.info(f"File '{output_file}' written.")


# Create a function to appropriately add a suffix to each corresponding file
def increment_suffix(suffix):
    alphabet = 'abcdefghijklmnopqrstuvwxyz'
    if suffix == 'z' * 6:
        raise ValueError('All possible files generated.')
    else:
        index = alphabet.index(suffix[0])
        next_char = alphabet[(index + 1) % 26]
        return next_char + suffix[1:]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bed", help="Path to input bed file", required=True)
    parser.add_argument("--n", help="number of variants per file", required=True, type=int)
    parser.add_argument("--bca", default=False, help="Flag to set to True if the VCF contains BCAs",
                        action='store_true')
    parser.add_argument("--log-level", required=False, default="INFO", help="Specify level of logging information")
    args = parser.parse_args()

    log_level = args.log_level
    numeric_level = getattr(logging, log_level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % log_level)
    logging.basicConfig(level=numeric_level, format='%(levelname)s: %(message)s')
    process_bed_file(args.bed, args.n, args.bca)


if __name__ == '__main__':
    main()
