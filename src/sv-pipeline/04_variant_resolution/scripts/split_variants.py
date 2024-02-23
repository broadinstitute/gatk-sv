#!/bin/python
import argparse
import logging


# Function to process the bed file by checking for conditions
def process_bed_file(input_bed, n_per_split, bca=True):
    SVTYPE_FIELD = 4
    END_FIELD = 2
    START_FIELD = 1

    # Dictionary to store the conditions to be checked with matching prefixes
    condition_prefixes = {
        'gt5kb': {
            'condition': lambda line: (line[SVTYPE_FIELD] == 'DEL' or line[SVTYPE_FIELD] == 'DUP') and (int(line[END_FIELD]) - int(line[START_FIELD]) >= 5000)},
        'lt5kb': {
            'condition': lambda line: (line[SVTYPE_FIELD] == 'DEL' or line[SVTYPE_FIELD] == 'DUP') and (int(line[END_FIELD]) - int(line[START_FIELD]) < 5000)},
        'bca': {'condition': lambda line: bca and (line[SVTYPE_FIELD] != 'DEL' and line[SVTYPE_FIELD] != 'DUP' and line[SVTYPE_FIELD] != 'INS')},
        'ins': {'condition': lambda line: bca and line[SVTYPE_FIELD] == 'INS'}
    }

    current_lines = {prefix: [] for prefix in condition_prefixes.keys()}
    current_counts = {prefix: 0 for prefix in condition_prefixes.keys()}
    current_suffixes = {prefix: 'a' for prefix in condition_prefixes.keys()}

    # Open the bed file and process
    with open(input_bed, 'r') as infile:
        for line in infile:
            # process bed file line by line
            line = line.strip('\n').split('\t')

            # Checks which condition and prefix the current line matches and appends it to the corresponding
            # array and increments the counter for that array
            for prefix, conditions in condition_prefixes.items():
                if conditions['condition'](line):
                    line[4], line[5] = line[5], line[4]
                    current_lines[prefix].append('\t'.join(line))
                    current_counts[prefix] += 1

                    # If the current array has the maximum allowed lines added to it create a new array
                    # with the preceding suffix and write the current array to a file
                    if current_counts[prefix] == n_per_split:
                        output_suffix = current_suffixes[prefix].rjust(6, 'a')
                        output_file = f"{prefix}.{output_suffix}.bed"
                        with open(output_file, 'w') as outfile:
                            outfile.write('\n'.join(current_lines[prefix]))

                        logging.info(f"File '{output_file}' written.")
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


# Function to generate the pattern for suffixes
def increment_suffix(suffix):
    # define the alphabet and ending
    alphabet = 'abcdefghijklmnopqrstuvwxyz'
    if suffix == 'z' * 6:
        raise ValueError('All possible files generated.')
    else:
        # if there are available suffixes increment to next available suffix
        index = alphabet.index(suffix[0])
        next_char = alphabet[(index + 1) % 26]
        return next_char + suffix[1:]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bed", help="Path to input bed file", required=True)
    parser.add_argument("--n", help="number of variants per file", required=True)
    parser.add_argument("--bca", default=False, help="Flag to set to True if the VCF contains BCAs",
                        action='store_true')
    parser.add_argument("--log-level", required=False, default="INFO", help="Specify level of logging information")
    args = parser.parse_args()

    # Set logging level from --log-level input
    log_level = args.log_level
    numeric_level = getattr(logging, log_level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % log_level)
    logging.basicConfig(level=numeric_level, format='%(levelname)s: %(message)s')
    process_bed_file(args.bed, args.n, args.bca)


if __name__ == '__main__':
    main()
