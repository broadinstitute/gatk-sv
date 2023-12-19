#!/bin/python
import argparse
import logging


# Function to process the bed file by checking for conditions
def process_bed_file(input_bed, n_per_split, bca=True):
    svtype_field = 4
    end_pos = 2
    start_pos = 1

    # Dictionary to store the conditions to be checked with matching prefixes
    condition_prefixes = {
        'gt5kb': {
            'condition': lambda curr_1: (curr_1[svtype_field] == 'DEL' or curr_1[svtype_field] == 'DUP') and (int(curr_1[end_pos]) - int(curr_1[start_pos]) >= 5000)},
        'lt5kb': {
            'condition': lambda curr_2: (curr_2[svtype_field] == 'DEL' or curr_2[svtype_field] == 'DUP') and (int(curr_2[end_pos]) - int(curr_2[start_pos]) < 5000)},
        'bca': {'condition': lambda curr_3: bca and (
                curr_3[svtype_field] != 'DEL' and curr_3[svtype_field] != 'DUP' and curr_3[svtype_field] != 'INS')},
        'ins': {'condition': lambda curr_4: bca and curr_4[svtype_field] == 'INS'}
    }

    current_lines = {prefix: [] for prefix in condition_prefixes.keys()}
    current_counts = {prefix: 0 for prefix in condition_prefixes.keys()}
    current_suffixes = {prefix: 'a' for prefix in condition_prefixes.keys()}

    # Open the bed file and process
    with open(input_bed, 'r') as infile:
        for line in infile:
            # process bed file line by line
            line = line.strip().split('\t')

            # Checks which condition and prefix the current line matches and appends it to the corresponding
            # array and increments the counter for that array
            for prefix, conditions in condition_prefixes.items():
                if conditions['condition'](line):
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
        # if there are available suffixes, increment with appropriate number
        # of padded zeroes
        index = alphabet.index(suffix[0])
        next_char = alphabet[(index + 1) % 26]
        return next_char + suffix[1:]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bed", help="Path to input bed file", required=True)
    parser.add_argument("--n", help="number of variants per file", required=True)
    parser.add_argument("--bca", default=False, help="If there are ", action='store_true')
    args = parser.parse_args()
    process_bed_file(args.bed, args.n, args.bca)


if __name__ == '__main__':
    main()
