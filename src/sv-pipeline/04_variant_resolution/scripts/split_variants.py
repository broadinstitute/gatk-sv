#!/bin/python
import argparse
import logging


def process_bed_file(input_bed, n_per_split, bca=True, digits=9):
    SVTYPE_FIELD = 5
    END_FIELD = 2
    START_FIELD = 1

    # Conditions for each category of variants
    condition_prefixes = {
        'gt5kb': {'condition': lambda line: (line[SVTYPE_FIELD] == 'DEL' or line[SVTYPE_FIELD] == 'DUP') and
                                            (int(line[END_FIELD]) - int(line[START_FIELD]) >= 5000)},
        'lt5kb': {'condition': lambda line: (line[SVTYPE_FIELD] == 'DEL' or line[SVTYPE_FIELD] == 'DUP') and
                                            (int(line[END_FIELD]) - int(line[START_FIELD]) < 5000)},
        'bca': {'condition': lambda line: bca and line[SVTYPE_FIELD] not in ['DEL', 'DUP'] and not line[SVTYPE_FIELD].startswith('INS')},
        'ins': {'condition': lambda line: bca and line[SVTYPE_FIELD].startswith('INS')}
    }

    # Create trackers for the current file information
    current_lines = {prefix: [] for prefix in condition_prefixes.keys()}
    current_counts = {prefix: 0 for prefix in condition_prefixes.keys()}
    current_suffixes = {prefix: 0 for prefix in condition_prefixes.keys()}

    with open(input_bed, 'r') as infile:
        for line in infile:
            line = line.strip('\n').split('\t')
            # This line swaps the last two columns so the sample names are in the fifth column and SV type in the last
            line[4], line[5] = line[5], line[4]
            for prefix, conditions in condition_prefixes.items():
                # If a line matches a condition add it to the appropriate category
                if conditions['condition'](line):
                    current_lines[prefix].append('\t'.join(line))
                    current_counts[prefix] += 1
                    # If a category has the specified number of records, create a new file and write the current records
                    if current_counts[prefix] == n_per_split:
                        output_file = get_file_name(prefix, current_suffixes[prefix], digits)
                        with open(output_file, 'w') as outfile:
                            outfile.write('\n'.join(current_lines[prefix]))
                        # Log the file name that was created
                        logging.info(f"File '{output_file}' written.")
                        # Update the tracking information
                        current_lines[prefix] = []
                        current_counts[prefix] = 0
                        current_suffixes[prefix] = current_suffixes[prefix] + 1
    # Handle the remaining records
    for prefix, lines in current_lines.items():
        if lines:
            output_file = get_file_name(prefix, current_suffixes[prefix], digits)
            with open(output_file, 'w') as outfile:
                outfile.write('\n'.join(lines))
            logging.info(f"File '{output_file}' written.")


def get_file_name(prefix, suffix, digits):
    if len(str(suffix)) > digits:
        raise ValueError('No more files can be generated with the current naming scheme. '
                         'Increase the digits parameter or the n parameter to proceed.')
    return f"{prefix}.{str(suffix).zfill(digits)}.bed"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bed", help="Path to input bed file", required=True)
    parser.add_argument("--n", help="number of variants per output file", required=True, type=int)
    parser.add_argument("--bca", default=False, help="Flag to set to True if the VCF contains BCAs",
                        action='store_true')
    parser.add_argument("--digits", "-d", default=9, type=int, help="Number of digits in filename suffix")
    parser.add_argument("--log-level", required=False, default="INFO", help="Specify level of logging information")
    args = parser.parse_args()

    log_level = args.log_level
    numeric_level = getattr(logging, log_level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % log_level)
    logging.basicConfig(level=numeric_level, format='%(levelname)s: %(message)s')
    process_bed_file(args.bed, args.n, args.bca, args.digits)


if __name__ == '__main__':
    main()
