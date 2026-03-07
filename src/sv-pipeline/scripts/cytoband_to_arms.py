import gzip
import argparse
import re
from collections import defaultdict
from dataclasses import dataclass


@dataclass
class ArmData:
    p_end: int = 0
    q_start: int = None
    q_end: int = 0


def natural_sort_key(s):
    return [int(text) if text.isdigit() else text.lower() for text in re.split("([0-9]+)", s)]


def main(input_file, output_file):
    arms = defaultdict(ArmData)

    with gzip.open(input_file, "rt") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.strip().split("\t")
            chrom, start, end, band = parts[0], int(parts[1]), int(parts[2]), parts[3]
            arm_type = band[0].lower()

            if arm_type == "p":
                if end > arms[chrom].p_end:
                    arms[chrom].p_end = end
            elif arm_type == "q":
                # Convert 0-based BED start to 1-based GATK start
                one_based_start = start + 1
                if arms[chrom].q_start is None:
                    arms[chrom].q_start = one_based_start
                arms[chrom].q_end = end

    with open(output_file, "w") as out:
        for chrom in sorted(arms.keys(), key=natural_sort_key):
            a = arms[chrom]

            # P-arm starts at 1, Q-arm starts at the first q-band + 1
            p_range = f"1-{a.p_end}" if a.p_end > 0 else "0"
            q_range = f"{a.q_start}-{a.q_end}" if a.q_start is not None else "0"

            # Output format: chrom \t p-interval \t q-interval
            out.write(f"{chrom}\t{p_range}\t{q_range}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()

    main(args.input, args.output)
