#!/usr/bin/env python

import argparse


DELIMITER = "\t"


def merge(genotypes_filename, gq_filename, merged_filename):
    with open(genotypes_filename, "r") as genotypes, open(gq_filename, "r") as gq, open(merged_filename, "w") as merged:

        # Integrity check: do the files have same columns?
        genotypes_header = genotypes.readline().rstrip().split(DELIMITER)
        gq_header = gq.readline().rstrip().split(DELIMITER)
        if not genotypes_header == gq_header:
            raise ValueError("The files do not have same number/order of columns")

        n_cols = len(gq_header)
        for genotypes_line, gq_line in zip(genotypes, gq):
            x = genotypes_line.rstrip().split(DELIMITER)
            y = gq_line.rstrip().split(DELIMITER)

            # Check if lines in the files are in the correct order.
            if not x[0:4] == y[0:4]:
                raise ValueError(f"The lines in the files are not in the same order; "
                                 f"expected the following lines to match.\n{x[0:4]}\n{y[0:4]}")

            h = DELIMITER.join(x[0:4])
            for i in range(4, n_cols):
                merged.write(DELIMITER.join([h, gq_header[i], x[i], y[i]]) + "\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('genotypes')
    parser.add_argument('GQ')
    parser.add_argument('fout')
    args = parser.parse_args()

    merge(args.genotypes, args.GQ, args.fout)
