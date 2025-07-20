#!/usr/bin/env python

"""
Preprocess GENCODE basic GTF to extract canonical protein-coding transcripts for functional consequence annotation.
"""

import argparse
import gzip


CHROM_FIELD = 0
ELEMENT_FIELD = 2
ATTRIBUTES_FIELD = 8
TRANSCRIPT_TYPES = {"protein_coding", "nonsense_mediated_decay"}
CANONICAL = {"MANE_Plus_Clinical", "MANE_Select", "Ensembl_canonical"}


# Flexibly open .gz or uncompressed file to read
def _open(filename):
    if filename.endswith(".gz"):
        return gzip.open(filename, 'rt')
    else:
        return open(filename, 'r')


# Extract transcript type and canonical status
def parse_attributes(field):
    # format: key1 "value1"; key2 "value2";
    # keys may be repeated so cannot convert directly to dictionary
    attributes_list = [tuple(x.replace('"', '').split(' ')) for x in field.rstrip(";").split("; ")]
    protein = False
    canonical = False
    for key, val in attributes_list:
        if key == "tag" and val in CANONICAL:
            canonical = True
        elif key == "transcript_type" and val in TRANSCRIPT_TYPES:
            protein = True
    return protein, canonical


def process(gtf, outfile):
    with _open(gtf) as inp, open(outfile, 'w') as out:
        gene_line = ""
        for line in inp:
            if line.startswith("#"):
                continue
            fields = line.rstrip('\n').split('\t')

            # Drop mitochondria
            if fields[CHROM_FIELD] == 'chrM':
                continue

            # Store gene line to print if transcript is eligible
            if fields[ELEMENT_FIELD] == "gene":
                gene_line = line
                continue

            # Select protein-coding and canonical transcripts only
            protein, canonical = parse_attributes(fields[ATTRIBUTES_FIELD])
            if protein and canonical:
                out.write(gene_line + line)
                gene_line = ""  # only print gene line before first transcript line


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('gtf', help="Input GTF from GENCODE")
    parser.add_argument('outfile', help="Output filename")
    args = parser.parse_args()

    process(args.gtf, args.outfile)


if __name__ == '__main__':
    main()
