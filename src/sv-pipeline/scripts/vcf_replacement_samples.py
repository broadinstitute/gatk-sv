#!/bin/python

"""
Writes list of new sample ids given a replacement dictionary and vcf
"""

import argparse
from pysam import VariantFile


def get_id_dictionary(path):
    with open(path, 'r') as f:
        return {tokens[0]: tokens[1] for tokens in [line.strip().split('\t') for line in f]}


def get_new_ids(vcf, id_dict):
    for sample in vcf.header.samples:
        if sample not in id_dict:
            raise Exception(
                "Header sample not found in dictionary: \"{}\"".format(sample))
        yield id_dict[sample]


def print_ids(sample_ids):
    print("\n".join(sample_ids))


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--vcf', help='Input vcf', required=True)
    parser.add_argument(
        '--dict', help='Tab-delimited sample id conversion table', required=True)

    args = parser.parse_args()
    vcf = VariantFile(args.vcf)
    id_dict = get_id_dictionary(args.dict)
    new_ids = get_new_ids(vcf, id_dict)
    print_ids(new_ids)
    vcf.close()


if __name__ == '__main__':
    main()
