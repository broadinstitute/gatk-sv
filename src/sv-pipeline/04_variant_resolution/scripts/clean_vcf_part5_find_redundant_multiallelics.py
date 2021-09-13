#!/usr/bin/env python

import argparse
import svtk.utils as svu
import sys


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('multiallelic_filename')
    parser.add_argument('fout')
    args = parser.parse_args()

    print("finding redundant overlapping sites", file=sys.stderr)
    multiallelic_bed = svu.vcf2bedtool(args.multiallelic_filename, include_filters=True)

    redundant_multiallelics = set()
    self_inter = multiallelic_bed.intersect(multiallelic_bed, wo=True)\
        .filter(lambda feature: feature[3] != feature[10])
    for feature in self_inter:
        a_len = int(feature.fields[2]) - int(feature.fields[1])
        b_len = int(feature.fields[9]) - int(feature.fields[8])
        overlap = int(feature.fields[14])
        small_coverage = overlap / min(a_len, b_len)
        if small_coverage > 0.50:
            if a_len < b_len:
                redundant_multiallelics.add(feature.fields[3])
            else:
                redundant_multiallelics.add(feature.fields[10])
    print("identified {} redundant multiallelic sites".format(len(redundant_multiallelics)), file=sys.stderr)
    with open(args.fout, "w") as list_file:
        for vid in redundant_multiallelics:
            print(vid, file=list_file)


if __name__ == '__main__':
    main()
