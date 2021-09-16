#!/usr/bin/env python

import argparse
import sys
import svtk.utils as svu


def process_features_for_size1(features_for_size1, redundant_multiallelics):
    for intersection in sorted(features_for_size1, key=lambda x: int(x[9]) - int(x[8]), reverse=True):
        b_len = int(intersection.fields[9]) - int(intersection.fields[8])
        overlap = int(intersection.fields[14])
        small_coverage = overlap / b_len
        if small_coverage > 0.50:
            if intersection.fields[3] not in redundant_multiallelics:
                redundant_multiallelics.add(intersection.fields[10])


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
        .filter(lambda feature: feature[3] != feature[10]) \
        .filter(lambda feature: (int(feature[2]) - int(feature[1])) >= (int(feature[9]) - int(feature[8]))) \
        .sort(sizeD=True)
    current_size1 = -1
    features_for_size1 = []
    for feature in self_inter:
        size1 = int(feature[2]) - int(feature[1])
        if size1 != current_size1:
            process_features_for_size1(features_for_size1, redundant_multiallelics)
            features_for_size1 = []

        current_size1 = size1
        features_for_size1.append(feature)

    process_features_for_size1(features_for_size1, redundant_multiallelics)
    print("identified {} redundant multiallelic sites".format(len(redundant_multiallelics)), file=sys.stderr)
    with open(args.fout, "w") as list_file:
        for vid in redundant_multiallelics:
            print(vid, file=list_file)


if __name__ == '__main__':
    main()
