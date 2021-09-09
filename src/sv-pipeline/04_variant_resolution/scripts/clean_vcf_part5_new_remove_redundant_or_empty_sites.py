#!/usr/bin/env python

import argparse
from collections import Counter
import gzip
import pysam
import svtk.utils as svu
import sys



def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('cleangq_filename')
    parser.add_argument('fout')
    args = parser.parse_args()

    print("finding redundant overlapping sites", file=sys.stderr)
    cleangq_bed = svu.vcf2bedtool(cleangq_filename, include_filters=True)

    multiallelic_bed = cleangq_bed.filter(lambda feature: 'MULTIALLELIC' in feature.fields[6].split(',')).saveas('multiallelics.bed')

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

    # one more pass through the VCF to remove variants with no called samples and the redundant multiallelics
    cleangq_vcf = pysam.VariantFile(cleangq_filename)

    if args.fout in '- stdout'.split():
        fout = pysam.VariantFile(sys.stdout, 'w', header=cleangq_vcf.header)
    else:
        fout = pysam.VariantFile(args.fout, 'w', header=cleangq_vcf.header)

    print("removing redundant overlapping sites", file=sys.stderr)
    for idx, record in enumerate(cleangq_vcf):
        if (idx - 1) % 1000 == 0:
            print("processed {} records".format(idx), file=sys.stderr)
        if record.id in redundant_multiallelics or len(svu.get_called_samples(record)) == 0:
            continue
        fout.write(record)
    fout.close()
    cleangq_vcf.close()
    print("done", file=sys.stderr)

if __name__ == '__main__':
    main()
