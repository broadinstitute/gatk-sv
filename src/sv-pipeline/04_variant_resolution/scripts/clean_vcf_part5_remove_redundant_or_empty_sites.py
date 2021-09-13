#!/usr/bin/env python

import argparse
import pysam
import svtk.utils as svu
import sys


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('cleangq_filename')
    parser.add_argument('redundant_multiallelics_file')
    parser.add_argument('fout')
    args = parser.parse_args()

    redundant_multiallelics = set([line.rstrip() for line in open(args.redundant_multiallelics_file, 'rt')])

    # one more pass through the VCF to remove variants with no called samples and the redundant multiallelics
    cleangq_vcf = pysam.VariantFile(args.cleangq_filename)

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
