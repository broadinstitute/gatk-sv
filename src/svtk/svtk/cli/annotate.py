#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Annotate resolved SV with genic effects and noncoding hits.

This tool permits the user to provide a GTF of Gencode gene annotations, a BED
of noncoding elements, or both. The BED of noncoding elements must contain four
columns: (chrom, start, end, element_class).

The following classes of genic effects are annotated as new VCF INFO fields if
the SV meets the defined criteria:
    1) LOF (and DUP_LOF) - Loss of function.
        * Deletions are annotated LOF if they overlap any exon.
        * Duplications are annotated DUP_LOF if they reside entirely within
        a gene boundary and overlap any exon.
        * Inversions are annotated LOF if reside entirely within an exon, if
        one breakpoint falls within an exon, if they reside entirely within a
        gene boundary and overlap an exon, or if only one breakpoint falls
        within a gene boundary.
        * Translocations are annotated LOF If they fall within a gene boundary.
    2) COPY_GAIN
        * Duplications are annotated COPY_GAIN if they span the entirety of a
        gene boundary.
    3) INTRONIC
        * Deletions, duplications, and inversions are annotated INTRONIC if
        they are localized to an intron.
    4) DUP_PARTIAL
        * Duplications are annotated DUP_PARTIAL if they overlap the start or
        end of a gene boundary but not its entirety, such that a whole copy of
        the gene is preserved.
    5) INV_SPAN
        * Inversions are annotated INV_SPAN if they overlap the entirety of a
        gene boundary without disrupting it.
    6) NEAREST_TSS
        * Intragenic events are annotated with the nearest transcription start
        site.

An SV is annotated with a new NONCODING INFO field containing all classes of
noncoding elements which the variant overlaps.
"""

import argparse
import sys
import pysam
import pybedtools as pbt
import svtk.annotation as anno


def main(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('vcf', help='Structural variants.')
    parser.add_argument('--gencode', help='Gencode gene annotations (GTF).')
    parser.add_argument('--noncoding', help='Noncoding elements (bed). '
                        'Columns = chr,start,end,element_class,element_name')
    parser.add_argument('annotated_vcf', help='Annotated variants.')

    # Print help if no arguments specified
    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    if args.gencode is None and args.noncoding is None:
        sys.stderr.write('ERROR: Neither Gencode annotations nor noncoding '
                         'elements provided. Must specify at least one to '
                         'annotate.\n\n')
        parser.print_help()
        sys.exit(1)

    vcf = pysam.VariantFile(args.vcf)

    gencode = None if args.gencode is None else pbt.BedTool(args.gencode)
    noncoding = None if args.noncoding is None else pbt.BedTool(args.noncoding)

    anno.annotate_vcf(vcf, gencode, noncoding, args.annotated_vcf)


if __name__ == '__main__':
    main(sys.argv[1:])
