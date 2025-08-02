#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Perform preliminary aesthetic cleanup to a VCF after svtk resolve
"""

import argparse
import sys
import pysam
from svtk.utils import parse_bnd_pos

INS_ALT_INFO = [
    '##ALT=<ID=INS:ME,Description="Mobile element insertion of unspecified ME class">',
    '##ALT=<ID=INS:ME:ALU,Description="Alu element insertion">',
    '##ALT=<ID=INS:ME:SVA,Description="SVA element insertion">',
    '##ALT=<ID=INS:ME:LINE1,Description="LINE1 element insertion">',
    '##ALT=<ID=INS:UNK,Description="Sequence insertion of unspecified origin">',
]

# info keys we need to make sure we have
INFO_KEYS = {
    'END2': '##INFO=<ID=END2,Number=1,Type=Integer,Description="Position of breakpoint on CHR2">',
    'CHR2': '##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END2 coordinate">',
    'UNRESOLVED': '##INFO=<ID=UNRESOLVED,Number=0,Type=Flag,Description="Variant is unresolved.">',
    'UNRESOLVED_TYPE': '##INFO=<ID=UNRESOLVED_TYPE,Number=1,Type=String,Description="Class of unresolved variant.">'
}

# Cleanup function


def cleanup(vcf, fout):

    # Iterate over records
    for record in vcf:

        # Skip records with SVLEN < 50bp
        if 0 <= record.info['SVLEN'] < 50:
            continue

        # Get basic info about record
        # chrom = record.chrom
        svtype = record.info['SVTYPE']
        # end = record.stop
        record.ref = 'N'
        # alts = record.alts

        # Clean up insertions
        if svtype == 'INS':

            # Set alt allele to specify kind of insertion, if known
            metypes = 'INS:ME INS:MEI INS:ME:ALU INS:ME:SVA INS:ME:LINE1'.split()
            if record.alts[0] not in metypes:

                # If ME specified in CPX_TYPE, store that as alt allele
                if 'CPX_TYPE' in record.info.keys():
                    if record.info['CPX_TYPE'] in metypes:
                        record.alts = ('<{0}>'.format(
                            record.info['CPX_TYPE']), )
                    elif record.info['CPX_TYPE'] in 'MEI_INS5 MEI_INS3'.split():
                        record.alts = ('<INS:ME>', )
                    else:
                        record.alts = ('<INS:UNK>', )

                # Otherwise, set alt allele to insertion of unknown origin
                else:
                    record.alts = ('<INS:UNK>', )
            else:
                record.alts = ('<INS:UNK>', )

            # Clear CPX_TYPE=INS (but leave others, as they're necessary for CPX GT)
            if 'CPX_TYPE' in record.info.keys():
                if record.info['CPX_TYPE'] in 'INS MEI_INS5 MEI_INS3'.split() \
                        or record.info['CPX_TYPE'] in metypes:
                    record.info.pop('CPX_TYPE')

            # Scrub all unresolved info from all insertions
            for info in 'UNRESOLVED UNRESOLVED_TYPE EVENT STRANDS'.split():
                if info in record.info.keys():
                    record.info.pop(info)

        # Clean up BNDs
        elif svtype == 'BND':

            # copy bnd data to CHR2, END2
            if "[" in record.alts[0] or "]" in record.alts[0]:
                chr2, end2 = parse_bnd_pos(record.alts[0])
                record.info['CHR2'] = chr2
                record.info['END2'] = end2
            else:
                record.info['END2'] = record.stop

            # Correct alt syntax
            record.alts = ('<BND>', )
            record.stop = record.pos + 1

            # All BNDs are unresolved by definition
            record.info['UNRESOLVED'] = True

            # Clean UNRESOLVED_TYPE
            trigger = 0
            if 'UNRESOLVED_TYPE' not in record.info.keys():
                trigger = 1
            elif record.info['UNRESOLVED_TYPE'] in 'UNK NOT_SPECIFIED'.split():
                trigger = 1

            if trigger == 0:

                # Add strand to SINGLE_ENDER notation, if possible
                unres = record.info['UNRESOLVED_TYPE']
                if unres == 'SINGLE_ENDER':
                    if 'STRANDS' in record.info.keys():
                        unres = 'SINGLE_ENDER_' + record.info['STRANDS']
                    else:
                        unres = 'SINGLE_ENDER_UNK'

            # Otherwise, assign UNRESOLVED_TYPE based on number of MEMBERS
            else:
                if 'MEMBERS' in record.info.keys():
                    if len(record.info['MEMBERS']) > 1:
                        unres = 'MIXED_BREAKENDS'
                    else:
                        if 'STRANDS' in record.info.keys():
                            unres = 'SINGLE_ENDER_' + record.info['STRANDS']
                        else:
                            unres = 'SINGLE_ENDER_UNK'

                # If record missing MEMBERS tag, assume single ender
                else:
                    if 'STRANDS' in record.info.keys():
                        unres = 'SINGLE_ENDER_' + record.info['STRANDS']
                    else:
                        unres = 'SINGLE_ENDER_UNK'

            # Set final UNRESOLVED_TYPE
            record.info['UNRESOLVED_TYPE'] = unres

        # Clean up INVs
        elif svtype == 'INV':

            # Make UNRESOLVED inversions SVTYPE=BNDs (retain INV alt allele)
            if 'UNRESOLVED' in record.info.keys():

                # All unresolved inversions become BNDs
                record.info['SVTYPE'] = 'BND'

            # Ensure all inversion single-enders are marked as UNRESOLVED
            else:

                # INV records without CPX_TYPE are unmarked unresolved single enders
                if 'CPX_TYPE' not in record.info.keys():
                    record.info['SVTYPE'] = 'BND'
                    record.info['UNRESOLVED'] = True
                    record.info['UNRESOLVED_TYPE'] = 'INVERSION_SINGLE_ENDER_' + \
                        record.info['STRANDS']

                # INV records with CPX_TYPE are resolved simple inversions
                # In these cases, strip CPX_TYPE and CPX_INTERVALS
                else:
                    for info in 'CPX_TYPE CPX_INTERVALS'.split():
                        if info in record.info.keys():
                            record.info.pop(info)

        # Remove various tags from CNVs
        elif svtype in 'DEL DUP'.split():
            for info in 'CPX_TYPE CPX_INTERVALS UNRESOLVED UNRESOLVED_TYPE'.split():
                if info in record.info.keys():
                    record.info.pop(info)

        # Clean some tags from all variants
        for info in 'EVENT'.split():
            if info in record.info.keys():
                record.info.pop(info)

        # Make sure all UNRESOLVED variants also have an UNRESOLVED_TYPE
        if 'UNRESOLVED' in record.info.keys():
            if 'UNRESOLVED_TYPE' not in record.info.keys():
                record.info['UNRESOLVED_TYPE'] = 'NOT_SPECIFIED'
            elif 'UNRESOLVED_TYPE' == 'UNK':
                record.info['UNRESOLVED_TYPE'] = 'NOT_SPECIFIED'

        # Write record to file
        fout.write(record)


# Main block
def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf')
    parser.add_argument('fout')

    args = parser.parse_args()

    if args.vcf in '- stdin'.split():
        vcf = pysam.VariantFile(sys.stdin)
    else:
        vcf = pysam.VariantFile(args.vcf)
    for line in INS_ALT_INFO:
        vcf.header.add_line(line)

    for key in INFO_KEYS.keys():
        if key not in vcf.header.info.keys():
            vcf.header.add_line(INFO_KEYS[key])

    if args.fout in '- stdout'.split():
        fout = pysam.VariantFile(sys.stdout, 'w', header=vcf.header)
    else:
        fout = pysam.VariantFile(args.fout, 'w', header=vcf.header)

    cleanup(vcf, fout)


if __name__ == '__main__':
    main()
