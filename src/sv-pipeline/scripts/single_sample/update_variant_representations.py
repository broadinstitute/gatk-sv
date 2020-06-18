#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Update variant VCF representations to be compatible with the letter and spirit of the VCF spec when possible.
Currently the major change made are to:
 -populate the reference allele with the actual reference base indicated by POS
 -re-write records with SVTYPE=BND, ALT=<BND>, and CHR2 and END2 INFO fields as two records with SVTYPE BND located at
  CHROM:POS and CHR2:END2 with alternate alleles in the BND alt format.
"""

import argparse
import pysam
import sys
import svtk.utils as svu


def get_ref_base(contig, pos, ref):
    return ref.fetch(contig, pos-1, pos)


def has_enough_info_to_make_bnd_mates(record):
    info_keys = record.info.keys()
    return 'STRANDS' in info_keys and 'CHR2' in info_keys and \
           ('END2' in info_keys or (record.info['CHR2'] == record.contig and 'SVLEN' in info_keys))


def make_paired_bnd_records(record, ref_fasta):
    chr2 = record.info['CHR2']
    if 'END2' in record.info.keys():
        end2 = record.info['END2']
    elif record.contig == record.info['CHR2']:
        end2 = record.pos + record.info['SVLEN']

    orig_id = record.id
    new_id = orig_id + "_M1"
    mate_id = orig_id + "_M2"
    strands = record.info['STRANDS']

    alt1 = svu.make_bnd_alt(chr2, end2, strands, ref_base=record.ref)
    record.alts = (alt1, )
    record.id = new_id
    record.info['MATEID'] = mate_id
    record.info.pop('CHR2')
    if 'END2' in record.info.keys():
        record.info.pop('END2')
    record.info.pop('SVLEN')

    mate_record = record.copy()
    mate_record.id = mate_id
    mate_record.contig = chr2
    mate_record.pos = int(end2)
    mate_record.ref = get_ref_base(chr2, int(end2), ref_fasta)
    mate_record.info['MATEID'] = new_id
    # strands are reversed for the mate BND
    mate_strands = strands[::-1]
    mate_alt = svu.make_bnd_alt(record.contig, record.pos, mate_strands, ref_base=mate_record.ref)
    mate_record.alts = (mate_alt, )
    return record, mate_record


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf', help='Input VCF')
    parser.add_argument('ref', help='Reference fasta')
    parser.add_argument('--ref_idx', help='Reference fasta index')
    parser.add_argument('-o', '--outfile', help='Output file [default: stdout]')

    args = parser.parse_args()

    if args.vcf in '- stdin'.split():
        vcf = pysam.VariantFile(sys.stdin)
    else:
        vcf = pysam.VariantFile(args.vcf)
    header = vcf.header

    header.add_line('##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">')

    # Open connection to outfile
    if args.outfile is None:
        fout = pysam.VariantFile(sys.stdout, 'w', header=header)
    else:
        out = args.outfile
        if '.gz' in out or '.bgz' in out:
            out = path.splitext(out)[0]
        fout = pysam.VariantFile(out, 'w', header=header)

    if args.ref_idx is None:
        ref_fasta = pysam.FastaFile(args.ref)
    else:
        ref_fasta = pysam.FastaFile(args.ref, args.ref_idx)

    for record in vcf:
        svtype = record.info['SVTYPE']
        ref_base = get_ref_base(record.contig, record.pos, ref_fasta)
        record.ref = ref_base
        alt = record.alts[0]
        if svtype == 'BND' and alt == '<BND>' and has_enough_info_to_make_bnd_mates(record):
            record, mate_record = make_paired_bnd_records(record, ref_fasta)
            fout.write(record)
            fout.write(mate_record)
        else:
            fout.write(record)

    fout.close()


if __name__ == '__main__':
    main()
