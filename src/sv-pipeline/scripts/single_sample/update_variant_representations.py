#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Update variant VCF representations to be compatible with the letter and spirit of the VCF spec when possible.
Currently the major change made are to:
 -populate the reference allele with the actual reference base indicated by POS
 -re-write records with SVTYPE=BND, ALT=<BND>, and CHR2 and END2 INFO fields as two records with SVTYPE BND located at
  CHROM:POS and CHR2:END2 with alternate alleles in the BND alt format.
 -re-write records with SVTYPE=CTX as two pairs of BND records
"""

import argparse
import os
import pysam
import sys
import svtk.utils as svu


def get_ref_base(contig, pos, ref):
    return ref.fetch(contig, pos - 1, pos)


def has_enough_info_to_make_bnd_mates(record):
    info_keys = record.info.keys()
    return 'STRANDS' in info_keys and 'CHR2' in info_keys and \
           ('END2' in info_keys or (
               record.info['CHR2'] == record.contig and 'SVLEN' in info_keys))


def has_enough_info_to_make_ctx_bnds(record):
    info_keys = record.info.keys()
    return 'CPX_TYPE' in info_keys \
           and (record.info['CPX_TYPE'] == "CTX_PQ/QP" or record.info['CPX_TYPE'] == "CTX_PP/QQ") \
           and 'CHR2' in info_keys


def make_paired_bnd_records(record, ref_fasta):
    chr2 = record.info['CHR2']
    if 'END2' in record.info.keys():
        end2 = record.info['END2']
    elif record.contig == record.info['CHR2'] and 'SVLEN' in record.info.keys():
        end2 = record.pos + record.info['SVLEN']
    else:
        raise ValueError("BND record must contain END2 or if it is intrachromosomal then SVLEN, "
                         f"but {record.id} contains neither")

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
    if 'SVLEN' in record.info:
        record.info.pop('SVLEN')

    mate_record = record.copy()
    mate_record.id = mate_id
    mate_record.contig = chr2
    mate_record.pos = int(end2)
    mate_record.ref = get_ref_base(chr2, int(end2), ref_fasta)
    mate_record.info['MATEID'] = new_id
    # strands are reversed for the mate BND
    mate_strands = strands[::-1]
    mate_alt = svu.make_bnd_alt(
        record.contig, record.pos, mate_strands, ref_base=mate_record.ref)
    mate_record.alts = (mate_alt, )
    return record, mate_record


def make_reciprocal_translocation_bnds(record, ref_fasta):
    chr1 = record.contig
    chr1_pos1 = record.pos
    chr1_pos1_ref = get_ref_base(chr1, chr1_pos1, ref_fasta)
    chr1_pos2 = chr1_pos1 + 1
    chr1_pos2_ref = get_ref_base(chr1, chr1_pos2, ref_fasta)
    chr2 = record.info['CHR2']
    if 'END2' in record.info.keys():
        chr2_pos1 = int(record.info['END2'])
    else:
        chr2_pos1 = record.stop
    chr2_pos1_ref = get_ref_base(chr2, chr2_pos1, ref_fasta)
    chr2_pos2 = chr2_pos1 + 1
    chr2_pos2_ref = get_ref_base(chr2, chr2_pos2, ref_fasta)

    orig_id = record.id
    m1_id = orig_id + "_M1"  # M1: chr1P->chr2 on chr1
    m2_id = orig_id + "_M2"  # M2: chr1P->chr2 on chr2
    m3_id = orig_id + "_M3"  # M3: chr1Q->chr2 on chr1
    m4_id = orig_id + "_M4"  # M3: chr1Q->chr2 on chr2

    event_id = orig_id

    if record.info['CPX_TYPE'] == "CTX_PQ/QP":
        strands_m1 = '+-'
        strands_m2 = '-+'
        strands_m3 = '-+'
        strands_m4 = '+-'
        m1_pos1 = chr1_pos1
        m1_ref = chr1_pos1_ref
        m1_pos2 = chr2_pos2
        m2_pos1 = chr2_pos2
        m2_ref = chr2_pos2_ref
        m2_pos2 = chr1_pos1
        m3_pos1 = chr1_pos2
        m3_ref = chr1_pos2_ref
        m3_pos2 = chr2_pos2
        m4_pos1 = chr2_pos1
        m4_ref = chr1_pos1_ref
        m4_pos2 = chr1_pos1
    elif record.info['CPX_TYPE'] == "CTX_PP/QQ":
        strands_m1 = '++'
        strands_m2 = '++'
        strands_m3 = '--'
        strands_m4 = '--'
        m1_pos1 = chr1_pos1
        m1_ref = chr1_pos1_ref
        m1_pos2 = chr2_pos1
        m2_pos1 = chr2_pos1
        m2_ref = chr2_pos1_ref
        m2_pos2 = chr1_pos1
        m3_pos1 = chr1_pos2
        m3_ref = chr1_pos2_ref
        m3_pos2 = chr2_pos2
        m4_pos1 = chr2_pos2
        m4_ref = chr2_pos2_ref
        m4_pos2 = chr1_pos2

    record.info['EVENT'] = event_id
    record.info['SVTYPE'] = 'BND'
    record.info.pop('CHR2')
    if 'SVLEN' in record.info:
        record.info.pop('SVLEN')

    m1 = record.copy()
    m1.id = m1_id
    m1.contig = chr1
    m1.pos = m1_pos1
    m1.stop = m1_pos1 + 1
    m1.ref = m1_ref
    m1.alts = (svu.make_bnd_alt(chr2, m1_pos2, strands_m1, ref_base=m1.ref), )
    m1.info['MATEID'] = m2_id

    m2 = record.copy()
    m2.id = m2_id
    m2.contig = chr2
    m2.pos = m2_pos1
    m2.stop = m2_pos1 + 1
    m2.ref = m2_ref
    m2.alts = (svu.make_bnd_alt(chr1, m2_pos2, strands_m2, ref_base=m2.ref), )
    m2.info['MATEID'] = m2_id

    m3 = record.copy()
    m3.id = m3_id
    m3.contig = chr1
    m3.pos = m3_pos1
    m3.stop = m3_pos1 + 1
    m3.ref = m3_ref
    m3.alts = (svu.make_bnd_alt(chr2, m3_pos2, strands_m3, ref_base=m3.ref), )
    m3.info['MATEID'] = m4_id

    m4 = record.copy()
    m4.id = m4_id
    m4.contig = chr2
    m4.pos = m4_pos1
    m4.stop = m4_pos1 + 1
    m4.ref = m4_ref
    m4.alts = (svu.make_bnd_alt(chr1, m4_pos2, strands_m4, ref_base=m4.ref), )
    m4.info['MATEID'] = m3_id

    return m1, m2, m3, m4


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf', help='Input VCF')
    parser.add_argument('ref', help='Reference fasta')
    parser.add_argument('--ref_idx', help='Reference fasta index')
    parser.add_argument('-o', '--outfile',
                        help='Output file [default: stdout]')

    args = parser.parse_args()

    if args.vcf in '- stdin'.split():
        vcf = pysam.VariantFile(sys.stdin)
    else:
        vcf = pysam.VariantFile(args.vcf)
    header = vcf.header

    header.add_line(
        '##INFO=<ID=MATEID,Number=.,Type=String,Description="ID of mate breakends">')
    header.add_line(
        '##INFO=<ID=EVENT,Number=1,Type=String,Description="ID of event associated to breakend">')

    # Open connection to outfile
    if args.outfile is None:
        fout = pysam.VariantFile(sys.stdout, 'w', header=header)
    else:
        out = args.outfile
        if '.gz' in out or '.bgz' in out:
            out = os.path.splitext(out)[0]
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
        elif svtype == 'CTX' and has_enough_info_to_make_ctx_bnds(record):
            m1, m2, m3, m4 = make_reciprocal_translocation_bnds(
                record, ref_fasta)
            fout.write(m1)
            fout.write(m3)
            fout.write(m2)
            fout.write(m4)
        else:
            fout.write(record)

    fout.close()


if __name__ == '__main__':
    main()
