#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
"""

import argparse
import os


def vcf_readin(vcf_name):
    out = []
    fin = os.popen(r'''zcat %s''' % (vcf_name))
    for line in fin:
        pin = line.strip().split()
        if not pin[0][0] == '#':
            out.append(pin)
    fin.close()
    return out


def vcf_name_readin(vcf_name):
    out = []
    fin = os.popen(r'''zcat %s''' % (vcf_name))
    for line in fin:
        pin = line.strip().split()
        if pin[0][0] == '#' and not pin[0][:2] == '##':
            out = pin[8:]
            break
    fin.close()
    return out


def info_cha_extract(pin, character='END'):
    out = ''
    for i in pin[7].split(';'):
        if i.split('=')[0] == character:
            out = i.split('=')[1]
    return out


def unify_list(list):
    out = []
    for i in list:
        if i not in out:
            out.append(i)
    return out


def extract_non_ref_samples(pin, sample_name):
    GT_POS = pin[8].split(':').index('GT')
    sample_GT = [i.split(':')[GT_POS] for i in pin[9:]]
    sample_out = [i for i in sample_name[1:]
                  if not sample_GT[sample_name.index(i) - 1] == '0/0']
    return ','.join(sample_out)


def write_output(final_out, fileout):
    fo = open(fileout, 'w')
    for i in final_out:
        i_new = [i[0].split('_')[1].split(':')[0]] + \
            i[0].split('_')[1].split(':')[1].split('-') + i[1:]
        print('\t'.join(i_new), file=fo)
    fo.close()


def extract_CNV_from_cpx(vcf_name, bed_name):
    # vcf_info=vcf_readin(vcf_name)
    sample_name = vcf_name_readin(vcf_name)
    final_out = []
    # take out duplicated region in DUP5/INS3
    fin = os.popen(r'''zcat %s''' % (vcf_name))
    for line in fin:
        pin = line.strip().split()
        if not pin[0][0] == '#':
            i = pin
            print(i[:8])
            if info_cha_extract(i, 'CPX_TYPE') == 'DUP5/INS3':
                # take out duplicated region in DUP5/INS3
                final_out.append([info_cha_extract(
                    i, 'SOURCE'), i[2], extract_non_ref_samples(i, sample_name), 'DUP'])
            if info_cha_extract(i, 'CPX_TYPE') == 'DUP3/INS5':
                # take out duplicated region in DUP3/INS5
                final_out.append([info_cha_extract(
                    i, 'SOURCE'), i[2], extract_non_ref_samples(i, sample_name), 'DUP'])
            if info_cha_extract(i, 'CPX_TYPE') in ['delINV', 'delINVdel', 'INVdel', 'dupINV', 'INVdup', 'dupINVdup', 'delINVdup']:
                # take out deleted and inverted region in delINV,delINVdel,INVdel, dupINVdel
                intervals = info_cha_extract(i, 'CPX_INTERVALS').split(',')
                for j in intervals:
                    final_out.append([j, i[2], extract_non_ref_samples(
                        i, sample_name), j.split('_')[0]])
            if info_cha_extract(i, 'CPX_TYPE') in ['INS_B2A', 'INS_A2B', 'CTX_INV_INS_B2A', 'CTX_INS_B2A', 'CTX_INS_A2B', 'CTX_INV_INS_A2B']:
                # take out deleted and inverted region in INS_B2A,INS_A2B
                final_out.append([info_cha_extract(
                    i, 'SOURCE'), i[2], extract_non_ref_samples(i, sample_name), 'DUP'])
    fin.close()
    write_output(final_out, bed_name)


def main():
    parser = argparse.ArgumentParser(
        description='Script to extract CNVs involved in complex SVs.')
    parser.add_argument('vcf', help='name of input vcf file in gzipped format')
    parser.add_argument('bed_name', help='name of output bed file')
    args = parser.parse_args()
    vcf_name = args.vcf
    # bed_name=vcf_name.replace('.vcf.gz','.CSV_CNV.bed')
    extract_CNV_from_cpx(vcf_name, args.bed_name)


if __name__ == '__main__':
    main()
