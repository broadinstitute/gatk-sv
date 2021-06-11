#!/usr/bin/env python
import os
import argparse
parser = argparse.ArgumentParser("merge_allosomes.py")
parser.add_argument("batch", type=str, help="")
parser.add_argument("source", type=str, help="")
parser.add_argument("chrom", type=str, help="")


def file_readin(file):
    info = []
    fin = open(file)
    for line in fin:
        pin = line.strip().split()
        info.append(pin)
    fin.close()
    return info


def unify_list(list):
    out = []
    for i in list:
        if i not in out:
            out.append(i)
    return out


def flag_record(record):
    # eg of record: ['MERGED_X_8601', '0.0', '0.0', '0.0']
    # return 'False' if 'NA' or '0.0' in record; else 'True'
    out = 'True'
    if record[2] in ['NA', '0.0'] and record[3] in ['NA', '0.0'] and record[4] in ['NA', '0.0'] and record[5] in ['NA', '0.0']:
        out = 'False'
    return out


def merge_info(fe_info, ma_info):
    out = [fe_info[0]]
    out_key = unify_list([i[0] for i in fe_info[1:]] + [i[0]
                                                        for i in ma_info[1:]])
    out_hash = {}
    for i in fe_info:
        out_hash[i[0]] = []
        out_hash[i[0]].append(i)
    for i in ma_info:
        if not i[0] in out_hash.keys():
            out_hash[i[0]] = []
        out_hash[i[0]].append(i)
    for i in out_hash.keys():
        if len(out_hash[i]) == 1:
            continue
        elif flag_record(out_hash[i][0]) == 'False':
            out_hash[i] = [out_hash[i][1]]
        else:
            out_hash[i] = [out_hash[i][0]]
    for i in out_key:
        out += out_hash[i]
    return out


args = parser.parse_args()
if args.chrom == 'Y':
    file_in = 'srtest_allosomes/' + \
        '.'.join([args.batch, args.source, args.chrom, 'males.stats'])
    file_to = 'srtest/' + \
        '.'.join([args.batch, args.source, args.chrom, 'stats'])
    os.system(r'''cp %s %s''' % (file_in, file_to))
elif args.chrom == 'X':
    female_in = 'srtest_allosomes/' + \
        '.'.join([args.batch, args.source, args.chrom, 'females.stats'])
    male_in = 'srtest_allosomes/' + \
        '.'.join([args.batch, args.source, args.chrom, 'males.stats'])
    file_to = 'srtest/' + \
        '.'.join([args.batch, args.source, args.chrom, 'stats'])
    fe_info = file_readin(female_in)
    ma_info = file_readin(male_in)
    male_femal_merge = merge_info(fe_info, ma_info)
    fo = open(file_to, 'w')
    for i in male_femal_merge:
        print('\t'.join(i), file=fo)
    fo.close()
