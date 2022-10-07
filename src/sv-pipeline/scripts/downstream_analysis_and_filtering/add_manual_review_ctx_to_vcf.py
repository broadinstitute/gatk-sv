#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Adds records representing manually reviewed translocation events to VCF.
Any events that match exactly on CHR1, POS, CHR2, END2, and CPX_TYPE will be merged into a single
record with the correct sample genotypes set.
The new record will have 'tloc' as a member of the ALGORITHMS INFO field
Genotypes for the new record are given GQ and PE_GQ values of 999.
Updates CTX records already present in the to store the location on CHR2 in the END2 info field and
sets their END info field to be CHR1 position + 1
"""

import argparse
import sys
import pandas as pd
import pysam


# removed type hints for backward compatibility with earlier python versions:
# loc1: tuple[str, int], loc2: tuple[str, int], contig_list: list[str]
def before() -> bool:
    return (contig_list.index(loc1[0]) < contig_list.index(loc2[0])) or \
        (contig_list.index(loc1[0]) == contig_list.index(loc2[0]) and
         loc1[1] < loc2[1])


def parse_info(data: pd.DataFrame):
    info_mapping = [{"info_" + token.split("=")[0]: token.split("=")[1]
                     for token in info_str.split(";")} for info_str in data['info']]
    info_df = pd.DataFrame(info_mapping, index=data.index)
    return pd.concat([data, info_df], axis=1)


# removed type hints for backward compatibility with earlier python versions:
# -> tuple[pysam.VariantRecord, int]
def create_record(new_idx: int,
                  data: pd.DataFrame,
                  fout: pysam.VariantFile,
                  ctx_number: int,
                  cohort_name: str,
                  batch_name: str):
    row = data.iloc[new_idx]
    contig = row["contig"]
    pos1 = row["pos"]
    id = new_record_id(cohort_name, batch_name, contig, ctx_number)
    ref = row["ref"]
    alt = row["alt"]
    filter = row["filter"]
    info_cols = [colname for colname in data.columns if colname.startswith("info_")]

    samples = {row['sample']: row['gt']}
    info_sets = {info_col: {row[info_col]} for info_col in info_cols}
    if 'info_ALGORITHMS' in info_sets:
        info_sets['info_ALGORITHMS'].add('tloc')
    else:
        info_sets['info_ALGORITHMS'] = {'tloc'}

    new_idx = new_idx + 1
    while new_idx < len(data.index) and \
            data.iloc[new_idx]['contig'] == contig and \
            data.iloc[new_idx]['pos'] == pos1 and \
            data.iloc[new_idx]['info_CHR2'] == row['info_CHR2'] and \
            data.iloc[new_idx]['info_END'] == row['info_END'] and \
            data.iloc[new_idx]['info_CPX_TYPE'] == row['info_CPX_TYPE']:
        info_sets['info_MEMBERS'].add(data.iloc[new_idx]['info_MEMBERS'])
        samples[data.iloc[new_idx]['sample']] = data.iloc[new_idx]['gt']
        new_idx = new_idx + 1

    new_record = fout.new_record(contig=contig,
                                 start=pos1,
                                 id=id,
                                 stop=pos1 + 1,
                                 alleles=(ref, alt),
                                 filter=filter)

    for info_key in info_sets:
        new_info_key = info_key.split("info_")[1]
        if new_info_key == 'END':
            continue
        if len(info_sets[info_key]) > 1:
            new_record.info[new_info_key] = ",".join(info_sets[info_key])
        else:
            elem = next(iter(info_sets[info_key]))
            if elem.lstrip("-+").isdigit():
                elem = int(elem)
            new_record.info[new_info_key] = elem

    for sample in samples:
        new_record.samples[sample]['GT'] = tuple(map(int, samples[sample].split(':')[0].split('/')))
        new_record.samples[sample]['CN'] = None
        new_record.samples[sample]['CNQ'] = None
        new_record.samples[sample]['EV'] = row['info_EVIDENCE']
        new_record.samples[sample]['GQ'] = 999
        new_record.samples[sample]['PE_GQ'] = 999
        new_record.samples[sample]['PE_GT'] = sum(map(int, samples[sample].split(':')[0].split('/')))
        new_record.samples[sample]['RD_CN'] = None
        new_record.samples[sample]['RD_GQ'] = None
        new_record.samples[sample]['SR_GQ'] = None
        new_record.samples[sample]['SR_GT'] = None

    new_record.info['END2'] = int(row["info_END"])
    return new_record, new_idx


def new_record_id(cohort_name: str, batch_name: str, contig: str, ctx_num: int) -> str:
    return "{}.{}_CTX_{}_{}".format(cohort_name, batch_name, contig, ctx_num)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--vcf', '-V', help='Input vcf (supports "stdin").', required=True)
    parser.add_argument('--reviewed-ctx-file', help='File with data on reviewed translocation events', required=True)
    parser.add_argument('--cohort-name', help='Name of the cohort (to be used in new record IDs)', required=True)
    parser.add_argument('--batch-name', help='Name of the batch run (to be used in new record IDs)', required=True)
    parser.add_argument('--out', '-O', help='Output file (supports "stdout").', required=True)

    args = parser.parse_args()

    if args.vcf in '- stdin'.split():
        vcf = pysam.VariantFile(sys.stdin)
    else:
        vcf = pysam.VariantFile(args.vcf)

    header = vcf.header

    contig_list = [contig.name for contig in vcf.header.contigs.values()]

    if args.out in '- stdout'.split():
        fout = pysam.VariantFile(sys.stdout, 'w', header=header)
    else:
        fout = pysam.VariantFile(args.out, 'w', header=header)

    data = pd.read_table(args.reviewed_ctx_file,
                         header=None,
                         names=['contig', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', 'gt', 'sample'])

    data = parse_info(data)
    data = data.sort_values(by=['contig', 'pos', 'info_CHR2', 'info_END', 'info_CPX_TYPE'])

    ctx_number = 0
    new_idx = 0
    prev_contig = None
    prev_pos = None

    new_loc = (data.loc[new_idx, 'contig'], data.loc[new_idx, 'pos'])

    for record in vcf:
        if record.contig != prev_contig:
            ctx_number = 0
        while new_idx < len(data.index) and\
                before(new_loc, (record.contig, record.pos), contig_list=contig_list) and \
                not before(new_loc, (prev_contig, prev_pos), contig_list=contig_list):
            ctx_number = ctx_number + 1
            new_record, new_idx = create_record(new_idx, data, fout, ctx_number, args.cohort_name, args.batch_name)
            fout.write(new_record)
            if new_idx < len(data.index):
                new_loc = (data.loc[new_idx, 'contig'], data.loc[new_idx, 'pos'])
        # adjust CTX records already in the VCF to store the location on CHR2 in the END2 INFO field
        if record.info['SVTYPE'] == "CTX":
            ctx_number = ctx_number + 1
            record.id = new_record_id(args.cohort_name, args.batch_name, record.contig, ctx_number)
            if not 'END2' in record.info:
                raise Exception("Record {}:{} {} does not have END2 set".format(record.contig, record.pos, record.id))
            else:
                record.stop = record.pos + 1
        fout.write(record)
        prev_contig = record.contig
        prev_pos = record.pos


if __name__ == '__main__':
    main()
