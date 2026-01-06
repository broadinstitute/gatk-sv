#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""

"""

import argparse
import numpy as np
import pandas as pd
import pysam


def make_evidence_int(ev):
    ev = ev.split(',')

    evidence = 0
    if 'RD' in ev:
        evidence += 1
    if 'PE' in ev:
        evidence += 2
    if 'SR' in ev:
        evidence += 4

    return evidence


def add_genotypes(record, genotypes, varGQ):
    """
    Update formats and add genotypes
    """

    # tag varGQ
    record.info['varGQ'] = int(np.round(varGQ['GQ']))

    # clear out unused formats (e.g. algorithm keys)
    for fmt in record.format.keys():
        if fmt not in 'GT GQ PE_GT PE_GQ SR_GT SR_GQ RD_CN RD_GQ EV'.split():
            del record.format[fmt]

    max_GT = genotypes['GT'].max()

    if max_GT > 2:
        record.alts = ('<CNV>',)
        if record.info['SVTYPE'] != 'DUP':
            msg = 'Invalid SVTYPE {0} for multiallelic record {1}'
            msg = msg.format(record.info['SVTYPE'], record.id)
            raise Exception(msg)
        record.info['SVTYPE'] = 'CNV'

    cols = 'name sample GT GQ RD_CN RD_GQ PE_GT PE_GQ SR_GT SR_GQ EV'.split()
    gt_matrix = genotypes.reset_index()[cols].to_numpy()

    # update genotype and other data for each sample
    for j, sample in enumerate(record.samples):
        #  data = genotypes.iloc[j]
        data = gt_matrix[j]

        #  if record.id != data.name[0] or sample != data.name[1]:
        if record.id != data[0] or sample != data[1]:
            msg = 'iloc failed. expected ({0}, {1}) found ({2}, {3}).'
            msg = msg.format(record.id, sample, data[0], data[1])
            raise Exception(msg)

        if max_GT > 2:
            record.samples[sample]['GT'] = (None, None)
        elif data[2] == 0:
            record.samples[sample]['GT'] = (0, 0)
        elif data[2] == 1:
            record.samples[sample]['GT'] = (0, 1)
        elif data[2] == 2:
            record.samples[sample]['GT'] = (1, 1)
        else:
            raise Exception('Invalid GT')

        for i, fmt in enumerate('GQ RD_CN RD_GQ PE_GT PE_GQ SR_GT SR_GQ'.split()):
            if data[i + 3] == '.':
                record.samples[sample][fmt] = None
            else:
                record.samples[sample][fmt] = int(data[i + 3])

        record.samples[sample]['EV'] = make_evidence_int(data[10])

    # Bug in pysam - sometimes gets rid of the first EV field
    data = gt_matrix[0]
    record.samples[0]['EV'] = make_evidence_int(data[10])


def update_vcf_header(header):
    """
    Add genotype FORMAT keys to header
    """
    for key in header.formats.keys():
        if key not in 'GT GQ'.split():
            header.formats.remove_header(key)

    if 'GT' not in header.formats.keys():
        header.add_line(
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    if 'GQ' not in header.formats.keys():
        header.add_line(
            '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">')

    header.add_line(
        '##FORMAT=<ID=RD_CN,Number=1,Type=Integer,Description="Predicted copy state">')
    header.add_line(
        '##FORMAT=<ID=RD_GQ,Number=1,Type=Integer,Description="Read-depth genotype quality">')
    header.add_line(
        '##FORMAT=<ID=PE_GT,Number=1,Type=Integer,Description="Paired-end genotype">')
    header.add_line(
        '##FORMAT=<ID=PE_GQ,Number=1,Type=Integer,Description="Paired-end genotype quality">')
    header.add_line(
        '##FORMAT=<ID=SR_GT,Number=1,Type=Integer,Description="Split-read genotype">')
    header.add_line(
        '##FORMAT=<ID=SR_GQ,Number=1,Type=Integer,Description="Split read genotype quality">')
    header.add_line(
        '##FORMAT=<ID=EV,Number=1,Type=Integer,Description="Classes of evidence supporting final genotype">')

    header.add_line(
        '##INFO=<ID=MULTIALLELIC,Number=0,Type=Flag,Description="Multiallelic site">')
    header.add_line(
        '##INFO=<ID=varGQ,Number=1,Type=Integer,Description="Variant genotype quality">')


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf')
    parser.add_argument('genotypes')
    parser.add_argument('varGQ')
    parser.add_argument('fout')
    args = parser.parse_args()

    vcf = pysam.VariantFile(args.vcf)
    update_vcf_header(vcf.header)

    # for whatever reason updating the records works fine but writing them
    # out through the pysam API produces a segfault, so just do it as strings
    fout = open(args.fout, 'w')
    fout.write(str(vcf.header))
    #  fout = pysam.VariantFile(args.fout, 'w', header=vcf.header)

    names = ('key name sample RD_CN RD_GQ PE_GT PE_GQ SR_GT SR_GQ GT GQ EV').split()
    genotypes = pd.read_table(args.genotypes, names=names, dtype={
        'sample': str}, sep='\s+')
    genotypes = genotypes.set_index('name sample'.split())

    # Catch empty shard case
    if genotypes.shape[0] == 0:
        return

    # sort genotype table for faster indexing
    ids = [record.id for record in vcf]
    vcf.reset()
    samples = list(vcf.header.samples)
    idx = pd.MultiIndex.from_product(iterables=[ids, samples],
                                     names=['name', 'sample'])
    genotypes = genotypes.reindex(idx)

    # TODO: catch missing samples after reindexing
    if genotypes['GT'].isnull().any():
        missing = genotypes.loc[genotypes['GT'].isnull()]
        missing = missing.reset_index()
        cols = ['name', 'sample']
        missing[cols].to_csv('missing_samples.txt', index=False, sep='\t')
        raise Exception("Missing sample genotypes. See missing_samples.txt")

    names = 'name RD_GQ PE_GQ SR_GQ GQ'.split()
    varGQ = pd.read_table(args.varGQ, names=names, sep='\s+')
    varGQ = varGQ.set_index('name').reindex(ids)

    n_samples = len(samples)
    for i, record in enumerate(vcf):
        try:
            start = i * n_samples
            end = (i + 1) * n_samples
            gt = genotypes.iloc[range(start, end)]
        except KeyError:
            msg = 'Record {0} not present in genotypes file'.format(record.id)
            raise Exception(msg)

        try:
            gq = varGQ.loc[record.id]
        except KeyError:
            msg = 'Record {0} not present in varGQ file'.format(record.id)
            raise Exception(msg)

        add_genotypes(record, gt, gq)
        fout.write(str(record))


if __name__ == '__main__':
    main()
