#!usr/bin/env python
from BAFpysam import *
import argparse
from collections import deque
import numpy as np
import pandas as pd
import pysam
import sys
##########


def preprocess(chrom, start, end, tbx, samples, window=None):
    """
    Report normalized BAFs in a set of samples across a desired region.

    Parameters
    ----------
    chrom : str
    start : int
    end : int
    vcf : pysam.VariantFile
        GATK VCF
    samples : list of str, optional
        Samples to consider (Defaults to all in VCF)
    window : int, optional
        Window around CNV to consider (Defaults to CNV length)

    Returns
    -------
    het_counts : pd.DataFrame
        Per-sample counts of heterozygous SNPs in windows before, inside, and
        after CNV
    called_bafs : pd.DataFrame
        BAF of each called SNP within CNV
    """

    # Load and filter SNP sites
    if window is None:
        window = end - start
    # records = tbx.fetch(chrom, start - window, end + window,parser=pysam.asTuple())
    # sites = deque()
    bafs = deque()
    if window < 1000000:
        for record in tbx.fetch(chrom, max(1, start - window), end + window):
            bafs.append(np.array(record.split('\t')))
    else:

        for record in tbx.fetch(chrom, max(1, start - 1000000), start):
            bafs.append(np.array(record.split('\t')))
        for record in tbx.fetch(chrom, (start + end) // 2 - 500000, (start + end) // 2 + 500000):
            bafs.append(np.array(record.split('\t')))
        for record in tbx.fetch(chrom, end, end + 1000000):
            bafs.append(np.array(record.split('\t')))
    bafs = np.array(bafs)
    # if bafs.shape[0] == 0:
    # return 0,0
    bafs = pd.DataFrame(bafs)
    if bafs.empty:
        return bafs, bafs
    bafs.columns = ['chr', 'pos', 'baf', 'sample']
    bafs = bafs[bafs['sample'].isin(samples)]
    # print(bafs)
    if bafs.empty:
        return bafs, bafs
    bafs['pos'] = bafs.pos.astype(int)
    bafs['baf'] = bafs.baf.astype(float)
    # print(bafqs)
    bafs.loc[bafs.pos <= start, 'region'] = 'before'
    bafs.loc[bafs.pos >= end, 'region'] = 'after'
    bafs.loc[(bafs.pos > start) & (bafs.pos < end), 'region'] = 'inside'
    het_counts = bafs.groupby('sample region'.split()).size(
    ).reset_index().rename(columns={0: 'count'})
    het_counts = het_counts.pivot_table(
        values='count', index='sample', columns='region', fill_value=0)
    het_counts = het_counts.reindex(samples).fillna(0).astype(int)
    cols = 'before inside after sample'.split()
    het_counts = het_counts.reset_index()
    for colname in cols:
        if colname not in het_counts.columns:
            het_counts[colname] = 0
    # het_counts = het_counts.reset_index()[cols]
    # Report BAF for variants inside CNV
    called_bafs = bafs.loc[bafs.region == 'inside'].copy()
    return het_counts, called_bafs


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('bed', help='GATK VCF.')
    parser.add_argument('file', help='Compiled snp file')
    parser.add_argument('-b', '--batch',)
    # help='Samples')
    args = parser.parse_args()
    # fi = args.file
    # if args.vcf.startswith('s3://'):
    # vcf_path = args.vcf[5:]
    # bucket = vcf_path.split('/')[0]
    # vcf_path = '/'.join(vcf_path.split('/')[1:])

    # vcf = load_s3vcf(bucket, vcf_path, args.tbi)

    # else:
    # vcf = pysam.VariantFile(args.vcf)
    splist = []
    with open(args.batch, 'r') as f:
        f.readline()
        for line in f:
            splist.append(line.split('\t')[0])
    total_sample = sum(1 for line in open(args.batch)) - 1
    tbx = pysam.TabixFile(args.file)
    with open(args.bed, 'r') as f:
        for line in f:
            if line[0] != "#":
                dat = line.rstrip().split('\t')
                chrom = dat[0]
                start = int(dat[1])
                end = int(dat[2])
                id = dat[3]
                samples = dat[4]
                samplelist = samples.split(',')
                type = dat[5]
                het_counts, called_bafs = preprocess(
                    chrom, start, end, tbx, samples=splist)
                # Running BAF testing
                if not het_counts.empty:
                    Del = DeletionTest(het_counts, samplelist,
                                       min(end - start, 1000000))
                    KS = KS2sample(called_bafs, samplelist)
                    ks, ksp = KS.test(samplelist)
                    mean, delp = Del.Ttest(samplelist)
                    statis = Del.stats(samplelist)
                    line = chrom + '\t' + str(start) + '\t' + str(end) + '\t' + id + '\t' + samples + '\t' + type + '\t' + str(
                        mean) + ',' + str(delp) + "\t" + str(ks) + ',' + str(ksp) + '\t' + statis
                else:
                    line = chrom + '\t' + str(start) + '\t' + str(end) + '\t' + id + '\t' + samples + '\t' + \
                        type + '\t' + 'NoSNP' + ',' + 'NoSNP' + "\t" + \
                        'NoSNP' + ',' + 'NoSNP' + '\t' + 'NoSNP'
                line = line.rstrip()
                dat = line.split('\t')
                if len(dat) > 11:
                    dat = dat[0:11]
                    dat[-1] = dat[-1][0:-len(dat[0])]
                    line = '\t'.join(dat)
                if len(dat) < 11:
                    dat = line.split('\t')
                    # samp = dat[4]
                    nsamp = len(dat[4].split(','))
                    ncontrol = total_sample - nsamp
                    dat = dat[0:-1] + ['0,0\t0,' +
                                       str(nsamp) + '\t0,0,' + str(ncontrol)]
                    line = '\t'.join(dat)
                dat = line.split('\t')
                delp = dat[6]
                dupp = dat[7]
                try:
                    dp = float(delp.split(',')[0])
                    dstat = float(delp.split(',')[1])
                    derr = 'delOK'
                except:
                    dstat = 'NA'
                    dp = 'NA'
                    derr = delp.split(',')[1]
                try:
                    up = float(dupp.split(',')[0])
                    ustat = float(dupp.split(',')[1])
                    uerr = 'dupOK'
                except:
                    ustat = 'NA'
                    up = 'NA'
                    uerr = dupp.split(',')[1]
                dat = dat[0:6] + [derr, str(dp), str(dstat), uerr, str(up), str(
                    ustat)] + dat[8].split(',') + dat[9].split(',') + dat[10].split(',')
                print('\t'.join(dat))
                sys.stdout.flush()


if __name__ == '__main__':
    main()
