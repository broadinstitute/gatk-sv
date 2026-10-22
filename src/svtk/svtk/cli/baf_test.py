#!usr/bin/env python
from svtk.baf.BAFpysam import *
import argparse
from collections import deque
import numpy as np
import pandas as pd
import pysam
import sys
##########


def preprocess(chrom, start, end, tbx, samples, window=None, called_samples=None, outlier_sample_ids=None):
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
    if window is None:
        window = end - start

    bafs = deque()
    if window < 1000000:
        for record in tbx.fetch(chrom, max(1, start - window), end + window, parser=pysam.asTuple()):
            bafs.append(np.array(record))
    else:
        for record in tbx.fetch(chrom, max(1, start - 1000000), start, parser=pysam.asTuple()):
            bafs.append(np.array(record))
        for record in tbx.fetch(chrom, (start + end) // 2 - 500000, (start + end) // 2 + 500000, parser=pysam.asTuple()):
            bafs.append(np.array(record))
        for record in tbx.fetch(chrom, end, end + 1000000, parser=pysam.asTuple()):
            bafs.append(np.array(record))

    bafs = pd.DataFrame(np.array(bafs))
    if bafs.empty:
        return bafs, bafs, called_samples
    bafs.columns = ['chr', 'pos', 'baf', 'sample']

    if outlier_sample_ids and len(outlier_sample_ids) > 0:
        # Create non-outlier sample lists
        background_samples = list(set(samples) - set(called_samples))
        non_outlier_called = [s for s in called_samples if s not in outlier_sample_ids]
        non_outlier_background = [s for s in background_samples if s not in outlier_sample_ids]

        # Exclude outlier samples only if non-outlier samples exist
        if len(non_outlier_called) > 0:
            called_samples = non_outlier_called
        if len(non_outlier_background) > 0:
            background_samples = non_outlier_background
    
        # Prune samples list
        samples = list(set(called_samples) | set(background_samples))
    
    bafs = bafs[bafs['sample'].isin(samples)]
    
    if bafs.empty:
        return bafs, bafs, called_samples
    
    bafs['pos'] = bafs.pos.astype(int)
    bafs['baf'] = bafs.baf.astype(float)
    bafs.loc[bafs.pos <= start, 'region'] = 'before'
    bafs.loc[bafs.pos >= end, 'region'] = 'after'
    bafs.loc[(bafs.pos > start) & (bafs.pos < end), 'region'] = 'inside'
    het_counts = bafs.groupby('sample region'.split()).size().reset_index().rename(columns={0: 'count'})
    het_counts = het_counts.pivot_table(values='count', index='sample', columns='region', fill_value=0)
    het_counts = het_counts.reindex(samples).fillna(0).astype(int)
    cols = 'before inside after sample'.split()
    het_counts = het_counts.reset_index()
    for colname in cols:
        if colname not in het_counts.columns:
            het_counts[colname] = 0

    # Report BAF for variants inside CNV
    called_bafs = bafs.loc[bafs.region == 'inside'].copy()
    return het_counts, called_bafs, called_samples


def main(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        prog='svtk baf-test',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('bed', help='BAF bed.')
    parser.add_argument('file', help='Compiled snp file')
    parser.add_argument('-b', '--batch',)
    parser.add_argument('--index', help='Tabix index for remote bed')
    parser.add_argument('--outlier-sample-ids', help='Path to file containing outlier sample IDs')

    # Print help if no arguments specified
    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

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

    if args.index is not None:
        tbx = pysam.TabixFile(args.file, index=args.index)
    else:
        if args.file.startswith('http'):
            raise Exception('Must provide tabix index with remote URL')
        tbx = pysam.TabixFile(args.file)

    outlier_samples = set()
    if args.outlier_sample_ids:
        with open(args.outlier_sample_ids, 'r') as f:
            outlier_samples = set(line.strip() for line in f)

    # this is necessary to avoid stochasticity in calculation of KS statistic
    np.random.seed(0)
    random_state = np.random.RandomState(0)

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
                try:
                    het_counts, called_bafs, samplelist = preprocess(
                        chrom, start, end, tbx, samples=splist,
                        outlier_sample_ids=outlier_samples,
                        called_samples=samplelist
                    )
                except ValueError:
                    het_counts = pd.DataFrame()
                    called_bafs = pd.DataFrame()
                # Running BAF testing
                if not het_counts.empty:
                    Del = DeletionTest(het_counts, samplelist,
                                       min(end - start, 1000000), random_state=random_state)
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
