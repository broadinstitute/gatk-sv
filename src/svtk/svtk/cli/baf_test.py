#!usr/bin/env python
from svtk.baf.BAFpysam import *
import argparse
from collections import deque
import numpy as np
import pandas as pd
import pysam
import sys

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
        return bafs, bafs

    bafs.columns = ['chr', 'pos', 'baf', 'sample']
    bafs = bafs[bafs['sample'].isin(samples)]

    if bafs.empty:
        return bafs, bafs

    bafs['pos'] = bafs.pos.astype(int)
    bafs['baf'] = bafs.baf.astype(float)

    bafs.loc[bafs.pos <= start, 'region'] = 'before'
    bafs.loc[bafs.pos >= end, 'region'] = 'after'
    bafs.loc[(bafs.pos > start) & (bafs.pos < end), 'region'] = 'inside'

    het_counts = bafs.groupby(['sample', 'region']).size().unstack(fill_value=0)
    het_counts = het_counts.reindex(samples, fill_value=0).astype(int).reset_index()

    for colname in ['before', 'inside', 'after']:
        if colname not in het_counts.columns:
            het_counts[colname] = 0

    called_bafs = bafs.loc[bafs.region == 'inside'].copy()
    return het_counts, called_bafs

def main(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        prog='svtk baf-test',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('bed', help='BAF bed.')
    parser.add_argument('file', help='Compiled snp file')
    parser.add_argument('-b', '--batch',)
    parser.add_argument('--index', help='Tabix index for remote bed')

    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    splist = []
    with open(args.batch, 'r') as f:
        f.readline()
        for line in f:
            splist.append(line.split('\t')[0])
    total_sample = len(splist)

    if args.index is not None:
        tbx = pysam.TabixFile(args.file, index=args.index)
    else:
        if args.file.startswith('http'):
            raise Exception('Must provide tabix index with remote URL')
        tbx = pysam.TabixFile(args.file)

    np.random.seed(0)
    random_state = np.random.RandomState(0)

    with open(args.bed, 'r') as f:
        for line in f:
            if line[0] != "#":
                dat = line.rstrip().split('\t')
                chrom, start, end, id, samples, type = dat[:6]
                start, end = int(start), int(end)
                samplelist = samples.split(',')

                try:
                    het_counts, called_bafs = preprocess(
                        chrom, start, end, tbx, samples=splist)
                except ValueError:
                    het_counts, called_bafs = pd.DataFrame(), pd.DataFrame()

                if not het_counts.empty:
                    Del = DeletionTest(het_counts, samplelist,
                                       min(end - start, 1000000), random_state=random_state)
                    KS = KS2sample(called_bafs, samplelist)
                    ks, ksp = KS.test(samplelist)
                    mean, delp = Del.Ttest(samplelist)
                    statis = Del.stats(samplelist)
                    line = f"{chrom}\t{start}\t{end}\t{id}\t{samples}\t{type}\t{mean},{delp}\t{ks},{ksp}\t{statis}"
                else:
                    line = f"{chrom}\t{start}\t{end}\t{id}\t{samples}\t{type}\tNoSNP,NoSNP\tNoSNP,NoSNP\tNoSNP"

                dat = line.rstrip().split('\t')

                if len(dat) > 11:
                    dat = dat[:11]
                    dat[-1] = dat[-1][0:-len(dat[0])]
                    line = '\t'.join(dat)
                elif len(dat) < 11:
                    nsamp = len(dat[4].split(','))
                    ncontrol = total_sample - nsamp
                    dat = dat[:-1] + [f"0,0\t0,{nsamp}\t0,0,{ncontrol}"]
                    line = '\t'.join(dat)

                dat = line.split('\t')
                delp, dupp = dat[6], dat[7]

                try:
                    dp, dstat = map(float, delp.split(','))
                    derr = 'delOK'
                except:
                    dstat, dp, derr = 'NA', 'NA', delp.split(',')[1]

                try:
                    up, ustat = map(float, dupp.split(','))
                    uerr = 'dupOK'
                except:
                    ustat, up, uerr = 'NA', 'NA', dupp.split(',')[1]

                dat = dat[:6] + [derr, str(dp), str(dstat), uerr, str(up), str(
                    ustat)] + dat[8].split(',') + dat[9].split(',') + dat[10].split(',')

                print('\t'.join(dat))
                sys.stdout.flush()

if __name__ == '__main__':
    main(sys.argv[1:])
