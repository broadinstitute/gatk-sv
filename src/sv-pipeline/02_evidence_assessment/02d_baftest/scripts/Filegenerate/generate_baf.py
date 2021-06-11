#!/usr/bin/env python

import argparse
import numpy as np
import pysam
import sys


def filter_record(record, unfiltered=False):
    """
    Filter VCF records to those informative for BAF genotyping.

    Returns only records which match all of the following criteria:
    1) Biallelic
    2) SNP
    3) FILTER == PASS (if filtered)

    Parameters
    ----------
    records : iterator of pysam.VariantRecords

    Returns
    ------
    record : true if filtered
    """

    # for record in records:
    # Restrict to biallelic sites
    if len(record.alleles) > 2:
        return True

    # If filtered, restrict to variants which PASS
    if (not unfiltered) and (record.filter.keys() != ['PASS']):
        return True

    # Restrict to SNPs
    ref, alt = record.alleles
    if len(ref) > 1 or len(alt) > 1:
        return True

    return False


def calc_BAF(record, samples=None):
    """

    Parameters
    ----------
    record : pysam.VariantRecord
    samples : list of str, optional
        Subset of samples in record to consider

    Returns
    -------
    bafs : np.ndarray of np.float
        BAF at site for each sample
    """

    def _is_het(sample):
        return record.samples[sample]['GT'] == (0, 1)

    def _calc_BAF(sample):
        if not _is_het(sample):
            return np.nan

        DP = record.samples[sample]['DP']
        AD = record.samples[sample]['AD']

        if (DP is not None) and (DP > 10):  # SNP sites with >10 DP are included in BAF profile
            return AD[0] / float(DP)
        else:
            return np.nan

    if samples is None:
        samples = record.samples.keys()

    bafs = np.array([_calc_BAF(sample) for sample in samples], dtype=np.float)

    return bafs, samples


def normalize_bafs(bafs, samples, max_std=0.2):
    """
    Normalize BAFs and exclude outlying sites
    Normalize so per variant median BAF==0.5. Ignore sites with more than 0.2 standard deviation across samples.

    Parameters
    ----------
    bafs : np.ndarray (n_sites x n_samples)
    max_std : float, optional
        Maximium standard deviation permitted at a site

    Returns
    -------
    normalized_bafs : np.ndarray
    """
    # Center each site's median BAF at 0.5
    nan_bafs = np.isnan(bafs)
    bafs = bafs[~nan_bafs]
    if len(bafs) == 0 or bafs.std(ddof=1) > max_std:
        return [], []
    samples = np.asarray(samples)[~nan_bafs]
    bafs = bafs - np.median(bafs) + 0.5
    return bafs, samples


def read_samples_list(path):
    with open(path, 'r') as f:
        samples = [line for line in f.read().splitlines() if line]
    if len(samples) == 0:
        raise ValueError("Samples list empty")
    return samples


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--samples-list', default=None)
    parser.add_argument('--unfiltered', action='store_true')
    args = parser.parse_args()
    vcf = pysam.VariantFile(sys.stdin)
    if args.samples_list is not None:
        samples_list = read_samples_list(args.samples_list)
    else:
        samples_list = None
    # While loop to iterate over all records, then break if reach the end
    for record in vcf:
        if not filter_record(record, args.unfiltered):
            baf, record_samples = calc_BAF(record, samples=samples_list)
            baf, record_samples = normalize_bafs(baf, record_samples)
            for i in range(len(record_samples)):
                print(str(record.chrom) + "\t" + str(record.pos) + "\t" +
                      "{0:.2f}".format(baf[i]) + "\t" + record_samples[i])


if __name__ == '__main__':
    main()
