#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2022 The Broad Institute of M.I.T. and Harvard
# Distributed under terms of the MIT license.
# Contact: Ryan Collins <rlcollins@g.harvard.edu>

"""
Average results from multi-filter benchmarking across multiple samples
"""


import pandas as pd
import argparse


def summarize_group(subdf):
    """
    Summarize stats for a single category across all samples
    """

    # Take median across all samples
    n_samples = len(subdf)
    vals = subdf.loc[:, 'N_all N_supported pct_supported'.split()].median()
    
    # Expand filtering criteria
    out_res = pd.Series([subdf.SVTYPE.values[0], subdf.SVLEN.values[0], 
                         subdf.SVLEN.values[0]+1, subdf.AF.values[0], 
                         subdf.AF.values[0]+1, subdf.EV.values[0], 
                         subdf.COMBO.values[0], n_samples, *vals.tolist()],
                         index=['SVTYPE', 'min_log10_SVLEN', 'max_log10_SVLEN',
                                'min_log10_AF', 'max_log10_AF', 'EV', 'passing_filters',
                                'n_samples_evaluated', 'median_SVs_per_sample',
                                'median_confirmed_SVs_per_sample', 
                                'confirmation_rate'])
    
    return out_res


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('inputs', help='List of input .tsvs; one per sample.')
    parser.add_argument('output', help='Path to output file.')
    args = parser.parse_args()

    # Loads each sample's input file as a pd.DataFrame
    with open(args.inputs) as inlist:
        sample_res = [pd.read_csv(file.rstrip(), sep='\t') for file in inlist.readlines()]

    # Pool results across all samples
    cohort_res = pd.concat(sample_res, axis=0).\
                    rename(columns={'#SVTYPE' : 'SVTYPE'}).\
                    reset_index(drop=True)

    # Group results by category and take median across all samples
    avg_res = cohort_res.groupby(['SVTYPE', 'SVLEN', 'AF', 'EV', 'COMBO']).\
                         apply(summarize_group).\
                         reset_index(drop=True)

    # Write results to outfile
    avg_res.rename(columns={'SVTYPE' : '#SVTYPE'}).\
            to_csv(args.output, sep='\t', header=True, index=False)


if __name__ == '__main__':
    main()

