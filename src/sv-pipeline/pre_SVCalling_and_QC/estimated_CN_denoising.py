#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def calculate_median_mad(group):
    median = np.median(group)
    mad = np.median(np.abs(group - median))
    return median, mad


def main():
    parser = argparse.ArgumentParser(description="Generate MAD/median stats and plots from copy number data.")
    parser.add_argument('--binwise-copy-number', required=True, help="TSV file with binwise copy number (wide format)")
    parser.add_argument('--estimated-copy-number', required=True, help="TSV file with estimated copy number per sample (wide format)")
    parser.add_argument('--output-stats', required=True, help="Output TSV file with median and MAD per sample and chromosome")
    parser.add_argument('--output-pdf', required=True, help="Output PDF file with error bar plots")

    args = parser.parse_args()

    binwise_cn = pd.read_csv(args.binwise_copy_number, sep='\t')
    estimated_cn = pd.read_csv(args.estimated_copy_number, sep='\t')

    binwise_melted = pd.melt(
        binwise_cn,
        id_vars=["#Chr", "Start", "End"],
        var_name="SampleID",
        value_name="BinwiseCopyNumber"
    )

    stats_results = binwise_melted.groupby(["#Chr", "SampleID"]).agg(
        BinwiseCopyNumber_median=('BinwiseCopyNumber', lambda x: calculate_median_mad(x)[0]),
        BinwiseCopyNumber_mad=('BinwiseCopyNumber', lambda x: calculate_median_mad(x)[1])
    ).reset_index()

    cn_melted = pd.melt(
        estimated_cn,
        id_vars=["sample_id"],
        var_name="Chr",
        value_name="CopyNumber"
    )
    cn_melted.rename(columns={"sample_id": "SampleID"}, inplace=True)
    cn_melted['Chr'] = cn_melted['Chr'].str.replace('_CopyNumber$', '', regex=True)

    merged_df = pd.merge(cn_melted, stats_results, how="left", left_on=["Chr", "SampleID"], right_on=["#Chr", "SampleID"])

    merged_df.to_csv(args.output_stats, sep='\t', index=False)

    unique_chromosomes = merged_df['Chr'].unique()
    with PdfPages(args.output_pdf) as pdf:
        for chrom in unique_chromosomes:
            chrom_data = merged_df[merged_df['Chr'] == chrom].copy()
            chrom_data = chrom_data.sort_values('CopyNumber')

            plt.figure(figsize=(14, 6))
            plt.errorbar(
                x=chrom_data['SampleID'],
                y=chrom_data['CopyNumber'],
                yerr=chrom_data['BinwiseCopyNumber_mad'],
                fmt='o',
                markerfacecolor='blue',
                markeredgecolor='blue',
                ecolor='black',
                capsize=2,
                linestyle='none'
            )

            plt.xlabel('SampleID')
            plt.ylabel('Copy Number')
            plt.title(f'Copy Number for {chrom} with MAD as Error Bars')
            plt.xticks(rotation=90, ticks=range(len(chrom_data['SampleID'])), labels=chrom_data['SampleID'])
            plt.tight_layout()
            pdf.savefig()
            plt.close()


if __name__ == "__main__":
    main()
