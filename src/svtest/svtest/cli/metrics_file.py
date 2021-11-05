#!/usr/bin/env python

"""
Collect metric file metrics. Writes stats to stdout.

Metrics:
  metrics_<common>_num_records   : Number of records
  metrics_<common>_mean_<metric> : Mean metric value
  metrics_<common>_num_empty_<metric> : Number of metric N/A entries

"""

import argparse
import math
import sys

from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
import pandas as pd
import svtest.utils.TestUtils as tu
import svtest.utils.IOUtils as iou

KEY_PREFIX = "metrics_"

EXPECTED_COLUMNS = ["name", "chrom", "svtype", "svsize", "vf", "poor_region_cov", "rmsk", "is_outlier_specific",
                    "PE_log_pval", "PE_called_median", "PE_bg_median", "PE_bg_frac", "SR_posA_log_pval",
                    "SR_posB_log_pval", "SR_sum_log_pval", "SR_posA_called_median", "SR_posB_called_median",
                    "SR_sum_called_median", "SR_posA_bg_median", "SR_posB_bg_median", "SR_sum_bg_median",
                    "SR_posA_bg_frac", "SR_posB_bg_frac", "SR_sum_bg_frac", "SR_posA_pos", "SR_posB_pos",
                    "PESR_log_pval", "PESR_called_median", "PESR_bg_median", "PESR_bg_frac", "RD_Median_Power",
                    "RD_P", "RD_2ndMaxP", "RD_Model", "RD_Median_Rank", "RD_Median_Separation", "RD_log_pval",
                    "RD_log_2ndMaxP", "BAF_delstat", "BAF_snp_ratio", "BAF_del_loglik", "BAF_dupstat", "BAF_KS_stat",
                    "BAF_KS_log_pval", "BAF_total_case_snps", "BAF_total_snps", "BAF_n_nonROH_cases", "BAF_n_samples",
                    "BAF_mean_control_snps", "BAF_n_nonROH_controls", "BAF_n_controls"]

EXPECTED_TYPES = ["DEL", "DUP", "INS", "INV", "BND"]

# Plotting config
WIDTH = 10.24
HEIGHT = 16
MAX_ROWS_PER_PLOT = 7
BINS = 30


def main(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        prog='svtest metrics-file',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('metrics_file', type=str)
    parser.add_argument('contig_list', type=str)
    parser.add_argument('--common', action='store_true')
    parser.add_argument('--plot', type=str, help="If provided, plots histogram to specified pdf path")

    # Print help if no arguments specified
    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    contigs = iou.read_contig_list(args.contig_list)
    tu.test_is_not_empty(contigs, "contigs")

    # Columns used as features for adjudication in module 03
    feature_cols = ["poor_region_cov", "is_outlier_specific",
                    "BAF_snp_ratio", "BAF_del_loglik", "BAF_KS_stat", "BAF_KS_log_pval",
                    "SR_sum_log_pval", "SR_sum_bg_frac", "RD_Median_Separation", "RD_log_pval",
                    "RD_log_2ndMaxP", "PE_log_pval", "PE_bg_frac"]

    # Read file
    df = pd.read_csv(args.metrics_file, sep='\t')
    metrics = get_metrics(df, contigs, feature_cols, args.common)

    # Generate histograms
    if args.plot is not None:
        plot_nonempty_data(df, feature_cols, args.plot)

    # Write metrics
    write_metrics(metrics)


def write_metrics(metrics):
    for key in metrics:
        sys.stdout.write("%s\t%s\n" % (key, str(metrics[key])))


def get_metrics(df, contigs, feature_cols, common):
    tu.test_sets_equal(df["chrom"], contigs, item_str="contig",
                       name_a="metric file contigs", name_b="contigs list")
    tu.test_sets_equal(df.columns, EXPECTED_COLUMNS, item_str="column",
                       name_a="metric file header", name_b="expected columns")
    metric_means = get_column_means(df, feature_cols)
    metric_empty_counts = get_columns_num_empty(df, feature_cols)

    if common:
        prefix = "{}common_".format(KEY_PREFIX)
    else:
        prefix = KEY_PREFIX
    metrics = {
        prefix + "num_records": df.size
    }
    for key in metric_means:
        metrics[prefix + "mean_" + key] = metric_means[key]
    for key in metric_empty_counts:
        metrics[prefix + "num_empty_" + key] = metric_empty_counts[key]
    return metrics


def get_columns_num_empty(df, columns):
    counts = {}
    for col in columns:
        counts[col] = get_column_num_empty(df, col)
    return counts


def get_column_num_empty(df, column):
    col = df[column]
    return col[col.isna()].size


def get_column_means(df, columns):
    means = {}
    for col in columns:
        means[col] = get_column_mean(df, col)
    return means


def get_column_mean(df, column):
    return df[column].mean(axis=0)


def plot_rows(df, cols, pdf):
    for i in range(len(cols)):
        plt.subplot(len(cols), 1, i + 1)
        df[cols[i]].astype('float64').plot.hist(figsize=(WIDTH, HEIGHT), bins=BINS)
        plt.xlabel(cols[i])
    plt.tight_layout()
    pdf.savefig()
    plt.close()


def plot_nonempty_data(df, cols, out_path):
    num_plots = math.ceil(len(cols) / float(MAX_ROWS_PER_PLOT))
    with PdfPages(out_path) as pdf:
        for i in range(num_plots):
            start = i * MAX_ROWS_PER_PLOT
            end = min((i + 1) * MAX_ROWS_PER_PLOT + 1, len(cols))
            plot_rows(df, cols[start:end], pdf)


if __name__ == '__main__':
    main()
