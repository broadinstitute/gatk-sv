#!/usr/bin/env python

"""
Compare two metrics files.
"""

import sys
import argparse
import math
import re
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import svtest.utils.IOUtils as iou

WIDTH = 10.24
HEIGHT_SCALE = 0.25
MAX_ROWS_PER_PLOT = 500


def main(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        prog='svtest plot-metrics',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('metrics_a', type=str)
    parser.add_argument('metrics_b', type=str)
    parser.add_argument('pdf_out', type=str)
    parser.add_argument('--a-name', type=str, default="metrics_a")
    parser.add_argument('--b-name', type=str, default="metrics_b")
    parser.add_argument('--sample-list', type=str, default=None)
    parser.add_argument('--changes-only', action='store_true',
                        help='Only plot values that are different')
    parser.add_argument('--linear', action='store_true',
                        help='Plot linear scale [default log]')
    parser.add_argument('--metrics-out', type=str,
                        help='Write plotted metrics to tsv', default=None)

    # Print help if no arguments specified
    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    # Read metric tables and join
    df_a = get_metrics(args.metrics_a)
    df_b = get_metrics(args.metrics_b)
    df = df_a.join(df_b, how='outer', lsuffix='_a', rsuffix='_b', sort=True)\
        .rename(columns={"value_a": args.a_name, "value_b": args.b_name})\
        .sort_values(by='name')

    # If sample ids are provided, consolidate sample-specific metrics
    if args.sample_list is not None:
        samples = iou.read_samples_list(args.sample_list)
        df = consolidate_sample_metrics(df, samples)

    # Only plot changed metrics
    if args.changes_only:
        df = df[df["value_a"] != df["value_b"]]

    # Write raw data to file
    if args.metrics_out is not None:
        df.to_csv(args.metrics_out, sep='\t')

    # Plot
    plot_data(df, args.pdf_out, args.linear)


def consolidate_sample_metrics(df, samples):
    sample_metric_rows = get_sample_metric_rows(df, samples)
    sample_df = df.loc[sample_metric_rows]
    df = df.loc[~sample_metric_rows]
    tags = get_sample_tag_groups(sample_df, samples)
    tags_df = pd.DataFrame(tags, index=sample_df.index, columns=["group"])
    sample_df = sample_df.join(tags_df)
    grouped_metrics = {}
    for tag in set(tags):
        compute_metrics(grouped_metrics, sample_df, tag)
    grouped_df = pd.DataFrame(grouped_metrics).T.sort_index()
    return df.append(grouped_df)


def compute_metrics(metrics_dict, metrics_df, group):
    metrics_dict[group + "<MEAN>"] = metrics_df[metrics_df["group"] == group].mean(axis=0)


def get_sample_metric_rows(df, samples):
    pat = '|'.join(samples)
    return df.index.str.contains(pat)


def get_sample_tag_groups(sample_df, samples):
    pat = r'({})'.format('|'.join(samples))
    pat_double_underscore = r'__'
    pat_underscore = r'(__|^_|_$)'
    tags = []
    for row in sample_df.index:
        row2 = re.sub(pat, "", row)
        row3 = re.sub(pat_double_underscore, "_", row2)
        row4 = re.sub(pat_underscore, "", row3)
        tags.append(row4)
    return tags


def get_metrics(path):
    df = pd.read_csv(path, sep='\t', names=["name", "value"])
    return df.set_index("name")


def plot_data(df, out_path, linear):
    if df.size == 0:
        plot_empty_data()
    else:
        plot_nonempty_data(df, out_path, linear)


def plot_empty_data():
    fig = plt.figure()
    ax = fig.add_axes([0, 0, 1, 1])
    left, width = .25, .5
    bottom, height = .25, .5
    right = left + width
    top = bottom + height
    ax.text(0.5 * (left + right), 0.5 * (bottom + top), 'No data',
            horizontalalignment='center',
            verticalalignment='center',
            fontsize=20, color='red',
            transform=ax.transAxes)


def plot_nonempty_data(df, out_path, linear):
    num_rows = df.index.size
    num_plots = math.ceil(num_rows / float(MAX_ROWS_PER_PLOT))
    with PdfPages(out_path) as pdf:
        for i in range(num_plots):
            start = i * MAX_ROWS_PER_PLOT
            end = min((i + 1) * MAX_ROWS_PER_PLOT + 1, num_rows)
            plot_rows(df, start, end, pdf, linear)


def plot_rows(df, start, end, pdf, linear):
    df.iloc[start:end].iloc[::-1].plot.barh(figsize=(WIDTH, HEIGHT_SCALE * (max(end - start, 15))))
    if not linear:
        plt.xscale('log')
    plt.legend(bbox_to_anchor=(0, 1.02, 1., 0.102), loc='lower left', ncol=2, borderaxespad=0.)
    plt.tight_layout()
    pdf.savefig()
    plt.close()


if __name__ == '__main__':
    main()
