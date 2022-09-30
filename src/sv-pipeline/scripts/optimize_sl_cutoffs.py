#!/bin/env python

import argparse
import os
import sys
from typing import Any, List, Text, Set, Dict, Optional

from matplotlib import pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn import metrics


def load_df(path):
    df = pd.read_csv(path, sep='\t', compression='gzip')
    df['LABEL'] = df['LABEL'].astype(int)
    for k in ['SVLEN', 'NON_REF_GENOTYPE_CONCORDANCE', 'GQ', 'OGQ', 'RD_GQ', 'PE_GQ', 'SR_GQ', 'AF']:
        if k not in df.columns:
            continue
        df.loc[df[k] == 'None', k] = 'nan'
        df[k] = df[k].astype(float)
    df['VID_SAMPLE'] = df['VID'] + "_" + df['SAMPLE']
    df = df.set_index('VID_SAMPLE')
    return df


def plot_hist(data, x, out_dir, out_name, hue=None, row=None, col=None, bins=20, stat='count'):
    sns.set_palette("muted")
    g = sns.displot(
        data, x=x,
        hue=hue,
        row=row,
        col=col,
        stat=stat,
        multiple="stack",
        edgecolor=".3",
        linewidth=.5,
        height=3, aspect=1.7,
        bins=bins
    )
    if str(data[x].dtype) == 'object':
        rotation = 90
        for i, ax in enumerate(g.fig.axes):
            ticks_loc = ax.get_xticks()
            ax.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
            ax.set_xticklabels([x for x in ax.get_xticklabels()], rotation = rotation)
    plt.savefig(fname=os.path.join(out_dir, f"{out_name}.histogram.{x}.png"), format="png",
                dpi=300, bbox_inches="tight")


# General QC:
# Check that there are sufficient labeled data across SV classes and sizes
# and that the SL scores are sensible (more negative for FP, more positive for TP)
def generate_qc_plots(df, out_dir, out_name):

    # Copy and perform any modifications to the dataset here
    df['LOG_SVLEN'] = np.log10(df['SVLEN'])
    df['TRUTH_LABEL'] = 'UNLABELED'
    df.loc[df['LABEL'] == 0, 'TRUTH_LABEL'] = 'FP'
    df.loc[df['LABEL'] == 1, 'TRUTH_LABEL'] = 'TP'
    df['TRUTH_LABEL'] = pd.Categorical(df['TRUTH_LABEL'], ['TP', 'FP', 'UNLABELED'])

    n_total = df.shape[0]
    n_true = df[df['LABEL'] == 1].shape[0]
    n_false = df[df['LABEL'] == 0].shape[0]
    n_unl = df[df['LABEL'] == -1].shape[0]
    print(f"True: {n_true} ({100*n_true/float(n_total)}%)")
    print(f"False: {n_false} ({100*n_false/float(n_total)}%)")
    print(f"Unlabeled: {n_unl} ({100*n_unl/float(n_total)}%)")
    print(f"Total: {n_total}")

    # Plot
    plot_hist(df, 'NON_REF_GENOTYPE_CONCORDANCE', hue='TRUTH_LABEL', bins=30, out_dir=out_dir, out_name=out_name)
    plot_hist(df, 'TRUTH_LABEL', hue='TRUTH_LABEL', out_dir=out_dir, out_name=out_name)
    plot_hist(df, 'LOG_SVLEN', hue='TRUTH_LABEL', out_dir=out_dir, out_name=out_name)
    plot_hist(df, 'SVTYPE', hue='TRUTH_LABEL', out_dir=out_dir, out_name=out_name)
    plot_hist(df, 'FILTER_CLASS', hue='TRUTH_LABEL', out_dir=out_dir, out_name=out_name)
    plot_hist(df, 'AF', hue='TRUTH_LABEL', out_dir=out_dir, out_name=out_name)
    plot_hist(df, 'SL', hue='TRUTH_LABEL', out_dir=out_dir, out_name=out_name)
    plot_hist(df, 'STATUS', hue='TRUTH_LABEL', out_dir=out_dir, out_name=out_name)
    plot_hist(df, 'NON_REF_GENOTYPE_CONCORDANCE', hue='TRUTH_LABEL', out_dir=out_dir, out_name=out_name)
    plot_hist(df, 'ALGORITHMS', hue='TRUTH_LABEL', col='SVTYPE', out_dir=out_dir, out_name=out_name)
    plot_hist(df, 'OGQ', hue='TRUTH_LABEL', col='SVTYPE', out_dir=out_dir, out_name=out_name)
    plot_hist(df, 'RD_GQ', hue='TRUTH_LABEL', col='SVTYPE', out_dir=out_dir, out_name=out_name)
    plot_hist(df, 'PE_GQ', hue='TRUTH_LABEL', col='SVTYPE', out_dir=out_dir, out_name=out_name)
    plot_hist(df, 'SR_GQ', hue='TRUTH_LABEL', col='SVTYPE', out_dir=out_dir, out_name=out_name)


# Plots precision / recall curves and performs cutoff optimization based on max F-score
def plot_precision_recall(data, name, out_dir, out_name, reverse=False, beta=1, plot=True, n_samples=None):
    results = dict()
    keys = [
        'ALL',
        'DEL_s',
        'DEL_m',
        'DEL_l',
        'DUP_s',
        'DUP_m',
        'DUP_l',
        'INS',
        'INV',
        'CNV',
        'BND',
        'CTX',
        'CPX']
    for key in keys:
        if key == 'ALL':
            df2 = data[data['FILTER_CLASS'] != 'BND']
        else:
            df2 = data[(data['FILTER_CLASS'] == key)]
        y_scores = df2[name]
        not_nan = ~y_scores.isnull()
        y_scores = y_scores[not_nan]
        results[key] = dict()
        results[key]['ppv'] = []
        results[key]['tpr'] = []
        results[key]['thresholds'] = []
        if len(y_scores) <= 1:
            continue
        if reverse:
            y_scores = -y_scores
        y_labels = df2.loc[not_nan, 'LABEL']
        ppv_i, tpr_i, thresholds_i = metrics.precision_recall_curve(y_labels, y_scores, pos_label=1)

        pos_counts, pos_bins = np.histogram(y_scores, bins=100)
        if reverse:
            pos_counts = np.cumsum(pos_counts)
        else:
            pos_counts = np.sum(pos_counts) - np.cumsum(pos_counts)
        if n_samples is not None:
            pos_counts = pos_counts / n_samples

        if len(thresholds_i) <= 1:
            continue
        fbeta = (1+beta*beta)*tpr_i*ppv_i/(tpr_i+(beta*beta*ppv_i))
        fmax = np.nanmax(fbeta)
        fmax_index = np.nanargmax(fbeta)
        fmax_threshold = thresholds_i[fmax_index]
        fmax_tpr = tpr_i[fmax_index]
        fmax_ppv = ppv_i[fmax_index]
        fmax_pos = np.interp(fmax_threshold, pos_bins[:-1], pos_counts)

        n = y_labels.shape[0]
        if reverse:
            fmax_threshold = -fmax_threshold
        results[key]['fmax'] = fmax
        results[key]['fmax_thresh'] = fmax_threshold
        results[key]['fmax_rec'] = fmax_tpr
        results[key]['fmax_prec'] = fmax_ppv
        results[key]['fmax_pos'] = fmax_pos

        results[key]['ppv'] = ppv_i
        results[key]['tpr'] = tpr_i
        results[key]['thresh'] = thresholds_i
        results[key]['pos_counts'] = pos_counts
        results[key]['pos_bins'] = pos_bins
        results[key]['n'] = n

    if not plot:
        return
    keys = [k for k in results.keys() if len(results[k]['tpr']) > 0]
    plt.figure(figsize=(16, 3))
    for k in keys:
        r = results[k]
        plt.subplot(1, 4, 1)
        plt.step(r['tpr'], r['ppv'])
        plt.xlabel('recall')
        plt.ylabel('precision')
        plt.legend(keys)
        t = np.asarray(r['thresh'])
        t_pos = np.asarray(r['pos_bins'])
        if reverse:
            t = -t
            t_pos = -t_pos
        plt.subplot(1, 4, 2)
        plt.step(t, r['tpr'][:-1], where='pre')
        plt.xlabel(name)
        plt.ylabel('recall')
        plt.subplot(1, 4, 3)
        plt.step(t, r['ppv'][:-1], where='pre')
        plt.xlabel(name)
        plt.ylabel('precision')

        plt.subplot(1, 4, 4)
        plt.step(t_pos[:-1], r['pos_counts'], where='pre')
        plt.xlabel(name)
        if n_samples is None:
            plt.ylabel('count')
        else:
            plt.ylabel('n_per_genome')

    plt.savefig(fname=os.path.join(out_dir, f"{out_name}.precision_recall.png"), format="png",
                dpi=300, bbox_inches="tight")
    return results


def write_stats(out_path, stats, stat_name):
    with open(out_path, 'w') as f:
        f.write("\t".join(['class', 'n', 'f_max', 'rec', 'prec', 'count', stat_name + '_thresh']) + "\n")
        for key in stats:
            if 'n' not in stats[key]:
                continue
            f.write("\t".join([key, str(stats[key]['n'])] +
                            ["{:.3f}".format(x) for x in [stats[key]['fmax'], stats[key]['fmax_rec'],
                                                          stats[key]['fmax_prec']]] +
                            [f"{stats[key]['fmax_pos']:.0f}", f"{stats[key]['fmax_thresh']:.0f}"]) + "\n")
        total = sum([stats[key]['fmax_pos'] for key in stats if 'fmax_pos' in stats[key] and key != 'ALL'])
        f.write(f"Total count across subclasses: {total:0.0f}\n")


# The next step is to filter genotypes using the FilterRecalibratedVcf workflow
# Use fmax thresholds to set the SL cutoffs
def write_sl_filter_args(out_path, stats):
    key_to_arg_map = {
        'DEL_s': 'small-del-threshold',
        'DEL_m': 'medium-del-threshold',
        'DEL_l': 'large-del-threshold',
        'DUP_s': 'small-dup-threshold',
        'DUP_m': 'medium-dup-threshold',
        'DUP_l': 'large-dup-threshold',
        'INS': 'ins-threshold',
        'INV': 'inv-threshold',
        'BND': 'bnd-threshold',
        'CPX': 'cpx-threshold',
        'CTX': 'ctx-threshold'
    }
    suggested_filter_args = []
    for key in stats:
        if key == 'ALL' or 'fmax_thresh' not in stats[key]:
            continue
        suggested_filter_args.append(f"--{key_to_arg_map[key]} {stats[key]['fmax_thresh']}")
    with open(out_path, 'w') as f:
        f.write(" ".join(suggested_filter_args) + "\n")


def _parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Optimizes cutoffs for SL filtering",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--table', type=str, help='Gzipped input vcf with SL annotations')
    parser.add_argument('--out-dir', type=str, help='Output directory', default="./")
    parser.add_argument('--out-name', type=str, help='Output filename base', default="filter_qc")
    parser.add_argument("--beta", type=float, default=1.0,
                        help="Beta parameter for F score (higher values will weight recall more over precision, "
                             "see https://en.wikipedia.org/wiki/F-score)")
    if len(argv) <= 1:
        parser.parse_args(["--help"])
        sys.exit(0)
    parsed_arguments = parser.parse_args(argv[1:])
    return parsed_arguments


def main(argv: Optional[List[Text]] = None):
    if argv is None:
        argv = sys.argv
    args = _parse_arguments(argv)

    # Read table into memory
    print("Loading data table...")
    df = load_df(args.table)

    # QC plotting
    print("Creating QC plots...")
    generate_qc_plots(df=df, out_dir=args.out_dir, out_name=args.out_name)

    # Perform filtering optimizations on SL score
    # Use plots and/or fmax to choose SL thresholds for each class

    # Remove unlabeled calls
    plot_df = df[df['LABEL'] != -1]

    # Calculate fmax and produce precision-recall plots across SL cutoffs
    stat_name = "SL"
    print("Optimizing cutoffs...")
    stats = plot_precision_recall(plot_df, stat_name, out_dir=args.out_dir, out_name=args.out_name, beta=args.beta)
    print("Writing stats file...")
    write_stats(out_path=os.path.join(args.out_dir, f"{args.out_name}.stats.txt"), stats=stats, stat_name=stat_name)
    print("Writing arguments file...")
    write_sl_filter_args(out_path=os.path.join(args.out_dir, f"{args.out_name}.filter_args.txt"), stats=stats)
    print("Done!")


if __name__ == "__main__":
    main()
