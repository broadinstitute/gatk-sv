#!/bin/python

import pandas as pd
import argparse
import numpy as np
import matplotlib.pyplot as plt
from os.path import basename, isfile, getsize
import logging

# Synopsis:
#  Generates summary statistics on Cromwell monitoring log summary table generated using the following command:
#       get_cromwell_resource_usage2.sh -u -r workflow_id > table.tsv
#  Cost estimates assume all machines are preemptible and have a fixed boot time. Resource
#  costs are estimated for requesting optimal resources (equal to the max observed) uniformly across all shards
#  ("static") and individually for each shard ("dynamic").
#
# Usage:
#   python analyze_monitoring_logs2.py /path/to/log_summary_table /path/to/output_base [optional parameters]
#
# Required parameters:
#   /path/to/logs : Path containing monitoring script log summary TSV from get_cromwell_resource_usage2.sh -u -r
#   /path/to/output_base : Base output path, to which extensions will be appended for each output file
# Optional parameters:
#   --overhead: Localization overhead in minutes
#   --semilog: Plot semilog y
#   --plot-norm: Specify number of samples to normalize plots to per sample
#   --log-level LEVEL (specify level of logging information to print, ie. INFO, WARNING, ERROR - not case-sensitive)
#
# Author: Emma Pierce-Hoffman (epierceh@broadinstitute.org)
# Modified from analyze_monitoring_logs.py by Mark Walker

COST_PER_GB_MEM_HR = 0.000892
COST_CPU_HR = 0.006655
COST_PER_GB_DISK_HR = 0.00005555555

MIN_CPU = 1
MIN_MEM_GB = 0.9
MIN_DISK_GB = 1

BOOT_DISK_GB = 10
DEFAULT_OVERHEAD_MIN = 5.


def check_table_columns(columns):
    column_set = set(columns)
    required_input_columns = ['ElapsedTime', 'nCPU', 'CPU', 'TotMem', 'Mem', 'MemPct', 'TotDisk', 'Disk',
                              'DiskPct', 'task']
    missing_cols = []
    missing = False
    for col in required_input_columns:
        if col not in column_set:
            missing = True
            missing_cols.append(col)
    if missing:
        raise RuntimeError(
            "Malformed input table; missing column(s): %s. Use TSV from get_cromwell_resource_usage2.sh -u -r" % ", ".join(missing_cols))


def load_data(log_file, overhead_mins):
    # columns in input:
    # ['ElapsedTime', 'nCPU', 'CPU', 'TotMem', 'Mem', 'MemPct', 'TotDisk', 'Disk', 'DiskPct', 'IORead', 'IOWrite', 'task']
    data = pd.read_table(
        log_file, usecols=lambda x: x not in ('IORead', 'IOWrite'))
    check_table_columns(data.columns)
    # rename some columns for consistency, clarity
    data.rename({'task': 'Task', 'ElapsedTime': 'Hours', 'Mem': 'MaxMem', 'CPU': 'PctCPU',
                 'Disk': 'MaxDisk', 'MemPct': 'PctMem', 'DiskPct': 'PctDisk'}, axis='columns', inplace=True)
    # add MaxCPU column
    data['MaxCPU'] = (data['PctCPU'] / 100) * data['nCPU']
    # reorder so Task column is first, MaxCPU is after nCPU (without assuming input order of columns)
    cols = data.columns.tolist()
    cols = [col for col in cols if col not in ('Task', 'MaxCPU')]
    cpu_ind = cols.index('nCPU')
    cols = ['Task'] + cols[:cpu_ind + 1] + ['MaxCPU'] + cols[cpu_ind + 1:]
    data = data[cols]
    # modify formats
    data['Hours'] = pd.to_timedelta(data['Hours']).dt.total_seconds(
    ) / 3600.0  # convert ElapsedTime to hours (float)
    data['Hours'] += overhead_mins / 60.0
    # keep last (most specific) task name, attempt number, and shard number, if present
    data['Task'] = data['Task'].str.replace('/shard', '.shard', regex=False) \
                               .str.replace('/attempt', '.attempt', regex=False) \
                               .str.rsplit('/', n=1).str[-1]

    return data


def estimate_costs_per_task(data):
    # columns after load_data():
    # ['Hours', 'nCPU', 'MaxCPU', 'PctCPU', 'TotMem', 'MaxMem', 'PctMem', 'TotDisk', 'MaxDisk', 'PctDisk', 'Task']
    # compute resource-hours : actual and with optimal settings based on maximum usage
    data['TotCPUHour'] = data['nCPU'] * data['Hours']
    data['MaxCPUHour'] = data['MaxCPU'] * data['Hours']
    data['TotMemHour'] = data['TotMem'] * data['Hours']
    data['MaxMemHour'] = data['MaxMem'] * data['Hours']
    data['TotDiskHour'] = data['TotDisk'] * data['Hours']
    data['MaxDiskHour'] = data['MaxDisk'] * data['Hours']

    # compute cost estimates : actual and with optimal resource settings based on maximum usage (per-task, so dynamic)
    data['TotCPUCost'] = data['TotCPUHour'] * COST_CPU_HR
    data['OptCPUCost'] = np.multiply(
        np.fmax(data['MaxCPU'], MIN_CPU), data['Hours']) * COST_CPU_HR
    data['TotMemCost'] = data['TotMemHour'] * COST_PER_GB_MEM_HR
    data['OptMemCost'] = np.multiply(
        np.fmax(data['MaxMem'], MIN_MEM_GB), data['Hours']) * COST_PER_GB_MEM_HR
    data['TotDiskCost'] = np.multiply(
        (data['TotDisk'] + BOOT_DISK_GB), data['Hours']) * COST_PER_GB_DISK_HR
    data['OptDiskCost'] = np.multiply((np.fmax(
        data['MaxDisk'], MIN_DISK_GB) + BOOT_DISK_GB), data['Hours']) * COST_PER_GB_DISK_HR
    data['TotTaskCost'] = data['TotCPUCost'] + \
        data['TotMemCost'] + data['TotDiskCost']
    data['OptTaskCost'] = data['OptCPUCost'] + \
        data['OptMemCost'] + data['OptDiskCost']

    data.sort_values(by='TotTaskCost', inplace=True, ascending=False)
    return data


def estimate_costs_per_group(data):
    # remove shard number, attempt number if present
    data['TaskGroup'] = data['Task'].str.split('.').str[0]
    groups = data['TaskGroup'].unique()
    data_grouped = pd.DataFrame(columns=['Task', 'Hours', 'AvgCPU', 'MaxCPU', 'PctCPU', 'AvgMem', 'MaxMem', 'PctMem',
                                         'AvgDisk', 'MaxDisk', 'PctDisk', 'TotCPUHour', 'PeakCPUHour', 'TotMemHour', 'PeakMemHour',
                                         'TotDiskHour', 'PeakDiskHour', 'TotCPUCost', 'StaticCPUCost', 'DynCPUCost', 'TotMemCost',
                                         'StaticMemCost', 'DynMemCost', 'TotDiskCost', 'StaticDiskCost', 'DynDiskCost', 'TotCost',
                                         'StaticCost', 'DynCost'])

    for group in groups:
        """
        columns of d: ['Task', 'Hours', 'nCPU', 'MaxCPU', 'PctCPU', 'TotMem', 'MaxMem',
           'PctMem', 'TotDisk', 'MaxDisk', 'PctDisk', 'TotCPUHour', 'MaxCPUHour',
           'TotMemHour', 'MaxMemHour', 'TotDiskHour', 'MaxDiskHour', 'TotCPUCost',
           'OptCPUCost', 'TotMemCost', 'OptMemCost', 'TotDiskCost', 'OptDiskCost',
           'TotTaskCost', 'OptTaskCost']
        """
        d = data.loc[data['TaskGroup'] == group]
        hours = np.sum(d['Hours'])
        max_cpu = np.nan if np.isnan(
            d['MaxCPU']).all() else np.max(d['MaxCPU'])
        max_mem = np.nan if np.isnan(
            d['MaxMem']).all() else np.max(d['MaxMem'])
        max_disk = np.nan if np.isnan(
            d['MaxDisk']).all() else np.max(d['MaxDisk'])
        group_data = {
            'Task': group,
            'Hours': hours,
            'AvgCPU': np.mean(d['nCPU']),
            'AvgMem': np.mean(d['TotMem']),
            'AvgDisk': np.mean(d['TotDisk']),
            'MaxCPU': max_cpu,
            'MaxMem': max_mem,
            'MaxDisk': max_disk,
            'PctCPU': np.nan if np.isnan(d['PctCPU']).all() else np.nanmax(d['PctCPU']),
            'PctMem': np.nan if np.isnan(d['PctMem']).all() else np.nanmax(d['PctMem']),
            'PctDisk': np.nan if np.isnan(d['PctDisk']).all() else np.nanmax(d['PctDisk']),
            'TotCPUHour': np.sum(d['TotCPUHour']),
            'TotMemHour': np.sum(d['TotMemHour']),
            'TotDiskHour': np.sum(d['TotDiskHour']),
            'PeakCPUHour': np.nan if np.isnan(d['MaxCPUHour']).all() else np.nanmax(d['MaxCPUHour']),
            'PeakMemHour': np.nan if np.isnan(d['MaxMemHour']).all() else np.nanmax(d['MaxMemHour']),
            'PeakDiskHour': np.nan if np.isnan(d['MaxDiskHour']).all() else np.nanmax(d['MaxDiskHour']),
            'TotCPUCost': np.sum(d['TotCPUCost']),
            'TotMemCost': np.sum(d['TotMemCost']),
            'TotDiskCost': np.sum(d['TotDiskCost']),
            'DynCPUCost': np.sum(d['OptCPUCost']),
            'DynMemCost': np.sum(d['OptMemCost']),
            'DynDiskCost': np.sum(d['OptDiskCost']),
            'StaticCPUCost': COST_CPU_HR * np.nanmax((max_cpu, MIN_CPU)) * hours,
            'StaticMemCost': COST_PER_GB_MEM_HR * np.nanmax((max_mem, MIN_MEM_GB)) * hours,
            'StaticDiskCost': COST_PER_GB_DISK_HR * (np.nanmax((max_disk, MIN_DISK_GB)) + BOOT_DISK_GB) * hours
        }
        group_data['TotCost'] = sum(
            (group_data['TotCPUCost'], group_data['TotMemCost'], group_data['TotDiskCost']))
        group_data['StaticCost'] = sum(
            (group_data['StaticCPUCost'], group_data['StaticMemCost'], group_data['StaticDiskCost']))
        group_data['DynCost'] = sum(
            (group_data['DynCPUCost'], group_data['DynMemCost'], group_data['DynDiskCost']))

        data_grouped = data_grouped.append(group_data, ignore_index=True)

    data_grouped.sort_values(by='TotCost', inplace=True, ascending=False)
    return data_grouped


def get_out_file_path(output_base, output_end):
    sep = "."
    if basename(output_base) == "":
        sep = ""
    out_file = output_base + sep + output_end
    return out_file


def write_data(data, out_file):
    logging.info("Writing %s" % out_file)
    data.to_csv(out_file, sep='\t', na_rep='NaN', index=False)


def do_simple_bar(data, xticks, path, bar_width=0.35, height=12, width=12,
                  xtitle='', ytitle='', title='', bottom_adjust=0, legend=[],
                  yscale='linear', sort_values=None):
    num_groups = max([d.shape[0] for d in data])
    if sort_values is not None:
        sort_indexes = np.flip(np.argsort(sort_values))
    else:
        sort_indexes = np.arange(num_groups)
    plt.figure(num=None, figsize=(width, height),
               dpi=100, facecolor='w', edgecolor='k')
    for i in range(len(data)):
        if i < len(legend):
            label = legend[i]
        else:
            label = "data" + str(i)
        x = (np.arange(num_groups) * len(data) + i) * bar_width
        plt.bar(x, data[i][sort_indexes], label=label)
    x = (np.arange(num_groups) * len(data)) * bar_width
    plt.xticks(x, [xticks[i] for i in sort_indexes], rotation='vertical')
    plt.xlabel(xtitle)
    plt.ylabel(ytitle)
    plt.title(title)
    plt.subplots_adjust(bottom=bottom_adjust)
    plt.yscale(yscale)
    plt.legend()
    plt.savefig(path)


def create_graphs(data, out_file, semilog=False, num_samples=None):
    logging.info("Writing %s" % out_file)
    # drop rows with any NA values before making plot
    data = data.loc[data.notna().all(axis=1)]
    data.reset_index(drop=True, inplace=True)
    if num_samples is not None:
        data = data / num_samples
        ytitle = "Cost ($/sample)"
        title = "Estimated Cost Per Sample"
    else:
        ytitle = "Cost ($)"
        title = "Estimated Total Cost"

    if semilog:
        yscale = "log"
    else:
        yscale = "linear"

    do_simple_bar(data=[data["TotCost"], data["StaticCost"], data["DynCost"]],
                  xticks=data['Task'],
                  path=out_file,
                  bar_width=1,
                  height=8,
                  width=12,
                  xtitle="Task",
                  ytitle=ytitle,
                  title=title,
                  bottom_adjust=0.35,
                  legend=["Current", "Uniform", "Dynamic"],
                  yscale=yscale,
                  sort_values=data["TotCost"])


def check_file_nonempty(f):
    if not isfile(f):
        raise RuntimeError("Required input file %s does not exist." % f)
    elif getsize(f) == 0:
        raise RuntimeError("Required input file %s is empty." % f)


# Main function
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "log_summary_file", help="Path to log summary TSV from get_cromwell_resource_usage2.sh -u -r")
    parser.add_argument("output_base", help="Output tsv file base path")
    parser.add_argument("--overhead", help="Localization overhead in minutes")
    parser.add_argument("--semilog", help="Plot semilog y",
                        action="store_true")
    parser.add_argument(
        "--plot-norm", help="Specify number of samples to normalize plots to per sample")
    parser.add_argument("--log-level",
                        help="Specify level of logging information, ie. info, warning, error (not case-sensitive)",
                        required=False, default="INFO")
    args = parser.parse_args()

    if not args.overhead:
        overhead = DEFAULT_OVERHEAD_MIN
    else:
        overhead = float(args.overhead)

    if args.plot_norm:
        plot_norm = int(args.plot_norm)
    else:
        plot_norm = None

    log_level = args.log_level
    numeric_level = getattr(logging, log_level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % log_level)
    logging.basicConfig(level=numeric_level,
                        format='%(levelname)s: %(message)s')

    log_file, output_base = args.log_summary_file, args.output_base
    check_file_nonempty(log_file)

    data = load_data(log_file, overhead)
    data = estimate_costs_per_task(data)
    write_data(data, get_out_file_path(output_base, "all.tsv"))
    grouped_data = estimate_costs_per_group(data)
    write_data(grouped_data, get_out_file_path(output_base, "grouped.tsv"))
    create_graphs(grouped_data, get_out_file_path(
        output_base, "cost.png"), semilog=args.semilog, num_samples=plot_norm)


if __name__ == "__main__":
    main()
