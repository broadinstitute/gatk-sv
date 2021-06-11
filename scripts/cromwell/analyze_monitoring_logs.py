#!/bin/python

import pandas as pd
import glob
import argparse
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt

# Synopsis:
#  Generates summary statistics on Cromwell monitoring logs collected using download_monitoring_logs.py.
#  Cost estimates assume all machines are preemptible and have a fixed bootup time. Resource
#  usage and costs are for requesting optimal resource (equal to the max observed) uniformly across all shards ("static")
#  and individually for each shard ("dynamic").
#
# Usage:
#   python analyze_monitoring_logs.py /path/to/logs /path/to/output_base
#
# Parameters:
#   /path/to/logs : Path containing monitoring script logs ending in ".monitoring.log"
#   /path/to/output_base : Base output path, to which extensions will be appended for each output file
#
# Author: Mark Walker (markw@broadinstitute.org)

TIME_FORMAT = "%a %b %d %H:%M:%S %Z %Y"
ALL_HEADER = '#job\ttask\thr\tmem_total\tmem_gb_max\tmem_pct_max\tdisk_total\tdisk_gb_max\tdisk_pct_max\tmem_gb_hr\tdisk_gb_hr\tmax_mem_gb_hr\tmax_disk_gb_hr\tcost_mem\tcost_mem_dyn\tcost_disk\tcost_disk_dyn\n'
GROUP_HEADER = '#task\thr\tmem_avg\tmem_gb_max\tmem_pct_max\tdisk_avg\tdisk_gb_max\tdisk_pct_max\tmem_gb_hr\tdisk_gb_hr\tmax_mem_gb_hr\tmax_disk_gb_hr\tcost_mem\tcost_mem_static\tcost_mem_dyn\tcost_disk\tcost_disk_static\tcost_disk_dyn\n'

COST_PER_GB_MEM_HR = 0.000892
COST_CPU_HR = 0.006655
COST_PER_GB_DISK_HR = 0.00005555555

MIN_CPU = 1
MIN_MEM_GB = 0.9
MIN_DISK_GB = 1

BOOT_DISK_GB = 10
DEFAULT_OVERHEAD_MIN = 5.


def write_data(data, file_path, header):
    with open(file_path, 'w') as f:
        f.write(header)
        for key in data.index:
            f.write(key + '\t' + '\t'.join([str(x)
                                            for x in data.loc(key)]) + '\n')


def read_data(dir, overhead_min=0):
    data = {}
    for filepath in glob.glob(dir + '/*.monitoring.log'):
        with open(filepath, 'r') as f:
            mem_gb_data_f = []
            disk_gb_data_f = []
            mem_pct_data_f = []
            disk_pct_data_f = []
            cpu_pct_data_f = []
            total_mem = 0
            total_disk = 0
            total_cpu = 0
            start_time = None
            end_time = None
            for line in f:
                tokens = line.strip().split(' ')
                if start_time is None and line.startswith('['):
                    start_time = datetime.strptime(
                        line.strip()[1:-1], TIME_FORMAT)
                if line.startswith('['):
                    end_time = datetime.strptime(
                        line.strip()[1:-1], TIME_FORMAT)
                if line.startswith('Total Memory:'):
                    total_mem = float(tokens[2])
                elif line.startswith('#CPU:'):
                    total_cpu = float(tokens[1])
                elif line.startswith('Total Disk space:'):
                    total_disk = float(tokens[3])
                elif line.startswith('* Memory usage:'):
                    mem_gb = float(tokens[3])
                    mem_pct = float(tokens[5][:-1]) / 100.0
                    mem_gb_data_f.append(mem_gb)
                    mem_pct_data_f.append(mem_pct)
                elif line.startswith('* Disk usage:'):
                    disk_gb = float(tokens[3])
                    disk_pct = float(tokens[5][:-1]) / 100.0
                    disk_gb_data_f.append(disk_gb)
                    disk_pct_data_f.append(disk_pct)
                elif line.startswith('* CPU usage:'):
                    if len(tokens) == 4:
                        cpu_pct = float(tokens[3].replace("%", "")) / 100.0
                    else:
                        cpu_pct = 1
                    cpu_pct_data_f.append(cpu_pct)
            if len(mem_gb_data_f) > 0 and len(disk_gb_data_f) > 0:
                filename = filepath.split('/')[-1]
                entry = filename.replace(".monitoring.log", "")
                task = entry.split('.')[0]

                max_mem_gb = max(mem_gb_data_f)
                max_mem_pct = max(mem_pct_data_f)
                max_disk_gb = max(disk_gb_data_f)
                max_disk_pct = max(disk_pct_data_f)
                max_cpu_pct = max(cpu_pct_data_f)
                max_cpu = max_cpu_pct * total_cpu

                delta_time = end_time - start_time
                delta_hour = (delta_time.total_seconds() /
                              3600.) + (overhead_min / 60.0)
                cpu_hour = total_cpu * delta_hour
                mem_hour = total_mem * delta_hour
                disk_hour = total_disk * delta_hour
                max_cpu_hour = max_cpu_pct * total_cpu * delta_hour
                max_mem_hour = max_mem_gb * delta_hour
                max_disk_hour = max_disk_gb * delta_hour

                cost_mem = COST_PER_GB_MEM_HR * mem_hour
                cost_mem_opt = COST_PER_GB_MEM_HR * \
                    max(max_mem_gb, MIN_MEM_GB) * delta_hour

                cost_disk = COST_PER_GB_DISK_HR * \
                    (total_disk + BOOT_DISK_GB) * delta_hour
                cost_disk_opt = COST_PER_GB_DISK_HR * \
                    (max(max_disk_gb, MIN_DISK_GB) + BOOT_DISK_GB) * delta_hour

                cost_cpu = COST_CPU_HR * total_cpu * delta_hour
                cost_cpu_opt = COST_CPU_HR * \
                    max(max_cpu, MIN_MEM_GB) * delta_hour

                data[entry] = {
                    "task": task,
                    "delta_hour": delta_hour,
                    "total_cpu": total_cpu,
                    "total_mem": total_mem,
                    "total_disk": total_disk,
                    "max_cpu": max_cpu,
                    "max_cpu_pct": max_cpu_pct,
                    "max_mem_gb": max_mem_gb,
                    "max_mem_pct": max_mem_pct,
                    "max_disk_gb": max_disk_gb,
                    "max_disk_pct": max_disk_pct,
                    "cpu_hour": cpu_hour,
                    "mem_hour": mem_hour,
                    "disk_hour": disk_hour,
                    "max_cpu_hour": max_cpu_hour,
                    "max_mem_hour": max_mem_hour,
                    "max_disk_hour": max_disk_hour,
                    "cost_cpu": cost_cpu,
                    "cost_cpu_opt": cost_cpu_opt,
                    "cost_mem": cost_mem,
                    "cost_mem_opt": cost_mem_opt,
                    "cost_disk": cost_disk,
                    "cost_disk_opt": cost_disk_opt
                }
    return data


def get_data_field(name, data):
    return [x[name] for x in data]


def calc_group(data):
    task_names = data.task.unique()
    group_data = {}
    for task in task_names:
        d = data.loc[data['task'] == task]
        hours = np.sum(d["delta_hour"])
        avg_cpu = np.mean(d["total_cpu"])
        avg_mem = np.mean(d["total_mem"])
        max_mem = np.max(d["max_mem_gb"])
        max_cpu = np.max(d["max_cpu"])
        max_cpu_pct = np.max(d["max_cpu_pct"])
        max_mem_pct = np.max(d["max_mem_pct"])
        avg_disk = np.mean(d["total_disk"])
        max_disk = np.max(d["max_disk_gb"])
        max_disk_pct = np.max(d["max_disk_pct"])
        cpu_hour = np.sum(d["cpu_hour"])
        mem_hour = np.sum(d["mem_hour"])
        disk_hour = np.sum(d["disk_hour"])
        max_cpu_hour = np.max(d["max_cpu_hour"])
        max_mem_hour = np.max(d["max_mem_hour"])
        max_disk_hour = np.max(d["max_disk_hour"])
        cost_cpu = np.sum(d["cost_cpu"])
        cost_cpu_dyn = np.sum(d["cost_cpu_opt"])
        cost_mem = np.sum(d["cost_mem"])
        cost_mem_dyn = np.sum(d["cost_mem_opt"])
        cost_disk = np.sum(d["cost_disk"])
        cost_disk_dyn = np.sum(d["cost_disk_opt"])

        cost_cpu_static = COST_CPU_HR * max(max_cpu, MIN_CPU) * hours
        cost_mem_static = COST_PER_GB_MEM_HR * max(max_mem, MIN_MEM_GB) * hours
        cost_disk_static = COST_PER_GB_DISK_HR * \
            (max(max_disk, MIN_DISK_GB) + BOOT_DISK_GB) * hours

        group_data[task] = {
            "hours": hours,
            "avg_cpu": avg_cpu,
            "avg_mem": avg_mem,
            "avg_disk": avg_disk,
            "max_cpu": max_cpu,
            "max_cpu_pct": max_cpu_pct,
            "max_mem": max_mem,
            "max_mem_pct": max_mem_pct,
            "max_disk": max_disk,
            "max_disk_pct": max_disk_pct,
            "cpu_hour": cpu_hour,
            "mem_hour": mem_hour,
            "disk_hour": disk_hour,
            "max_cpu_hour": max_cpu_hour,
            "max_mem_hour": max_mem_hour,
            "max_disk_hour": max_disk_hour,
            "cost_cpu": cost_cpu,
            "cost_cpu_static": cost_cpu_static,
            "cost_cpu_dyn": cost_cpu_dyn,
            "cost_mem": cost_mem,
            "cost_mem_static": cost_mem_static,
            "cost_mem_dyn": cost_mem_dyn,
            "cost_disk": cost_disk,
            "cost_disk_static": cost_disk_static,
            "cost_disk_dyn": cost_disk_dyn,
            "total_cost": cost_cpu + cost_mem + cost_disk,
            "total_cost_static": cost_cpu_static + cost_mem_static + cost_disk_static,
            "total_cost_dyn": cost_cpu_dyn + cost_mem_dyn + cost_disk_dyn
        }
    return group_data


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


def create_graphs(data, out_files_base, semilog=False, num_samples=None):
    tasks = data.index
    if num_samples is not None:
        data = data / num_samples
        ytitle = "Cost, $/sample"
        title = "Estimated Cost Per Sample"
    else:
        ytitle = "Cost, $"
        title = "Estimated Total Cost"

    if semilog:
        yscale = "log"
    else:
        yscale = "linear"

    do_simple_bar(data=[data["total_cost"], data["total_cost_static"], data["total_cost_dyn"]],
                  xticks=tasks,
                  path=out_files_base + ".cost.png",
                  bar_width=1,
                  height=8,
                  width=12,
                  xtitle="Task",
                  ytitle=ytitle,
                  title=title,
                  bottom_adjust=0.35,
                  legend=["Current", "Unif", "Pred"],
                  yscale=yscale,
                  sort_values=data["total_cost"])


# Main function
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "log_dir", help="Path containing monitoring script logs ending in \".monitoring.log\"")
    parser.add_argument("output_file", help="Output tsv file base path")
    parser.add_argument("--overhead", help="Localization overhead in minutes")
    parser.add_argument("--semilog", help="Plot semilog y",
                        action="store_true")
    parser.add_argument(
        "--plot-norm", help="Specify number of samples to normalize plots to per sample")
    args = parser.parse_args()

    if not args.overhead:
        overhead = DEFAULT_OVERHEAD_MIN
    else:
        overhead = float(args.overhead)

    if args.plot_norm:
        plot_norm = int(args.plot_norm)
    else:
        plot_norm = None

    log_dir = args.log_dir
    out_file = args.output_file
    data = read_data(log_dir, overhead_min=overhead)
    df = pd.DataFrame(data).T
    group_data = calc_group(df)
    group_df = pd.DataFrame(group_data).T
    df.to_csv(path_or_buf=out_file + ".all.tsv", sep="\t")
    group_df.to_csv(path_or_buf=out_file + ".grouped.tsv", sep="\t")
    create_graphs(group_df, out_file, semilog=args.semilog,
                  num_samples=plot_norm)


if __name__ == "__main__":
    main()
