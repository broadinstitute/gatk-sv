#!/bin/python

import json
import argparse
import numpy as np
import dateutil.parser
from dateutil.tz import tzutc
import datetime
import time
import pandas as pd
import math
import random
import matplotlib.pyplot as plt
import matplotlib.ticker
import sys
SEED = 42

"""
Summary: Scrapes workflow metadata to analyze resource acquisition (VMs, CPU, RAM, disk memory), for the purpose
	of understanding resource peaks and timing jobs to avoid hitting quotas.

Usage: 
	python analyze_resource_acquisition workflow_metadata.json /path/to/output_basename

Parameters:
	workflow_metadata.json: path to Cromwell metadata file for workflow of interest
	/path/to/output_basename: base output path, to which extensions will be appended for each output file, 
		ie. /output_dir/basename will yield /output_dir/basename_plot.png, /output_dir/basename_table.tsv, etc

Outputs:
	- plot (png) of VMs, CPUs (total, preemptible, and nonpreemptible), RAM, and disk (HDD, SSD) acquisitioned over time 
	- table of resource acquisition over time for each type of resource listed above
	- text file containing peak resource acuqisition for each type of resource listed above
	- [preemptible things]
	- text file containing list of tasks that used call-caching
"""

def get_seconds_from_epoch(time_string):
	return (dateutil.parser.parse(time_string) - datetime.datetime(1900,1,1,tzinfo=tzutc())).total_seconds()

def get_disk_info(metadata):
	"""
	Modified from: https://github.com/broadinstitute/dsde-pipelines/blob/develop/scripts/calculate_cost.py
	Modified to return (hdd_size, ssd_size)
	"""
	if "runtimeAttributes" in metadata and "disks" in metadata['runtimeAttributes']:
		bootDiskSizeGb = 0.0
		if "bootDiskSizeGb" in metadata['runtimeAttributes']:
			bootDiskSizeGb = float(metadata['runtimeAttributes']['bootDiskSizeGb'])
		# Note - am lumping boot disk in with requested disk.  Assuming boot disk is same type as requested.
		# i.e. is it possible that boot disk is HDD when requested is SDD.
		(name, disk_size, disk_type) = metadata['runtimeAttributes']["disks"].split()
		if disk_type == "HDD":
			return float(disk_size), float(0)
		elif disk_type == "SSD":
			return float(0), float(disk_size)
		else:
			return float(0), float(0)
	else:
		# we can't tell disk size in this case so just return nothing
		return float(0),float(0)
		
def was_preemptible_vm(metadata):
	if "runtimeAttributes" in metadata and "preemptible" in metadata['runtimeAttributes']:
		pe_count = int(metadata['runtimeAttributes']["preemptible"])
		attempt = int(metadata['attempt'])

		return attempt <= pe_count
	else:
		# we can't tell (older metadata) so conservatively return false
		return False

def used_cached_results(metadata):
	"""
	Source: https://github.com/broadinstitute/dsde-pipelines/blob/develop/scripts/calculate_cost.py
	"""
	return "callCaching" in metadata and metadata["callCaching"]["hit"]

def calculate_start_end(call_info):
	"""
	Modified from: https://github.com/broadinstitute/dsde-pipelines/blob/develop/scripts/calculate_cost.py
	"""
	# get start (start time of VM start) & end time (end time of 'ok') according to metadata
	start = None
	end = None

	if 'executionEvents' in call_info:
		for x in call_info['executionEvents']:
			y = x['description']

			if 'backend' in call_info and call_info['backend'] == 'PAPIv2':
				if y.startswith("PreparingJob"):
					start = dateutil.parser.parse(x['startTime'])
				if y.startswith("Worker released"):
					end = dateutil.parser.parse(x['endTime'])
			else:
				if y.startswith("start"):
					start = dateutil.parser.parse(x['startTime'])
				if y.startswith("ok"):
					end = dateutil.parser.parse(x['endTime'])

	# if we are preempted or if cromwell used previously cached results, we don't even get a start time from JES.
	# if cromwell was restarted, the start time from JES might not have been written to the metadata.
	# in either case, use the Cromwell start time which is earlier but not wrong.
	if start is None:
		start = dateutil.parser.parse(call_info['start'])

	# if we are preempted or if cromwell used previously cached results, we don't get an endtime from JES right now.
	# if cromwell was restarted, the start time from JES might not have been written to the metadata.
	# in either case, use the Cromwell end time which is later but not wrong
	if end is None:
		end = dateutil.parser.parse(call_info['end'])

	# start = get_seconds_from_epoch(start)
	# end = get_seconds_from_epoch(end)
	return start, end

def was_preempted(call_info):
	# We treat Preempted and RetryableFailure the same.  The latter is a general case of the former
	return call_info['executionStatus'] in ['Preempted', 'RetryableFailure']

def calculate_start_end_simple(m):
	start = 'na'
	if 'start' in m:
		start = get_seconds_from_epoch(m['start'])

	end = 'na'
	if 'end' in m:
		end = get_seconds_from_epoch(m['end'])

	return start, end

def getCalls(m, alias=None):
	"""
	Modified from download_monitoring_logs.py script by Mark Walker
	https://github.com/broadinstitute/gatk-sv/blob/master/scripts/cromwell/download_monitoring_logs.py
	"""
	if isinstance(m, list):
		call_metadata = []
		for m_shard in m:
			call_metadata.extend(getCalls(m_shard, alias=alias))
		return call_metadata

	if 'labels' in m:
		if 'wdl-call-alias' in m['labels']:
			alias = m['labels']['wdl-call-alias']
		elif 'wdl-task-name' in m['labels']:
			alias = m['labels']['wdl-task-name']

	shard_index = '-2'
	if 'shardIndex' in m:
		shard_index = m['shardIndex']

	attempt = '0'
	if 'attempt' in m:
		attempt = m['attempt']

	job_id = 'na'
	if 'jobId' in m:
		job_id = m['jobId'].split('/')[-1]

	start, end = calculate_start_end(m)

	cpu = 'na'
	memory = 'na'
	if 'runtimeAttributes' in m:
		if 'cpu' in m['runtimeAttributes']:
			cpu = int(m['runtimeAttributes']['cpu'])
		if 'memory' in m['runtimeAttributes']:
			mem_str = m['runtimeAttributes']['memory']
			memory = float(mem_str[:mem_str.index(" ")])

	cached = used_cached_results(m)

	callroot = 'na'
	if 'callRoot' in m:
		callroot = m['callRoot']

	preemptible = was_preemptible_vm(m)
	preemptible_cpu = 0
	nonpreemptible_cpu = 0
	if preemptible:
		preemptible_cpu = cpu
	else:
		nonpreemptible_cpu = cpu

	hdd_size, ssd_size = get_disk_info(m)

	call_metadata = []
	if 'calls' in m:
		for call in m['calls']:
			# Skips scatters that don't contain calls
			if '.' not in call:
				continue
			call_alias = call.split('.')[1]
			call_metadata.extend(getCalls(m['calls'][call], alias=call_alias))

	if 'subWorkflowMetadata' in m:
		call_metadata.extend(getCalls(m['subWorkflowMetadata'], alias=alias))

	# in a call
	if alias and ('monitoringLog' in m):
		#call_metadata.append((alias, start,end,shard_index, attempt, job_id, cpu, cached, callroot, memory, preemptible, hdd_size, ssd_size))
		call_metadata.append((start, 1, cpu, preemptible_cpu, nonpreemptible_cpu, memory, hdd_size, ssd_size))
		call_metadata.append((end, -1, -1*cpu, -1*preemptible_cpu, -1*nonpreemptible_cpu, -1*memory, -1*hdd_size, -1*ssd_size))

	return call_metadata

def get_data_table(metadata_file):
	"""
	Based on: 
		https://github.com/broadinstitute/dsde-pipelines/blob/master/scripts/quota_usage.py
		https://github.com/broadinstitute/gatk-sv/blob/master/scripts/cromwell/download_monitoring_logs.py
	"""
	metadata = json.load(open(metadata_file, 'r'))
	colnames = ['timestamp', 'vm_delta', 'cpu_all_delta', 'cpu_preemptible_delta', 'cpu_nonpreemptible_delta', 'memory_delta', 'hdd_delta', 'ssd_delta']

	call_metadata = getCalls(metadata, metadata['workflowName'])
	call_metadata = pd.DataFrame(call_metadata, columns=colnames)
	#call_metadata[['timestamp']] = call_metadata[['timestamp']].astype('datetime64')

	call_metadata = call_metadata.sort_values(by='timestamp')
	call_metadata['timestamp'] -= call_metadata.timestamp.iloc[0] # make timestamps start from 0 by subtracting minimum (at index 0 after sorting)
	# print(call_metadata[['timestamp']])
	call_metadata['seconds'] = call_metadata['timestamp'].dt.total_seconds()
	print(call_metadata.dtypes)
	
	call_metadata['vm'] = call_metadata.vm_delta.cumsum()
	call_metadata['cpu_all'] = call_metadata.cpu_all_delta.cumsum()
	call_metadata['cpu_preemptible'] = call_metadata.cpu_preemptible_delta.cumsum()
	call_metadata['cpu_nonpreemptible'] = call_metadata.cpu_nonpreemptible_delta.cumsum()
	call_metadata['memory'] = call_metadata.memory_delta.cumsum()
	call_metadata['ssd'] = call_metadata.ssd_delta.cumsum()
	call_metadata['hdd'] = call_metadata.hdd_delta.cumsum()

	return call_metadata

def plot(df, title_name, output_name):
	"""
	Modified from: https://github.com/broadinstitute/dsde-pipelines/blob/master/scripts/quota_usage.py
	"""
	print("Writing " + output_name, file = sys.stderr)
	colors = {
		"vm": "#006FA6", #blue (dark blue, could use if change to turquoise/cyan: "#0000A6")
		"cpu_all": "black",
		"cpu_preemptible": "#10a197", #turquoise 
		"cpu_nonpreemptible": "#A30059", #dark pink
		"memory": "#FF4A46", #coral red
		"hdd": "#72418F", #purple
		"ssd": "#008941", #green
		#  "#FF913F", # orange # "gray"
	}
	LABEL_SIZE = 17
	TITLE_SIZE = 20
	TICK_SIZE = 15

	fig,ax = plt.subplots(4,1, figsize = (14,26), sharex = True)
	ax[0].set_title(title_name + "Resource Acquisition Over Time", fontsize=TITLE_SIZE)

	ax[0].plot(df['seconds'], df['vm'], color = colors["vm"])
	ax[0].set_ylabel("VMs", fontsize = LABEL_SIZE)
	plt.setp(ax[0].get_yticklabels(), fontsize = TICK_SIZE)

	ax[1].plot(df['seconds'], df['cpu_all'], color = colors["cpu_all"], linewidth = 2, label = "All")
	ax[1].plot(df['seconds'], df['cpu_preemptible'], color = colors["cpu_preemptible"], linestyle = "dashed", label = "Preemptible")
	ax[1].plot(df['seconds'], df['cpu_nonpreemptible'], color = colors["cpu_nonpreemptible"], linestyle = "dashed", label = "Non-preemptible")
	ax[1].set_ylabel("CPU Cores", fontsize = LABEL_SIZE)
	plt.setp(ax[1].get_yticklabels(), fontsize = TICK_SIZE)
	ax[1].legend(loc = "upper right", title = "CPU Types", fontsize = TICK_SIZE, title_fontsize = TICK_SIZE)

	ax[2].plot(df['seconds'], df['memory'], color = colors["memory"])
	ax[2].set_ylabel("RAM (GiB)", fontsize = LABEL_SIZE)
	plt.setp(ax[2].get_yticklabels(), fontsize = TICK_SIZE)

	ax[3].plot(df['seconds'], df['hdd'], color = colors["hdd"], label = "HDD")
	ax[3].plot(df['seconds'], df['ssd'], color = colors["ssd"], label = "SSD")
	ax[3].set_ylabel("Disk Memory (GiB)", fontsize = LABEL_SIZE)
	plt.setp(ax[3].get_yticklabels(), fontsize = TICK_SIZE)
	ax[3].legend(loc = "upper right", title = "Disk Types", fontsize = TICK_SIZE, title_fontsize = TICK_SIZE)

	formatter = matplotlib.ticker.FuncFormatter(lambda x,pos: str(datetime.timedelta(seconds = x)))
	ax[3].xaxis.set_major_formatter(formatter) 
	plt.setp(ax[3].get_xticklabels(), rotation = 15, fontsize = TICK_SIZE)
	ax[3].set_xlabel("Time", fontsize = LABEL_SIZE)
	
	fig.savefig(output_name, bbox_inches='tight')

def print_peak_usage(m):
	print("Peak VMs: " + str(max(m['vm'])))
	print("Peak CPUs (all): " + str(max(m['cpu_all'])))
	print("Peak Preemptible CPUs: " + str(max(m['cpu_preemptible'])))
	print("Peak Non-preemptible CPUs: " + str(max(m['cpu_nonpreemptible'])))
	print("Peak RAM: " + "{:.2f}".format(max(m['memory'])) + " GiB")
	print("Peak HDD: " + str(max(m['hdd'])) + " GiB")
	print("Peak SSD: " + str(max(m['ssd'])) + " GiB")

# Main function
def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("workflow_metadata", help="Workflow metadata JSON file")
	parser.add_argument("output_base", help="Output directory + basename")
	parser.add_argument("--plot_title", help="Workflow name for plot title: <name> Resource Acquisition Over Time", required=False, default = "")
	args = parser.parse_args()
	random.seed(SEED)

	metadata_file = args.workflow_metadata
	output_base = args.output_base
	plt_title = args.plot_title
	if plt_title != "":
		plt_title += " "

	call_metadata = get_data_table(metadata_file)
	
	plot_file = output_base + ".plot.png"
	plot(call_metadata, plt_title, plot_file)

	table_file = output_base + ".table.tsv"
	print("Writing " + table_file)
	call_metadata.to_csv(table_file, sep='\t', index = False)

	print_peak_usage(call_metadata)


if __name__== "__main__":
	main()










