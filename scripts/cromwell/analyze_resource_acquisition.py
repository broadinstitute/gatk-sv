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
from os.path import basename, isfile, getsize
SEED = 42

"""
Summary: Scrapes workflow metadata to analyze resource acquisition (VMs, CPU, RAM, disk memory), for the purpose
	of understanding resource peaks and timing jobs to avoid hitting quotas.

Usage: 
	python analyze_resource_acquisition workflow_metadata.json /path/to/output_basename [optional flags]
	Optional flags:
	--plot_title WorkflowName (WorkflowName will be prepended to plot title)
	--save_table (to save a TSV file of the table used to create the plot)
	--override_warning (to attempt to run script on workflow that raised a warning)

Parameters:
	workflow_metadata.json: path to Cromwell metadata file for workflow of interest
	/path/to/output_basename: base output path, to which extensions will be appended for each output file, 
		ie. /output_dir/basename will yield /output_dir/basename.plot.png, /output_dir/basename.table.tsv, etc
		or, /output_dir/ will yield plot.png, table.tsv, peaks.txt, etc

Outputs:
	- plot (png) of VMs, CPUs (total, preemptible, and nonpreemptible), RAM, and disk (HDD, SSD) acquisitioned over time 
	- (optional, with --save_table flag) table of resource acquisition over time for each type of resource listed above
	- peak resource acquisition for each type of resource listed above (TSV file)
	- number of non-preemptible VMs used (prints to stdout) and the names of the tasks that used them (TSV file, if any)
	- number of tasks call-cached (prints to stdout), and their task names (TSV file, if any)
"""
NUM_CACHED = 0
CACHED = dict()
NUM_NONPREEMPTIBLE = 0
NONPREEMBTIBLE_TASKS = dict()

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
	"""
	Source: https://github.com/broadinstitute/dsde-pipelines/blob/develop/scripts/calculate_cost.py
	"""
	if "runtimeAttributes" in metadata and "preemptible" in metadata['runtimeAttributes']:
		pe_count = int(metadata['runtimeAttributes']["preemptible"])
		attempt = int(metadata['attempt'])

		return attempt <= pe_count
	else:
		# we can't tell (older metadata) so conservatively return false
		return False

def used_cached_results(metadata):
	"""
	Modified from: https://github.com/broadinstitute/dsde-pipelines/blob/develop/scripts/calculate_cost.py
	"""
	return "callCaching" in metadata and "hit" in metadata["callCaching"] and metadata["callCaching"]["hit"]

def calculate_start_end(call_info):
	"""
	Modified from: https://github.com/broadinstitute/dsde-pipelines/blob/develop/scripts/calculate_cost.py
	"""
	# get start (start time of VM start) & end time (end time of 'ok') according to metadata
	start = None
	end = None

	if 'executionEvents' in call_info:
		for x in call_info['executionEvents']:
			# ignore incomplete executionEvents (could be due to server restart or similar)
			if 'description' not in x:
				continue
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

	return start, end

def get_mem_cpu(m):
	"""
	Modified from: https://github.com/broadinstitute/dsde-pipelines/blob/develop/scripts/calculate_cost.py
	"""
	cpu = 'na'
	memory = 'na'
	if 'runtimeAttributes' in m:
		if 'cpu' in m['runtimeAttributes']:
			cpu = int(m['runtimeAttributes']['cpu'])
		if 'memory' in m['runtimeAttributes']:
			mem_str = m['runtimeAttributes']['memory']
			memory = float(mem_str[:mem_str.index(" ")])
	return cpu, memory

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
		# In alias, track hierarchy of workflow/task up to current task nicely without repetition
		if alias is None:
			alias = ""
		to_add = ""
		if 'wdl-call-alias' in m['labels']:
			to_add = m['labels']['wdl-call-alias']
		elif 'wdl-task-name' in m['labels']:
			to_add = m['labels']['wdl-task-name']
		if to_add != "" and not alias.endswith(to_add):
			if alias != "" and alias[-1] != ".":
				alias += "."
			alias += to_add

	

	call_metadata = []
	if 'calls' in m:
		for call in m['calls']:
			# Skips scatters that don't contain calls
			if '.' not in call:
				continue
			# In call_alias, track hierarchy of workflow/task up to current call nicely without repetition
			if alias is None:
				alias = ""
			call_split = call.split('.')
			call_name = call
			if alias.endswith(call_split[0]):
				call_name = call_split[1]
			call_alias = alias
			if call_alias != "" and call_alias[-1] != ".":
				call_alias += "."
			call_alias += call_name
			# recursively get metadata
			call_metadata.extend(getCalls(m['calls'][call], alias=call_alias))

	if 'subWorkflowMetadata' in m:
		call_metadata.extend(getCalls(m['subWorkflowMetadata'], alias=alias))

	# in a call
	if alias and ('monitoringLog' in m):
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

		cpu, memory = get_mem_cpu(m)

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

		call_metadata.append((start, 1, cpu, preemptible_cpu, nonpreemptible_cpu, memory, hdd_size, ssd_size))
		call_metadata.append((end, -1, -1*cpu, -1*preemptible_cpu, -1*nonpreemptible_cpu, -1*memory, -1*hdd_size, -1*ssd_size))
		if not preemptible:
			global NUM_NONPREEMPTIBLE
			global NONPREEMBTIBLE_TASKS
			NUM_NONPREEMPTIBLE += 1
			if alias in NONPREEMBTIBLE_TASKS:
				NONPREEMBTIBLE_TASKS[alias] += 1
			else:
				NONPREEMBTIBLE_TASKS[alias] = 1

		if cached:
			global CACHED
			global NUM_CACHED
			NUM_CACHED += 1
			if alias in CACHED:
				CACHED[alias] += 1
			else:
				CACHED[alias] = 1

	return call_metadata

def check_workflow_valid(metadata, metadata_file, override_warning):
	# these errors cannot be overcome
	if 'status' not in metadata:
		print("ERROR: Incomplete metadata input file %s. File lacks workflow status field." % metadata_file, file=sys.stderr)
		exit(1)
	if metadata['status'] == "fail": # Unrecognized workflow ID failure - unable to download metadata
		err_msg = "ERROR: Workflow metadata download failure."
		if 'message' in metadata:
			err_msg += " Message: " + metadata['message']
		print(err_msg, file=sys.stderr)
		exit(1)
	# these errors may be able to be overcome for partial output
	found_retryable_error = False
	if metadata['status'] == "Failed":
		print("Warning: Workflow failed.", file=sys.stderr)
		found_retryable_error = True
	for event in metadata['workflowProcessingEvents']:
		if event['description'] == "Released":
			print("Warning: Server was interrupted during workflow execution, which is likely to impact plot accuracy.", file=sys.stderr)
			found_retryable_error = True
			break
	if found_retryable_error:
		if override_warning:
			print("Override_warning=TRUE. Proceeding with caution.", file=sys.stderr)
		else:
			print("To attempt to proceed anyway, re-run the script with the --override_warning flag.", file=sys.stderr)
			exit(1)

def get_call_metadata(metadata_file, override_warning=False):
	"""
	Based on: https://github.com/broadinstitute/gatk-sv/blob/master/scripts/cromwell/download_monitoring_logs.py
	"""
	metadata = json.load(open(metadata_file, 'r'))
	check_workflow_valid(metadata, metadata_file, override_warning)
	colnames = ['timestamp', 'vm_delta', 'cpu_all_delta', 'cpu_preemptible_delta', 'cpu_nonpreemptible_delta', 'memory_delta', 'hdd_delta', 'ssd_delta']

	call_metadata = getCalls(metadata)
	if len(call_metadata) == 0:
		print("ERROR: No calls in workflow metadata.", file=sys.stderr)
		exit(1)
	call_metadata = pd.DataFrame(call_metadata, columns=colnames)

	return call_metadata

def transform_call_metadata(call_metadata):
	"""
	Based on: https://github.com/broadinstitute/dsde-pipelines/blob/master/scripts/quota_usage.py
	"""
	call_metadata = call_metadata.sort_values(by='timestamp')
	call_metadata['timestamp'] -= call_metadata.timestamp.iloc[0] # make timestamps start from 0 by subtracting minimum (at index 0 after sorting)
	call_metadata['seconds'] = call_metadata['timestamp'].dt.total_seconds() # get timedelta in seconds because plot labels won't format correctly otherwise
	
	call_metadata['vm'] = call_metadata.vm_delta.cumsum()
	call_metadata['cpu_all'] = call_metadata.cpu_all_delta.cumsum()
	call_metadata['cpu_preemptible'] = call_metadata.cpu_preemptible_delta.cumsum()
	call_metadata['cpu_nonpreemptible'] = call_metadata.cpu_nonpreemptible_delta.cumsum()
	call_metadata['memory'] = call_metadata.memory_delta.cumsum()
	call_metadata['ssd'] = call_metadata.ssd_delta.cumsum()
	call_metadata['hdd'] = call_metadata.hdd_delta.cumsum()

	return call_metadata

def plot_resources_time(df, title_name, output_name):
	"""
	Modified from: https://github.com/broadinstitute/dsde-pipelines/blob/master/scripts/quota_usage.py
	"""
	print("Writing " + output_name, file = sys.stderr)
	colors = {
		"vm": "#006FA6", #blue 
		"cpu_all": "black",
		"cpu_preemptible": "#10a197", #turquoise 
		"cpu_nonpreemptible": "#A30059", #dark pink
		"memory": "#FF4A46", #coral red
		"hdd": "#72418F", #purple
		"ssd": "#008941", #green
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

def write_resources_time_table(call_metadata, table_file):
	print("Writing " + table_file, file = sys.stderr)
	call_metadata.to_csv(
		table_file, 
		columns = ["timestamp", "seconds", "vm", "cpu_all", "cpu_preemptible", "cpu_nonpreemptible", "memory", "hdd", "ssd"],
		sep='\t', 
		index = False
	)

def write_peak_usage(m, peak_file):
	print("Writing " + peak_file, file = sys.stderr)
	with open(peak_file, 'w') as out:
		out.write("peak_vms\t" + str(max(m['vm'])) + "\n")
		out.write("peak_cpu_all\t" + str(max(m['cpu_all'])) + "\n")
		out.write("peak_cpu_preemptible\t" + str(max(m['cpu_preemptible'])) + "\n")
		out.write("peak_cpu_nonpreemptible\t" + str(max(m['cpu_nonpreemptible'])) + "\n")
		out.write("peak_ram_gib\t" + "{:.2f}".format(max(m['memory'])) + "\n")
		out.write("peak_disk_hdd_gib\t" + str(max(m['hdd'])) + "\n")
		out.write("peak_disk_ssd_gib\t" + str(max(m['ssd'])) + "\n")

def write_cached_warning(cached_file):
	global CACHED
	global NUM_CACHED
	if NUM_CACHED > 0:
		print("%d cached tasks found, writing tasks to %s." % (NUM_CACHED, cached_file))
		with open(cached_file, 'w') as cached_out:
			cached_out.write("#task_name\tnum_cached\n")
			cached_out.write("all_tasks\t%d\n" % NUM_CACHED)
			cached_out.write("\n".join([x + '\t' + str(CACHED[x]) for x in sorted(list(CACHED.keys()))]) + "\n")
	else:
		print("0 cached tasks found.")

def write_nonpreemptible_vms(vms_file):
	global NUM_NONPREEMPTIBLE
	global NONPREEMBTIBLE_TASKS
	if NUM_NONPREEMPTIBLE > 0:
		print("%d non-preemptible VMs found, writing tasks to %s." % (NUM_NONPREEMPTIBLE, vms_file))
		with open(vms_file, 'w') as vms_out:
			vms_out.write("#task_name\tnum_nonpreemptible\n")
			vms_out.write("all_tasks\t%d\n" % NUM_NONPREEMPTIBLE)
			vms_out.write("\n".join([x + '\t' + str(NONPREEMBTIBLE_TASKS[x]) for x in sorted(list(NONPREEMBTIBLE_TASKS.keys()))]) + "\n")
	else:
		print("0 non-preemptible VMs found.")
			
def check_file_nonempty(f):
	if not isfile(f): 
		print("ERROR: Required metadata input file %s does not exist." % f, file=sys.stderr)
		exit(1)
	elif getsize(f) == 0:
		print("ERROR: Required metadata input file %s is empty." % f, file=sys.stderr)
		exit(1)

# Main function
def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("workflow_metadata", help="Workflow metadata JSON file")
	parser.add_argument("output_base", help="Output directory + basename")
	parser.add_argument("--plot_title", help="Provide workflow name for plot title: <name> Resource Acquisition Over Time", required=False, default = "")
	parser.add_argument("--override_warning", help="Execute script workflow warning (server interrupted, workflow failed, etc.), which may impact plot accuracy", \
		required=False, default=False, action='store_true')
	parser.add_argument("--save_table", help="Save TSV copy of resources over time table used to make plot", required=False, default=False, action='store_true')
	args = parser.parse_args()
	random.seed(SEED)

	metadata_file, output_base, plt_title, override_warning, save_table = args.workflow_metadata, args.output_base, args.plot_title, args.override_warning, args.save_table
	if plt_title != "":
		plt_title += " "
	sep = "."
	if basename(output_base) == "":
		sep = ""

	check_file_nonempty(metadata_file)
	call_metadata = get_call_metadata(metadata_file, override_warning)
	call_metadata = transform_call_metadata(call_metadata)

	cached_file = output_base + sep + "cached.tsv"
	write_cached_warning(cached_file)

	vms_file = output_base + sep + "vms_file.tsv"
	write_nonpreemptible_vms(vms_file)
	
	plot_file = output_base + sep + "plot.png"
	plot_resources_time(call_metadata, plt_title, plot_file)

	if save_table:
		table_file = output_base + sep + "table.tsv"
		write_resources_time_table(call_metadata, table_file)

	peak_file = output_base + sep + "peaks.tsv"
	write_peak_usage(call_metadata, peak_file)


if __name__== "__main__":
	main()





