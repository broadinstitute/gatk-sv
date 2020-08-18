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
SEED = 42

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
		call_metadata.append((start, cpu, preemptible_cpu, nonpreemptible_cpu, memory, hdd_size, ssd_size))
		call_metadata.append((end, -1*cpu, -1*preemptible_cpu, -1*nonpreemptible_cpu, -1*memory, -1*hdd_size, -1*ssd_size))

	return call_metadata

def get_data_table(metadata_file):
	metadata = json.load(open(metadata_file, 'r'))
	#call_metadata = np.asarray(getCalls(metadata, metadata['workflowName']))
	#colnames = ['alias','start','end','shard_index','attempt','job_id','cpu', 'cached','callroot', 'memory', 'preemptible', 'hdd_size', 'ssd_size']
	colnames = ['timestamp','cpu_all', 'memory', 'preemptible', 'hdd_size', 'ssd_size']

	call_metadata = getCalls(metadata, metadata['workflowName'])
	print(len(call_metadata))
	print(len(set(call_metadata)))
	call_metadata = pd.DataFrame(call_metadata, columns=colnames)
	print(call_metadata)

	mintime = min(call_metadata['start'])
	call_metadata['start'] -= mintime
	call_metadata['end'] -= mintime

	return call_metadata

def get_data_table_cumulative(metadata_file):
	metadata = json.load(open(metadata_file, 'r'))
	colnames = ['timestamp','cpu_all_delta', 'cpu_preemptible_delta', 'cpu_nonpreemptible_delta','memory_delta', 'hdd_delta', 'ssd_delta']

	call_metadata = getCalls(metadata, metadata['workflowName'])
	call_metadata = pd.DataFrame(call_metadata, columns=colnames)
	#call_metadata[['timestamp']] = call_metadata[['timestamp']].astype('datetime64')

	call_metadata = call_metadata.sort_values(by='timestamp')
	#print(call_metadata[['timestamp']])

	call_metadata['cpu_all'] = call_metadata.cpu_all_delta.cumsum()
	call_metadata['cpu_preemptible'] = call_metadata.cpu_preemptible_delta.cumsum()
	call_metadata['cpu_nonpreemptible'] = call_metadata.cpu_nonpreemptible_delta.cumsum()
	call_metadata['memory'] = call_metadata.memory_delta.cumsum()
	call_metadata['ssd'] = call_metadata.ssd_delta.cumsum()
	call_metadata['hdd'] = call_metadata.hdd_delta.cumsum()

	return call_metadata

def make_time_bins(call_metadata, count_by):
	total_elapsed_seconds = max(call_metadata['end'])
	time_ranges = np.array(range(0,math.ceil(total_elapsed_seconds),count_by))
	bins = len(time_ranges)
	time_ranges = np.vstack((time_ranges, time_ranges + count_by)).T

	return bins, time_ranges

def get_binned_data(call_metadata):
	#wip 
	count_by = 5 # bin resolution, in seconds
	bins, time_ranges = make_time_bins(call_metadata, count_by)
	
	cols_list = ["vms", "all_cpus", "preemptibles", "non-preemptibles","memory", "hdd", "ssd"]
	cols = {x:i for i,x in enumerate(cols_list)}
	data = np.zeros((bins, len(cols_list)))

	cpus = np.zeros(bins)
	for rownum in range(len(call_metadata)):
		call = call_metadata.loc[rownum]
		#print("alias: %s, shard_index: %s, attempt: %s, id: %s, start: %.2f, end: %.2f"% (call['alias'], call['shard_index'], call['attempt'],call['job_id'],call['start'], call['end']))
		if call['cpu'] != 'na':
			cpus[np.intersect1d(np.where(call['start'] <= time_ranges[:,1])[0], np.where(call['end'] > time_ranges[:,0])[0])] += call['cpu']
	#print(cpus)

def plot(df, title_name, output_name):
	"""
	Modified from: https://github.com/broadinstitute/dsde-pipelines/blob/master/scripts/quota_usage.py
	"""
	ax = df[['timestamp','cpu_all','cpu_preemptible', 'cpu_nonpreemptible', 'ssd','hdd']].plot(x='timestamp', secondary_y=['ssd', 'hdd'], title=title_name, drawstyle="steps-post")
	ax.set_xlabel("Time")
	ax.set_ylabel("Cores")
	ax.legend(loc='center')
	ax.right_ax.set_ylabel("Disk (GiB)")
	
	h1, l1 = ax.get_legend_handles_labels()
	h2, l2 = ax.right_ax.get_legend_handles_labels()
	ax.legend(h1+h2, l1+l2, loc="upper center")
	
	ax.get_figure().savefig(output_name)

def print_peak_usage(m):
	print("Peak CPUs (all): " + str(max(m['cpu_all'])))
	print("Peak Preemptible CPUs: " + str(max(m['cpu_preemptible'])))
	print("Peak Non-preemptible CPUs: " + str(max(m['cpu_nonpreemptible'])))
	print("Peak HDD: " + str(max(m['hdd'])) + " GiB")
	print("Peak SSD: " + str(max(m['ssd'])) + " GiB")

# Main function
def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("workflow_metadata", help="Workflow metadata JSON file")
	parser.add_argument("output_base", help="Output directory + basename")
	args = parser.parse_args()
	random.seed(SEED)

	metadata_file = args.workflow_metadata
	output_base = args.output_base

	call_metadata = get_data_table_cumulative(metadata_file)
	
	png_file = output_base + ".png"
	print("Writing " + png_file)
	plot(call_metadata, metadata_file, png_file)

	txt_file = output_base + ".txt"
	print("Writing " + txt_file)
	call_metadata.to_csv(txt_file, sep='\t')

	print_peak_usage(call_metadata)


if __name__== "__main__":
	main()










