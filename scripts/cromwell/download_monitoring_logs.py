#!/bin/python

import json
import os.path
import argparse
from google.cloud import storage
from joblib import Parallel, delayed
import random

# Synopsis:
#  Downloads all Cromwell monitoring logs from a given workflow. Monitoring logs should be generated using
#  cromwell_monitoring_script.sh. Log files are named with the format:
#    <task_name>.<shard#>.<attempt#>.<job_id>.monitoring.log
#  Note downloaded logs will not be overwritten if they already exist. If you interrupt a download, make sure to
#  delete any partially-downloaded files. Empty files are created for logs that could not be found in GCS.
#
# Usage:
#   python download_monitoring_logs.py workflow_metadata.json /output/dir
#
# Parameters:
#   workflow_metadata.json : Workflow metadata file
#   /output/dir : Directory to place logs
#
# Author: Mark Walker (markw@broadinstitute.org)

# Download threads
NUM_THREADS = 8
RAND_SEED = 7282993


def getCalls(m, alias=None):
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
        call_metadata.append((m, alias, shard_index, attempt, job_id))

    return call_metadata


def download(data, output_dir):
    (m, alias, shard_index, attempt, job_id) = data
    if job_id != 'na':
        output_dest = output_dir + '/' + alias + '.' + \
            str(shard_index) + '.' + str(attempt) + \
            '.' + job_id + '.monitoring.log'
        log_url = m['monitoringLog']
        if os.path.isfile(output_dest):
            print("skipping " + log_url)
            return
        with open(output_dest, 'wb') as f:
            client = storage.Client()
            tokens = log_url.split('/')
            bucket_name = tokens[2]
            bucket_object = '/'.join(tokens[3:])
            bucket = client.get_bucket(bucket_name)
            blob = bucket.get_blob(bucket_object)
            if blob:
                print(log_url)
                blob.download_to_file(f)

# Main function


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("workflow_metadata",
                        help="Workflow metadata JSON file")
    parser.add_argument("output_dir", help="Output directory")
    args = parser.parse_args()
    random.seed(RAND_SEED)

    metadata_file = args.workflow_metadata
    output_dir = args.output_dir

    metadata = json.load(open(metadata_file, 'r'))
    call_metadata = getCalls(metadata, metadata['workflowName'])
    Parallel(n_jobs=NUM_THREADS)(delayed(download)(d, output_dir)
                                 for d in call_metadata)


if __name__ == "__main__":
    main()
