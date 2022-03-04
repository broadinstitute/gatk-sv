#!/bin/python

# Synopsis:
#  Creates input values for a new test batch from a GATKSVPipelineBatch run
#

import argparse
import json
import sys


INPUT_KEYS = set([
    "name",
    "samples",
    "bam_or_cram_files",
    "requester_pays_crams",
    "gvcfs",
    "ped_file",
    "contig_ploidy_model_tar",
    "gcnv_model_tars",
    "qc_definitions",
    "outlier_cutoff_table"
])


def replace_output_dir_maybe_list(value, execution_bucket, outputs_dir):
    if outputs_dir is None:
        return value
    if isinstance(value, list):
        return [replace_output_dir(v, execution_bucket, outputs_dir) for v in value]
    else:
        return replace_output_dir(value, execution_bucket, outputs_dir)


def replace_output_dir(value, execution_bucket, outputs_dir):
    if execution_bucket not in value:
        raise ValueError(f"Execution bucket {execution_bucket} not found in output: {value}")
    return value.replace(execution_bucket, outputs_dir)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("metadata", help="GATKSVPipelineBatch metadata JSON file")
    parser.add_argument("--final-workflow-outputs-dir", help="If used, final_workflow_outputs_dir "
                                                             "option from Cromwell config file")
    parser.add_argument("--execution-bucket", help="Cromwell execution bucket, required if "
                                                   "using --final-workflow-outputs-dir")
    args = parser.parse_args()

    execution_bucket = args.execution_bucket
    outputs_dir = args.final_workflow_outputs_dir
    if outputs_dir is not None:
        if execution_bucket is None:
            raise ValueError("Must supply --execution-bucket if using --final-workflow-outputs-dir")
        if not execution_bucket.startswith("gs://"):
            raise ValueError("--execution-bucket must start with gs://")
        if not outputs_dir.startswith("gs://"):
            raise ValueError("--final-workflow-outputs-dir must start with gs://")
        if execution_bucket.endswith('/'):
            execution_bucket = execution_bucket[:-1]
        if outputs_dir.endswith('/'):
            outputs_dir = outputs_dir[:-1]

    with open(args.metadata, 'r') as f:
        metadata = json.load(f)
    values = {key.replace("GATKSVPipelineBatch.", ""):
              replace_output_dir_maybe_list(value, execution_bucket, outputs_dir)
              for key, value in metadata["outputs"].items() if value is not None}
    inputs = metadata["inputs"]
    for raw_key in set(inputs.keys()).intersection(INPUT_KEYS):
        key = raw_key.split('.')[-1]
        values[key] = inputs[key]
    for key in INPUT_KEYS - set(values.keys()):
        sys.stderr.write(f"Warning: expected workflow input '{key}' not found in metadata. You will need to add "
                         f"this entry manually.\n")
        values[key] = None

    print(json.dumps(values, sort_keys=True, indent=4))


if __name__ == "__main__":
    main()
