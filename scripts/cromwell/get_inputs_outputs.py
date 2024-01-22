#!/bin/python

import json
import argparse
import os

# Synopsis:
#  Generates JSON files with the inputs and outputs of every (sub)workflow
#
# Usage:
#   python get_inputs_outputs.py workflow_metadata.json /output/dir
#
# Parameters:
#   workflow_metadata.json : Workflow metadata file
#   /output/dir : Directory to place logs
#
# Author: Mark Walker (markw@broadinstitute.org)


def get_subworkflows(m, alias):
    if isinstance(m, list):
        return get_subworkflows(m[0], alias)

    task = ''
    if 'workflowName' in m:
        task = m['workflowName']

    # in a call
    if not ('subWorkflowMetadata' in m or 'calls' in m):
        return []

    call_metadata = []
    if 'calls' in m:
        for call in m['calls']:
            call_metadata.extend(get_subworkflows(m['calls'][call], call))

    if 'subWorkflowMetadata' in m:
        call_metadata.extend(get_subworkflows(m['subWorkflowMetadata'], alias))

    if ('inputs' in m or 'outputs' in m) and task:
        call_metadata.append((m, task, alias))

    return call_metadata


def write_files(workflow_metadata, output_dir):
    for (m, task, alias) in workflow_metadata:
        inputs_path = os.path.join(output_dir, alias + '.inputs.json')
        outputs_path = os.path.join(output_dir, alias + '.outputs.json')
        if 'inputs' in m:
            inputs = {task + "." + key: val for key, val in m['inputs'].items() if val}
            with open(inputs_path, 'w') as f:
                f.write(json.dumps(inputs, sort_keys=True, indent=2))
        if 'outputs' in m:
            outputs = {key: val for key, val in m['outputs'].items() if val}
            with open(outputs_path, 'w') as f:
                f.write(json.dumps(outputs, sort_keys=True, indent=2))


# Main function
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("workflow_metadata",
                        help="Workflow metadata JSON file")
    parser.add_argument("output_dir", help="Output directory")
    args = parser.parse_args()

    metadata_file = args.workflow_metadata
    output_dir = args.output_dir

    metadata = json.load(open(metadata_file, 'r'))
    workflow_metadata = get_subworkflows(metadata, metadata['workflowName'])
    write_files(workflow_metadata, output_dir)


if __name__ == "__main__":
    main()
