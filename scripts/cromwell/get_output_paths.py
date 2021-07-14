#!/bin/python

import argparse
import json
import logging
import re
import os.path
from urllib.parse import urlparse

from google.cloud import storage

"""
Summary: Find GCS paths for specified workflow file outputs for multiple workflows at once without downloading metadata.

Caveats: Assumes cromwell file structure. Recommended for use with cromwell final_workflow_outputs_dir
    to reduce number of files to search. Requires file suffixes for each output file that are
    unique within the workflow directory.

For usage & parameters: Run python get_output_paths.py --help

Output: TSV file with columns for each output variable and a row for each
    batch (or entity, if providing --entities-file), containing GCS output paths

Author: Emma Pierce-Hoffman (epierceh@broadinstitute.org)
"""


def check_file_nonempty(f):
    # Validate existence of file and that it is > 0 bytes
    if not os.path.isfile(f):
        raise RuntimeError("Required input file %s does not exist." % f)
    elif os.path.getsize(f) == 0:
        raise RuntimeError("Required input file %s is empty." % f)


def read_entities_file(entities_file):
    # Get list of entities from -e entities file
    entities = []
    if entities_file is not None:
        # proceed with reading file - must not be None at this point
        check_file_nonempty(entities_file)
        with open(entities_file, 'r') as f:
            for line in f:
                entities.append(line.strip())
    return entities


def load_filenames(filenames):
    # Read -f filenames / output names JSON
    files_dict = json.load(open(filenames, 'r'))
    output_names = sorted(files_dict.keys())
    if len(output_names) == 0:
        raise ValueError("No output files to search for found in required -f/--filenames JSON %s." % filenames)
    return files_dict, output_names


def split_bucket_subdir(directory):
    # Parse -b URI input into top-level bucket name (no gs://) and subdirectory path
    uri = urlparse(directory)
    return uri.netloc, uri.path.lstrip("/")


def get_batch_dirs(workflows, workflow_id, directory):
    # Return list of (batch_name, batch_subdirectory) and top-level bucket parsed from -b URI input
    batches_dirs = []  # to hold tuples of (batch, dir) in order given in input
    bucket, subdir = split_bucket_subdir(directory)
    # If using -i input, just add workflow ID to subdirectory path and return
    if workflow_id is not None:
        return [("placeholder_batch", os.path.join(subdir, workflow_id))], bucket
    # If using -w input, read workflows file to get batch names and workflow IDs
    with open(workflows, 'r') as inp:
        for line in inp:
            if line.strip() == "":
                continue
            (batch, workflow) = line.strip().split('\t')
            batch_dir = os.path.join(subdir, workflow)
            batches_dirs.append((batch, batch_dir))
    return batches_dirs, bucket


def find_batch_output_files(batch, bucket, prefix, files_dict, output_names, num_outputs):
    # Search batch directory for files with specified prefixes

    # Get all objects in directory
    storage_client = storage.Client()
    blobs = storage_client.list_blobs(bucket, prefix=prefix,
                                      delimiter=None)  # only one workflow per batch - assumes caching if multiple

    # Go through each object in directory once, checking if it matches any filenames not yet found
    batch_outputs = {file: [] for file in output_names}
    names_left = list(output_names)
    num_found = 0
    for blob in blobs:
        blob_name = blob.name.strip()
        # in case multiple files, continue matching on suffixes even if already found file match(es)
        for name in output_names:
            if blob_name.endswith(files_dict[name]):
                blob_path = os.path.join("gs://", bucket, blob_name)  # reconstruct URI
                if len(batch_outputs[name]) == 0:
                    num_found += 1
                    names_left.remove(name)
                batch_outputs[name].append(blob_path)
                break

    # Warn if some outputs not found
    if num_found < num_outputs:
        for name in names_left:
            logging.warning(f"{batch} output file {name} not found in gs://{bucket}/{prefix}. Outputting empty string")

    return batch_outputs


def sort_files_by_shard(file_list):
    # Attempt to sort file list by shard number based on last occurrence of "shard-" in URI
    if len(file_list) < 2:
        return file_list
    regex = r'^(shard-)([0-9]+)(/.*)'  # extract shard number for sorting - group 2
    shard_numbers = []
    check_different_shard = None
    for file in file_list:
        index = file.rfind("shard-")  # find index of last occurrence of shard- substring in file path
        if index == -1:
            return file_list  # abandon sorting if no shard- substring
        shard = int(re.match(regex, file[index:]).group(2))
        # make sure first two shard numbers actually differ
        if check_different_shard is None:
            check_different_shard = shard
        elif check_different_shard != -1:
            if shard == check_different_shard:
                return file_list  # if first two shard numbers match, then abandon sorting by shard
            check_different_shard = -1
        shard_numbers.append(shard)
    return [x for _, x in sorted(zip(shard_numbers, file_list), key=lambda pair: pair[0])]


def format_batch_line(batch, output_names, batch_outputs):
    # Format line with batch and outputs (if not using entities option)
    batch_line = batch + "\t"
    batch_line += "\t".join(",".join(sort_files_by_shard(batch_outputs[name])) for name in output_names)
    batch_line += "\n"
    return batch_line


def update_entity_outputs(output_names, batch_outputs, entities, entity_outputs):
    # Edit entity_outputs dict in place: add new batch outputs to each corresponding entity
    for output_index, name in enumerate(output_names):
        filepaths = batch_outputs[name]
        filenames = [path.split("/")[-1] for path in filepaths]
        for entity in entities:  # not efficient but should be <500 entities and filenames to search
            for i, filename in enumerate(filenames):
                # cannot handle Array[File] output for one entity
                if entity in filename and entity_outputs[entity][output_index] == "":
                    entity_outputs[entity][output_index] = filepaths[i]
                    entity_outputs[entity].append(filepaths[i])
                    filenames.remove(filename)
                    filepaths.remove(filepaths[i])
                    break


def write_entity_outputs(entity_outputs, keep_all_entities, entities, output_stream):
    # Check, format, and write entity outputs
    # do write inside function to be able to print line-by-line
    for entity in entities:
        # check for blank entities
        if all(element == "" for element in entity_outputs[entity]):
            if keep_all_entities:
                logging.info(f"No output files found for entity '{entity}' in provided directories. "
                             f"Outputting blank entry. Remove -k argument to exclude empty entities.")
            else:
                logging.info(f"No output files found for entity '{entity}' in provided directories. "
                             f"Omitting from output. Use -k argument to include empty entities.")
                continue
        output_stream.write(entity + "\t" + "\t".join(entity_outputs[entity]) + "\n")


def retrieve_and_write_output_files(batches_dirs, bucket, files_dict, output_names, output_file,
                                    entities, entity_type, keep_all_entities):
    num_outputs = len(output_names)
    num_entities = len(entities)
    entity_outputs = {entity: [""] * num_outputs for entity in entities}  # empty if entities is empty
    logging.info("Writing %s" % output_file)
    with open(output_file, 'w') as out:
        out.write(entity_type + "\t" + "\t".join(output_names) + "\n")
        for batch, batch_dir in batches_dirs:
            logging.info("Searching for outputs for %s" % batch)
            batch_outputs = find_batch_output_files(batch, bucket, batch_dir, files_dict, output_names, num_outputs)
            if num_entities > 0:
                update_entity_outputs(output_names, batch_outputs, entities, entity_outputs)
            else:
                batch_line = format_batch_line(batch, output_names, batch_outputs)
                out.write(batch_line)
        if num_entities > 0:
            write_entity_outputs(entity_outputs, keep_all_entities, entities, out)
    logging.info("Done!")


# Main function
def main():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-w", "--workflows-file",
                       help="TSV file (no header) with batch (or sample) names and workflow IDs (one workflow "
                            "per batch). Either -i or -w required.")
    group.add_argument("-i", "--workflow-id",
                       help="Workflow ID provided directly on the command line; alternative to -w if only "
                            "one workflow. Either -i or -w required.")
    parser.add_argument("-f", "--filenames", required=True,
                        help="JSON file with workflow output file names (for column names in output TSV) and a "
                             "unique filename suffix expected for each workflow output. "
                             "Format is { \"output_file_name\": \"unique_file_suffix\" }.")
    parser.add_argument("-o", "--output-file", required=True, help="Output file path to create")
    parser.add_argument("-b", "--bucket", required=True,
                        help="Google bucket path to search for files - should include all subdirectories "
                             "preceding the workflow ID, including the workflow name.")
    parser.add_argument("-l", "--log-level", required=False, default="INFO",
                        help="Specify level of logging information, ie. info, warning, error (not case-sensitive). "
                             "Default: INFO")
    parser.add_argument("-e", "--entities-file", required=False,
                        help="Newline-separated text file of entity (ie. sample, batch) names (no header). "
                             "Entity here refers to units, like samples within a batch or batches within a cohort, "
                             "for which the workflow(s) produced outputs; the script expects one output per entity "
                             "for all outputs, with the filename containing the entity ID provided in the entities "
                             "file. Output will have one line per entity in the order provided. "
                             "If multiple batches, outputs will be concatenated and order may be affected.")
    parser.add_argument("-t", "--entity-type", required=False, default="batch",
                        help="Entity type (ie. sample, batch) of each line of output. If using -e, then define "
                             "what each entity name in the file is (ie. a sample, a batch). Otherwise, define "
                             "what each workflow corresponds to. This type will be the first column name. "
                             "Default: batch")
    parser.add_argument("-k", "--keep-all-entities", required=False, default=False, action='store_true',
                        help="With --entities-file, output a line for every entity, even if none of the "
                             "output files are found.")
    args = parser.parse_args()

    # Set logging level from -l input
    log_level = args.log_level
    numeric_level = getattr(logging, log_level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % log_level)
    logging.basicConfig(level=numeric_level, format='%(levelname)s: %(message)s')

    # Set required arguments. Validate existence of & read filenames JSON
    filenames, output_file, bucket = args.filenames, args.output_file, args.bucket  # required
    check_file_nonempty(filenames)
    files_dict, output_names = load_filenames(filenames)

    # Determine workflow IDs from -w or -i arguments. Get subdirectories
    workflows, workflow_id = args.workflows_file, args.workflow_id
    if workflows is not None:
        check_file_nonempty(workflows)
    batches_dirs, bucket = get_batch_dirs(workflows, workflow_id, bucket)

    # Set entity arguments and read entities file
    entity_type, entities_file, keep_all_entities = args.entity_type, args.entities_file, args.keep_all_entities
    entities = read_entities_file(entities_file)

    # Core functionality
    retrieve_and_write_output_files(batches_dirs, bucket, files_dict, output_names, output_file,
                                    entities, entity_type, keep_all_entities)


if __name__ == "__main__":
    main()
