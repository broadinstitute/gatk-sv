#!/bin/python

import argparse
from os.path import isfile, getsize
import logging
import json
from google.cloud import storage
import re


"""
Summary: find GCS paths for specified outputs for multiple batches without downloading metadata

Usage:
  python get_output_paths.py -w workflows.tsv -f filenames.json -o output.tsv -b gs://bucket/workflow-name [-l LEVEL]

Parameters:
  Required:
  -w,--workflows: TSV file (no header) with batch (or sample) names and workflow IDs (one per batch)
  -f,--filenames: JSON file with output names and filename suffixes (assumes ONE file per output, not an array, string, or int)
  -o,--output_file: Output file path (TSV)
  -b,--bucket: Google bucket path to search for files (common to all workflows, should include all subdirectories preceding workflow ID)
  Optional:
  -l,--log-level LEVEL (specify level of logging information to print, ie. INFO, WARNING, ERROR - not case-sensitive)

Outputs:
  - TSV file with columns for batch and each output variable, a row for each batch, containing GCS output paths

Author: Emma Pierce-Hoffman (epierceh@broadinstitute.org)
"""


def check_file_nonempty(f):
  if not isfile(f):
    raise RuntimeError("Required input file %s does not exist." % f)
  elif getsize(f) == 0:
    raise RuntimeError("Required input file %s is empty." % f)


def load_filenames(filenames):
  files_dict = json.load(open(filenames, 'r'))
  output_names = sorted(files_dict.keys())
  num_outputs = len(output_names)
  return files_dict, output_names, num_outputs


def split_bucket_subdir(directory):
  regex = r'^(gs://)?([^/]+)(/)?(.*)'
  bucket = re.match(regex, directory).group(2)
  subdir = re.match(regex, directory).group(4)
  if subdir[-1] != '/':
    subdir += "/"
  return bucket, subdir


def make_batch_dir_dict(workflows, directory):
  batch_dirs = {}
  batches = []  # to hold batches in order given in input
  bucket, subdir = split_bucket_subdir(directory)
  with open(workflows, 'r') as inp:
    for line in inp:
      (batch, workflow_id) = line.strip().split('\t')
      batch_dirs[batch] = subdir + workflow_id + "/"
      batches.append(batch)
  return batches, batch_dirs, bucket


def find_batch_output_files(batch, bucket, prefix, files_dict, output_names, num_outputs):
  batch_outputs = {file: None for file in output_names}
  num_found = 0
  storage_client = storage.Client()
  blobs = storage_client.list_blobs(bucket, prefix=prefix, delimiter=None)  # only one workflow per batch - assumes caching if multiple
  names_left = [name for name in output_names]
  # go through each object in bucket once, checking if it matches any filenames not yet found
  for blob in blobs:
    blob_name = blob.name.strip()
    for name in names_left:
      if batch_outputs[name] is None and blob_name.endswith(files_dict[name]):
        num_found += 1
        batch_outputs[name] = "gs://" + bucket + "/" + blob_name  # reconstruct URI
        names_left.remove(name)  # does not handle array outputs - don't search for this filename again
        break
    if num_found >= num_outputs:  # stop search early if all outputs found
      break
  # warn if some outputs not found
  if num_found < num_outputs:
    for name in output_names:
      if batch_outputs[name] is None:
        logging.warning(f"{batch} output file {name} not found in provided directories. Outputting empty string")
        batch_outputs[name] = ""
  return batch_outputs


def get_output_files(batches, bucket, batch_dirs, files_dict, output_names, num_outputs, output_file):
  logging.info("Writing %s" % output_file)
  with open(output_file, 'w') as out:
    out.write("batch\t" + "\t".join(output_names) + "\n")
    for batch in batches:
      logging.info("Searching for outputs for %s" % batch)
      prefix = batch_dirs[batch]
      batch_outputs = find_batch_output_files(batch, bucket, prefix, files_dict, output_names, num_outputs)
      out.write(batch + "\t" + "\t".join([batch_outputs[name] for name in output_names]) + "\n")
  logging.info("Done!")


# Main function
def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("-w", "--workflows", required=True, help="TSV file (no header) with batch (or sample) names and \
                                                                workflow IDs (one workflow per batch)")
  parser.add_argument("-f", "--filenames", required=True, help="JSON file with output names and filename suffixes \
                                                                (assumes ONE file per output, not an array)")  # TODO
  parser.add_argument("-o", "--output_file", required=True, help="Output file path")
  parser.add_argument("-b", "--bucket", required=True, help="Google bucket path to search for files - should include \
                                                            all subdirectories preceding workflow ID")
  parser.add_argument("-l", "--log-level",
                      help="Specify level of logging information, ie. info, warning, error (not case-sensitive)",
                      required=False, default="INFO")
  args = parser.parse_args()

  log_level = args.log_level
  numeric_level = getattr(logging, log_level.upper(), None)
  if not isinstance(numeric_level, int):
    raise ValueError('Invalid log level: %s' % log_level)
  logging.basicConfig(level=numeric_level, format='%(levelname)s: %(message)s')

  workflows, filenames, output_file, bucket = args.workflows, args.filenames, args.output_file, args.bucket
  check_file_nonempty(workflows)
  check_file_nonempty(filenames)

  batches, batch_dirs, bucket = make_batch_dir_dict(workflows, bucket)
  files_dict, output_names, num_outputs = load_filenames(filenames)

  get_output_files(batches, bucket, batch_dirs, files_dict, output_names, num_outputs, output_file)


if __name__ == "__main__":
  main()
