#!/bin/python

import argparse
from os.path import isfile, getsize
import logging
import subprocess
import json


"""
Summary: find GCS paths for specified outputs for multiple batches without downloading metadata

Usage:
  python get_output_paths.py -w workflows.tsv -f filenames.json -o output.tsv -b gs://bucket/ [-l LEVEL]

Parameters:
  Required:
  -w,--workflows: TSV file (no header) with batch (or sample) names and workflow IDs (comma-separated if multiple)
  -f,--filenames: JSON file with output names and filename suffixes (assumes ONE file per output, not an array, string, or int)
  -o,--output_file: Output file path (TSV)
  -b,--bucket: Google bucket path to search for files (common to all workflows)
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


def make_batch_dir_dict(workflows, bucket):
  batch_dirs = {}
  batches = [] # to hold batches in order given in input
  if bucket[-1] != '/':
    bucket += "/"
  with open(workflows, 'r') as inp:
    for line in inp:
      (batch, ids) = line.strip().split('\t')
      ids_list = ids.split(',')
      batch_dirs[batch] = [bucket + "*/" + workflow_id + "/**" for workflow_id in ids_list]
      batches.append(batch)
  return batches, batch_dirs


def find_output_file(batch, paths, filename):
  for path in paths:
    path_search = path + filename
    command = "gsutil ls " + path_search
    cmd_obj = subprocess.run(command, shell=True, capture_output=True, text=True)
    # return code is 0 if found, 1 if not
    rc = cmd_obj.returncode
    if rc == 0:
      file = cmd_obj.stdout.strip()  # will not be empty string because return code would be 1
      return file
    elif rc == 1:
      continue  # object not found, try next path
    else:  # do want to warn of other error besides not found
      logging.warning("Command gsutil ls " + path_search + " exited with code " + str(cmd_obj.returncode))
  # if get to this point, have not found file in any of the possible workflow directories for the batch
  logging.warning(f"{batch} output file {filename} not found in provided directories. Outputting empty string")
  return ""


def get_output_files(batches, batch_dirs, files_dict, output_names, num_outputs, output_file):
  logging.info("Writing %s" % output_file)
  with open(output_file, 'w') as out:
    out.write("batch\t" + "\t".join(output_names) + "\n")
    for batch in batches:
      logging.info("Searching for outputs for %s" % batch)
      paths = batch_dirs[batch]
      batch_line = [""] * (num_outputs + 1)
      batch_line[0] = batch
      for i, output_name in enumerate(output_names):
        batch_line[i + 1] = find_output_file(batch, paths, files_dict[output_name])
      out.write("\t".join(batch_line) + "\n")
  logging.info("Done!")


# Main function
def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("-w", "--workflows", required=True, help="TSV file (no header) with batch (or sample) names and \
                                                                workflow IDs")
  parser.add_argument("-f", "--filenames", required=True, help="JSON file with output names and filename suffixes \
                                                                (assumes ONE file per output, not an array)")  # TODO
  parser.add_argument("-o", "--output_file", required=True, help="Output file path")
  parser.add_argument("-b", "--bucket", required=True, help="Google bucket path to search for files")
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

  batches, batch_dirs = make_batch_dir_dict(workflows, bucket)
  files_dict = json.load(open(filenames, 'r'))
  output_names = sorted(files_dict.keys())
  num_outputs = len(output_names)

  get_output_files(batches, batch_dirs, files_dict, output_names, num_outputs, output_file)


if __name__ == "__main__":
  main()
