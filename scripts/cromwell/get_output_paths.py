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
  python get_output_paths.py [-w workflows.tsv | -i workflow-id] -f filenames.json -o output.tsv 
            -b gs://bucket/workflow-name [-l LEVEL] [-s] [-k] [-e entities.txt]

Parameters:
  -w WORKFLOWS, --workflows WORKFLOWS
                        TSV file (no header) with batch (or sample) names and
                        workflow IDs (one workflow per batch). Either -i or -w
                        required. -i takes precedence if both provided.
  -i ID, --id ID        Workflow ID (alternative to -w if only one workflow).
                        Either -i or -w required. -i takes precedence if both
                        provided.
  -f FILENAMES, --filenames FILENAMES
                        JSON file with output names and filename suffixes
                        (assumes ONE file per output, not an array)
  -o OUTPUT_FILE, --output-file OUTPUT_FILE
                        Output file path
  -b BUCKET, --bucket BUCKET
                        Google bucket path to search for files - should
                        include all subdirectories preceding workflow ID
  -l LOG_LEVEL, --log-level LOG_LEVEL
                        Specify level of logging information, ie. info,
                        warning, error (not case-sensitive)
  -s, --single-file     Assume only one file per output/filename (more
                        efficient if true)
  -e ENTITIES_FILE, --entities-file ENTITIES_FILE
                        Newline-separated text file of entity names. First
                        line is entity type (ie. sample, batch). Expect one
                        output per entity for all outputs, with filename
                        containing entity ID. Output will have one line per
                        entity. If multiple batches, outputs will be concatenated.
  -k, --keep-all-entities
                        With --entities-file, output a line for every entity,
                        even if none of the output files are found

Outputs:
  - TSV file with columns for each output variable and a row for each 
    batch (or entity, if providing --entities-file), containing GCS output paths

Author: Emma Pierce-Hoffman (epierceh@broadinstitute.org)
"""


def check_file_nonempty(f):
  if not isfile(f):
    raise RuntimeError("Required input file %s does not exist." % f)
  elif getsize(f) == 0:
    raise RuntimeError("Required input file %s is empty." % f)


def read_entities_file(entities_file, single_file):
  if entities_file is None:
    return "batch", 0, None
  elif single_file:
    logging.warning("--single-file overrides --entities")
    return "batch", 0, None
  # proceed with reading file - must not be None at this point
  check_file_nonempty(entities_file)
  first = True
  entity_type = None
  entities = []
  num_entities = 0
  with open(entities_file, 'r') as f:
    for line in f:
      item = line.strip()
      if first:
        entity_type = item
        first = False
      else:
        entities.append(item)
        num_entities += 1
  return entity_type, num_entities, entities


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


def make_batch_dir_dict(workflows, workflow_id, directory):
  single_workflow = True
  batch_dirs = {}
  batches = []  # to hold batches in order given in input
  bucket, subdir = split_bucket_subdir(directory)
  if workflow_id is not None:
    return ["placeholder_batch"], {"placeholder_batch": subdir + workflow_id + "/"}, bucket, single_workflow
  count = 0
  with open(workflows, 'r') as inp:
    for line in inp:
      count += 1
      (batch, workflow) = line.strip().split('\t')
      batch_dirs[batch] = subdir + workflow + "/"
      batches.append(batch)
  if count > 1:
    single_workflow = False
  return batches, batch_dirs, bucket, single_workflow


def file_search_single(blobs, names_left, batch_outputs, files_dict, num_found, num_outputs, bucket):
  # assumes one file per output/filename suffix --> more efficient search if true
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
  return batch_outputs, num_found, names_left


def file_search_multi(blobs, names_left, batch_outputs, files_dict, num_found, bucket, output_names):
  # go through each object in bucket once, checking if it matches any filenames not yet found
  for blob in blobs:
    blob_name = blob.name.strip()
    for name in output_names:  # for multi, continue matching on suffixes even if already found file match(es)
      if blob_name.endswith(files_dict[name]):
        blob_path = "gs://" + bucket + "/" + blob_name  # reconstruct URI
        if batch_outputs[name] is None:
          num_found += 1
          names_left.remove(name)
          batch_outputs[name] = [blob_path]
        else:
          batch_outputs[name].append(blob_path)
        break
  return batch_outputs, num_found, names_left


def find_batch_output_files(batch, bucket, prefix, files_dict, output_names, num_outputs, single_file):
  batch_outputs = {file: None for file in output_names}
  num_found = 0
  storage_client = storage.Client()
  blobs = storage_client.list_blobs(bucket, prefix=prefix, delimiter=None)  # only one workflow per batch - assumes caching if multiple
  names_left = [name for name in output_names]
  
  if single_file:
    batch_outputs, num_found, names_left = file_search_single(blobs, names_left, batch_outputs, files_dict, num_found, num_outputs, bucket)
  else:
    batch_outputs, num_found, names_left = file_search_multi(blobs, names_left, batch_outputs, files_dict, num_found, bucket, output_names)
  # warn if some outputs not found
  if num_found < num_outputs:
    for name in names_left:
      logging.warning(f"{batch} output file {name} not found in gs://{bucket}/{prefix}. Outputting empty string")
      batch_outputs[name] = ""
  return batch_outputs


def sort_files_by_shard(file_list):
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


def format_multifile_line(batch, output_names, batch_outputs):
  batch_line = batch
  for name in output_names:
    file_list = batch_outputs[name]
    if file_list == "":
      batch_line += "\t"
    elif len(file_list) == 1:
      batch_line += "\t" + file_list[0]
    else:
      batch_line += '\t["' + '", "'.join(sort_files_by_shard(file_list)) + '"]'
  batch_line += "\n"
  return batch_line


def get_entity_outputs(output_names, batch_outputs, entities, entity_outputs):
  # edit entity_outputs dict in place
  for name in output_names:
    filepaths = batch_outputs[name]
    filenames = [path.split("/")[-1] for path in filepaths]
    for entity in entities:  # not efficient but should be <500 entities and filenames to search
      found = False
      for i, filename in enumerate(filenames):
        if entity in filename:
          entity_outputs[entity].append(filepaths[i])
          filenames.remove(filename)
          filepaths.remove(filepaths[i])
          found = True
          break
      if not found:
        # logging.warning(f"{entity} output file {name} not found in provided directories. Outputting empty string")
        entity_outputs[entity].append("")
  entities_left = []
  for entity in entities:
    if all(element == "" for element in entity_outputs[entity]):
      entity_outputs[entity] = []
      entities_left.append(entity)
  return entities_left


def format_entity_outputs(entity_outputs, keep_all_entities, entities):
  output_string = ""
  for entity in entities:
    if not keep_all_entities and all(element == "" for element in entity_outputs[entity]):
      # logging.info(f"No output files found for {entity} in provided directories. Omitting from output")
      continue
    output_string += entity + "\t" + "\t".join(entity_outputs[entity]) + "\n"
  return output_string


def get_output_files(batches, bucket, batch_dirs, files_dict, output_names, num_outputs, output_file, single_file, num_entities, entities, entity_type, keep_all_entities):
  entity_outputs = None
  if num_entities > 0:
    entities_left = [x for x in entities]
    entity_outputs = {entity: [] for entity in entities}
  logging.info("Writing %s" % output_file)
  with open(output_file, 'w') as out:
    out.write(entity_type + "\t" + "\t".join(output_names) + "\n")
    for batch in batches:
      logging.info("Searching for outputs for %s" % batch)
      prefix = batch_dirs[batch]
      batch_outputs = find_batch_output_files(batch, bucket, prefix, files_dict, output_names, num_outputs, single_file)
      if single_file:
        out.write(batch + "\t" + "\t".join([batch_outputs[name] for name in output_names]) + "\n")
      elif num_entities > 0:
        entities_left = get_entity_outputs(output_names, batch_outputs, entities_left, entity_outputs)
      else:
        batch_line = format_multifile_line(batch, output_names, batch_outputs)
        out.write(batch_line)
    if num_entities > 0:
      out.write(format_entity_outputs(entity_outputs, keep_all_entities, entities))
  logging.info("Done!")


# Main function
def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("-w", "--workflows", required=False, help="TSV file (no header) with batch (or sample) names and \
                                                                workflow IDs (one workflow per batch). Either -i or -w \
                                                                required. -i takes precedence if both provided.",
                      default=None)
  parser.add_argument("-i", "--id", required=False, help="Workflow ID (alternative to -w if only one workflow). \
                                                        Either -i or -w required. -i takes precedence if both given.",
                      default=None)
  parser.add_argument("-f", "--filenames", required=True, help="JSON file with output names and filename suffixes \
                                                                (assumes ONE file per output, not an array)")  # TODO
  parser.add_argument("-o", "--output-file", required=True, help="Output file path")
  parser.add_argument("-b", "--bucket", required=True, help="Google bucket path to search for files - should include \
                                                            all subdirectories preceding workflow ID")
  parser.add_argument("-l", "--log-level",
                      help="Specify level of logging information, ie. info, warning, error (not case-sensitive)",
                      required=False, default="INFO")
  parser.add_argument("-s", "--single-file", help="Assume only one file per output/filename (more efficient if true)",
                      required=False, default=False, action='store_true')
  parser.add_argument("-e", "--entities-file", help="Newline-separated text file of entity names. First line is entity \
                                              type (ie. sample, batch). Expect one output per entity for all outputs, \
                                              with filename containing entity ID. Output will have one line per entity.\
                                               If multiple batches, outputs will be concatenated.",
                      required=False)
  parser.add_argument("-k", "--keep-all-entities", help="With --entities-file, output a line for every entity, \
                                                        even if none of the output files are found",
                      required=False, default=False, action='store_true')
  args = parser.parse_args()

  log_level = args.log_level
  numeric_level = getattr(logging, log_level.upper(), None)
  if not isinstance(numeric_level, int):
    raise ValueError('Invalid log level: %s' % log_level)
  logging.basicConfig(level=numeric_level, format='%(levelname)s: %(message)s')

  filenames, output_file, bucket = args.filenames, args.output_file, args.bucket  # required
  check_file_nonempty(filenames)
  files_dict, output_names, num_outputs = load_filenames(filenames)

  workflows, workflow_id = args.workflows, args.id
  if workflow_id is None:
    if workflows is None: 
      raise ValueError("Must provide either -w or -i input")
    else:
      check_file_nonempty(workflows)
  else:
    if workflows is not None:
      logging.warning("Using workflow ID provided by --id and ignoring --workflows input")
  batches, batch_dirs, bucket, single_workflow = make_batch_dir_dict(workflows, workflow_id, bucket)
  
  single_file, entities_file, keep_all_entities = args.single_file, args.entities_file, args.keep_all_entities
  entity_type, num_entities, entities = read_entities_file(entities_file, single_file)

  get_output_files(batches, bucket, batch_dirs, files_dict, output_names, num_outputs, output_file, single_file, num_entities, entities, entity_type, keep_all_entities)


if __name__ == "__main__":
  main()
