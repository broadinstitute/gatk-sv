#!/bin/python

import subprocess
import argparse
import json
from os import listdir
from os.path import basename, isfile
import logging
import atexit


"""
Summary: Cursory validation of Terra input JSONs for pipeline workflows and metrics workflows. 
  Checks if JSONs contain all required inputs and that they do not contain extraneous inputs.
  Does not perform type-checking. Also checks that all expected JSONs have been generated.

Usage:
  python scripts/test/terra_validation.py -d /path/to/base/dir -j /path/to/womtool/jar [optional flags]

Parameters:
  /path/to/base/dir: path to base directory of gatk-sv repo
  /path/to/womtool/jar: path to user's womtool jar file
  Optional flags:
  --log-level LEVEL (specify level of logging information to print, ie. INFO, WARNING, ERROR - not case-sensitive)

Outputs: If successful, last line of printout should read "19 of 19 Terra input JSONs exist and passed validation."
  Prior lines will detail the JSONs that were examined and any errors found.
"""


TMP_JSON = "womtool_inputs.tmp.json"
WDLS_PATH = "/wdl/"
TERRA_INPUTS_PATH = "/inputs/terra_workspaces/cohort_mode/workflow_configurations/"
METRICS_INPUTS_PATH = "metrics_workflows/"
NUM_TERRA_INPUT_JSONS = 12
NUM_TERRA_METRICS_JSONS = 7


def exit_handler():
  if isfile(TMP_JSON):
    subprocess.run("rm " + TMP_JSON, shell=True)


def list_jsons(inputs_path):
  jsons = [x for x in listdir(inputs_path) if x.endswith(".json")]
  num_input_jsons = len(jsons)
  if num_input_jsons < NUM_TERRA_INPUT_JSONS:
    logging.warning(f"Expected {NUM_TERRA_INPUT_JSONS} Terra input JSONs but found {num_input_jsons}.")
  jsons.sort()

  metrics_jsons = [METRICS_INPUTS_PATH + x for x in listdir(inputs_path + METRICS_INPUTS_PATH) if x.endswith(".json")]
  num_metrics_jsons = len(metrics_jsons)
  if num_metrics_jsons < NUM_TERRA_METRICS_JSONS:
    logging.warning(f"Expected {NUM_TERRA_METRICS_JSONS} Terra metrics input JSONs but found {num_metrics_jsons}.")
  metrics_jsons.sort()

  jsons.extend(metrics_jsons)
  return jsons


def get_wdl_json_pairs(wdl_path, terra_inputs_path):
  jsons = list_jsons(terra_inputs_path)

  for json_file in jsons:
    path_to_wdl = wdl_path + basename(json_file)[:-5] + ".wdl"
    if isfile(path_to_wdl):
      yield path_to_wdl, terra_inputs_path + json_file
    else:
      logging.warning(f"Can't find WDL corresponding to {basename(json_file)} at {path_to_wdl}.")


def validate_terra_json(wdl, terra_json, womtool_jar):
  womtool_input_generation_command = "java -jar " + womtool_jar + " inputs " + wdl + " > " + TMP_JSON
  subprocess.run(womtool_input_generation_command, shell=True)

  with open(TMP_JSON, 'r') as womtool_json_file, open(terra_json, 'r') as terra_json_file:
    womtool_inputs = json.load(womtool_json_file)
    terra_inputs = json.load(terra_json_file)
    wdl_name = basename(wdl)

    valid = True
    for inp in womtool_inputs:
      if "optional" in womtool_inputs[inp]:
        continue
      elif inp not in terra_inputs:
        logging.error(f"Missing input: Required input {inp} for {wdl_name} missing from {terra_json}")
        valid = False

    for inp in terra_inputs:
      if inp not in womtool_inputs:
        logging.error(f"Unexpected input: {terra_json} contains unexpected input {inp} for {wdl_name}")
        valid = False

  if valid:
    logging.info(f"PASS: {terra_json} is a valid Terra input JSON for {wdl_name}")

  subprocess.run("rm " + TMP_JSON, shell=True)
  return int(valid)  # return 1 if valid, 0 if not valid


def validate_all_terra_jsons(base_dir, womtool_jar):
  successes = 0
  for wdl, json_file in get_wdl_json_pairs(base_dir + WDLS_PATH, base_dir + TERRA_INPUTS_PATH):
    successes += validate_terra_json(wdl, json_file, womtool_jar)

  print("\n")
  total_jsons = NUM_TERRA_METRICS_JSONS + NUM_TERRA_INPUT_JSONS
  print(f"{successes} of {total_jsons} Terra input JSONs exist and passed validation.")


# Main function
def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("-d", "--base-dir", help="Relative path to base of gatk-sv repo", required=True)
  parser.add_argument("-j", "--womtool-jar", help="Path to womtool jar", required=True)
  parser.add_argument("--log-level",
                      help="Specify level of logging information, ie. info, warning, error (not case-sensitive)",
                      required=False, default="INFO")
  args = parser.parse_args()

  atexit.register(exit_handler)

  # get args as variables
  base_dir, womtool_jar, log_level = args.base_dir, args.womtool_jar, args.log_level

  numeric_level = getattr(logging, log_level.upper(), None)
  if not isinstance(numeric_level, int):
    raise ValueError('Invalid log level: %s' % log_level)
  logging.basicConfig(level=numeric_level, format='%(levelname)s: %(message)s')

  validate_all_terra_jsons(base_dir, womtool_jar)


if __name__ == "__main__":
  main()
