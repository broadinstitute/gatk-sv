#!/bin/python

import argparse
import json
import logging
import subprocess
from os import listdir
import os.path

"""
Summary: Cursory validation of Terra input JSONs for pipeline workflows.
    Checks if JSONs contain all required inputs and that
    they do not contain extraneous inputs. Does not perform type-checking.
    Also checks that all expected JSONs have been generated.

Usage:
    python scripts/test/terra_validation.py -d /path/to/base/dir -j /path/to/womtool/jar [optional inputs]

Parameters:
    /path/to/base/dir: path to base directory of gatk-sv repo
    /path/to/womtool/jar: path to user's womtool jar file
    Optional inputs:
    -n,--num-input-jsons INT: override default expected number of Terra input
        JSONs
    --log-level LEVEL: specify level of logging information to print,
        ie. INFO, WARNING, ERROR - not case-sensitive)

Outputs: If successful, last line of printout should read
    "<N> of <N> Terra input JSONs exist and passed validation."
    Prior lines will detail the JSONs that were examined and any errors found.
"""

WDLS_PATH = "wdl/"
TERRA_INPUTS_PATH = "inputs/terra_workspaces/cohort_mode/workflow_configurations/"


def list_jsons(inputs_path, expected_num_jsons, subdir="", description=""):
    jsons = [os.path.join(subdir, x) for x in listdir(os.path.join(inputs_path, subdir)) if x.endswith(".json")]
    num_input_jsons = len(jsons)
    if num_input_jsons < expected_num_jsons:
        raise Exception(f"Expected {expected_num_jsons} Terra {description}input JSONs but found {num_input_jsons}.")
    jsons.sort()
    return jsons


def get_wdl_json_pairs(wdl_path, terra_inputs_path, expected_num_inputs):
    jsons = list_jsons(terra_inputs_path, expected_num_inputs)

    for json_file in jsons:
        path_to_wdl = os.path.join(wdl_path, os.path.basename(json_file).split(".")[0] + ".wdl")
        if os.path.isfile(path_to_wdl):
            yield path_to_wdl, os.path.join(terra_inputs_path, json_file)
        else:
            logging.warning(f"Can't find WDL corresponding to {os.path.basename(json_file)} at {path_to_wdl}.")


def validate_terra_json(wdl, terra_json, womtool_jar):
    womtool_command = "java -jar " + womtool_jar + " inputs " + wdl
    womtool_json = subprocess.run(womtool_command, shell=True, stdout=subprocess.PIPE)

    with open(terra_json, 'r') as terra_json_file:
        womtool_inputs = json.loads(womtool_json.stdout)
        terra_inputs = json.load(terra_json_file)
        wdl_name = os.path.basename(wdl)

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

    return int(valid)  # return 1 if valid, 0 if not valid


def validate_all_terra_jsons(base_dir, womtool_jar, expected_num_inputs):
    successes = 0
    for wdl, json_file in get_wdl_json_pairs(os.path.join(base_dir, WDLS_PATH),
                                             os.path.join(base_dir, TERRA_INPUTS_PATH),
                                             expected_num_inputs):
        successes += validate_terra_json(wdl, json_file, womtool_jar)

    print("\n")
    if successes != expected_num_inputs:
        raise RuntimeError(f"{successes} of {expected_num_inputs} Terra input JSONs passed validation.")
    else:
        print(f"{successes} of {expected_num_inputs} Terra input JSONs passed validation.")


# Main function
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--base-dir", help="Relative path to base of gatk-sv repo", required=True)
    parser.add_argument("-j", "--womtool-jar", help="Path to womtool jar", required=True)
    parser.add_argument("-n", "--num-input-jsons",
                        help="Number of Terra input JSONs expected",
                        required=False, default=19, type=int)
    parser.add_argument("--log-level",
                        help="Specify level of logging information, ie. info, warning, error (not case-sensitive)",
                        required=False, default="INFO")
    args = parser.parse_args()

    # get args as variables
    base_dir, womtool_jar, log_level = args.base_dir, args.womtool_jar, args.log_level
    expected_num_inputs = args.num_input_jsons

    numeric_level = getattr(logging, log_level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % log_level)
    logging.basicConfig(level=numeric_level, format='%(levelname)s: %(message)s')

    validate_all_terra_jsons(base_dir, womtool_jar, expected_num_inputs)


if __name__ == "__main__":
    main()
