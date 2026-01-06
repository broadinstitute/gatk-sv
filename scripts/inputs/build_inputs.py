#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os
import os.path
import glob
import json
from jinja2 import Environment, FileSystemLoader, Undefined
import logging

# This script generates input files (for example JSON inputs to be passed to cromwell, but also
# potentially tsv files for Terra import) based on the input templates in the repository and the
# set of specified input values.
#
# Any json files in the values directory are read into maps corresponding to the file name. The input json files
# should consist of single dictionaries. If an input dictionary contains a "name" key, its value will be used to report
# missing value accesses.
#
# Input value maps can then be aliased in to resource bundles accessible in the template files.  For example, if the "resources_h38"
# is aliased to "reference_resources", a template accessing the "reference_resources" bundle will pull values from the
# "resources_hg38.json" input values file.
#
# Values can be referred to by their resource bundle name + "." + attribute. For example, if the values
# directory contains a file called dockers.json containing the map { "sv_pipeline_docker" : "gatksv/sv-pipeline:tag" },
# and the "dockers.json" input file has been aliased to the "dockers" resource bundle, then in a template
# the string {{ dockers.sv_pipeline_docker }} will be replaced with the string gatksv/sv-pipeline:tag.
#
# By default the following resource bundle aliases are applied:
#
#   dockers -> dockers
#   ref_panel -> ref_panel_1kg
#   reference_resources -> resources_hg38
#   test_batch -> empty
#
# Where the empty resource bundle is just an empty map.
#
# Resource bundles can be aliased in on the command line with the -a parameter, which takes a JSON dict string. For example,
# the parameter
#
#    -a '{"test_batch" : "test_batch_small"}'
#
# Will cause the "test_batch_small" input value set to be aliased to the "test_batch" resource bundle.
#
# If a template refers to missing property from a resource bundle, it will be skipped, with an info message listing which
# properties are missing. This feature can be used purposefully to generate different sets of input files from the same sets
# of templates depending on which properties are present in the input value files. For example, the build_default_inputs.sh
# script generates inputs three times from the test_input_templates directory, with the test_batch bundle aliased to the
# small, large, and single sample test batches, respectively. Each run will generate a different set of outputs -- for example
# the run based on test_batch_small does not currently produce results for module03 and beyond because those files depend
# resources not defined for the small batch (but defined for the large batch).
#
# Any template files with a filename that begins with the string 'workspace' will be transposed before the final inputs are written.
# This functionality allows us to simplify the editing, viewing of diffs,
# and tracking of changes in Git of workspace TSV files that Terra expects in a two-row format.
#
# Jinja2 filters can be applied. For example, to ensure that string values are quoted in json files, use the tojson
# filter: {{ dockers.sv_pipeline_docker | tojson }}
# this class drops logs undefined value references in the "undefined_names" list
undefined_names = []


class TrackMissingValuesUndefined(Undefined):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # dict.__init__(self, fname=fname)
        if 'name' in self._undefined_obj:
            undefined_names.append(
                self._undefined_obj['name'] + "." + self._undefined_name)
        else:
            undefined_names.append(self._undefined_name)
        # self._fail_with_undefined_error()

    def __str__(self):
        return ""


def to_json_custom(value, *args, **kwargs):
    if isinstance(value, Undefined):
        return "UNDEFINED"
    else:
        return json.dumps(value, *args, **kwargs)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("input_values_directory",
                        help="Directory containing input value map JSON files")
    parser.add_argument(
        "template_path", help="Path to template directory or file (directories will be processed recursively)")
    parser.add_argument("output_directory",
                        help="Directory to create output files in")
    parser.add_argument('-a', '--aliases', type=json.loads,
                        default={}, help="Aliases for input value bundles")
    parser.add_argument('--log-info', action='store_true',
                        help="Show INFO-level logging messages. Use for troubleshooting.")
    args = parser.parse_args()

    # Set logger
    logging_fmt = '%(levelname)s: %(message)s'
    if args.log_info:
        logging.basicConfig(level=logging.INFO, format=logging_fmt)
    else:
        logging.basicConfig(level=logging.WARNING, format=logging_fmt)

    # prepare input values and bundle aliases
    input_directory = args.input_values_directory
    input_files = glob.glob(input_directory + "/*.json")
    raw_input_bundles = {os.path.splitext(os.path.basename(input_file))[
        0]: json.load(open(input_file, "r")) for input_file in input_files}
    raw_input_bundles['ref_panel_empty'] = {}
    raw_input_bundles['ref_panel_empty']['name'] = 'ref_panel'
    raw_input_bundles['test_batch_empty'] = {}
    raw_input_bundles['test_batch_empty']['name'] = 'test_batch'
    raw_input_bundles['single_sample_none'] = {}
    raw_input_bundles['single_sample_none']['name'] = 'single_sample'

    default_aliases = {'dockers': 'dockers',
                       'ref_panel': 'ref_panel_empty',
                       'reference_resources': 'resources_hg38',
                       'test_batch': 'test_batch_empty',
                       'single_sample': 'single_sample_none'}

    # prepare the input_dict using default, document default, and user-specified aliases
    input_dict = {}
    for alias in default_aliases:
        input_dict[alias] = raw_input_bundles[default_aliases[alias]]

    user_aliases = args.aliases
    logging.info("Using user aliases: " + str(user_aliases))
    for alias in user_aliases:
        input_dict[alias] = raw_input_bundles[user_aliases[alias]]

    template_path = args.template_path
    target_directory = args.output_directory

    if os.path.isdir(template_path):
        process_directory(input_dict, template_path, target_directory)
    else:
        process_file(input_dict, os.path.dirname(template_path),
                     os.path.basename(template_path), target_directory)


def transpose_tsv(input_str):
    # Split input string into lines and remove trailing whitespace
    lines = [line.rstrip('\n') for line in input_str.split('\n')]

    # Group the lines by their column number
    groups = {}
    for line in lines:
        if line:
            columns = line.split('\t')
            for i, column in enumerate(columns):
                if i not in groups:
                    groups[i] = []
                groups[i].append(column)

    # Transpose the groups and write the result to a new string
    transposed_groups = []
    for i in sorted(groups.keys()):
        transposed_groups.append('\t'.join(groups[i]))
    transposed_str = '\n'.join(transposed_groups) + '\n'

    return transposed_str


def process_directory(input_dict, template_dir, target_directory):
    template_dir_split = template_dir.split(os.sep)
    template_root = template_dir_split[len(template_dir_split) - 1]
    template_base = os.sep.join(
        template_dir_split[0:len(template_dir_split) - 1])
    target_dir_split = target_directory.split(os.sep)
    target_root = target_dir_split[len(target_dir_split) - 1]
    target_base = os.sep.join(target_dir_split[0:len(target_dir_split) - 1])
    for subdir, subdirList, fileList in os.walk(template_dir):
        stripped_subdir = subdir[(len(template_base) + len(os.sep)):]
        stripped_subdir = stripped_subdir[(len(template_root) + len(os.sep)):]
        if len(stripped_subdir) > 0:
            target_subdir = os.sep.join(
                [target_base, target_root, stripped_subdir])
        else:
            target_subdir = os.sep.join([target_base, target_root])
        for file in fileList:
            undefined_names.clear()
            process_file(input_dict, subdir, file, target_subdir)


def process_file(input_dict, template_subdir, template_file, target_subdir):
    template_file_path = os.sep.join([template_subdir, template_file])

    # only process files that end with .tmpl
    if not template_file.endswith(".tmpl"):
        logging.warning("skipping file " + template_file_path +
                        " because it does not have .tmpl extension")
        return

    target_file = template_file.rsplit('.', 1)[0]
    target_file_path = os.sep.join([target_subdir, target_file])
    env = Environment(loader=FileSystemLoader(template_subdir),
                      undefined=TrackMissingValuesUndefined)
    env.policies['json.dumps_function'] = to_json_custom
    logging.info(template_file_path + " -> " + target_file_path)
    processed_content = env.get_template(template_file).render(input_dict)
    # Check if the template_file starts with "workspace"
    if template_file.startswith("workspace"):
        # Transpose the TSV data in processed_content
        processed_content = transpose_tsv(processed_content)
    if len(undefined_names) > 0:
        logging.info("skipping file " + template_file_path +
                     " due to missing values " + str(undefined_names))
    else:
        os.makedirs(target_subdir, exist_ok=True)
        target_file = open(target_file_path, "w")
        target_file.write(processed_content)
        target_file.close()


if __name__ == "__main__":
    main()
