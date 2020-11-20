#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse, os, os.path
import glob
import json
from jinja2 import Environment, Template, FileSystemLoader, meta, Undefined

# This script generates input files (for example JSON inputs to be passed to cromwell, but also
# potentially tsv files for Terra import) based on the input templates in the repository and the
# set of specified input values.
#
# Any json files in the input_values directory are read into maps corresponding to the file name. The input json files
# should consist of single dictionaries. If an input dictionary contains a "name" key, its value will be used to report
# missing value accesses.
#
# Input value maps can then be aliased in to resource bundles accessible in the template files.  For example, if the "resources_h38"
# is aliased to "reference_resources", a template accessing the "reference_resources" bundle will pull values from the
# "resources_hg38.json" input values file.
#
# Values can be referred to by their resource bundle name + "." + attribute. For example, if the input_values
# directory contains a file called dockers.json containing the map { "sv_pipeline_docker" : "gatksv/sv-pipeline:tag" },
# and the "dockers.json" input file has been aliased to the "dockers" resource bundle, then in a template
# the string {{ dockers.sv_pipeline_docker }} will be replaced with the string gatksv/sv-pipeline:tag.
#
# By default the following resource bundle aliases are applied:
#
#   dockers -> dockers
#   ref_panel -> ref_panel_v1b
#   reference_resources -> resources_hg38
#   test_batch -> empty
#   test_single_sample -> empty
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
# If a template refers to missing property from a resource bundle, it will be skipped, with a warning message listing which
# properties are missing. This feature can be used purposefully to generate different sets of input files from the same sets
# of templates depending on which properties are present in the input value files. For example, the build_default_inputs.sh
# script generates inputs three times from the test_input_templates directory, with the test_batch bundle aliased to the
# small, large, and single sample test batches, respectively. Each run will generate a different set of outputs -- for example
# the run based on test_batch_small does not currently produce results for module03 and beyond because those files depend
# resources not defined for the small batch (but defined for the large batch).
#
# Jinja2 filters can be applied. For example, to ensure that string values are quoted in json files, use the tojson
# filter: {{ dockers.sv_pipeline_docker | tojson }}

# this class drops logs undefined value references in the "undefined_names" list
undefined_names = []
class TrackMissingValuesUndefined(Undefined):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        #dict.__init__(self, fname=fname)
        if 'name' in self._undefined_obj:
            undefined_names.append(self._undefined_obj['name'] + "." + self._undefined_name)
        else:
            undefined_names.append(self._undefined_name)
        #self._fail_with_undefined_error()

    def __str__(self):
        return ""

def to_json_custom(value, *args, **kwargs):
    if isinstance(value, Undefined):
        return "UNDEFINED"
    else:
        return json.dumps(value, *args, **kwargs)

def make_target_subdir(split_path):
    for i in [i+1 for i in range(0,len(split_path))]:
        target_subdir = os.sep.join(split_path[0:i])
        if not os.path.isdir(target_subdir):
            os.mkdir(target_subdir)

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("template_directory")
    parser.add_argument("output_directory")
    parser.add_argument('-a', '--aliases', type=json.loads, default={})
    args = parser.parse_args()

    template_dir = args.template_directory
    target_directory = args.output_directory

    input_directory = 'input_values'

    input_files = glob.glob(input_directory + "/*.json")
    raw_input_bundles = {os.path.splitext(os.path.basename(input_file))[0]:json.load(open(input_file, "r")) for input_file in input_files}
    raw_input_bundles['empty'] = {}

    default_aliases = { 'dockers' : 'dockers',
                        'ref_panel' : 'ref_panel_v1b',
                        'reference_resources' : 'resources_hg38',
                        'test_batch' : 'empty',
                        'test_single_sample' : 'empty' }

    # prepare the input_dict using default, document default, and user-specified aliases
    input_dict = {}
    for alias in default_aliases:
        input_dict[alias] = raw_input_bundles[default_aliases[alias]]

    user_aliases = args.aliases
    print("Using user aliases: " + str(user_aliases))

    for alias in user_aliases:
        input_dict[alias] = raw_input_bundles[user_aliases[alias]]

    for subdir, subdirList, fileList in os.walk(template_dir):
        split_path = subdir.split(os.sep)
        split_path[0] = target_directory
        target_subdir = os.sep.join(split_path)
        for file in fileList:
            undefined_names.clear()

            # only process files that end with .tmpl
            if not file.endswith(".tmpl"):
                continue
            target_file = file.rsplit('.', 1)[0]
            env = Environment(loader=FileSystemLoader(subdir), undefined=TrackMissingValuesUndefined)
            env.policies['json.dumps_function'] = to_json_custom

            template_file_path = os.sep.join([subdir, file])
            target_file_path = os.sep.join([target_subdir, target_file])
            print(template_file_path + " -> " + target_file_path)

            processed_content = env.get_template(file).render(input_dict)
            if len(undefined_names) > 0:
                print("WARNING: skipping file " + template_file_path + " due to missing values " + str(undefined_names))
            else:
                make_target_subdir(split_path)
                target_file = open(target_file_path, "w")
                target_file.write(processed_content)
                target_file.close()


if __name__ == "__main__":
    main()