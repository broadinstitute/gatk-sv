#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse, os, os.path
import glob
import json
import jinja2

# This script generates input files (for example JSON inputs to be passed to cromwell, but also
# potentially tsv files for Terra import) based on the input templates in the repository and the
# set of specified input values.
#
# Any json files in the input_values directory are read into maps corresponding to the file name.
#
# The script then finds any template files in the directories listed in the keys to the template_targets map
# and renders them into output files using the jinja2 templating library.
#
# Values can be referred to by their input file name + "." + input file key. For example, if the input_values
# directory contains a file called dockers.json containing the map { "sv_pipeline_docker" : "gatksv/sv-pipeline:tag" },
# the string {{ dockers.sv_pipeline_docker }} will be replaced with the string gatksv/sv-pipeline:tag.
#
# Jinja2 filters can be applied. For example, to ensure that string values are quoted in json files, use the tojson
# filter: {{ dockers.sv_pipeline_docker | tojson }}

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    args = parser.parse_args()

    input_directory = 'input_values'
    template_targets = { 'input_templates' : 'inputs', 'test_input_templates' : 'test'}

    input_files = glob.glob(input_directory + "/*.json")
    input_dict = {os.path.splitext(os.path.basename(input_file))[0]:json.load(open(input_file, "r")) for input_file in input_files}

    for template_dir in template_targets.keys():
        for subdir, subdirList, fileList in os.walk(template_dir):
            target_root_dir = template_targets[template_dir]
            split_path = subdir.split(os.sep)
            split_path[0] = target_root_dir
            target_subdir = os.sep.join(split_path)
            if not os.path.isdir(target_subdir):
                os.mkdir(target_subdir)
            for file in fileList:
                template_file_path = os.sep.join([subdir, file])
                target_file_path = os.sep.join([target_subdir, file])
                print(template_file_path + " -> " + target_file_path)
                template_file = open(template_file_path, "r")
                target_file = open(target_file_path, "w")
                template_string = template_file.read()
                processed_content = jinja2.Template(template_string).render(input_dict)
                target_file.write(processed_content)
                target_file.close()
                template_file.close()

if __name__ == "__main__":
    main()