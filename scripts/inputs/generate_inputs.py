#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse, os, os.path
import glob
import json
from jinja2 import Template

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input_directory', help='Directory containing parameter JSON files')
    parser.add_argument('template', help='Template File', type=argparse.FileType('r'))
    parser.add_argument('fout', help='Output File', type=argparse.FileType('w'),
                        default=sys.stdout, nargs='?')
    args = parser.parse_args()

    input_files = glob.glob(args.input_directory + "/*.json")
    input_dict = {os.path.splitext(os.path.basename(input_file))[0]:json.load(open(input_file, "r")) for input_file in input_files}

    print(input_dict)

    template_string = args.template.read()

    print(Template(template_string).render(input_dict))

if __name__ == "__main__":
    main()