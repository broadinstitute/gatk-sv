#!/bin/python

import argparse
import gzip
import logging
import os
import sys

from collections import defaultdict
from typing import List, Text, Optional

import pysam

IMAGE_WIDTH = 500

def _parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Resets DEL/DUP genotypes to homozygous-reference",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--image-list', type=str, required=True, help='Input images')
    parser.add_argument('--region-bed', type=str, required=True, help='GDR bed')
    parser.add_argument('--out', type=str, required=True, help='Output pdf')
    if len(argv) <= 1:
        parser.parse_args(["--help"])
        sys.exit(0)
    parsed_arguments = parser.parse_args(argv[1:])
    return parsed_arguments


class ImageData:

    def __init__(self, kind, interval, path, name, region):
        self.kind = kind
        self.interval = interval
        self.path = path
        self.name = name
        self.region = region

    def title(self):
        if self.kind == "before_revise":
            title = f"Variant before revision"
        elif self.kind == "after_revise":
            title = f"Variant after revision"
        elif self.kind == "new":
            title = f"New variant resulting from another partially invalidated call"
        elif self.kind == "gdr":
            title = f"Raw region (random sample)"
        elif self.kind == "gdr2var":
            title = f"Region overlapping variant"
        elif self.kind == "var2gdr":
            title = f"Variant overlapping region"
        elif self.kind == "subdiv":
            title = f"Region subdivision (random sample)"
        return title

    def image_string(self):
        return f"<img src=\"{self.path}\" width={IMAGE_WIDTH}><br>\n\n"


def main(argv: Optional[List[Text]] = None):
    if argv is None:
        argv = sys.argv
    args = _parse_arguments(argv)
    logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

    with open(args.image_list) as f:
        image_paths = [line.strip() for line in f]

    with open(args.region_bed) as f:
        region_names = [line.strip().split("\t")[3] for line in f]

    image_list = list()
    for path in image_paths:
        #print(path)
        filename = os.path.basename(path)
        tokens = filename.split('_')
        chrom = tokens[0]
        pos = int(tokens[1])
        stop = int(tokens[2])
        sample = tokens[3]
        name = filename.replace(".jpg", "")
        region = None
        for r in region_names:
            if r in filename:
                region = r
                break
        if "rdtest_before_revise" in filename:
            kind = "before_revise"
        elif "rdtest_after_revise" in filename:
            kind = "after_revise"
        elif "rdtest_full" in filename:
            kind = "gdr"
        elif "rdtest_gdr2var" in filename:
            kind = "gdr2var"
        elif "rdtest_var2gdr" in filename:
            kind = "var2gdr"
        elif "rdtest_subdiv" in filename:
            kind = "subdiv"
        else:
            raise ValueError(f"Unknown type for {path}")
        # TODO create separate pdf for subdivisions
        if kind == "subdiv":
            continue
        interval_str = f"{chrom}:{pos}-{stop}"
        image_list.append(ImageData(kind=kind, interval=interval_str, path=path, name=name, region=region))

    region_to_image_dict = defaultdict(list)
    for image in image_list:
        if image.region is not None:
            region_to_image_dict[image.region].append(image)

    with open("temp.html", "w") as f:
        for region in region_names:
            f.write("<h2>" + region + "</h2>\n\n")
            for image in region_to_image_dict[region]:
                f.write(f"<h3>{image.name}<br>\n")
                f.write(f"{image.title()}</h3>\n")
                f.write(image.image_string())
            f.write("\n<hr>\n")


if __name__ == "__main__":
    main()
