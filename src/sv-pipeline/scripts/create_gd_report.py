#!/bin/python

import argparse
import textwrap
import logging
import os
import sys

from collections import defaultdict
from typing import List, Text, Optional

IMAGE_WIDTH = 500
TABLE_PADDING = 20
TEXT_WRAP_WIDTH = 50

KIND_ORDERING = ["gdr", "gdr2var", "var2gdr", "before_revise", "after_revise"]


def _parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Resets DEL/DUP genotypes to homozygous-reference",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--image-list', type=str, required=True, help='Input images')
    parser.add_argument('--region-bed', type=str, required=True, help='GDR bed')
    parser.add_argument('--out', type=str, default="gdr", help='Output prefix')
    if len(argv) <= 1:
        parser.parse_args(["--help"])
        sys.exit(0)
    parsed_arguments = parser.parse_args(argv[1:])
    return parsed_arguments


class ImageData:

    def __init__(self, kind, interval, path, name, region, batch):
        self.kind = kind
        self.interval = interval
        self.path = path
        self.name = name
        self.region = region
        self.batch = batch

    def title(self):
        if self.kind == "before_revise":
            title = f"Existing variant before revision"
        elif self.kind == "after_revise":
            if "new_rescue" in self.name:
                title = f"New variant"
            else:
                title = f"Existing variant after revision"
        elif self.kind == "new":
            title = f"New variant resulting from another partially invalidated call"
        elif self.kind == "gdr":
            title = f"Raw region (random sample)"
        elif self.kind == "gdr2var":
            title = f"Region overlapping existing variant"
        elif self.kind == "var2gdr":
            title = f"Existing variant overlapping region"
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
            batch_sep = "rdtest_before_revise"
        elif "rdtest_after_revise" in filename:
            kind = "after_revise"
            batch_sep = "rdtest_after_revise"
        elif "rdtest_full" in filename:
            kind = "gdr"
            batch_sep = "rdtest_full"
        elif "rdtest_gdr2var" in filename:
            kind = "gdr2var"
            batch_sep = "rdtest_gdr2var"
        elif "rdtest_var2gdr" in filename:
            kind = "var2gdr"
            batch_sep = "rdtest_var2gdr"
        elif "rdtest_subdiv" in filename:
            kind = "subdiv"
            batch_sep = "rdtest_subdiv"
        else:
            raise ValueError(f"Unknown type for {path}")
        # TODO create separate pdf for subdivisions
        if kind == "subdiv":
            continue
        batch = filename.split(batch_sep + "_")[-1].replace(".jpg", "")
        interval_str = f"{chrom}:{pos}-{stop}"
        image_list.append(ImageData(kind=kind, interval=interval_str, path=path, name=name, region=region, batch=batch))

    region_to_image_dict = dict()
    batches = set()
    for image in image_list:
        if image.region is not None:
            if image.region not in region_to_image_dict:
                region_to_image_dict[image.region] = dict()
            if image.batch not in region_to_image_dict[image.region]:
                region_to_image_dict[image.region][image.batch] = list()
            batches.add(image.batch)
            region_to_image_dict[image.region][image.batch].append(image)

    batches = sorted(list(batches))
    f_dict = {batch: open(f"{args.out}.{batch}.html", "w") for batch in batches}
    for region in region_names:
        for f in f_dict.values():
            f.write("<h2>" + region + "</h2>\n\n")
        for batch in batches:
            if region not in region_to_image_dict:
                continue
            f = f_dict[batch]
            f.write("<h3>" + batch + "</h3>\n\n")
            images = region_to_image_dict.get(region, dict).get(batch, list)
            num_images = len(images)
            table_width = num_images * (IMAGE_WIDTH + TABLE_PADDING)
            f.write(f"<table width={table_width} border=1>\n")
            f.write("<tr>\n")
            images = sorted(images, key=lambda x: (KIND_ORDERING.index(x.kind), x.name))
            for image in images:
                name = "<br>".join(textwrap.wrap(image.name, width=TEXT_WRAP_WIDTH))
                if image.kind != "gdr":
                    f.write("<td align=center bgcolor=\"#FFAAAA\">\n")
                else:
                    f.write("<td align=center>\n")
                f.write(image.image_string())
                f.write(f"{name}<br>\n")
                f.write(f"<b>{image.title()}</b><br>\n")
                f.write("</td>\n")
            f.write("</tr>\n")
            f.write(f"</table>\n")
        for f in f_dict.values():
            f.write("\n<hr>\n")
    for f in f_dict.values():
        f.close()


if __name__ == "__main__":
    main()
