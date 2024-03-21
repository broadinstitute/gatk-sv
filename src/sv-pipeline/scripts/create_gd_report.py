#!/bin/python

import argparse
import gzip
import logging
import os
import sys

from collections import defaultdict
from typing import List, Text, Optional

import pysam


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

    def __init__(self, kind, interval, path, vid, region):
        self.kind = kind
        self.interval = interval
        self.path = path
        self.vid = vid
        self.region = region

    def title(self):
        if self.kind == "orig_invalid":
            title = f"Flagged original variant {self.vid}"
        elif self.kind == "subtract_invalid":
            title = f"Flagged variant after correction {self.vid}"
        elif self.kind == "new":
            title = f"New variant {self.vid} resulting from another partially invalidated call"
        elif self.kind == "gdr":
            title = f"Raw region {self.region}"
        elif self.kind == "gdr2var":
            title = f"Region {self.region} overlapping variant {self.vid}"
        elif self.kind == "var2gdr":
            title = f"Variant {self.vid} overlapping region {self.region}"
        title = "## " + title + "\n\n"
        return title

    def image_string(self):
        return f"![Image]({self.path})\n\n"


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
        vid = "_".join(tokens[4:10])
        region = None
        for r in region_names:
            if r in filename:
                region = r
                break
        if "rdtest_orig_invalid" in filename:
            kind = "orig_invalid"
        elif "rdtest_subtract_invalid" in filename:
            kind = "subtract_invalid"
        elif "rdtest_new" in filename:
            kind = "new"
        elif "rdtest_full" in filename:
            kind = "gdr"
            vid = None
        elif "rdtest_gdr2var" in filename:
            kind = "gdr2var"
        elif "rdtest_var2gdr" in filename:
            kind = "var2gdr"
        else:
            raise ValueError(f"Unknown type for {path}")
        interval_str = f"{chrom}:{pos}-{stop}"
        image_list.append(ImageData(kind=kind, interval=interval_str, path=path, vid=vid, region=region))

    region_to_image_dict = defaultdict(list)
    vid_to_image_dict = defaultdict(list)
    for image in image_list:
        if image.region is not None:
            region_to_image_dict[image.region].append(image)
        if image.vid is not None:
            vid_to_image_dict[image.vid].append(image)

    with open("temp.md", "w") as f:
        for region in region_names:
            if len(region_to_image_dict[region]) > 1:
                f.write("# " + region + "\n\n")
                for image in region_to_image_dict[region]:
                    f.write(image.title())
                    f.write(image.image_string())
                f.write("\n---\n")

        for vid in sorted(list(vid_to_image_dict.keys())):
            f.write("# " + vid + "\n")
            for image in vid_to_image_dict[vid]:
                f.write(image.title())
                f.write(image.image_string())
            f.write("\n---\n")


if __name__ == "__main__":
    main()
