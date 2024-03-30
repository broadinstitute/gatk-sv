#!/bin/python

import argparse
import gzip
import textwrap
import logging
import os
import sys

from enum import Enum
from typing import List, Text, Optional

IMAGE_WIDTH = 256
TABLE_WIDTH = 700
TABLE_PADDING = 20
TEXT_WRAP_WIDTH = 100
HIGHLIGHT_BG_COLOR = "#FFAAAA"

PLOT_TYPE_ENUM = Enum("PlotType", ["GDR", "GDR2VAR", "VAR2GDR", "BEFORE_REVISE", "AFTER_REVISE", "SUBDIVISION"])
KIND_ORDERING = [PLOT_TYPE_ENUM.GDR, PLOT_TYPE_ENUM.GDR2VAR, PLOT_TYPE_ENUM.VAR2GDR, PLOT_TYPE_ENUM.BEFORE_REVISE,
                 PLOT_TYPE_ENUM.AFTER_REVISE]


class ImageData:
    def __init__(self, plot_type, interval, path, name, region, batch):
        self.plot_type = plot_type
        self.interval = interval
        self.path = path
        self.name = name
        self.region = region
        self.batch = batch

    def title(self):
        if self.plot_type == PLOT_TYPE_ENUM.BEFORE_REVISE:
            title = f"Existing variant before revision"
        elif self.plot_type == PLOT_TYPE_ENUM.AFTER_REVISE:
            if "new_rescue" in self.name:
                title = f"New variant"
            else:
                title = f"Existing variant after revision"
        elif self.plot_type == "new":
            title = f"New variant resulting from another partially invalidated call"
        elif self.plot_type == PLOT_TYPE_ENUM.GDR:
            title = f"Raw region (random sample)"
        elif self.plot_type == PLOT_TYPE_ENUM.GDR2VAR:
            title = f"Region overlapping existing variant"
        elif self.plot_type == PLOT_TYPE_ENUM.VAR2GDR:
            title = f"Existing variant overlapping region"
        elif self.plot_type == PLOT_TYPE_ENUM.SUBDIVISION:
            title = f"Region subdivision (random sample)"
        return title

    def image_string(self):
        return f"<img src=\"{self.path}\" width={IMAGE_WIDTH}><br>\n\n"


def load_images(paths, region_names):
    image_list = list()
    for path in paths:
        filename = os.path.basename(path)
        tokens = filename.split('_')
        chrom = tokens[0]
        pos = int(tokens[1])
        stop = int(tokens[2])
        name = filename.replace(".jpg", "")
        region = None
        for r in region_names:
            if r in filename:
                region = r
                break
        if "rdtest_before_revise" in filename:
            plot_type = PLOT_TYPE_ENUM.BEFORE_REVISE
            batch_sep = "rdtest_before_revise"
        elif "rdtest_after_revise" in filename:
            plot_type = PLOT_TYPE_ENUM.AFTER_REVISE
            batch_sep = "rdtest_after_revise"
        elif "rdtest_full" in filename:
            plot_type = PLOT_TYPE_ENUM.GDR
            batch_sep = "rdtest_full"
        elif "rdtest_gdr2var" in filename:
            plot_type = PLOT_TYPE_ENUM.GDR2VAR
            batch_sep = "rdtest_gdr2var"
        elif "rdtest_var2gdr" in filename:
            plot_type = PLOT_TYPE_ENUM.VAR2GDR
            batch_sep = "rdtest_var2gdr"
        elif "rdtest_subdiv" in filename:
            plot_type = PLOT_TYPE_ENUM.SUBDIVISION
            batch_sep = "rdtest_subdiv"
        else:
            raise ValueError(f"Unknown type for {path}")
        # TODO create separate report for subdivisions
        if plot_type == PLOT_TYPE_ENUM.SUBDIVISION:
            continue
        batch = filename.split(batch_sep + "_")[-1].replace(".jpg", "")
        interval_str = f"{chrom}:{pos}-{stop}"
        image_list.append(ImageData(plot_type=plot_type, interval=interval_str, path=path, name=name, region=region,
                                    batch=batch))
    return image_list


def create_image_dict(image_list):
    region_to_image_dict = dict()
    batches = set()
    no_region_images = list()
    for image in image_list:
        if image.region is None:
            no_region_images.append(image)
        else:
            if image.region not in region_to_image_dict:
                region_to_image_dict[image.region] = dict()
            if image.batch not in region_to_image_dict[image.region]:
                region_to_image_dict[image.region][image.batch] = list()
            batches.add(image.batch)
            region_to_image_dict[image.region][image.batch].append(image)
    batches = sorted(list(batches))
    return region_to_image_dict, no_region_images, batches


def show_image(f, image):
    name = "<br>".join(textwrap.wrap(image.name, width=TEXT_WRAP_WIDTH))
    f.write(f"<table width={TABLE_WIDTH} border=1 cellpadding=10>\n")
    f.write("<tr>\n")
    if image.plot_type != PLOT_TYPE_ENUM.GDR:
        f.write(f"<td align=center bgcolor=\"{HIGHLIGHT_BG_COLOR}\">\n")
    else:
        f.write("<td align=center>\n")
    f.write(image.image_string())
    f.write(f"{name}<br>\n")
    f.write(f"<b>{image.title()}</b><br>\n")
    f.write("</td>\n")
    f.write("</tr>\n")
    f.write(f"</table>\n")


def write_reports(base_path, region_to_image_dict, no_region_images, region_names, batches, genotypes):
    f_dict = {batch: open(f"{base_path}.{batch}.html", "w") for batch in batches}
    for batch in f_dict:
        f_dict[batch].write("<style>@page { size: letter portrait; margin: 2cm; } </style>\n")
    for region in region_names:
        for f in f_dict.values():
            f.write("<h2>" + region + "</h2>\n\n")
        for batch in batches:
            if region not in region_to_image_dict:
                continue
            f = f_dict[batch]
            f.write("<h3>" + batch + "</h3>\n\n")
            images = region_to_image_dict.get(region, dict).get(batch, list)
            images = sorted(images, key=lambda x: (KIND_ORDERING.index(x.plot_type), x.name))
            image_names = [image.name for image in images]
            for image in images:
                show_image(f, image)
            image_genotypes = [g for g in genotypes if any(g["vid"] + "_" in img for img in image_names)]
            if len(image_genotypes) > 0:
                f.write("<h3 color=#FF0000>Possibly relevant genotype revisions</h3>\n\n")
                f.write("<p>Samples may exist in other batches. "
                        "For large variants, these revisions may have been triggered by another region. "
                        "Also, in rare cases these may be unrelated if the variant ID is a substring of "
                        "another.</p>\n\n")
                f.write("<table border=1 cellpadding=5>\n")
                for g in image_genotypes:
                    f.write(f"\t<tr><td width=300>{g['vid']}</td><td width=100>{g['sample']}</td><td width=50>{g['genotype']}</td></tr>\n")
                f.write("</table>\n")
        for f in f_dict.values():
            f.write("\n<hr>\n")
    for f in f_dict.values():
        f.close()


def _parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Resets DEL/DUP genotypes to homozygous-reference",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--image-list', type=str, required=True, help='Input images')
    parser.add_argument('--genotypes', type=str, required=True, help='Genotype revisions file (.tsv.gz)')
    parser.add_argument('--region-bed', type=str, required=True, help='GDR bed')
    parser.add_argument('--out', type=str, default="gdr", help='Output prefix')
    if len(argv) <= 1:
        parser.parse_args(["--help"])
        sys.exit(0)
    parsed_arguments = parser.parse_args(argv[1:])
    return parsed_arguments


def main(argv: Optional[List[Text]] = None):
    if argv is None:
        argv = sys.argv
    args = _parse_arguments(argv)
    logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

    with open(args.image_list) as f:
        image_paths = [line.strip() for line in f]

    with open(args.region_bed) as f:
        region_names = [line.strip().split("\t")[3] for line in f]

    with gzip.open(args.genotypes, "rt") as f:
        genotypes = list()
        for line in f:
            tokens = line.strip().split("\t")
            vid = tokens[1]
            sample = tokens[2]
            genotype = tokens[3]
            if genotype == "0":
                genotype = "hom_ref"
            elif genotype == "1":
                genotype = "het"
            elif genotype == "2":
                genotype = "hom_var"
            else:
                raise ValueError(f"Unrecognized genotype code {genotype}, should be in " + '{0, 1, 2}')
            genotypes.append({"vid": vid, "sample": sample, "genotype": genotype})

    image_list = load_images(paths=image_paths, region_names=region_names)
    region_to_image_dict, no_region_images, batches = create_image_dict(image_list=image_list)
    write_reports(base_path=args.out, region_to_image_dict=region_to_image_dict, no_region_images=no_region_images,
                  region_names=region_names, batches=batches, genotypes=genotypes)


if __name__ == "__main__":
    main()
