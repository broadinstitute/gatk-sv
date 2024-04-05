#!/bin/python

import argparse
import gzip
import textwrap
import logging
import os
import sys

from enum import Enum
from itertools import groupby
from operator import itemgetter
from typing import List, Text, Optional

IMAGE_WIDTH = 300
TABLE_WIDTH = 600
TABLE_PADDING = 20
TEXT_WRAP_WIDTH = 100
CODE_TABLE_TEXT_WRAP_WIDTH = 15
HIGHLIGHT_BG_COLOR_VAR2GDR = "#FFAAFF"
HIGHLIGHT_BG_COLOR_GDR2VAR = "#AAAAFF"
HIGHLIGHT_BG_COLOR = "#FFAAAA"

PLOT_TYPE_ENUM = Enum("PlotType", ["GDR", "GDR2VAR", "VAR2GDR", "BEFORE_REVISE", "AFTER_REVISE", "NEW", "SUBDIVISION"])
KIND_ORDERING = [PLOT_TYPE_ENUM.GDR, PLOT_TYPE_ENUM.GDR2VAR, PLOT_TYPE_ENUM.VAR2GDR, PLOT_TYPE_ENUM.BEFORE_REVISE,
                 PLOT_TYPE_ENUM.AFTER_REVISE, PLOT_TYPE_ENUM.NEW]

RDTEST_GDR = "rdtest_full"
RDTEST_GDR2VAR = "rdtest_gdr2var"
RDTEST_VAR2GDR = "rdtest_var2gdr"
RDTEST_BEFORE_REVISE = "rdtest_before_revise"
RDTEST_AFTER_REVISE = "rdtest_after_revise"
RDTEST_SUBDIVISION = "rdtest_subdiv"
RDTEST_NEW = "rdtest_new"

MANIFEST_EMPTY_FIELD = "."

RDTEST_NAME_TO_TYPE = {
    RDTEST_GDR: PLOT_TYPE_ENUM.GDR,
    RDTEST_GDR2VAR: PLOT_TYPE_ENUM.GDR2VAR,
    RDTEST_VAR2GDR: PLOT_TYPE_ENUM.VAR2GDR,
    RDTEST_BEFORE_REVISE: PLOT_TYPE_ENUM.BEFORE_REVISE,
    RDTEST_AFTER_REVISE: PLOT_TYPE_ENUM.AFTER_REVISE,
    RDTEST_SUBDIVISION: PLOT_TYPE_ENUM.SUBDIVISION,
    RDTEST_NEW: PLOT_TYPE_ENUM.NEW
}

# Delimiter suffix appended to the end of interval IDs before the index, e.g. "intervalA__0", "intervalA__1", ...
INDEX_DELIMITER = "__"

CODE_FALSE_NEGATIVE_IN_EXISTING_VARIANT = "FALSE_NEGATIVE_IN_EXISTING_VARIANT"
CODE_FALSE_POSITIVE_IN_EXISTING_VARIANT = "FALSE_POSITIVE_IN_EXISTING_VARIANT"
CODE_REVISED_BREAKPOINTS_OF_EXISTING_VARIANT = "REVISED_BREAKPOINTS_OF_EXISTING_VARIANT"
CODE_NEW_VARIANT_FOR_FALSE_NEGATIVE_IN_REGION = "NEW_VARIANT_FOR_FALSE_NEGATIVE_IN_REGION"


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
            title = f"Existing variant after revision"
        elif self.plot_type == PLOT_TYPE_ENUM.NEW:
            title = f"Revised breakpoints or new variant"
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


def show_images_table(f, images):
    unique_images = list()
    unique_image_names = set()
    for image in images:
        if image.name not in unique_image_names:
            unique_image_names.add(image.name)
            unique_images.append(image)
    images = unique_images
    f.write(f"<table border=1 cellpadding=10 width={TABLE_WIDTH}>\n")
    if len(images) == 0:
        f.write("<tr><td><h3>Warning: no images found!</h3></td></tr>")
    else:
        for image in images:
            f.write("<tr>\n")
            if image.plot_type == PLOT_TYPE_ENUM.VAR2GDR:
                f.write(f"<td align=center bgcolor=\"{HIGHLIGHT_BG_COLOR_VAR2GDR}\">\n")
            elif image.plot_type == PLOT_TYPE_ENUM.GDR2VAR:
                f.write(f"<td align=center bgcolor=\"{HIGHLIGHT_BG_COLOR_GDR2VAR}\">\n")
            elif image.plot_type != PLOT_TYPE_ENUM.GDR:
                f.write(f"<td align=center bgcolor=\"{HIGHLIGHT_BG_COLOR}\">\n")
            else:
                f.write("<td align=center>\n")
            f.write(image.image_string())
            name = "<br>".join(textwrap.wrap(image.name, width=TEXT_WRAP_WIDTH))
            f.write(f"{name}<br>\n")
            f.write(f"<b>{image.title()}</b><br>\n")
            f.write("</td>\n")
            f.write("</tr>\n")
    f.write(f"</table>\n")


def wrap_text_table(text):
    if text is None:
        return None
    else:
        return "<br>".join(textwrap.wrap(text, width=CODE_TABLE_TEXT_WRAP_WIDTH))


def show_manifest_table(f, manifest_records, genotypes):
    f.write(f"<table border=1 cellpadding=5 width={TABLE_WIDTH}>\n")
    f.write("<tr bgcolor=#EEEEEE>\n")
    f.write(f"<td>new_vid</td>\n")
    f.write(f"<td>coords</td>\n")
    f.write(f"<td>old_vid</td>\n")
    f.write(f"<td>old_coords</td>\n")
    f.write(f"<td>code</td>\n")
    f.write(f"<td>geno</td>\n")
    f.write("</tr>\n")
    for m in manifest_records:
        new_vid = wrap_text_table(m['new_vid'])
        old_vid = wrap_text_table(m['old_vid'])
        coords = wrap_text_table(f"{m['chrom']}:{m['pos']}-{m['stop']}")
        if m['old_pos'] is not None and m['old_stop'] is not None:
            old_coords = wrap_text_table(f"{m['chrom']}:{m['old_pos']}-{m['old_stop']}")
        else:
            old_coords = wrap_text_table(None)
        code = wrap_text_table(m['code'])
        table_genotypes = [g["genotype"] for g in genotypes if g['vid'] == m['new_vid']]
        genotypes_col = wrap_text_table(', '.join(table_genotypes))
        f.write(f"<td>{new_vid}</td>\n")
        f.write(f"<td>{coords}</td>\n")
        f.write(f"<td>{old_vid}</td>\n")
        f.write(f"<td>{old_coords}</td>\n")
        f.write(f"<td>{code}</td>\n")
        f.write(f"<td>{genotypes_col}</td>\n")
        f.write("</tr>\n")
    f.write(f"</table>\n")


def get_candidate_tokens(rdtest_name, m):
    subdir = f"/{rdtest_name}_{m['batch']}/"
    if rdtest_name == RDTEST_GDR2VAR:
        middle = f"_{m['region']}__"
        suffix = f"_{rdtest_name}_{m['batch']}.jpg"
    elif rdtest_name == RDTEST_VAR2GDR:
        middle = None
        suffix = f"__{m['region']}_{rdtest_name}_{m['batch']}.jpg"
    elif rdtest_name == RDTEST_GDR:
        middle = None
        suffix = f"_{m['region']}_{rdtest_name}_{m['batch']}.jpg"
    elif rdtest_name == RDTEST_SUBDIVISION:
        middle = f"_{m['region']}{INDEX_DELIMITER}"
        suffix = f"_{rdtest_name}_{m['batch']}.jpg"
    elif rdtest_name == RDTEST_AFTER_REVISE:
        middle = None
        suffix = f"_{m['new_vid']}_{rdtest_name}_{m['batch']}.jpg"
    elif rdtest_name == RDTEST_BEFORE_REVISE:
        middle = None
        suffix = f"_{m['new_vid']}_{rdtest_name}_{m['batch']}.jpg"
    elif rdtest_name == RDTEST_NEW:
        middle = None
        suffix = f"_{m['new_vid']}_{rdtest_name}_{m['batch']}.jpg"
    else:
        raise ValueError(f"Unknown rdtest name: {rdtest_name}")
    return subdir, middle, suffix


def get_candidate_paths(rdtest_names, m, image_paths):
    found_paths = set()
    for name in rdtest_names:
        subdir, middle, suffix = get_candidate_tokens(rdtest_name=name, m=m)
        for p in image_paths:
            if p.endswith(suffix) and subdir in p and (middle is None or middle in p) and p not in found_paths:
                found_paths.add(p)
                yield name, p


def get_image_paths(m, image_paths):
    # m : manifest record
    if m["code"] == CODE_FALSE_NEGATIVE_IN_EXISTING_VARIANT \
            or m["code"] == CODE_FALSE_POSITIVE_IN_EXISTING_VARIANT:
        rdtest_names = [RDTEST_BEFORE_REVISE, RDTEST_AFTER_REVISE]
    elif m["code"] == CODE_REVISED_BREAKPOINTS_OF_EXISTING_VARIANT:
        rdtest_names = [RDTEST_BEFORE_REVISE, RDTEST_AFTER_REVISE, RDTEST_NEW]
    elif m["code"] == CODE_NEW_VARIANT_FOR_FALSE_NEGATIVE_IN_REGION:
        rdtest_names = [RDTEST_AFTER_REVISE, RDTEST_NEW]
    elif m["code"] is None:
        # special case for region-wide images
        rdtest_names = [RDTEST_GDR, RDTEST_VAR2GDR, RDTEST_GDR2VAR]
    candidate_paths = list(set(get_candidate_paths(rdtest_names=rdtest_names, m=m, image_paths=image_paths)))
    return candidate_paths


def get_images(manifest, image_paths, region, batch):
    images_out = list()
    for m in manifest:
        m_image_paths = get_image_paths(m=m, image_paths=image_paths)
        images = list()
        for rdtest_name, path in m_image_paths:
            name = os.path.basename(path).replace(".jpg", "")
            images.append(
                ImageData(plot_type=RDTEST_NAME_TO_TYPE[rdtest_name], interval=None, path=path, name=name,
                          region=region, batch=batch))
        images = sorted(images, key=lambda x: (KIND_ORDERING.index(x.plot_type), x.name))
        images_out.extend(images)
    return images_out


def write_reports(base_path, image_paths, region_names, batches, genotypes, manifest):
    f_dict = {batch: open(f"{base_path}.{batch}.html", "w") for batch in batches}
    for batch in f_dict:
        f_dict[batch].write("<style>@page { size: letter portrait; margin: 2cm; } </style>\n")
    for region in region_names:
        region_manifest = [m for m in manifest if m["region"] == region]
        for f in f_dict.values():
            f.write("<h2>" + region + "</h2>\n\n")
        for batch in batches:
            f = f_dict[batch]
            m_region = {"new_vid": None, "old_vid": None, "svtype": None, "region": region,
                        "sample": None, "batch": batch, "code": None}
            gdr_image_paths = [p for _, p in get_candidate_paths(rdtest_names=[RDTEST_GDR, RDTEST_GDR2VAR, RDTEST_VAR2GDR],
                                                                 m=m_region, image_paths=image_paths)]
            gdr_images = get_images(manifest=[m_region], image_paths=gdr_image_paths, region=region, batch=batch)
            show_images_table(f, gdr_images)
            # Need to sort for groupby
            batch_manifest = sorted([m for m in region_manifest if m["batch"] == batch], key=itemgetter("sample"))
            for sample, sample_manifest in groupby(batch_manifest, key=itemgetter("sample")):
                sample_genotypes = [g for g in genotypes if g["sample"] == sample]
                sample_manifest = list(sample_manifest)
                sample_images = list()
                for m in sample_manifest:
                    images = get_images(manifest=sample_manifest, image_paths=image_paths, region=region, batch=batch)
                    images = sorted(images, key=lambda x: (KIND_ORDERING.index(x.plot_type), x.name))
                    sample_images.extend(images)
                f.write("<h3>" + m['sample'] + "</h3>\n\n")
                show_manifest_table(f=f, manifest_records=sample_manifest, genotypes=sample_genotypes)
                show_images_table(f, sample_images)
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
    parser.add_argument('--genotypes', type=str, required=True, help='Genotype revision file (.tsv.gz)')
    parser.add_argument('--manifest', type=str, required=True, help='Revision manifest file (.tsv.gz)')
    parser.add_argument('--batch-list', type=str, required=True, help='List of batch names')
    parser.add_argument('--region-bed', type=str, required=True, help='Region bed')
    parser.add_argument('--out', type=str, default="gdr", help='Output prefix')
    if len(argv) <= 1:
        parser.parse_args(["--help"])
        sys.exit(0)
    parsed_arguments = parser.parse_args(argv[1:])
    return parsed_arguments


def read_manifest_entry(val, integer=False):
    if val == MANIFEST_EMPTY_FIELD:
        return None
    elif integer:
        return int(val)
    else:
        return val


def main(argv: Optional[List[Text]] = None):
    if argv is None:
        argv = sys.argv
    args = _parse_arguments(argv)
    logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

    with open(args.image_list) as f:
        image_paths = [line.strip() for line in f]

    with open(args.batch_list) as f:
        batches = [line.strip() for line in f]

    with open(args.region_bed) as f:
        region_names = [line.strip().split("\t")[3] for line in f]

    with gzip.open(args.manifest, "rt") as f:
        manifest = list()
        for line in f:
            tokens = line.strip().split("\t")
            chrom = tokens[0]
            pos = int(tokens[1])
            stop = int(tokens[2])
            new_vid = tokens[3]
            old_vid = read_manifest_entry(tokens[4])
            old_pos = read_manifest_entry(tokens[5], integer=True)
            old_stop = read_manifest_entry(tokens[6], integer=True)
            svtype = tokens[7]
            region = tokens[8]
            sample = tokens[9]
            batch = tokens[10]
            code = tokens[11]
            manifest.append({"chrom": chrom, "pos": pos, "stop": stop, "new_vid": new_vid,
                             "old_vid": old_vid, "old_pos": old_pos, "old_stop": old_stop,
                             "svtype": svtype, "region": region, "sample": sample, "batch": batch, "code": code})

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

    # image_list = load_images(paths=image_paths, region_names=region_names)
    # region_to_image_dict, no_region_images, batches = create_image_dict(image_list=image_list)
    write_reports(base_path=args.out, image_paths=image_paths,
                  region_names=region_names, batches=batches, genotypes=genotypes, manifest=manifest)


if __name__ == "__main__":
    main()
