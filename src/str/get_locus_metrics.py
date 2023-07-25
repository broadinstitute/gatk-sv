"""
This scripts implements a hacky solution to measure a some metrics on the STR expansion loci.
These metrics could be ideally extracted from the STR expansion caller, ExpansionHunter;
but the tool is not currently maintained. Hence, this script leverages parsing SVG images
that render a locus to measure the metrics.
"""

import argparse
import json
import numpy as np
import os
import pandas as pd
import re
import subprocess
import tempfile
import time

from dataclasses import dataclass
from typing import List, Tuple, Union


SAMPLE_ID_COL_NAME = "sample_id"
LOCUS_ID_COL_NAME = "locus_id"


@dataclass
class ExtendedMetrics:
    total_upper_interrupt = np.nan
    results_total_upper_count = np.nan
    content_c_y_axis_count = np.nan
    total_bottom_interrupt = np.nan
    results_total_bottom_count = np.nan
    gap_bottom_count = np.nan
    gap_upper_count = np.nan
    results_total_upper = []
    results_total_bottom = []
    start_upper_x_orange = np.nan
    start_upper_y = np.nan
    start_bottom_x_orange = np.nan
    end_bottom_x_orange = np.nan
    start_bottom_y = np.nan
    end_upper_x_orange = np.nan


def get_loci(variant_catalog: str) -> List[str]:
    with open(variant_catalog, "r") as var_catalog_file:
        catalog = json.load(var_catalog_file)
    return [x["LocusId"] for x in catalog]


def call_reviewer(
        sample_id: str,
        realigned_reads: str,
        reference: str,
        genotypes: str,
        catalog: str,
        locus_id: str) -> Tuple[int, Union[str, None], Union[pd.DataFrame, None], Union[pd.DataFrame]]:

    prefix = f"{sample_id}_{locus_id.replace('-', '_')}"

    start_time = time.time()
    cmd = f"reviewer " \
          f"--reads {realigned_reads} " \
          f"--vcf {genotypes} " \
          f"--reference {reference} " \
          f"--catalog {catalog} " \
          f"--locus {locus_id} " \
          f"--output-prefix {prefix}"
    result = subprocess.run(cmd, shell=True)
    print(f"REViewer runtime: {time.time() - start_time :.2f} seconds")

    # Expected filename of the reviewer outputs.
    ro_svg_filename = f"{prefix}.{locus_id}.svg"
    ro_metrics_filename = f"{prefix}.metrics.tsv"
    ro_phasing_filename = f"{prefix}.phasing.tsv"

    if any(not os.path.isfile(x) for x in [ro_svg_filename, ro_metrics_filename, ro_phasing_filename]):
        return result.returncode, None, None, None

    svg_filename = f"{prefix}.svg"
    metrics_filename = f"{prefix}_metrics.csv"
    phasing_filename = f"{prefix}_phasing.csv"
    os.rename(ro_svg_filename, svg_filename)
    os.rename(ro_metrics_filename, metrics_filename)
    os.rename(ro_phasing_filename, phasing_filename)

    metrics_df = pd.read_csv(metrics_filename, delimiter="\t")
    phasing_df = pd.read_csv(phasing_filename, delimiter="\t")

    # The two dfs may have different shapes, so row count is get for each df separately.
    metrics_df[SAMPLE_ID_COL_NAME] = [sample_id] * metrics_df.shape[0]
    phasing_df[SAMPLE_ID_COL_NAME] = [sample_id] * phasing_df.shape[0]
    metrics_df[LOCUS_ID_COL_NAME] = [locus_id] * metrics_df.shape[0]
    phasing_df[LOCUS_ID_COL_NAME] = [locus_id] * phasing_df.shape[0]

    return result.returncode, svg_filename, metrics_df, phasing_df


def get_metrics_from_image(svg_filename: str) -> ExtendedMetrics:
    """
    This method attempts to extract some additional metrics of STR expansion loci
    than those exported by ExpansionHunter (EH). It would be ideal if EH outputs
    these metrics, but since EH is not developed/extended at the moment, this
    method attempts to extract those metrics by parsing the SVG images the
    Illumina/REViewer tool generates.

    This method is based on an ad-hoc python Script that was
    developed for exploratory purposes.
    """

    with open(svg_filename, 'rt') as f:
        _ = f.readline()
        img_contents = f.read()

    img_contents_no_def = img_contents.replace("</svg>", "").split("</defs>")[-1]
    matches = list(re.finditer(r"<line[^>]+y1=\"(\d+)\"[^>]+#arrow[^>]+>", img_contents_no_def, re.DOTALL))

    ex = ExtendedMetrics()

    if len(matches) == 0 or len(matches) > 2:
        return ex

    # For STRs with a single haplotype:
    if len(matches) == 1:
        # T, A, G, C upper: how many times sequence is interrupted starting with nucleotide T, A, G, or C.
        # This section depends on how the lines are described in image it will count >T<, >A<, >G< or, >C< interruptions.
        # The output is not the total number of nucleotides.

        section_upper_contents = img_contents_no_def
        ex.t_upper = list(re.finditer(r"fill=\"#cdcdcd\"\sopacity=\"(1|0.7)\"\s/>\s<text\sx=\"(\d+)\"\sy=\"(\d+)\"\sdy=\"0.25em\"\stext-anchor=\"middle\"\sfont-family=\"monospace\"\sfont-size=\"11px\"\sstroke=\"#FCA100\">T<", section_upper_contents, re.DOTALL))
        ex.a_upper = list(re.finditer(r"fill=\"#cdcdcd\"\sopacity=\"(1|0.7)\"\s/>\s<text\sx=\"(\d+)\"\sy=\"(\d+)\"\sdy=\"0.25em\"\stext-anchor=\"middle\"\sfont-family=\"monospace\"\sfont-size=\"11px\"\sstroke=\"#FF6347\">A<", section_upper_contents, re.DOTALL))
        ex.g_upper = list(re.finditer(r"fill=\"#cdcdcd\"\sopacity=\"(1|0.7)\"\s/>\s<text\sx=\"(\d+)\"\sy=\"(\d+)\"\sdy=\"0.25em\"\stext-anchor=\"middle\"\sfont-family=\"monospace\"\sfont-size=\"11px\"\sstroke=\"#2F8734\">G<", section_upper_contents, re.DOTALL))
        ex.c_upper = list(re.finditer(r"fill=\"#cdcdcd\"\sopacity=\"(1|0.7)\"\s/>\s<text\sx=\"(\d+)\"\sy=\"(\d+)\"\sdy=\"0.25em\"\stext-anchor=\"middle\"\sfont-family=\"monospace\"\sfont-size=\"11px\"\sstroke=\"#393939\">C<", section_upper_contents, re.DOTALL))
        ex.t_upper_count = len(ex.t_upper)
        ex.a_upper_count = len(ex.a_upper)
        ex.g_upper_count = len(ex.g_upper)
        ex.c_upper_count = len(ex.c_upper)
        ex.total_upper_interrupt = ex.a_upper_count + ex.t_upper_count + ex.g_upper_count + ex.c_upper_count
        ex.results_total_upper = []

        upper_match = str(matches[0])
        start_upper_x_orange_match = list(re.finditer(r"x1=\"(\d+)\"", upper_match, re.DOTALL))
        end_upper_x_orange_match = list(re.finditer(r"x2=\"(\d+)\"", upper_match, re.DOTALL))
        start_upper_y_match = list(re.finditer(r"y1=\"(\d+)\"", upper_match, re.DOTALL))
        ex.start_upper_x_orange = start_upper_x_orange_match[0].group(1)
        ex.end_upper_x_orange = end_upper_x_orange_match[0].group(1)
        ex.start_upper_y = start_upper_y_match[0].group(1)

        for i in range(int(ex.t_upper_count)):
            content_t_y_axis = int(ex.t_upper[i].group(3))
            if content_t_y_axis not in ex.results_total_upper:
                ex.results_total_upper.append(content_t_y_axis)
        for i in range(int(ex.a_upper_count)):
            content_a_yaxis = int(ex.a_upper[i].group(3))
            if content_a_yaxis not in ex.results_total_upper:
                ex.results_total_upper.append(content_a_yaxis)
        for i in range(int(ex.g_upper_count)):
            content_g_yaxis = int(ex.g_upper[i].group(3))
            if content_g_yaxis not in ex.results_total_upper:
                ex.results_total_upper.append(content_g_yaxis)
        for i in range(int(ex.c_upper_count)):
            content_c_y_axis = int(ex.c_upper[i].group(3))
            if content_c_y_axis not in ex.results_total_upper:
                ex.results_total_upper.append(content_c_y_axis)
        gap_upper = list(re.finditer(r"</text><line\sx1=\"(\d+)\".y1=\"(\d+)\"\sx2=\"(\d+)\"\sy2=\"(\d+)\"\sstroke=\"black\"\s/>", section_upper_contents, re.DOTALL))
        ex.gap_upper_count = len(gap_upper)
        ex.gap_bottom_count = np.nan
        ex.results_total_bottom = []
        ex.total_bottom_interrupt = np.nan
        ex.start_bottom_x_orange = np.nan
        ex.end_bottom_x_orange = np.nan
        ex.start_bottom_y = np.nan
    if len(matches) == 2:
        section_upper_contents = img_contents_no_def[matches[0].start():matches[1].start()]
        ex.t_upper = list(re.finditer(r"fill=\"#cdcdcd\"\sopacity=\"(1|0.7)\"\s/>\s<text\sx=\"(\d+)\"\sy=\"(\d+)\"\sdy=\"0.25em\"\stext-anchor=\"middle\"\sfont-family=\"monospace\"\sfont-size=\"11px\"\sstroke=\"#FCA100\">T<", section_upper_contents, re.DOTALL))
        ex.a_upper = list(re.finditer(r"fill=\"#cdcdcd\"\sopacity=\"(1|0.7)\"\s/>\s<text\sx=\"(\d+)\"\sy=\"(\d+)\"\sdy=\"0.25em\"\stext-anchor=\"middle\"\sfont-family=\"monospace\"\sfont-size=\"11px\"\sstroke=\"#FF6347\">A<", section_upper_contents, re.DOTALL))
        ex.g_upper = list(re.finditer(r"fill=\"#cdcdcd\"\sopacity=\"(1|0.7)\"\s/>\s<text\sx=\"(\d+)\"\sy=\"(\d+)\"\sdy=\"0.25em\"\stext-anchor=\"middle\"\sfont-family=\"monospace\"\sfont-size=\"11px\"\sstroke=\"#2F8734\">G<", section_upper_contents, re.DOTALL))
        ex.c_upper = list(re.finditer(r"fill=\"#cdcdcd\"\sopacity=\"(1|0.7)\"\s/>\s<text\sx=\"(\d+)\"\sy=\"(\d+)\"\sdy=\"0.25em\"\stext-anchor=\"middle\"\sfont-family=\"monospace\"\sfont-size=\"11px\"\sstroke=\"#393939\">C<", section_upper_contents, re.DOTALL))
        ex.t_upper_count = len(ex.t_upper)
        ex.a_upper_count = len(ex.a_upper)
        ex.g_upper_count = len(ex.g_upper)
        ex.c_upper_count = len(ex.c_upper)
        ex.total_upper_interrupt = ex.a_upper_count + ex.t_upper_count + ex.g_upper_count + ex.c_upper_count
        ex.results_total_upper = []
        upper_match = str(matches[0])
        start_upper_x_orange_match = list(re.finditer(r"x1=\"(\d+)\"", upper_match, re.DOTALL))
        end_upper_x_orange_match = list(re.finditer(r"x2=\"(\d+)\"", upper_match, re.DOTALL))
        start_upper_y_match = list(re.finditer(r"y1=\"(\d+)\"", upper_match, re.DOTALL))
        ex.start_upper_x_orange = start_upper_x_orange_match[0].group(1)
        ex.end_upper_x_orange = end_upper_x_orange_match[0].group(1)
        ex.start_upper_y = start_upper_y_match[0].group(1)
        for i in range(int(ex.t_upper_count)):
            content_t_y_axis = int(ex.t_upper[i].group(3))
            if content_t_y_axis not in ex.results_total_upper:
                ex.results_total_upper.append(content_t_y_axis)
        for i in range(int(ex.a_upper_count)):
            content_a_yaxis = int(ex.a_upper[i].group(3))
            if content_a_yaxis not in ex.results_total_upper:
                ex.results_total_upper.append(content_a_yaxis)
        for i in range(int(ex.g_upper_count)):
            content_g_yaxis = int(ex.g_upper[i].group(3))
            if content_g_yaxis not in ex.results_total_upper:
                ex.results_total_upper.append(content_g_yaxis)
        for i in range(int(ex.c_upper_count)):
            content_c_y_axis = int(ex.c_upper[i].group(3))
            if content_c_y_axis not in ex.results_total_upper:
                ex.results_total_upper.append(content_c_y_axis)
        gap_upper = list(re.finditer(r"</text><line\sx1=\"(\d+)\".y1=\"(\d+)\"\sx2=\"(\d+)\"\sy2=\"(\d+)\"\sstroke=\"black\"\s/>", section_upper_contents, re.DOTALL))
        ex.gap_upper_count = len(gap_upper)
        section_bottom_contents = img_contents_no_def[matches[1].start():]
        ex.results_total_bottom = []
        t_bottom = list(re.finditer(r"fill=\"#cdcdcd\"\sopacity=\"(1|0.7)\"\s/>\s<text\sx=\"(\d+)\"\sy=\"(\d+)\"\sdy=\"0.25em\"\stext-anchor=\"middle\"\sfont-family=\"monospace\"\sfont-size=\"11px\"\sstroke=\"#FCA100\">T<", section_bottom_contents, re.DOTALL))
        a_bottom = list(re.finditer(r"fill=\"#cdcdcd\"\sopacity=\"(1|0.7)\"\s/>\s<text\sx=\"(\d+)\"\sy=\"(\d+)\"\sdy=\"0.25em\"\stext-anchor=\"middle\"\sfont-family=\"monospace\"\sfont-size=\"11px\"\sstroke=\"#FF6347\">A<", section_bottom_contents, re.DOTALL))
        g_bottom = list(re.finditer(r"fill=\"#cdcdcd\"\sopacity=\"(1|0.7)\"\s/>\s<text\sx=\"(\d+)\"\sy=\"(\d+)\"\sdy=\"0.25em\"\stext-anchor=\"middle\"\sfont-family=\"monospace\"\sfont-size=\"11px\"\sstroke=\"#2F8734\">G<", section_bottom_contents, re.DOTALL))
        c_bottom = list(re.finditer(r"fill=\"#cdcdcd\"\sopacity=\"(1|0.7)\"\s/>\s<text\sx=\"(\d+)\"\sy=\"(\d+)\"\sdy=\"0.25em\"\stext-anchor=\"middle\"\sfont-family=\"monospace\"\sfont-size=\"11px\"\sstroke=\"#393939\">C<", section_bottom_contents, re.DOTALL))
        t_bottom_count = len(t_bottom)
        a_bottom_count = len(a_bottom)
        g_bottom_count = len(g_bottom)
        c_bottom_count = len(c_bottom)
        ex.total_bottom_interrupt = a_bottom_count + t_bottom_count + g_bottom_count + c_bottom_count
        bottom_match = str(matches[1])
        start_bottom_x_orange_match = list(re.finditer(r"x1=\"(\d+)\"", bottom_match, re.DOTALL))
        end_bottom_x_orange_match = list(re.finditer(r"x2=\"(\d+)\"", bottom_match, re.DOTALL))
        start_bottom_y_match = list(re.finditer(r"y1=\"(\d+)\"", bottom_match, re.DOTALL))
        ex.start_bottom_x_orange = start_bottom_x_orange_match[0].group(1)
        ex.end_bottom_x_orange = end_bottom_x_orange_match[0].group(1)
        ex.start_bottom_y = start_bottom_y_match[0].group(1)
        for i in range(int(t_bottom_count)):
            content_t_y_axis = int(t_bottom[i].group(3))
            if content_t_y_axis not in ex.results_total_bottom:
                ex.results_total_bottom.append(content_t_y_axis)
        for i in range(int(a_bottom_count)):
            content_a_yaxis = int(a_bottom[i].group(3))
            if content_a_yaxis not in ex.results_total_bottom:
                ex.results_total_bottom.append(content_a_yaxis)
        for i in range(int(g_bottom_count)):
            content_g_yaxis = int(g_bottom[i].group(3))
            if content_g_yaxis not in ex.results_total_bottom:
                ex.results_total_bottom.append(content_g_yaxis)
        for i in range(int(c_bottom_count)):
            content_c_y_axis = int(c_bottom[i].group(3))
            if content_c_y_axis not in ex.results_total_bottom:
                ex.results_total_bottom.append(content_c_y_axis)
        gap_upper = list(re.finditer(r"</text><line\sx1=\"(\d+)\".y1=\"(\d+)\"\sx2=\"(\d+)\"\sy2=\"(\d+)\"\sstroke=\"black\"\s/>", section_upper_contents, re.DOTALL))
        ex.gap_upper_count = len(gap_upper)
        gap_bottom = list(re.finditer(r"</text><line\sx1=\"(\d+)\".y1=\"(\d+)\"\sx2=\"(\d+)\"\sy2=\"(\d+)\"\sstroke=\"black\"\s/>", section_bottom_contents, re.DOTALL))
        ex.gap_bottom_count = len(gap_bottom)
    return ex


def get_nucleotides_count(
        svg_filename: str,
        metrics: ExtendedMetrics,
        sample_id: str,
        locus_id: str) -> pd.DataFrame:

    tmp_info_header = "\t".join(["start_upper_X_orange", "end_upper_X_orange", "start_upper_Y", "start_bottom_X_orange", "end_bottom_X_orange", "start_bottom_Y"])
    tmp_info_line = "\t".join([str(metrics.start_upper_x_orange), str(metrics.end_upper_x_orange), str(metrics.start_upper_y), str(metrics.start_bottom_x_orange), str(metrics.end_bottom_x_orange), str(metrics.start_bottom_y)])
    tmp_info_filename = tempfile.NamedTemporaryFile(delete=False).name
    with open(tmp_info_filename, "w") as tmp_info_file:
        tmp_info_file.write(tmp_info_header + "\n")
        tmp_info_file.write(tmp_info_line + "\n")

    tmp_out_total_filename = tempfile.NamedTemporaryFile(delete=False).name
    cmd = f"bash /opt/str/count_nucleotides_from_reads_total.sh {svg_filename} {tmp_info_filename} > {tmp_out_total_filename}"
    result = subprocess.run(cmd, shell=True)
    if result.returncode != 0:
        raise

    tmp_out_orange_filename = tempfile.NamedTemporaryFile(delete=False).name
    cmd = f"bash /opt/str/count_nucleotides_from_reads_orange.sh {svg_filename} {tmp_info_filename} > {tmp_out_orange_filename}"
    result = subprocess.run(cmd, shell=True)
    if result.returncode != 0:
        raise

    total_counts_df = pd.read_csv(tmp_out_total_filename, delimiter="\t", index_col=False)
    total_counts_df[LOCUS_ID_COL_NAME] = [locus_id] * total_counts_df.shape[0]

    orange_counts_df = pd.read_csv(tmp_out_orange_filename, delimiter="\t", index_col=False)
    orange_counts_df[LOCUS_ID_COL_NAME] = [locus_id] * orange_counts_df.shape[0]

    nucleotides_count_df = pd.merge(total_counts_df, orange_counts_df, on=LOCUS_ID_COL_NAME)
    nucleotides_count_df[SAMPLE_ID_COL_NAME] = [sample_id] * nucleotides_count_df.shape[0]

    os.remove(tmp_info_filename)
    os.remove(tmp_out_total_filename)
    os.remove(tmp_out_orange_filename)
    return nucleotides_count_df


def extract_metrics(
        sample_id: str,
        realigned_reads: str,
        reference: str,
        genotypes: str,
        variant_catalog: str,
        output_filename: str) -> None:

    loci = get_loci(variant_catalog)
    write_header = True

    # This is added to address a corner case where the realigned bam
    # does not have reads overlapping with any of the given target loci.
    open(output_filename, "w").close()

    for locus_id in loci:
        _, svg_f, metrics_df, phasing_df = call_reviewer(
            sample_id, realigned_reads, reference, genotypes, variant_catalog, locus_id)

        if any(x is None for x in (svg_f, metrics_df, phasing_df)):
            continue

        start_time = time.time()
        ex = get_metrics_from_image(svg_f)
        print(f"Get metrics from image runtime: {time.time() - start_time:.2f} seconds")

        metrics_df["Times_Upper_Interruptions"] = ex.total_upper_interrupt
        metrics_df["Count_Upper_reads_with_interruptions"] = len(ex.results_total_upper)
        metrics_df["Times_Upper_Gaps"] = ex.gap_upper_count
        metrics_df["Times_Bottom_Interruptions"] = ex.total_bottom_interrupt
        metrics_df["Count_Bottom_reads_with_interruptions"] = len(ex.results_total_bottom)
        metrics_df["Times_Bottom_Gaps"] = ex.gap_bottom_count

        start_time = time.time()
        nucleotides_count_df = get_nucleotides_count(svg_f, ex, sample_id, locus_id)
        print(f"Get nucleotides count from image runtime: {time.time() - start_time:.2f} seconds")
        metrics_df = pd.merge(metrics_df, nucleotides_count_df, on=[SAMPLE_ID_COL_NAME, LOCUS_ID_COL_NAME])

        # It could be ideal to collect all the metrics of all the loci in one dataframe;
        # however, that may lead to out-of-memory issues if processing many loci.
        # Hence, the measurements of each locus is appended to a file.
        metrics_df.to_csv(output_filename, mode="a", header=write_header, sep="\t", index=False)
        write_header = False


def main():
    parser = argparse.ArgumentParser(
        description="Extracts STR-expansion related metrics of given locus. "
                    "This script has the following requirements."
                    "(1) It requires Illumina/REViewer tool is "
                    "installed and can be invoked using `reviewer`."
                    "(2) It needs the utility scripts `/opt/str/count_nucleotides_from_reads_total.sh`"
                    "(3) It needs the utility scripts `/opt/str/count_nucleotides_from_reads_orange.sh`",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "-s", "--sample-id",
        help="The ID of the sample."
    )

    parser.add_argument(
        "-i", "--realigned-bam",
        help="The file containing realigned read in the BAM format. "
             "This full MUST be the output of ExpansionHunter, "
             "other BAM files are not supported. "
    )

    parser.add_argument(
        "-r", "--reference",
        help="The filename of the reference sample."
    )

    parser.add_argument(
        "-g", "--genotypes",
        help="The name of the file containing genotypes in the VCF format."
    )

    parser.add_argument(
        "-c", "--variant-catalog",
        help="The filename of a file containing variant catalogs in JSON format. "
             "Metrics will be generated for the locus in this file. The file is expected "
             "to contain a list of dictionaries, each dictionary must have a value for the key `LocusId`."
    )

    parser.add_argument(
        "-o", "--output",
        default="output.csv",
        help="The output CSV filename."
    )

    args = parser.parse_args()
    start_time = time.time()
    extract_metrics(
        args.sample_id,
        args.realigned_bam,
        args.reference,
        args.genotypes,
        args.variant_catalog,
        args.output)
    print(f"Total runtime: {time.time() - start_time:.2f} seconds")


if __name__ == '__main__':
    main()
