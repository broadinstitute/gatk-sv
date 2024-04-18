#!/bin/python

import argparse
import logging
import sys

from collections import defaultdict
from operator import itemgetter
from typing import List, Text, Optional

from intervaltree import Interval, IntervalTree
import pysam

GENOMIC_DISORDER_KEY = "GENOMIC_DISORDER"
GENOMIC_DISORDER_HEADER = \
    f"##INFO=<ID={GENOMIC_DISORDER_KEY},Number=1,Type=String,Description=\"Genomic disorder region\">"

_gt_sum_map = dict()


def _cache_gt_sum(gt):
    if gt is None:
        return 0
    s = _gt_sum_map.get(gt, None)
    if s is None:
        s = sum([1 for a in gt if a is not None and a > 0])
        _gt_sum_map[gt] = s
    return s


# Creates dictionary of trees[sv_type][contig] from bed file
def create_trees_from_bed_records_by_type(bed_path):
    trees = defaultdict(dict)
    variant_types = ["DEL", "DUP"]
    for svtype in variant_types:
        trees[svtype] = defaultdict(IntervalTree)
    del_ids = set()
    dup_ids = set()
    with open(bed_path) as f:
        for line in f:
            record = line.strip().split('\t')
            contig = record[0]
            start = int(record[1])
            stop = int(record[2])
            name = record[3]
            svtype = record[4]
            if svtype not in trees:
                raise ValueError("Unexpected SVTYPE in bed file: %s" % type)
            if svtype == "DEL":
                del_ids.add(name)
            elif svtype == "DUP":
                dup_ids.add(name)
            trees[svtype][contig].addi(start, stop, name)
    return trees, del_ids, dup_ids


def overlap_size(interval_a, interval_b):
    return max(0, min(interval_a[1], interval_b[1]) - max(interval_a[0], interval_b[0]))


def reciprocal_overlap(interval_a, interval_b):
    return overlap_size(interval_a, interval_b) / max(interval_a[1] - interval_a[0], interval_b[1] - interval_b[0])


def _parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description=f"Annotates variants with {GENOMIC_DISORDER_KEY} INFO field and produces bed files for "
                    f"RDTest plotting",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--vcf', type=str, required=True, help='Input vcf, GATK-formatted with ECN format fields')
    parser.add_argument('--region-bed', type=str, required=True, help='Preprocessed genomic disorder regions bed')
    parser.add_argument('--out', type=str, required=True, help='Output base path')
    parser.add_argument('--overlap', type=float, default=0.5, help='Reciprocal overlap cutoff')
    parser.add_argument('--plot-padding', type=float, default=0.5, help='Padding for RDTest bed as a fraction of '
                                                                        'the variant length')
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
    gdr_trees, gdr_del_ids, gdr_dup_ids = create_trees_from_bed_records_by_type(bed_path=args.region_bed)
    vcf_out_path = f"{args.out}.vcf.gz"
    variant_manifest_path = f"{args.out}.manifest.tsv"
    variant_records_bed_path = f"{args.out}.padded_variants.bed"
    region_records_bed_path = f"{args.out}.padded_regions.bed"
    region_data = {interval.data: (svtype, chrom, interval.begin, interval.end, set())
                   for svtype in gdr_trees for chrom in gdr_trees[svtype] for interval in gdr_trees[svtype][chrom]}
    with pysam.VariantFile(args.vcf) as vcf_in:
        header = vcf_in.header
        header.add_line(GENOMIC_DISORDER_HEADER)
        with pysam.VariantFile(vcf_out_path, mode="w", header=vcf_in.header) as vcf_out, \
                open(variant_manifest_path, "w") as manifest_out, \
                open(variant_records_bed_path, "w") as bed_out:
            manifest_out.write(f"#CHROM\tPOS\tEND\tVARIANT_ID\tSVTYPE\tREGION_ID\tREGION_POS\tREGION_END\tSAMPLES\n")
            for record in vcf_in:
                svtype = record.info.get("SVTYPE", "")
                if svtype in gdr_trees and record.chrom in gdr_trees[svtype]:
                    record_interval = Interval(record.pos, record.stop)
                    tree = gdr_trees[svtype][record.chrom]
                    overlappers = sorted([(ov, reciprocal_overlap(ov, record_interval))
                                           for ov in tree.overlap(record_interval)], key=itemgetter(1))
                    if len(overlappers) > 0 and overlappers[-1][1] >= args.overlap:
                        region_start = overlappers[-1][0].begin
                        region_stop = overlappers[-1][0].end
                        region_name = overlappers[-1][0].data
                        record.info[GENOMIC_DISORDER_KEY] = region_name
                        logging.info(f"{record.id} : {region_name}")
                        carriers = sorted([s for s, gt in record.samples.items() if _cache_gt_sum(gt["GT"]) > 0])
                        if len(carriers) > 0:
                            # Exclude variants with no carriers
                            sample_col = ",".join(carriers)
                            manifest_out.write(f"{record.chrom}\t{record.pos}\t{record.stop}\t{record.id}\t{svtype}\t"
                                               f"{region_name}\t{region_start}\t{region_stop}\t{sample_col}\n")
                            padding = int(args.plot_padding * (record.stop - record.pos))
                            padded_pos = max(0, record.pos - padding)
                            padded_stop = record.stop + padding
                            bed_out.write(
                                f"{record.chrom}\t{padded_pos}\t{padded_stop}\t{record.id}\t{svtype}\t{sample_col}\n")
                            for c in carriers:
                                region_data[region_name][4].add(c)
                vcf_out.write(record)
    pysam.tabix_index(vcf_out_path, preset="vcf", force=True)
    # Write GDRs
    with open(region_records_bed_path, "w") as bed_out:
        for region_name in region_data:
            svtype, chrom, pos, stop, carriers = region_data[region_name]
            if len(carriers) > 0:
                padding = int(args.plot_padding * (stop - pos))
                padded_pos = max(0, pos - padding)
                padded_stop = stop + padding
                sample_col = ",".join(carriers)
                bed_out.write(f"{chrom}\t{padded_pos}\t{padded_stop}\t{region_name}\t{svtype}\t{sample_col}\n")


if __name__ == "__main__":
    main()
