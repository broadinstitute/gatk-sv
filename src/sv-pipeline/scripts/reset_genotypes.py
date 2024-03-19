#!/bin/python

import argparse
import gzip
import logging
import sys

from collections import defaultdict
from typing import List, Text, Optional

import pysam


RESET_PESR_FORMATS_DICT = {
    "SR_GT": 0,
    "SR_GQ": 99,
    "PE_GT": 0,
    "PE_GQ": 99
}

RESET_RD_GQ_VALUE = 99

_gt_set_hom_ref_map = dict()


def _cache_gt_set_hom_ref(gt):
    s = _gt_set_hom_ref_map.get(gt, None)
    if s is None:
        if all(a is None for a in gt):
            s = tuple(None for _ in gt)
        else:
            s = tuple(0 for _ in gt)
        _gt_set_hom_ref_map[gt] = s
    return s


def reset_format_if_exists(gt, key, value):
    if key in gt:
        gt[key] = value


def read_gzip_text_mode(path):
    return gzip.open(path, mode="rt")


def read_bed(path):
    bed_func = read_gzip_text_mode if path.endswith(".gz") else open
    with bed_func(path) as f:
        bed_data_sets = defaultdict(set)
        for line in f:
            tokens = line.strip().split("\t")
            vid = tokens[0]
            sample = tokens[1]
            bed_data_sets[vid].add(sample)
        return {key: list(val) for key, val in bed_data_sets.items()}


def process_vcf(in_path, out_path, bed_data):
    with pysam.VariantFile(in_path) as fin, pysam.VariantFile(out_path, mode="w", header=fin.header) as fout:
        current_chrom = None
        for record in fin:
            if record.chrom != current_chrom:
                current_chrom = record.chrom
                logging.info(f"  {record.chrom}")
            if record.id in bed_data:
                for sample in bed_data[record.id]:
                    if sample not in record.samples:
                        raise ValueError(f"Sample {sample} not found in the vcf")
                    gt = record.samples[sample]
                    gt["GT"] = _cache_gt_set_hom_ref(gt["GT"])
                    # Reset RD genotyping fields but not PESR since we did not re-examine that evidence
                    # Note we do not take PAR into account here, in order to match the rest of the pipeline
                    reset_format_if_exists(gt, "RD_CN", len(gt["GT"]))
                    reset_format_if_exists(gt, "RD_GQ", RESET_RD_GQ_VALUE)
                    reset_format_if_exists(gt, "CN", len(gt["GT"]))
                    reset_format_if_exists(gt, "CNQ", RESET_RD_GQ_VALUE)
            fout.write(record)


def _parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Resets genotypes to homozygous-reference",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--vcf', type=str, required=True, help='Input vcf')
    parser.add_argument('--genotype-bed', type=str, required=True,
                        help='Bed file of genotypes to reset, with variant and sample ID columns (.bed or .bed.gz)')
    parser.add_argument('--out', type=str, required=True, help='Output vcf')
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

    logging.info("Reading bed file...")
    bed_data = read_bed(args.genotype_bed)

    logging.info("Processing vcf...")
    process_vcf(args.vcf, args.out, bed_data)
    pysam.tabix_index(args.out, preset="vcf", force=True)


if __name__ == "__main__":
    main()
