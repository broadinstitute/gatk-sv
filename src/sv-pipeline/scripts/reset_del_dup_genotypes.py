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
RESET_GQ_VALUE = 99

_gt_set_hom_ref_map = dict()
_gt_set_het_map = dict()
_gt_set_hom_var_map = dict()


def _cache_gt_set_hom_ref(gt):
    s = _gt_set_hom_ref_map.get(gt, None)
    if s is None:
        s = tuple(0 for _ in gt)
        _gt_set_hom_ref_map[gt] = s
    return s


def _cache_gt_set_het(gt):
    s = _gt_set_het_map.get(gt, None)
    if s is None:
        s = list(0 for _ in gt)
        if len(s) > 0:
            s[-1] = 1
        _gt_set_het_map[gt] = s
    return s


def _cache_gt_set_hom_var(gt):
    s = _gt_set_hom_var_map.get(gt, None)
    if s is None:
        s = tuple(1 for _ in gt)
        _gt_set_hom_var_map[gt] = s
    return s


def reset_format_if_exists(gt, key, value):
    if key in gt:
        gt[key] = value


def reset_format_if_not_exists(gt, key, value):
    if gt.get(key, None) is None:
        gt[key] = value


def read_gzip_text_mode(path):
    return gzip.open(path, mode="rt")


def read_tsv(path):
    if path is None:
        return dict()
    open_func = read_gzip_text_mode if path.endswith(".gz") else open
    with open_func(path) as f:
        data_sets = defaultdict(set)
        for line in f:
            tokens = line.strip().split("\t")
            if len(tokens) != 4:
                raise ValueError(f"Encountered record without 4 columns: {tokens}")
            vid = tokens[1]
            sample = tokens[2]
            genotype = int(tokens[3])
            data_sets[vid].add((sample, genotype))
        return {key: list(val) for key, val in data_sets.items()}


def reset_format_if_exists(gt, key, value):
    if key in gt:
        gt[key] = value


def get_ecn(gt):
    ecn = gt.get("ECN", None)
    if ecn is None:
        raise ValueError("Missing ECN format field")
    return ecn


def reset_format_fields(gt, n_alt_alleles):
    if any(a is None for a in gt["GT"]):
        raise ValueError(f"Attempted to reset a no-call genotype")
    # Note we do not take PAR into account here to match the rest of the pipeline
    ecn = get_ecn(gt)
    if ecn == 0:
        # Should already be empty
        return
    gt["GQ"] = RESET_GQ_VALUE
    gt["RD_CN"] = ecn
    gt["RD_GQ"] = RESET_RD_GQ_VALUE
    if n_alt_alleles == 0:
        gt["GT"] = _cache_gt_set_hom_ref(gt["GT"])
    elif n_alt_alleles == 1:
        gt["GT"] = _cache_gt_set_het(gt["GT"])
    elif n_alt_alleles == 2:
        gt["GT"] = _cache_gt_set_hom_var(gt["GT"])
    else:
        raise ValueError("Unsupported genotype code " + n_alt_alleles + ", must be in {0, 1, 2}")


def process_vcf(in_path, out_path, genotype_data):
    with pysam.VariantFile(in_path) as fin, pysam.VariantFile(out_path, mode="w", header=fin.header) as fout:
        current_chrom = None
        for record in fin:
            if record.chrom != current_chrom:
                current_chrom = record.chrom
                logging.info(f"  {record.chrom}")
            if record.id in genotype_data:
                svtype = record.info.get("SVTYPE", "None")
                if svtype != "DEL" and svtype != "DUP":
                    raise ValueError(f"Record {record.id} has SVTYPE {svtype} but must be DEL or DUP")
                for sample, n_alt_alleles in genotype_data[record.id]:
                    if sample not in record.samples:
                        raise ValueError(f"Sample {sample} not found in the vcf")
                    gt = record.samples[sample]
                    reset_format_fields(gt, n_alt_alleles)
            fout.write(record)


def _parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Resets DEL/DUP genotypes to homozygous-reference",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--vcf', type=str, required=True, help='Input vcf')
    parser.add_argument('--genotype-tsv', type=str, required=False,
                        help='If provided, genotypes to reset. Headerless, with chrom, variant, and sample ID columns '
                             '(.tsv or .tsv.gz)')
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
    genotype_data = read_tsv(args.genotype_tsv)

    logging.info("Processing vcf...")
    process_vcf(in_path=args.vcf, out_path=args.out, genotype_data=genotype_data)
    pysam.tabix_index(args.out, preset="vcf", force=True)


if __name__ == "__main__":
    main()
