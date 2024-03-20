#!/bin/python

import argparse
import gzip
import logging
import sys

from collections import defaultdict
from typing import List, Text, Optional

import pysam


SEX_MALE = "M"
SEX_FEMALE = "F"
SEX_UNKNOWN = "U"

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


def reset_format_if_not_exists(gt, key, value):
    if gt.get(key, None) is None:
        gt[key] = value


def read_gzip_text_mode(path):
    return gzip.open(path, mode="rt")


def read_tsv(path):
    open_func = read_gzip_text_mode if path.endswith(".gz") else open
    with open_func(path) as f:
        data_sets = defaultdict(set)
        for line in f:
            tokens = line.strip().split("\t")
            if len(tokens) != 2:
                raise ValueError(f"Encountered record with more than 2 columns: {tokens}")
            vid = tokens[0]
            sample = tokens[1]
            data_sets[vid].add(sample)
        return {key: list(val) for key, val in data_sets.items()}


def get_rd_cn(chrom, sample_sex, chr_x, chr_y):
    if sample_sex == SEX_UNKNOWN:
        return None
    elif chrom != chr_x and chrom != chr_y:
        return 2
    elif sample_sex == SEX_MALE:
        return 1
    elif sample_sex == SEX_FEMALE:
        if chrom == chr_x:
            return 2
        else:
            # chrY
            return None
    else:
        raise ValueError(f"Unknown sex assignment {sample_sex} (bug)")


def reset_rd_format_fields(gt, sample, chrom, sample_sex_dict, chr_x, chr_y):
    # Reset RD genotyping fields but not PESR since we did not re-examine that evidence
    # Note we do not take PAR into account here, in order to match the rest of the pipeline
    sample_sex = sample_sex_dict.get(sample, None)
    if sample_sex is None:
        raise ValueError(f"No sex assignment for sample {sample}, check ped file")
    gt["RD_CN"] = get_rd_cn(chrom, sample_sex, chr_x, chr_y)
    gt["RD_GQ"] = RESET_RD_GQ_VALUE
    print(f"{sample} {sample_sex} {gt['GT']} {gt['RD_CN']} {gt['RD_GQ']}")


def process_vcf(in_path, out_path, genotype_data, reset_rd_genotype, sample_sex_dict, chr_x, chr_y):
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
                for sample in genotype_data[record.id]:
                    if sample not in record.samples:
                        raise ValueError(f"Sample {sample} not found in the vcf")
                    gt = record.samples[sample]
                    gt["GT"] = _cache_gt_set_hom_ref(gt["GT"])
                    reset_rd_format_fields(gt, sample, record.chrom, sample_sex_dict, chr_x, chr_y)
                if reset_rd_genotype:
                    for sample, gt in record.samples.items():
                        rd_cn = gt.get("RD_CN", None)
                        rd_gq = gt.get("RD_GQ", None)
                        if rd_cn is None or rd_gq is None:
                            reset_rd_format_fields(gt, sample, record.chrom, sample_sex_dict, chr_x, chr_y)
            fout.write(record)


def read_ped_file(path):
    with open(path) as f:
        data = dict()
        for line in f:
            record = line.strip().split('\t')
            sample = record[1]
            x_ploidy = record[4]
            if x_ploidy == "1":
                data[sample] = SEX_MALE
            elif x_ploidy == "2":
                data[sample] = SEX_FEMALE
            else:
                data[sample] = SEX_UNKNOWN
    return data


def _parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Resets DEL/DUP genotypes to homozygous-reference",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--vcf', type=str, required=True, help='Input vcf')
    parser.add_argument('--genotype-tsv', type=str, required=False,
                        help='If provided, genotypes to reset. Headerless, with variant and sample ID columns '
                             '(.tsv or .tsv.gz)')
    parser.add_argument('--ped-file', type=str, required=True, help='Ped file')
    parser.add_argument('--reset-rd-genotype', action='store_true', help='Reset RD_CN/RD_GQ FORMAT fields if empty')
    parser.add_argument('--out', type=str, required=True, help='Output vcf')
    parser.add_argument('--chr-x', type=str, default="chrX", help='Chromosome X identifier')
    parser.add_argument('--chr-y', type=str, default="chrY", help='Chromosome Y identifier')
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

    logging.info("Reading ped file...")
    sample_sex_dict = read_ped_file(args.ped_file)

    logging.info("Reading bed file...")
    genotype_data = read_tsv(args.genotype_tsv)

    logging.info("Processing vcf...")
    process_vcf(in_path=args.vcf, out_path=args.out, genotype_data=genotype_data,
                reset_rd_genotype=args.reset_rd_genotype, sample_sex_dict=sample_sex_dict,
                chr_x=args.chr_x, chr_y=args.chr_y)
    pysam.tabix_index(args.out, preset="vcf", force=True)


if __name__ == "__main__":
    main()
