#!/bin/python

import argparse
import sys
import pysam
import gzip
from collections import defaultdict
from typing import List, Text, Optional

from intervaltree import IntervalTree, Interval

SEX_MALE = "M"
SEX_FEMALE = "F"
SEX_UNKNOWN = "U"

# Creates dictionary of trees[sv_type][contig] from bed file
def create_trees_from_bed_records(bed_path):
    trees = {}
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


class GDRMatch:

    def __init__(self, name, chrom, pos, end, sample, frac, overlaps):
        self.name = name
        self.chrom = chrom
        self.pos = pos
        self.end = end
        self.sample = sample
        self.frac = frac
        self.overlaps = overlaps

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return "(" + ", ".join(str(x) for x in [self.name, self.chrom, self.pos, self.end, self.sample, self.frac, self.overlaps]) + ")"


def get_overlapping_samples_bed(bed_path, intervals, medians, cutoff):

    bed_overlappers = dict()
    with gzip.open(bed_path) as f:
        for line in f:
            record = line.decode().strip().split('\t')
            if record[0].startswith("#"):
                continue
            chrom = record[0]
            pos = int(record[1])
            stop = int(record[2])
            name = record[3]
            sample = record[4]
            if chrom not in intervals:
                continue
            if sample not in bed_overlappers:
                bed_overlappers[sample] = dict()
            intersect = intervals[chrom].overlap(pos, stop)
            for interval in intersect:
                interval_name = interval.data
                # Only load potential false positives
                if sample in medians and interval_name in medians[sample]:
                    continue
                if interval_name not in bed_overlappers[sample]:
                    bed_overlappers[sample][interval_name] = list()
                bed_overlappers[sample][interval_name].append((pos, stop, name))

    intervals_dict = {a.data: (a.begin, a.end, chrom) for chrom in intervals for a in intervals[chrom]}
    bed_matches = list()
    for sample in bed_overlappers:
        for interval_name in bed_overlappers[sample]:
            tree = IntervalTree()
            for call in bed_overlappers[sample][interval_name]:
                tree.addi(call[0], call[1])
            # Merge so we don't double-count
            tree.merge_overlaps()
            lookup_interval = intervals_dict[interval_name]
            overlap = sum([overlap_size(call, lookup_interval) for call in tree])
            overlap_frac = overlap / (lookup_interval[1] - lookup_interval[0])
            if overlap_frac > cutoff:
                bed_matches.append(GDRMatch(name=interval_name, chrom=lookup_interval[2], pos=lookup_interval[0],
                                        end=lookup_interval[1], sample=sample, frac=overlap_frac,
                                        overlaps=bed_overlappers[sample][interval_name]))
    return bed_matches

def get_overlapping_samples_vcf(vcf_path, intervals, medians, cutoff, vcf_min_size):
    _gt_sum_map = dict()

    def _cache_gt_sum(gt):
        if gt is None:
            return 0
        s = _gt_sum_map.get(gt, None)
        if s is None:
            s = sum([1 for a in gt if a is not None and a > 0])
            _gt_sum_map[gt] = s
        return s

    # TODO abstract above blocks and reuse here
    vcf_overlappers = dict()
    with pysam.VariantFile(vcf_path) as vcf:
        for record in vcf:
            chrom = record.chrom
            pos = record.pos
            stop = record.stop
            name = record.id
            svtype = record.info.get("SVTYPE", "")
            if svtype not in intervals.keys():
                continue
            if record.info.get("SVLEN", 0) < vcf_min_size:
                continue
            if record.chrom not in intervals[svtype]:
                continue
            carriers = [s for s in record.samples if _cache_gt_sum(record.samples[s]["GT"]) > 0]
            intersect = intervals[svtype][chrom].overlap(pos, stop)
            for sample in carriers:
                if sample not in vcf_overlappers:
                    vcf_overlappers[sample] = dict()
                for interval in intersect:
                    interval_name = interval.data
                    # Only load potential false positives
                    if sample in medians and interval_name in medians[sample]:
                        continue
                    if interval_name not in vcf_overlappers[sample]:
                        vcf_overlappers[sample][interval_name] = list()
                    vcf_overlappers[sample][interval_name].append((pos, stop, name))

    for svtype in ["DEL", "DUP"]:
        intervals_dict = {a.data: (a.begin, a.end, chrom) for chrom in intervals[svtype] for a in intervals[svtype][chrom]}
        vcf_matches = list()
        for sample in vcf_overlappers:
            for interval_name in vcf_overlappers[sample]:
                tree = IntervalTree()
                for call in vcf_overlappers[sample][interval_name]:
                    tree.addi(call[0], call[1])
                # Merge so we don't double-count
                tree.merge_overlaps()
                if interval_name not in intervals_dict:
                    continue
                lookup_interval = intervals_dict[interval_name]
                overlap = sum([overlap_size(call, lookup_interval) for call in tree])
                overlap_frac = overlap / (lookup_interval[1] - lookup_interval[0])
                if overlap_frac > cutoff:
                    vcf_matches.append(GDRMatch(name=interval_name, chrom=lookup_interval[2], pos=lookup_interval[0],
                                                end=lookup_interval[1], sample=sample, frac=overlap_frac,
                                                overlaps=vcf_overlappers[sample][interval_name]))
        yield vcf_matches


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


def read_median_geno(path, del_ids, dup_ids, del_cutoff, dup_cutoff, sample_sex_dict, chr_x, chr_y):
    with open(path) as f:
        header = f.readline().strip().split('\t')
        samples = header[4:]
        data = dict()
        for line in f:
            record = line.strip().split('\t')
            chrom = record[0]
            vid = record[3]
            if vid in del_ids:
                isDel = True
            elif vid in dup_ids:
                isDel = False
            else:
                continue
            for i, median in enumerate(record[4:]):
                median = float(median)
                if isDel:
                    cutoff = del_cutoff
                else:
                    cutoff = dup_cutoff
                sample = samples[i]
                if chrom == chr_x or chrom == chr_y:
                    sex = sample_sex_dict[sample]
                    if sex == SEX_MALE:
                        cutoff -= 0.5
                    elif sex == SEX_FEMALE and chrom == chr_y:
                        continue
                    elif sex == SEX_UNKNOWN:
                        continue
                if isDel and median > cutoff:
                    continue
                elif (not isDel) and median < cutoff:
                    continue
                if sample not in data:
                    data[sample] = dict()
                data[sample][vid] = median
    return data


def _parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Cleans up CNVs in genomic disorder regions",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--del-bed', type=str, help='Raw DEL bed file (gzipped)')
    parser.add_argument('--dup-bed', type=str, help='Raw DUP bed file (gzipped)')
    parser.add_argument('--vcf', type=str, help='Final vcf')
    parser.add_argument('--region-bed', type=str, help='Genomic disorder regions bed')
    parser.add_argument('--median-geno', type=str, help='RDTest ".median_geno" file')
    parser.add_argument('--ped-file', type=str, help='Ped file')
    parser.add_argument('--out', type=str, help='Output base path')
    parser.add_argument('--overlap', type=float, default=0.3, help='Overlap fraction cutoff')
    parser.add_argument('--del-median', type=float, default=0.6, help='DEL median cutoff')
    parser.add_argument('--dup-median', type=float, default=1.4, help='DUP median cutoff')
    parser.add_argument('--vcf-min-size', type=int, default=1000, help='Min size for vcf variants')
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
    sample_sex_dict = read_ped_file(args.ped_file)
    gdr_trees, gdr_del_ids, gdr_dup_ids = create_trees_from_bed_records(args.region_bed)
    medians = read_median_geno(args.median_geno, del_ids=gdr_del_ids, dup_ids=gdr_dup_ids,
                               del_cutoff=args.del_median, dup_cutoff=args.dup_median,
                               sample_sex_dict=sample_sex_dict, chr_x=args.chr_x, chr_y=args.chr_y)
    print(medians)
    matches_bed_del = get_overlapping_samples_bed(bed_path=args.del_bed, intervals=gdr_trees["DEL"],
                                                  medians=medians, cutoff=args.overlap)
    matches_bed_dup = get_overlapping_samples_bed(bed_path=args.dup_bed, intervals=gdr_trees["DUP"],
                                                  medians=medians, cutoff=args.overlap)
    matches_vcf_del, matches_vcf_dup = get_overlapping_samples_vcf(vcf_path=args.vcf, intervals=gdr_trees,
                                                                   medians=medians, cutoff=args.overlap,
                                                                   vcf_min_size=args.vcf_min_size)

    vcf_fp_del = defaultdict(set)
    for record in matches_vcf_del:
        vcf_fp_del[record.sample].add(record.name)
    vcf_fp_dup = defaultdict(set)
    for record in matches_vcf_dup:
        vcf_fp_dup[record.sample].add(record.name)

    # Raw bed calls corresponding to VCF calls that are FPs against RDTest
    # These need to be checked against the VCF call breakpoints
    matches_bed_del_vcf_fp = [record for record in matches_bed_del if record.sample in vcf_fp_del and record.name in vcf_fp_del[record.sample]]
    matches_bed_dup_vcf_fp = [record for record in matches_bed_dup if record.sample in vcf_fp_dup and record.name in vcf_fp_dup[record.sample]]

    with open(args.out + ".vcf_del.bed", "w") as f:
        for record in matches_vcf_del:
            for ov in record.overlaps:
                f.write(f"{record.chrom}\t{ov[0]}\t{ov[1]}\t{record.sample}\n")
    with open(args.out + ".vcf_dup.bed", "w") as f:
        for record in matches_vcf_dup:
            for ov in record.overlaps:
                f.write(f"{record.chrom}\t{ov[0]}\t{ov[1]}\t{record.sample}\n")
    with open(args.out + ".bed_del.bed", "w") as f:
        for record in matches_bed_del_vcf_fp:
            for ov in record.overlaps:
                f.write(f"{record.chrom}\t{ov[0]}\t{ov[1]}\t{record.sample}\n")
    with open(args.out + ".bed_dup.bed", "w") as f:
        for record in matches_bed_dup_vcf_fp:
            for ov in record.overlaps:
                f.write(f"{record.chrom}\t{ov[0]}\t{ov[1]}\t{record.sample}\n")


if __name__ == "__main__":
    main()
