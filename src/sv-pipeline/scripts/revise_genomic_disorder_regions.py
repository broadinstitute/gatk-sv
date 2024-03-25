#!/bin/python

import argparse
import gzip
import logging
import math
import subprocess
import sys
import tempfile

from collections import defaultdict
from collections.abc import Iterable
from itertools import groupby
from operator import attrgetter
from typing import List, Text, Optional

from intervaltree import IntervalTree
import pysam

SEX_MALE = "M"
SEX_FEMALE = "F"
SEX_UNKNOWN = "U"

# Delimiter suffix appended to the end of interval IDs before the index, e.g. "intervalA__0", "intervalA__1", ...
INDEX_DELIMITER = "__"

RESET_PESR_FORMATS_DICT = {
    "SR_GT": 0,
    "SR_GQ": 99,
    "PE_GT": 0,
    "PE_GQ": 99
}

RESET_RD_GQ_VALUE = 99
RESET_GQ_VALUE = 99

RESCUE_RD_GQ_VALUE = 99
RESCUE_GQ_VALUE = 99

_gt_sum_map = dict()
_gt_set_hom_ref_map = dict()
_gt_set_het_map = dict()


def _cache_gt_sum(gt):
    if gt is None:
        return 0
    s = _gt_sum_map.get(gt, None)
    if s is None:
        s = sum([1 for a in gt if a is not None and a > 0])
        _gt_sum_map[gt] = s
    return s


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


def create_region_intervals_dict(trees):
    regions_dict = dict()
    # Get interval sets for each region
    for svtype, svtype_dict in trees.items():
        if svtype not in regions_dict:
            regions_dict[svtype] = dict()
        for contig, contig_tree in svtype_dict.items():
            for interval in contig_tree:
                region_subinterval = interval.data
                region, _ = get_base_name_and_index(region_subinterval)
                if region not in regions_dict[svtype]:
                    regions_dict[svtype][region] = list()
                regions_dict[svtype][region].append((contig, interval.begin, interval.end))
    # Merge intervals using the min/max endpoints, assuming they are all contiguous
    for svtype, svtype_dict in regions_dict.items():
        for region, interval_list in svtype_dict.items():
            contigs = list(set(interval[0] for interval in interval_list))
            if len(contigs) > 1:
                raise ValueError(f"Region {region} is on multiple contigs: {contigs}")
            svtype_dict[region] = (contigs[0], min(interval[1] for interval in interval_list),
                                   max(interval[2] for interval in interval_list))
    return regions_dict


def create_trees_from_bed_records(bed_path):
    trees = defaultdict(IntervalTree)
    with open(bed_path) as f:
        for line in f:
            record = line.strip().split('\t')
            contig = record[0]
            start = int(record[1])
            stop = int(record[2])
            trees[contig].addi(start, stop)
    return trees


def overlap_size(interval_a, interval_b):
    return max(0, min(interval_a[1], interval_b[1]) - max(interval_a[0], interval_b[0]))


class GDRMatch:

    def __init__(self, name, chrom, pos, end, sample, overlapping_region_intervals, valid_region_intervals,
                 supporting_intervals, region, n_region_subdivisions, genotype_medians):
        self.name = name
        self.chrom = chrom
        self.pos = pos
        self.end = end
        self.sample = sample
        self.overlapping_region_intervals = overlapping_region_intervals
        self.valid_region_intervals = valid_region_intervals
        self.supporting_intervals = supporting_intervals
        self.region = region
        self.n_region_subdivisions = n_region_subdivisions,
        self.genotype_medians = genotype_medians

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return "(" + ", ".join(str(x) for x in [self.name, self.chrom, self.pos, self.end, self.sample, self.region,
                                                len(self.supporting_intervals), len(self.valid_region_intervals),
                                                len(self.overlapping_region_intervals),
                                                self.n_region_subdivisions]) + ")"


class GenotypeMedian:
    def __init__(self, median, threshold, hom_ref_cn):
        self.median = median
        self.threshold = threshold
        self.hom_ref_cn = hom_ref_cn


def get_base_name_and_index(interval_id):
    tokens = interval_id.split(INDEX_DELIMITER)
    if len(tokens) != 2:
        raise ValueError(f"Expected to find \"{INDEX_DELIMITER}\" in interval name: {interval_id}")
    return tokens[0], int(tokens[1])


def add_match(matches_dict, svtype, name, chrom, pos, end, sample, overlapping_region_intervals, valid_region_intervals,
              supporting_intervals, region, n_region_subdivisions, genotype_medians):
    if svtype not in matches_dict:
        matches_dict[svtype] = dict()
    if sample not in matches_dict[svtype]:
        matches_dict[svtype][sample] = dict()
    if name not in matches_dict[svtype][sample]:
        matches_dict[svtype][sample][name] = list()
    matches_dict[svtype][sample][name].append(
        GDRMatch(name=name, chrom=chrom, pos=pos, end=end, sample=sample,
                 overlapping_region_intervals=overlapping_region_intervals,
                 valid_region_intervals=valid_region_intervals,
                 supporting_intervals=supporting_intervals,
                 region=region,
                 n_region_subdivisions=n_region_subdivisions,
                 genotype_medians=genotype_medians)
    )


def get_n_alt_alleles(genotype_median):
    # Returns either 1 or 2
    return min(2, max(1, math.ceil(math.fabs(genotype_median.median - genotype_median.threshold) / 0.5)))


def get_overlapping_samples_vcf(vcf_path, intervals, medians, median_vids, cutoff, vcf_min_size, min_supported_valid_overlapping_intervals_frac,
                                min_region_overlap, min_supported_region_intervals_frac):

    interval_subdivision_counts = dict()
    unique_subdivision_names = set()
    for svtype in intervals:
        for chrom in intervals[svtype]:
            for interval in intervals[svtype][chrom]:
                if interval.data in unique_subdivision_names:
                    raise ValueError(f"Encountered duplicate interval id: {interval.data}")
                unique_subdivision_names.add(interval.data)
                name, _ = get_base_name_and_index(interval.data)
                if name not in interval_subdivision_counts:
                    interval_subdivision_counts[name] = 0
                interval_subdivision_counts[name] += 1

    false_positive_matches = dict()
    false_negative_matches = dict()
    vcf_vids = set()
    with pysam.VariantFile(vcf_path) as vcf:
        current_chrom = None
        for record in vcf:
            chrom = record.chrom
            if chrom != current_chrom:
                current_chrom = chrom
                logging.info(f"  Processing {chrom}")
            pos = record.pos
            stop = record.stop
            name = record.id
            if name in vcf_vids:
                raise ValueError(f"All variant ids in the input vcf must be unique; encountered duplicate: {name}")
            vcf_vids.add(name)
            svtype = record.info.get("SVTYPE", "")
            svlen = record.info.get("SVLEN", 0)
            if svtype not in intervals.keys():
                continue
            if svlen < vcf_min_size:
                continue
            if record.chrom not in intervals[svtype]:
                continue
            intersect = intervals[svtype][chrom].overlap(pos, stop)
            valid_intersect = [interval for interval in intersect
                               if overlap_size((pos, stop), interval) / (interval[1] - interval[0]) >= cutoff]
            if len(valid_intersect) == 0:
                continue
            distinct_regions = set()
            for interval in valid_intersect:
                base_name, _ = get_base_name_and_index(interval.data)
                distinct_regions.add(base_name)
            carriers = None
            for region in distinct_regions:
                overlapping_region_intervals = [interval for interval in intersect
                                                if get_base_name_and_index(interval.data)[0] == region]
                valid_region_intersect = [interval for interval in valid_intersect
                                          if get_base_name_and_index(interval.data)[0] == region]
                num_region_subdivisions = interval_subdivision_counts[region]
                # Require sufficient fraction of subdivisions in the region to overlap
                if len(valid_region_intersect) / num_region_subdivisions < min_region_overlap:
                    continue
                if carriers is None:
                    carriers = set(s for s in record.samples if _cache_gt_sum(record.samples[s]["GT"]) > 0)
                for sample in record.samples:
                    supported_intersect = list()
                    genotype_medians = list()
                    for interval in valid_region_intersect:
                        interval_name = interval.data
                        if interval_name not in median_vids:
                            raise ValueError(f"Interval {interval_name} not found in median_geno file")
                        if sample in medians and interval_name in medians[sample]:
                            supported_intersect.append(interval)
                            genotype_medians.append(medians[sample][interval_name])
                    n_support = len(supported_intersect)
                    n_valid_overlap = len(valid_region_intersect)
                    frac_valid_overlapping_intervals_supported = n_support / n_valid_overlap
                    frac_region_intervals_supported = n_support / num_region_subdivisions
                    if frac_valid_overlapping_intervals_supported >= min_supported_valid_overlapping_intervals_frac \
                            and frac_region_intervals_supported >= min_supported_region_intervals_frac \
                            and sample not in carriers:
                        # False negative
                        # We will try to add this sample to the variant
                        add_match(false_negative_matches, svtype, name=name, chrom=chrom, pos=pos, end=stop,
                                  sample=sample,
                                  overlapping_region_intervals=overlapping_region_intervals,
                                  valid_region_intervals=valid_region_intersect,
                                  supporting_intervals=supported_intersect,
                                  region=region,
                                  n_region_subdivisions=num_region_subdivisions,
                                  genotype_medians=genotype_medians)
                    elif frac_valid_overlapping_intervals_supported < min_supported_valid_overlapping_intervals_frac \
                            and frac_region_intervals_supported < min_supported_region_intervals_frac \
                            and sample in carriers:
                        # False positive
                        # We will subtract this sample from the variant
                        add_match(false_positive_matches, svtype, name=name, chrom=chrom, pos=pos, end=stop,
                                  sample=sample,
                                  overlapping_region_intervals=overlapping_region_intervals,
                                  valid_region_intervals=valid_region_intersect,
                                  supporting_intervals=supported_intersect,
                                  region=region,
                                  n_region_subdivisions=num_region_subdivisions,
                                  genotype_medians=genotype_medians)
    return false_positive_matches, false_negative_matches


def read_ploidy_table(path):
    """
    Parses tsv of sample ploidy values.

    Parameters
    ----------
    path: Text
        table path

    Returns
    -------
    header: Dict[Text, Dict[Text, int]]
        map of sample to contig to ploidy, i.e. Dict[sample][contig] = ploidy
    """
    ploidy_dict = dict()
    with open(path, 'r') as f:
        header = f.readline().strip().split('\t')
        for line in f:
            tokens = line.strip().split('\t')
            sample = tokens[0]
            ploidy_dict[sample] = {header[i]: int(tokens[i]) for i in range(1, len(header))}
    return ploidy_dict


def is_in_par_region(chrom, pos, stop, par_trees, cutoff):
    length = stop - pos
    if chrom in par_trees:
        par_overlaps = par_trees[chrom].overlap(pos, stop)
        for overlap in par_overlaps:
            if overlap_size(overlap, (pos, stop)) / length >= cutoff:
                return True
    return False


def get_expected_cn(chrom, sample, ploidy_table_dict, is_par=False):
    if is_par:
        return 2
    else:
        if sample not in ploidy_table_dict:
            raise ValueError(f"Sample {sample} not defined in ploidy table")
        if chrom not in ploidy_table_dict[sample]:
            raise ValueError(f"Contig {chrom} not defined in ploidy table")
        return ploidy_table_dict[sample][chrom]


def read_median_geno(list_path, del_ids, dup_ids, del_cutoff, dup_cutoff, ploidy_table_dict, par_trees,
                     min_par_overlap):
    with open(list_path) as f:
        file_paths = [line.strip() for line in f]
    data = defaultdict(dict)
    vids = set()
    for path in file_paths:
        with open(path) as f:
            header = f.readline().strip().split('\t')
            samples_list = header[4:]
            for line in f:
                record = line.strip().split('\t')
                chrom = record[0]
                pos = int(record[1])
                stop = int(record[2])
                vid = record[3]
                vids.add(vid)
                if vid in del_ids:
                    is_del = True
                elif vid in dup_ids:
                    is_del = False
                else:
                    continue
                is_par = is_in_par_region(chrom=chrom, pos=pos, stop=stop, par_trees=par_trees, cutoff=min_par_overlap)
                medians = [float(m) for m in record[4:]]
                for i, median in enumerate(medians):
                    if is_del:
                        cutoff = del_cutoff
                    else:
                        cutoff = dup_cutoff
                    sample = samples_list[i]
                    expected_cn = get_expected_cn(chrom=chrom, sample=sample, ploidy_table_dict=ploidy_table_dict,
                                                  is_par=is_par)
                    if expected_cn == 0:
                        continue
                    cutoff = max(cutoff - 0.5 * (2 - expected_cn), 0)
                    if is_del and median >= cutoff:
                        continue
                    elif (not is_del) and median <= cutoff:
                        continue
                    data[sample][vid] = GenotypeMedian(median=median, threshold=cutoff, hom_ref_cn=expected_cn)
    return data, vids


def get_interval_indices_dict_by_region(intervals):
    result = defaultdict(set)
    for interval in intervals:
        region, index = get_base_name_and_index(interval.data)
        result[region].add(index)
    return result


def get_flanking_intervals_to_remove(sample_overlappers, invalidated_intervals, region):
    region_sample_overlappers = [overlapper for overlapper in sample_overlappers if overlapper.region == region]
    invalidated_interval_ids = set(interval.data for interval in invalidated_intervals)
    sufficiently_overlapping_interval_ids = set(interval.data for ov in region_sample_overlappers
                                                for interval in ov.valid_region_intervals)
    region_overlapping_intervals = [interval for ov in region_sample_overlappers
                                    for interval in ov.overlapping_region_intervals]
    if len(region_overlapping_intervals) < 2:
        return list()
    region_overlapping_intervals = sorted(region_overlapping_intervals, key=attrgetter('begin'))
    flanking_intervals_to_remove = list()
    if region_overlapping_intervals[1].data in invalidated_interval_ids and \
            region_overlapping_intervals[0].data not in sufficiently_overlapping_interval_ids:
        flanking_intervals_to_remove.append(region_overlapping_intervals[0])
    if region_overlapping_intervals[-2].data in invalidated_interval_ids and \
            region_overlapping_intervals[-1].data not in sufficiently_overlapping_interval_ids:
        flanking_intervals_to_remove.append(region_overlapping_intervals[-1])
    return flanking_intervals_to_remove


def revise_genotypes(f_in, f_geno, f_before, f_after, false_positives_dict, revise_false_negative_dict):
    # Find genotypes to reset from existing records
    current_chrom = None
    for record in f_in:
        if record.chrom != current_chrom:
            current_chrom = record.chrom
            logging.info(f"  Subtracting {record.chrom}")
        if record.id in false_positives_dict or record.id in revise_false_negative_dict:
            svtype = record.info.get("SVTYPE", "")
            # Write original record for review
            f_before.write(record)
            if record.id in false_positives_dict:
                gdr_records = false_positives_dict[record.id]
                samples = [r.sample for r in gdr_records]
                # Reset genotypes to hom-ref for invalidated samples
                for s in samples:
                    f_geno.write(f"{record.chrom}\t{record.id}\t{s}\t0\n")
                    # Reset RD genotyping fields but not PESR since we did not re-examine that evidence
                    reset_format_fields(record.samples[s], reset_genotype=True, reset_pesr=False)
            if record.id in revise_false_negative_dict:
                gdr_records = revise_false_negative_dict[record.id]
                # Set genotypes to het for supported samples
                for gdr_match in gdr_records:
                    f_geno.write(f"{record.chrom}\t{record.id}\t{gdr_match.sample}\t1\n")
                    # Set RD genotyping fields but not PESR since we did not re-examine that evidence
                    rescue_format_fields(gt=record.samples[gdr_match.sample],
                                         svtype=svtype,
                                         genotype_medians=gdr_match.genotype_medians,
                                         rescue_genotype=True, reset_pesr=False)
            # Write revised record for review
            f_after.write(record)


def get_revised_intervals(sample_overlappers, regions, pos, stop, dangling_fraction):
    tree = IntervalTree()
    tree.addi(pos, stop)
    # Remove overlapping intervals
    supported_interval_ids = set(interval.data for ov in sample_overlappers for interval in ov.supporting_intervals)
    # Remove intervals that failed genotyping. Also removes flanking GDR subdivisions that partially
    # overlap the variant and are adjacent to a failed interval.
    for region in regions:
        invalidated_intervals = [interval for ov in sample_overlappers
                                 for interval in ov.valid_region_intervals
                                 if interval.data not in supported_interval_ids and ov.region == region]
        flanking_intervals_to_remove = get_flanking_intervals_to_remove(sample_overlappers=sample_overlappers,
                                                                        invalidated_intervals=invalidated_intervals,
                                                                        region=region)
        for ov in invalidated_intervals + flanking_intervals_to_remove:
            tree.chop(ov.begin, ov.end)

    # Clean up small leftover fragments
    for interval in list(tree):
        if (interval.end - interval.begin) / (stop - pos) < dangling_fraction:
            tree.remove(interval)
    return sorted(tree)


def retain_values_if_present(data_dict, key, values_list):
    value = data_dict.get(key, list())
    if isinstance(value, Iterable):
        return [ev for ev in data_dict.get(key, list()) if ev in values_list]
    else:
        return value


def write_revised_variant_record(f_revise_after, interval, index, base_record, vid, samples, original_gt_dict):
    new_record = base_record.copy()
    new_record.id = f"{vid}_revised_{index}"
    new_record.pos = interval[0]
    new_record.stop = interval[1]
    new_record.info["SVLEN"] = new_record.stop - new_record.pos
    # Restore carriers' original GT and RD genotype data
    for sample in samples:
        sample_gt = new_record.samples[sample]
        original_gt_tuple = original_gt_dict[sample]
        sample_gt["GT"] = original_gt_tuple[0]
        sample_gt["RD_CN"] = original_gt_tuple[1]
        sample_gt["RD_GQ"] = original_gt_tuple[2]
    f_revise_after.write(new_record)


def get_ecn(gt):
    ecn = gt.get("ECN", None)
    if ecn is None:
        raise ValueError("Missing ECN format field")
    return ecn


def reset_format_fields(gt, reset_genotype, reset_pesr):
    # Note we do not take PAR into account here to match the rest of the pipeline
    ecn = get_ecn(gt)
    if ecn == 0:
        # Should already be empty
        return
    gt["RD_CN"] = ecn
    gt["RD_GQ"] = RESET_RD_GQ_VALUE
    gt["GQ"] = RESET_GQ_VALUE
    if reset_genotype:
        gt["GT"] = _cache_gt_set_hom_ref(gt["GT"])
    if reset_pesr:
        gt["EV"] = ("RD",)
        for key, val in RESET_PESR_FORMATS_DICT.items():
            gt[key] = val


def rescue_format_fields(gt, svtype, genotype_medians, rescue_genotype, reset_pesr):
    ecn = get_ecn(gt)
    if ecn == 0:
        # Do not rescue on ploidy 0
        return
    # Use median alt allele count for new genotype
    genotype_alt_counts = sorted(genotype_medians, key=get_n_alt_alleles)
    median_genotype = genotype_alt_counts[len(genotype_alt_counts) // 2]
    alt_count = get_n_alt_alleles(median_genotype)
    if svtype == "DEL":
        rd_cn = max(median_genotype.hom_ref_cn - alt_count, 0)
    elif svtype == "DUP":
        rd_cn = min(median_genotype.hom_ref_cn + alt_count, 4)
    else:
        raise ValueError(f"Attempted to rescue unsupported SVTYPE: {svtype}")
    gt["RD_CN"] = rd_cn
    gt["RD_GQ"] = RESCUE_RD_GQ_VALUE
    gt["GQ"] = RESCUE_GQ_VALUE
    if rescue_genotype:
        gt["GT"] = _cache_gt_set_het(gt["GT"])
    if reset_pesr:
        gt["EV"] = ("RD",)
        for key, val in RESET_PESR_FORMATS_DICT.items():
            gt[key] = val


def reset_record(record):
    new_record = record.copy()
    new_record.info["EVIDENCE"] = retain_values_if_present(new_record.info, "EVIDENCE", ["RD", "BAF"])
    # Reset all other samples to hom-ref
    for s, gt in new_record.samples.items():
        reset_format_fields(gt, reset_genotype=True, reset_pesr=True)
    return new_record


def get_interval_key(interval):
    return interval.begin, interval.end


def revise_partially_supported_variants(f_before, f_revise_after, false_positives_dict, dangling_fraction):
    # Create partial events from existing variants
    current_chrom = None
    for vcf_record in f_before:
        if vcf_record.chrom != current_chrom:
            current_chrom = vcf_record.chrom
            logging.info(f"  Revising {vcf_record.chrom}")
        vid = vcf_record.id
        pos = vcf_record.pos
        stop = vcf_record.stop
        overlappers = false_positives_dict[vid]
        regions = list(set(overlapper.region for overlapper in overlappers))
        if len(overlappers) == 0:
            raise ValueError(f"Should have found variant {vid} in overlap record set (bug)")
        # Store original GT and RD genotype data, since we'll be modifying the vcf record in place
        overlapper_samples = set(overlapper.sample for overlapper in overlappers)
        original_gt_dict = {sample: (vcf_record.samples[sample]["GT"], vcf_record.samples[sample]["RD_CN"],
                                     vcf_record.samples[sample]["RD_GQ"]) for sample in overlapper_samples}
        # Copy of record with all genotypes reset to hom ref
        wiped_record = reset_record(vcf_record)
        # Dictionary keyed on interval mapping to a list of carrier samples
        intervals_dict = defaultdict(list)
        # Iterate over samples, subtracting invalidated intervals and creating new records from the result
        for sample, sample_overlappers in groupby(overlappers, key=attrgetter("sample")):
            # copy iterable to permanent list
            sample_overlappers = list(sample_overlappers)
            new_intervals = get_revised_intervals(sample_overlappers=sample_overlappers, regions=regions,
                                                  pos=pos, stop=stop, dangling_fraction=dangling_fraction)
            for interval in new_intervals:
                intervals_dict[get_interval_key(interval)].append(sample)
        sorted_intervals = sorted(intervals_dict.keys())
        for i, interval in enumerate(sorted_intervals):
            # Write out new revised record. Note that these are not globally sorted.
            write_revised_variant_record(f_revise_after=f_revise_after, interval=interval, index=i,
                                         base_record=wiped_record, vid=vid, samples=intervals_dict[interval],
                                         original_gt_dict=original_gt_dict)


def create_new_variants(f_revise_after, new_records_dict, region_intervals_dict, ploidy_table_dict):
    # Create new rescued calls not derived from existing variants
    # note new_records_dict structure: svtype -> sample -> region -> record corresponding to new variants to create
    #      region_intervals_dict structure: svtype -> region -> interval
    record_index = 0
    for svtype, svtype_dict in new_records_dict.items():
        for sample, sample_dict in svtype_dict.items():
            for region, match in sample_dict.items():
                if svtype not in region_intervals_dict:
                    raise ValueError(f"Attempted to rescue variant with SVTYPE {svtype} "
                                     f"but no regions of that type exist")
                if region not in region_intervals_dict[svtype]:
                    raise ValueError(f"Region {region} interval undefined (bug)")
                region_contig, region_start, region_end = region_intervals_dict[svtype][region]
                # TODO: get reference allele
                alleles = ("N", f"<{svtype}>")
                vid = f"{match.name}_rescue_{record_index}"
                record = f_revise_after.new_record(contig=region_contig, start=region_start, stop=region_end, alleles=alleles,
                                         id=vid, filter=None)
                record.info["SVTYPE"] = svtype
                record.info["SVLEN"] = region_end - region_start
                record.info["EVIDENCE"] = ("RD",)
                record.info["ALGORITHMS"] = ("depth",)
                for s, gt in record.samples.items():
                    # Set GT and ECN so we can use the rescue/reset format field methods
                    # TODO : we use diploid genotypes to follow gatk-sv format, but ideally we'd
                    #  detect it from an existing record
                    gt["GT"] = (None, None)
                    gt["ECN"] = get_expected_cn(chrom=region_contig, sample=sample, ploidy_table_dict=ploidy_table_dict)
                    if s == sample:
                        rescue_format_fields(gt=gt, svtype=svtype, genotype_medians=match.genotype_medians,
                                             rescue_genotype=True, reset_pesr=True)
                    else:
                        reset_format_fields(gt, reset_genotype=True, reset_pesr=True)
                f_revise_after.write(record)
                record_index += 1


def remap_overlapper_dict_by_vid(overlappers_dict):
    # make dictionary: vid -> list of GDR overlap records, for fast lookup as we traverse the vcf
    result = defaultdict(list)
    for _, svtype_dict in overlappers_dict.items():
        for _, sample_dict in svtype_dict.items():
            for record_name, record_list in sample_dict.items():
                result[record_name].extend(record_list)
    return result


def remap_overlapper_dict_by_sample_and_region(overlappers_dict):
    # make dictionary: svtype -> sample -> region -> list of GDR overlap records, for fast lookup
    result = dict()
    for svtype, svtype_dict in overlappers_dict.items():
        for sample, sample_dict in svtype_dict.items():
            for _, record_list in sample_dict.items():
                for record in record_list:
                    if svtype not in result:
                        result[svtype] = dict()
                    if sample not in result[svtype]:
                        result[svtype][sample] = defaultdict(list)
                    result[svtype][sample][record.region].append(record)
    return result


def fraction_supporting_intervals_in_variant(record):
    return sum([interval.end - interval.begin for interval in record.supporting_intervals]) / (record.end - record.pos)


def get_maximal_record(record_list, key):
    max_value = 0
    max_records = list()
    for record in record_list:
        value = key(record)
        if value > max_value:
            max_records.clear()
            max_records.append(record)
            max_value = value
        elif value == max_value:
            max_records.append(record)
    # Should be very rare, but use variant ID as tiebreaker
    max_record = sorted(max_records, key=attrgetter("name"))[0]
    return max_value, max_record


def adjudicate_false_negatives(false_negative_dict, min_false_negative_rescue_overlap):
    # make dictionary: vid -> record list for false negatives to rescue
    revise_false_negative_dict = defaultdict(list)
    # make dictionary: svtype -> sample -> region -> record corresponding to new variants to create
    new_records_dict = dict()
    for svtype, svtype_dict in false_negative_dict.items():
        for sample, sample_dict in svtype_dict.items():
            for region, record_list in sample_dict.items():
                if len(record_list) == 0:
                    # Should not happen but catch it in case
                    continue
                max_overlap, max_record = get_maximal_record(record_list=record_list,
                                                             key=fraction_supporting_intervals_in_variant)
                if max_overlap >= min_false_negative_rescue_overlap:
                    # Enforce sufficient overlap of supporting intervals on the variant
                    # We will rescue the call by adding it to the existing record
                    revise_false_negative_dict[max_record.name].append(max_record)
                else:
                    # We will rescue the call by creating a new record
                    # Note there is an edge case where overlapping GDRs could result in duplicated calls, but this
                    # should be handled by downstream reclustering/deduplication
                    if svtype not in new_records_dict:
                        new_records_dict[svtype] = dict()
                    if sample not in new_records_dict[svtype]:
                        new_records_dict[svtype][sample] = dict()
                    new_records_dict[svtype][sample][region] = max_record
    return revise_false_negative_dict, new_records_dict


def revise_genotypes_and_create_records(input_vcf_path, revised_genotypes_tsv_path, revised_records_before_update_path,
                                        unsorted_revised_records_after_update_path, new_revised_records_vcf_path,
                                        region_intervals_dict, false_positive_matches, false_negative_matches,
                                        dangling_fraction, min_false_negative_rescue_overlap, ploidy_table_dict):
    false_positives_dict = remap_overlapper_dict_by_vid(false_positive_matches)
    false_negative_dict = remap_overlapper_dict_by_sample_and_region(false_negative_matches)
    revise_false_negative_dict, new_records_dict = \
        adjudicate_false_negatives(false_negative_dict=false_negative_dict,
                                   min_false_negative_rescue_overlap=min_false_negative_rescue_overlap)
    # Find invalidated records and reset their genotypes
    with pysam.VariantFile(input_vcf_path) as f_in, \
            gzip.open(revised_genotypes_tsv_path, mode="wt") as f_geno, \
            pysam.VariantFile(revised_records_before_update_path, mode="w", header=f_in.header) as f_before, \
            pysam.VariantFile(unsorted_revised_records_after_update_path, mode="w", header=f_in.header) as f_after:
        revise_genotypes(f_in=f_in, f_geno=f_geno, f_before=f_before, f_after=f_after,
                         false_positives_dict=false_positives_dict,
                         revise_false_negative_dict=revise_false_negative_dict)
    # Revise invalidated records
    with pysam.VariantFile(revised_records_before_update_path) as f_before, \
            pysam.VariantFile(new_revised_records_vcf_path, mode="w", header=f_before.header) as f_revise_after:
        create_new_variants(f_revise_after=f_revise_after, new_records_dict=new_records_dict,
                            region_intervals_dict=region_intervals_dict, ploidy_table_dict=ploidy_table_dict)
        revise_partially_supported_variants(f_before=f_before,
                                            f_revise_after=f_revise_after,
                                            false_positives_dict=false_positives_dict,
                                            dangling_fraction=dangling_fraction)


def sort_vcf(vcf_path, out_path, temp_dir):
    proc = subprocess.Popen(
        ['bcftools', 'sort', '-T', temp_dir, vcf_path, '-O', 'z', '-o', out_path]
    )
    logging.info(f"Executing bcftools subprocess: {' '.join(proc.args)}")
    proc_out, proc_err = proc.communicate()
    return_code = proc.returncode
    if return_code != 0:
        sys.stderr.write(proc_err)
        raise Exception('bcftools sort returned non-zero exit code: {}'.format(return_code))
    logging.info(f"Subprocess completed successfully")


def _parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Cleans up biallelic DEL/DUP variants in genomic disorder regions",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--vcf', type=str, required=True, help='Input vcf, GATK-formatted with ECN format fields')
    parser.add_argument('--region-bed', type=str, required=True, help='Preprocessed genomic disorder regions bed')
    parser.add_argument('--par-bed', type=str, required=True, help='PAR region bed')
    parser.add_argument('--median-geno-list', type=str, required=True,
                        help='List of paths to RDTest ".median_geno" files')
    parser.add_argument('--ploidy-table', type=str, required=True, help='Sample ploidy table')
    parser.add_argument('--out', type=str, required=True, help='Output base path')
    parser.add_argument('--overlap', type=float, default=0.9, help='Subdivision overlap fraction cutoff')
    parser.add_argument('--del-median', type=float, default=0.6, help='DEL RDTest median cutoff')
    parser.add_argument('--dup-median', type=float, default=1.4, help='DUP RDTest median cutoff')
    parser.add_argument('--vcf-min-size', type=int, default=1000, help='Min size of vcf variants')
    parser.add_argument('--min-region-overlap', type=float, default=0.3,
                        help='Min fraction of sufficiently overlapping genotyped intervals in the GD region')
    parser.add_argument('--min-frac-supporting-genotypes', type=float, default=0.7,
                        help='Min fraction of sufficiently overlapping genotyped intervals that must support a call')
    parser.add_argument('--min-supported-region-intervals-frac', type=float, default=0.5,
                        help='Min fraction of region intervals that must be supported to accept a call')
    parser.add_argument('--min-dangling-frac', type=float, default=0.3,
                        help='Min interval fraction of variant size for removing leftover fragments')
    parser.add_argument('--min-par-overlap', type=float, default=0.5,
                        help='Min overlap fraction for an interval to be treated as PAR')
    parser.add_argument('--min-false-negative-rescue-overlap', type=float, default=0.8,
                        help='Min overlap fraction of a variant to rescue false negatives in an existing call. '
                             'If not met, a new variant will be created')
    parser.add_argument('--temp', type=str, default="./", help='Temporary directory path')
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
    ploidy_table_dict = read_ploidy_table(args.ploidy_table)
    logging.info("Reading genomic disorder regions bed file...")
    gdr_trees, gdr_del_ids, gdr_dup_ids = create_trees_from_bed_records_by_type(args.region_bed)
    region_intervals_dict = create_region_intervals_dict(gdr_trees)
    logging.info("Reading PAR bed file...")
    par_trees = create_trees_from_bed_records(args.par_bed)
    logging.info("Reading \"median_geno\" file...")
    medians, median_vids = read_median_geno(list_path=args.median_geno_list, del_ids=gdr_del_ids, dup_ids=gdr_dup_ids,
                                            del_cutoff=args.del_median, dup_cutoff=args.dup_median,
                                            ploidy_table_dict=ploidy_table_dict, par_trees=par_trees,
                                            min_par_overlap=args.min_par_overlap)
    logging.info("Parsing VCF for variants overlapping genomic disorder regions...")
    false_positive_matches, false_negative_matches = \
        get_overlapping_samples_vcf(vcf_path=args.vcf, intervals=gdr_trees,
                                    medians=medians,
                                    median_vids=median_vids,
                                    cutoff=args.overlap,
                                    vcf_min_size=args.vcf_min_size,
                                    min_supported_valid_overlapping_intervals_frac=args.min_frac_supporting_genotypes,
                                    min_region_overlap=args.min_region_overlap,
                                    min_supported_region_intervals_frac=args.min_supported_region_intervals_frac)
    revised_genotypes_tsv_path = f"{args.out}.revised_genotypes.tsv.gz"
    revised_records_before_update_path = f"{args.out}.revised_before_update.vcf.gz"
    unsorted_revised_records_after_update_path = f"{args.out}.revised_after_update.unsorted.vcf.gz"
    sorted_revised_records_after_update_path = f"{args.out}.revised_after_update.vcf.gz"
    logging.info("Revising genotypes and variants...")
    with tempfile.NamedTemporaryFile(dir=args.temp, suffix=".vcf.gz") as temp_vcf:
        revise_genotypes_and_create_records(input_vcf_path=args.vcf,
                                            revised_genotypes_tsv_path=revised_genotypes_tsv_path,
                                            revised_records_before_update_path=revised_records_before_update_path,
                                            unsorted_revised_records_after_update_path=unsorted_revised_records_after_update_path,
                                            new_revised_records_vcf_path=temp_vcf.name,
                                            region_intervals_dict=region_intervals_dict,
                                            false_positive_matches=false_positive_matches,
                                            false_negative_matches=false_negative_matches,
                                            dangling_fraction=args.min_dangling_frac,
                                            min_false_negative_rescue_overlap=args.min_false_negative_rescue_overlap,
                                            ploidy_table_dict=ploidy_table_dict)
        pysam.tabix_index(revised_records_before_update_path, preset="vcf", force=True)
        pysam.tabix_index(unsorted_revised_records_after_update_path, preset="vcf", force=True)
        with tempfile.TemporaryDirectory(dir=args.temp) as temp_dir:
            sort_vcf(vcf_path=temp_vcf.name, out_path=sorted_revised_records_after_update_path, temp_dir=temp_dir)
            pysam.tabix_index(sorted_revised_records_after_update_path, preset="vcf", force=True)


if __name__ == "__main__":
    main()
