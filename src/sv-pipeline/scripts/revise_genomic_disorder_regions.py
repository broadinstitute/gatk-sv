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

from intervaltree import Interval, IntervalTree
import pysam

# Delimiter suffix appended to the end of interval IDs before the index, e.g. "intervalA__0", "intervalA__1", ...
INDEX_DELIMITER = "__"
LEFT_INDEX_PREFIX = "L"
MIDDLE_INDEX_PREFIX = "M"
RIGHT_INDEX_PREFIX = "R"

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

CODE_FALSE_NEGATIVE_IN_EXISTING_VARIANT = "FALSE_NEGATIVE_IN_EXISTING_VARIANT"
CODE_FALSE_POSITIVE_IN_EXISTING_VARIANT = "FALSE_POSITIVE_IN_EXISTING_VARIANT"
CODE_REVISED_BREAKPOINTS_OF_EXISTING_VARIANT = "REVISED_BREAKPOINTS_OF_EXISTING_VARIANT"
CODE_NEW_VARIANT_FOR_FALSE_NEGATIVE_IN_REGION = "NEW_VARIANT_FOR_FALSE_NEGATIVE_IN_REGION"

_gt_sum_map = dict()
_gt_set_hom_ref_map = dict()
_gt_set_het_map = dict()
_gt_set_hom_var_map = dict()


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


def _cache_gt_set_hom_var(gt):
    s = _gt_set_hom_var_map.get(gt, None)
    if s is None:
        s = tuple(1 for _ in gt)
        _gt_set_hom_var_map[gt] = s
    return s


# Creates dictionary of trees[sv_type][contig] from bed file
def create_trees_from_bed_records_by_type(bed_path, min_region_size):
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
            if stop - start < min_region_size:
                logging.warning(f"Region {name} smaller than minimum region size ({min_region_size})")
                continue
            if svtype not in trees:
                raise ValueError("Unexpected SVTYPE in bed file: %s" % type)
            if svtype == "DEL":
                del_ids.add(name)
            elif svtype == "DUP":
                dup_ids.add(name)
            trees[svtype][contig].addi(start, stop, name)
    return trees, del_ids, dup_ids


# Region intervals dictionary, only over the middle section
def create_region_intervals_dict(trees):
    regions_dict = dict()
    # Get interval sets for each region
    for svtype, svtype_dict in trees.items():
        if svtype not in regions_dict:
            regions_dict[svtype] = dict()
        for contig, contig_tree in svtype_dict.items():
            for interval in contig_tree:
                region_subinterval = interval.data
                region, index_prefix, _ = get_base_name_and_index(region_subinterval)
                if index_prefix == MIDDLE_INDEX_PREFIX:
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


def reciprocal_overlap(interval_a, interval_b):
    return overlap_size(interval_a, interval_b) / max(interval_a[1] - interval_a[0], interval_b[1] - interval_b[0])


class GDRMatch:
    def __init__(self, name, chrom, pos, end, sample, svtype, overlapping_region_intervals, valid_region_intervals,
                 left_supporting_intervals, middle_supporting_intervals, right_supporting_intervals, region, n_region_subdivisions, genotype_medians):
        self.name = name
        self.chrom = chrom
        self.pos = pos
        self.end = end
        self.sample = sample
        self.svtype = svtype
        self.overlapping_region_intervals = overlapping_region_intervals
        self.valid_region_intervals = valid_region_intervals
        self.left_supporting_intervals = left_supporting_intervals,
        self.middle_supporting_intervals = middle_supporting_intervals,
        self.right_supporting_intervals = right_supporting_intervals,
        self.region = region
        self.n_region_subdivisions = n_region_subdivisions
        self.genotype_medians = genotype_medians

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return "(" + ", ".join(str(x) for x in [self.name, self.chrom, self.pos, self.end, self.sample, self.svtype,
                                                self.region,
                                                len(self.middle_supporting_intervals), len(self.valid_region_intervals),
                                                len(self.overlapping_region_intervals),
                                                self.n_region_subdivisions]) + ")"


class GenotypeMedian:
    def __init__(self, median, threshold, hom_ref_cn, pos, stop):
        self.median = median
        self.threshold = threshold
        self.hom_ref_cn = hom_ref_cn
        self.pos = pos
        self.stop = stop


def get_base_name_and_index(interval_id):
    tokens = interval_id.split(INDEX_DELIMITER)
    if len(tokens) != 2:
        raise ValueError(f"Expected to find \"{INDEX_DELIMITER}\" in interval name: {interval_id}")
    region_name = tokens[0]
    index_prefix = None
    for prefix in [LEFT_INDEX_PREFIX, MIDDLE_INDEX_PREFIX, RIGHT_INDEX_PREFIX]:
        if tokens[1].startswith(prefix):
            index_prefix = prefix
            break
    if index_prefix is None:
        raise ValueError(f"Could find valid index prefix for interval {interval_id}. Valid prefixes are: "
                         f"{set([LEFT_INDEX_PREFIX, MIDDLE_INDEX_PREFIX, RIGHT_INDEX_PREFIX])}")
    index = int(tokens[1].replace(index_prefix, ""))
    return region_name, index_prefix, index


def add_match(matches_dict, svtype, name, chrom, pos, end, sample, overlapping_region_intervals, valid_region_intervals,
              left_supporting_intervals, middle_supporting_intervals, right_supporting_intervals,
              region, n_region_subdivisions, genotype_medians):
    if svtype not in matches_dict:
        matches_dict[svtype] = dict()
    if sample not in matches_dict[svtype]:
        matches_dict[svtype][sample] = dict()
    if name not in matches_dict[svtype][sample]:
        matches_dict[svtype][sample][name] = list()
    matches_dict[svtype][sample][name].append(
        GDRMatch(name=name, chrom=chrom, pos=pos, end=end, sample=sample, svtype=svtype,
                 overlapping_region_intervals=overlapping_region_intervals,
                 valid_region_intervals=valid_region_intervals,
                 left_supporting_intervals=left_supporting_intervals,
                 middle_supporting_intervals=middle_supporting_intervals,
                 right_supporting_intervals=right_supporting_intervals,
                 region=region,
                 n_region_subdivisions=n_region_subdivisions,
                 genotype_medians=genotype_medians)
    )


def get_n_alt_alleles(genotype_median):
    # Returns either 1 or 2
    return min(2, max(1, math.ceil(math.fabs(genotype_median.median - genotype_median.threshold) / 0.5)))


def get_overlapping_samples_vcf(vcf_path, gdr_trees, region_intervals_dict,
                                medians, median_vids, cutoff, vcf_min_size,
                                min_supported_valid_overlapping_intervals_frac,
                                min_region_overlap):
    middle_interval_subdivision_counts = dict()
    unique_subdivision_names = set()
    for svtype in gdr_trees:
        for chrom in gdr_trees[svtype]:
            for interval in gdr_trees[svtype][chrom]:
                if interval.data in unique_subdivision_names:
                    raise ValueError(f"Encountered duplicate interval id: {interval.data}")
                unique_subdivision_names.add(interval.data)
                name, index_prefix, _ = get_base_name_and_index(interval.data)
                if index_prefix == MIDDLE_INDEX_PREFIX:
                    if name not in middle_interval_subdivision_counts:
                        middle_interval_subdivision_counts[name] = 0
                    middle_interval_subdivision_counts[name] += 1

    false_positive_matches = dict()
    false_negative_matches = dict()
    true_positives = dict()
    vcf_vids = set()
    with pysam.VariantFile(vcf_path) as vcf:
        samples_list = list(vcf.header.samples)
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
            if svtype not in gdr_trees.keys():
                continue
            if svlen < vcf_min_size:
                continue
            if record.chrom not in gdr_trees[svtype]:
                continue
            get_revisions_over_regions(
                regions=None,
                false_negative_matches=false_negative_matches, false_positive_matches=false_positive_matches,
                true_positives=true_positives, name=record.id, pos=pos, stop=stop, svtype=svtype, chrom=chrom,
                samples=record.samples, gdr_trees=gdr_trees, medians=medians, median_vids=median_vids, cutoff=cutoff,
                middle_interval_subdivision_counts=middle_interval_subdivision_counts,
                min_supported_valid_overlapping_intervals_frac=min_supported_valid_overlapping_intervals_frac,
                min_region_overlap=min_region_overlap
            )

    logging.info(f"  Finding raw variants in full regions")
    # Find false negatives over full GDRs
    raw_region_false_negatives = dict()
    # Set all samples to hom-ref and reuse get_revisions_over_regions() method to find false negatives
    gdr_samples = {s: {"GT": (0, 0)} for s in samples_list}
    for svtype in region_intervals_dict:
        for region_name in region_intervals_dict[svtype]:
            region_chrom, region_start, region_end = region_intervals_dict[svtype][region_name]
            name = f"{region_name}_new"
            get_revisions_over_regions(
                regions=[region_name],
                false_negative_matches=raw_region_false_negatives, false_positive_matches=dict(), true_positives=dict(),
                name=name, pos=region_start, stop=region_end, svtype=svtype, chrom=region_chrom, samples=gdr_samples,
                gdr_trees=gdr_trees, medians=medians, median_vids=median_vids, cutoff=cutoff,
                middle_interval_subdivision_counts=middle_interval_subdivision_counts,
                min_supported_valid_overlapping_intervals_frac=min_supported_valid_overlapping_intervals_frac,
                min_region_overlap=min_region_overlap
            )
    return false_positive_matches, false_negative_matches, raw_region_false_negatives, true_positives


def get_revisions_over_regions(regions, false_negative_matches, false_positive_matches, true_positives,
                               name, pos, stop, svtype, chrom, samples, gdr_trees, medians, median_vids, cutoff,
                               middle_interval_subdivision_counts, min_supported_valid_overlapping_intervals_frac,
                               min_region_overlap):
    intersect = gdr_trees[svtype][chrom].overlap(pos, stop)
    valid_intersect = [interval for interval in intersect
                       if overlap_size((pos, stop), interval) / (interval[1] - interval[0]) >= cutoff]
    if len(valid_intersect) == 0:
        return
    if regions is None:
        regions = set()
        for interval in valid_intersect:
            base_name, _, _ = get_base_name_and_index(interval.data)
            regions.add(base_name)
    carriers = None
    for region in regions:
        overlapping_region_intervals = [interval for interval in intersect
                                        if get_base_name_and_index(interval.data)[0] == region]
        padded_region_intersect = [interval for interval in valid_intersect
                                   if get_base_name_and_index(interval.data)[0] == region]
        left_padding_names = list()
        middle_names = list()
        right_padding_names = list()
        for interval in padded_region_intersect:
            index_prefix = get_base_name_and_index(interval.data)[1]
            if index_prefix == LEFT_INDEX_PREFIX:
                left_padding_names.append(interval.data)
            elif index_prefix == MIDDLE_INDEX_PREFIX:
                middle_names.append(interval.data)
            elif index_prefix == RIGHT_INDEX_PREFIX:
                right_padding_names.append(interval.data)
            else:
                raise ValueError(f"Unexpected index prefix: {index_prefix}")
        num_region_subdivisions = middle_interval_subdivision_counts[region]
        # Require sufficient fraction of subdivisions in the region to overlap
        middle_intersect = [interval for interval in valid_intersect if interval.data in middle_names]
        if len(middle_intersect) / num_region_subdivisions < min_region_overlap:
            continue
        # Lazily load carriers
        if carriers is None:
            carriers = set(s for s, gt in samples.items() if _cache_gt_sum(gt["GT"]) > 0)
        for sample in samples:
            left_supported_intersect = list()
            right_supported_intersect = list()
            middle_supported_intersect = list()
            genotype_medians = list()
            for interval in padded_region_intersect:
                interval_name = interval.data
                # if interval_name not in median_vids:
                # TODO
                # raise ValueError(f"Interval {interval_name} not found in median_geno file")
                # logging.warning(f"Interval {interval_name} not found in median_geno file")
                if sample in medians and interval_name in medians[sample]:
                    if interval_name in left_padding_names:
                        left_supported_intersect.append(interval)
                    elif interval_name in middle_names:
                        middle_supported_intersect.append(interval)
                    elif interval_name in right_padding_names:
                        right_supported_intersect.append(interval)
                    else:
                        raise ValueError(f"Left/middle/right status of region not found: {interval_name} (bug)")
                    genotype_medians.append(medians[sample][interval_name])
            n_support = len(middle_supported_intersect)
            n_valid_overlap = len(middle_intersect)
            # Fraction of supported intervals in the actual region (middle section)
            frac_valid_overlapping_intervals_supported = n_support / n_valid_overlap
            if frac_valid_overlapping_intervals_supported >= min_supported_valid_overlapping_intervals_frac:
                if sample not in carriers:
                    # False negative
                    # We will try to add this sample to the variant
                    add_match(false_negative_matches, svtype, name=name, chrom=chrom, pos=pos, end=stop,
                              sample=sample,
                              overlapping_region_intervals=overlapping_region_intervals,
                              valid_region_intervals=padded_region_intersect,
                              left_supporting_intervals=left_supported_intersect,
                              middle_supporting_intervals=middle_supported_intersect,
                              right_supporting_intervals=right_supported_intersect,
                              region=region,
                              n_region_subdivisions=num_region_subdivisions,
                              genotype_medians=genotype_medians)
                else:
                    # True positive
                    if svtype not in true_positives:
                        true_positives[svtype] = dict()
                    if sample not in true_positives[svtype]:
                        true_positives[svtype][sample] = dict()
                    if region not in true_positives[svtype][sample]:
                        true_positives[svtype][sample][region] = list()
                    true_positives[svtype][sample][region].append((pos, stop))
            elif sample in carriers:
                # False positive
                # We will subtract this sample from the variant
                add_match(false_positive_matches, svtype, name=name, chrom=chrom, pos=pos, end=stop,
                          sample=sample,
                          overlapping_region_intervals=overlapping_region_intervals,
                          valid_region_intervals=padded_region_intersect,
                          left_supporting_intervals=left_supported_intersect,
                          middle_supporting_intervals=middle_supported_intersect,
                          right_supporting_intervals=right_supported_intersect,
                          region=region,
                          n_region_subdivisions=num_region_subdivisions,
                          genotype_medians=genotype_medians)


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
                    data[sample][vid] = GenotypeMedian(median=median, threshold=cutoff, hom_ref_cn=expected_cn,
                                                       pos=pos, stop=stop)
    return data, vids


def get_fragment_intervals_to_remove(sample_overlappers, invalidated_intervals, region):
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


def revise_genotypes(f_in, f_geno, f_before, f_after, f_manifest, batch,
                     false_positives_dict, revise_false_negative_dict):
    # Find genotypes to reset from existing records
    current_chrom = None
    for record in f_in:
        if record.chrom != current_chrom:
            current_chrom = record.chrom
            logging.info(f"  Revising genotypes on {record.chrom}")
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
                for gdr_match in gdr_records:
                    f_manifest.write(manifest_record(new_id=record.id, old_id=record.id, svtype=svtype,
                                                     region=gdr_match.region, sample=gdr_match.sample,
                                                     batch=batch,
                                                     code=CODE_FALSE_POSITIVE_IN_EXISTING_VARIANT))
            if record.id in revise_false_negative_dict:
                gdr_records = revise_false_negative_dict[record.id]
                # Set genotypes to het/hom-var for supported samples
                for gdr_match in gdr_records:
                    f_geno.write(f"{record.chrom}\t{record.id}\t{gdr_match.sample}\t1\n")
                    # Set RD genotyping fields but not PESR since we did not re-examine that evidence
                    rescue_format_fields(gt=record.samples[gdr_match.sample],
                                         svtype=svtype,
                                         genotype_medians=gdr_match.genotype_medians,
                                         rescue_genotype=True, reset_pesr=False)
                    f_manifest.write(manifest_record(new_id=record.id, old_id=record.id, svtype=svtype,
                                                     region=gdr_match.region, sample=gdr_match.sample,
                                                     batch=batch,
                                                     code=CODE_FALSE_NEGATIVE_IN_EXISTING_VARIANT))
            # Write revised record for review
            f_after.write(record)


def get_revised_intervals(sample_overlappers, regions, pos, stop, dangling_fraction):
    tree = IntervalTree()
    tree.addi(pos, stop)
    # Get supported interval names
    supported_interval_ids = set(
        interval.data for ov in sample_overlappers
        for interval in ov.left_supporting_intervals + ov.middle_supporting_intervals + ov.right_supporting_intervals
    )
    # Remove intervals that failed genotyping. Also removes leftover fragments.
    for region in regions:
        invalidated_intervals = [interval for ov in sample_overlappers
                                 for interval in ov.valid_region_intervals
                                 if interval.data not in supported_interval_ids and ov.region == region]
        fragment_intervals_to_remove = get_fragment_intervals_to_remove(sample_overlappers=sample_overlappers,
                                                                        invalidated_intervals=invalidated_intervals,
                                                                        region=region)
        for ov in invalidated_intervals + fragment_intervals_to_remove:
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


def manifest_record(new_id, old_id, svtype, region, sample, batch, code):
    return f"{new_id}\t{old_id}\t{svtype}\t{region}\t{sample}\t{batch}\t{code}\n"


def write_revised_variant_record(f_revise_after, f_manifest, interval, index, base_record, vid, samples,
                                 original_gt_dict, regions, svtype, batch, code):
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
    for sample in samples:
        for region in regions:
            f_manifest.write(manifest_record(new_id=new_record.id, old_id=vid, svtype=svtype, region=region,
                                             sample=sample, batch=batch, code=code))


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
        if alt_count == 1:
            gt["GT"] = _cache_gt_set_het(gt["GT"])
        elif alt_count == 2:
            gt["GT"] = _cache_gt_set_hom_var(gt["GT"])
        else:
            raise ValueError(f"Cannot rescue genotype with alt count {alt_count}")
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


def revise_partially_supported_variants(f_before, f_new, f_manifest, batch,
                                        false_positives_dict, dangling_fraction):
    # Create partial events from existing variants
    # (svtype, sample, chrom) -> list(intervals)
    new_partial_events_tree_dict = dict()
    current_chrom = None
    for vcf_record in f_before:
        vid = vcf_record.id
        if vid not in false_positives_dict:
            continue
        if vcf_record.chrom != current_chrom:
            current_chrom = vcf_record.chrom
            logging.info(f"  Revising partial variants on {vcf_record.chrom}")
        pos = vcf_record.pos
        stop = vcf_record.stop
        svtype = vcf_record.info.get("SVTYPE", "")
        overlappers = false_positives_dict[vid]
        regions = list(set(overlapper.region for overlapper in overlappers))
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
        if len(sorted_intervals) > 0 and (svtype, sample, vcf_record.chrom) not in new_partial_events_tree_dict:
            new_partial_events_tree_dict[(svtype, sample, vcf_record.chrom)] = IntervalTree()
        for i, interval in enumerate(sorted_intervals):
            new_partial_events_tree_dict[(svtype, sample, vcf_record.chrom)].addi(interval[0], interval[1])
            # Write out new revised record. Note that these are not globally sorted.
            write_revised_variant_record(f_revise_after=f_new, f_manifest=f_manifest, interval=interval,
                                         index=i, base_record=wiped_record, vid=vid, samples=intervals_dict[interval],
                                         original_gt_dict=original_gt_dict, regions=regions,
                                         svtype=svtype, batch=batch, code=CODE_REVISED_BREAKPOINTS_OF_EXISTING_VARIANT)
    return new_partial_events_tree_dict


def create_new_variants(f_new, f_manifest, new_records_dict, region_intervals_dict, batch, ploidy_table_dict):
    # Create new rescued calls not derived from existing variants
    # note new_records_dict structure: svtype -> sample -> region -> record corresponding to new variants to create
    #      region_intervals_dict structure: svtype -> region -> interval
    logging.info(f"  Creating new variants for novel false negatives")
    record_index = 0
    for svtype, svtype_dict in new_records_dict.items():
        for sample, sample_records_list in svtype_dict.items():
            for match in sample_records_list:
                if svtype not in region_intervals_dict:
                    raise ValueError(f"Attempted to rescue variant with SVTYPE {svtype} "
                                     f"but no regions of that type exist")
                if match.region not in region_intervals_dict[svtype]:
                    raise ValueError(f"Region {match.region} interval undefined (bug)")
                region_contig, region_start, region_end = region_intervals_dict[svtype][match.region]
                # TODO: get reference allele
                alleles = ("N", f"<{svtype}>")
                vid = f"{match.name}_rescue_{record_index}"
                record = f_new.new_record(contig=region_contig, start=match.pos, stop=match.end,
                                                   alleles=alleles, id=vid, filter=None)
                record.info["SVTYPE"] = svtype
                record.info["SVLEN"] = record.stop - record.pos
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
                f_new.write(record)
                f_manifest.write(manifest_record(new_id=vid, old_id=".", svtype=svtype, region=match.region,
                                                 sample=sample, batch=batch,
                                                 code=CODE_NEW_VARIANT_FOR_FALSE_NEGATIVE_IN_REGION))
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


def fraction_supporting_middle_intervals_in_variant(record):
    return sum([interval.end - interval.begin
                for interval in record.middle_supporting_intervals]) / (record.end - record.pos)


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
    for svtype, svtype_dict in false_negative_dict.items():
        for sample, sample_dict in svtype_dict.items():
            for _, record_list in sample_dict.items():
                if len(record_list) == 0:
                    # Should not happen but catch it in case
                    continue
                # Look over candidate variants to use in the FN rescue and pick one with most support
                max_overlap, max_record = get_maximal_record(record_list=record_list,
                                                             key=fraction_supporting_middle_intervals_in_variant)
                if max_overlap >= min_false_negative_rescue_overlap:
                    # Enforce sufficient overlap of supporting intervals on the variant
                    # We will rescue the call by adding it to the existing record
                    revise_false_negative_dict[max_record.name].append(max_record)
    return revise_false_negative_dict


def defragment_region(intervals, pos, stop, dangling_fraction):
    tree = IntervalTree()
    for interval in intervals:
        tree.addi(interval.begin, interval.end)
    tree.merge_overlaps(strict=False)
    # Remove fragments
    for interval in list(tree):
        if (interval.end - interval.begin) / (stop - pos) < dangling_fraction:
            tree.remove(interval)
    # Edge case: if defragmentation deletes the event, just use original intervals
    if len(tree) == 0:
        tree.addi(pos, stop)
    return tree


def get_supported_intervals_on_raw_region_rescue(record, dangling_fraction):
    # Get cleaned up support over the middle region
    tree = defragment_region(
        intervals=record.middle_supporting_intervals,
        pos=record.pos, stop=record.end,
        dangling_fraction=dangling_fraction
    )
    # Add flanking supporting intervals
    for interval in record.left_supporting_intervals + record.right_supporting_intervals:
        tree.add(interval)
    # Merge contiguous intervals
    tree.merge_overlaps(strict=False)
    # Return resulting intervals that overlap with the middle
    region_interval = Interval(record.pos, record.end)
    intervals = list(tree.overlap(region_interval))
    # Should have at least one interval here
    if len(intervals) == 0:
        raise ValueError(f"No intervals left for raw region rescue for {record.region} / {record.sample} (bug)")
    # Take the one with the greatest reciprocal overlap
    largest_overlap_interval = sorted(intervals, key=lambda x: reciprocal_overlap(x, region_interval))[-1]
    return largest_overlap_interval


def adjudicate_raw_region_false_negatives(revise_false_negative_dict, raw_region_false_negatives, true_positives,
                                          new_partial_events_tree_dict, min_supported_valid_overlapping_intervals_frac,
                                          dangling_fraction, min_true_positive_reciprocal_overlap):
    # Re-key for fast lookup
    false_negative_svtype_sample_region_dict = defaultdict(list)
    for record_list in revise_false_negative_dict.values():
        for record in record_list:
            false_negative_svtype_sample_region_dict[(record.svtype, record.sample, record.region)].append(record)
    # make dictionary: svtype -> sample -> region -> record corresponding to new variants to create
    new_records_dict = dict()
    # Fill in missing false negatives not already rescued in an existing variant or are existing true positives
    for svtype, svtype_dict in raw_region_false_negatives.items():
        for sample, sample_dict in svtype_dict.items():
            for _, record_list in sample_dict.items():
                if len(record_list) > 1:
                    raise ValueError(f"Raw region false negative should only have one call but "
                                     f"found {len(record_list)} (bug)")
                record = record_list[0]
                region = record.region
                chrom = record.chrom
                region_interval = Interval(record.pos, record.end)
                is_true_positive = svtype in true_positives and sample in true_positives[svtype] \
                    and region in true_positives[svtype][sample] and any(
                    reciprocal_overlap(region_interval, interval) >= min_true_positive_reciprocal_overlap
                    for interval in true_positives[svtype][sample][region])
                is_rescued = False
                # Check that a genotype revision rescued this region
                record_size = record.end - record.pos
                if (not is_true_positive) and (svtype, sample, region) in false_negative_svtype_sample_region_dict:
                    false_negative_record_list = false_negative_svtype_sample_region_dict[(svtype, sample, region)]
                    total_overlap = 0
                    for false_negative_record in false_negative_record_list:
                        interval = Interval(false_negative_record.pos, false_negative_record.end)
                        total_overlap += overlap_size(interval, region_interval)
                    if total_overlap / record_size >= min_supported_valid_overlapping_intervals_frac:
                        is_rescued = True
                # Check that a partial event rescued this region
                if (not is_true_positive) and (not is_rescued) \
                        and (svtype, sample, chrom) in new_partial_events_tree_dict:
                    # Since records could span multiple regions, we must intersect with this region
                    tree = new_partial_events_tree_dict[(svtype, sample, chrom)]
                    overlappers = tree.overlap(region_interval)
                    total_overlap = sum(overlap_size(interval, region_interval) for interval in overlappers)
                    if total_overlap / record_size >= min_supported_valid_overlapping_intervals_frac:
                        is_rescued = True
                if (not is_true_positive) and (not is_rescued):
                    if svtype not in new_records_dict:
                        new_records_dict[svtype] = dict()
                    if sample not in new_records_dict[svtype]:
                        new_records_dict[svtype][sample] = dict()
                    # We will rescue the call by creating a new record
                    # Get supported intervals only
                    new_intervals = get_revised_intervals([record], [region], record.pos, record.end, dangling_fraction)
                    new_records = list()
                    for interval in new_intervals:
                        # Subset genotypes to the interval
                        genotype_medians = [g for g in record.genotype_medians
                                            if overlap_size((g.pos, g.stop), interval) > 0]
                        # Copy record but with revised interval coordinates
                        new_records.append(GDRMatch(name=record.name, chrom=record.chrom,
                                                    pos=interval.begin, end=interval.end, sample=record.sample,
                                                    svtype=record.svtype,
                                                    overlapping_region_intervals=record.overlapping_region_intervals,
                                                    valid_region_intervals=record.valid_region_intervals,
                                                    left_supporting_intervals=record.left_supporting_intervals,
                                                    middle_supporting_intervals=record.middle_supporting_intervals,
                                                    right_supporting_intervals=record.right_supporting_intervals,
                                                    region=record.region,
                                                    n_region_subdivisions=record.n_region_subdivisions,
                                                    genotype_medians=genotype_medians))
                    new_records_dict[svtype][sample][region] = new_records
    # Merge overlapping calls
    merged_new_records_dict = dict()
    for svtype, svtype_dict in new_records_dict.items():
        for sample, sample_dict in svtype_dict.items():
            tree = IntervalTree()
            for region, record_list in sample_dict.items():
                for record in record_list:
                    tree.addi(record.pos, record.end, record)
            # Arbitrarily use first record when merging, which shouldn't affect the end result
            tree.merge_overlaps(strict=False, data_reducer=lambda x, y: x)
            if svtype not in merged_new_records_dict:
                merged_new_records_dict[svtype] = dict()
            if sample not in merged_new_records_dict[svtype]:
                merged_new_records_dict[svtype][sample] = list()
            for interval in tree:
                record = interval.data
                record.pos = interval.begin
                record.end = interval.end
                merged_new_records_dict[svtype][sample].append(record)
    return merged_new_records_dict


def revise_genotypes_and_create_records(input_vcf_path, revised_genotypes_tsv_path, revised_records_before_update_path,
                                        revised_records_after_update_path, unsorted_new_records_path,
                                        revision_manifest_path,
                                        region_intervals_dict, false_positive_matches, false_negative_matches,
                                        raw_region_false_negatives, true_positives, dangling_fraction,
                                        min_false_negative_rescue_overlap,
                                        min_supported_valid_overlapping_intervals_frac,
                                        min_true_positive_reciprocal_overlap,
                                        ploidy_table_dict, batch):
    false_positives_dict = remap_overlapper_dict_by_vid(false_positive_matches)
    false_negative_dict = remap_overlapper_dict_by_sample_and_region(false_negative_matches)
    revise_false_negative_dict = \
        adjudicate_false_negatives(
            false_negative_dict=false_negative_dict,
            min_false_negative_rescue_overlap=min_false_negative_rescue_overlap
        )
    with gzip.open(revision_manifest_path, "wt") as f_manifest, \
            pysam.VariantFile(input_vcf_path) as f_in:
        # Find invalidated records and reset their genotypes
        with gzip.open(revised_genotypes_tsv_path, mode="wt") as f_geno, \
                pysam.VariantFile(revised_records_before_update_path, mode="w", header=f_in.header) as f_before, \
                pysam.VariantFile(revised_records_after_update_path, mode="w", header=f_in.header) as f_after:
            revise_genotypes(f_in=f_in, f_geno=f_geno, f_before=f_before, f_after=f_after, f_manifest=f_manifest,
                             batch=batch, false_positives_dict=false_positives_dict,
                             revise_false_negative_dict=revise_false_negative_dict)
        # Revise invalidated records
        with pysam.VariantFile(revised_records_before_update_path) as f_before, \
                pysam.VariantFile(unsorted_new_records_path, mode="w", header=f_in.header) as f_new:
            new_partial_events_tree_dict = revise_partially_supported_variants(
                f_before=f_before,
                f_new=f_new,
                f_manifest=f_manifest,
                batch=batch,
                false_positives_dict=false_positives_dict,
                dangling_fraction=dangling_fraction
            )
            new_records_dict = \
                adjudicate_raw_region_false_negatives(
                    revise_false_negative_dict=revise_false_negative_dict,
                    raw_region_false_negatives=raw_region_false_negatives,
                    true_positives=true_positives,
                    new_partial_events_tree_dict=new_partial_events_tree_dict,
                    min_supported_valid_overlapping_intervals_frac=min_supported_valid_overlapping_intervals_frac,
                    dangling_fraction=dangling_fraction,
                    min_true_positive_reciprocal_overlap=min_true_positive_reciprocal_overlap
                )
            create_new_variants(f_new=f_new, f_manifest=f_manifest, new_records_dict=new_records_dict,
                                region_intervals_dict=region_intervals_dict, batch=batch,
                                ploidy_table_dict=ploidy_table_dict)


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
    parser.add_argument('--batch', type=str, required=True, help='Batch name')
    parser.add_argument('--out', type=str, required=True, help='Output base path')
    parser.add_argument('--overlap', type=float, default=0.9, help='Subdivision overlap fraction cutoff')
    parser.add_argument('--del-median', type=float, default=0.6, help='DEL RDTest median cutoff')
    parser.add_argument('--dup-median', type=float, default=1.4, help='DUP RDTest median cutoff')
    parser.add_argument('--vcf-min-size', type=int, default=1000, help='Min size of vcf variants')
    parser.add_argument('--min-region-overlap', type=float, default=0.3,
                        help='Min fraction of sufficiently overlapping genotyped intervals in the GD region')
    parser.add_argument('--min-frac-supporting-genotypes', type=float, default=0.5,
                        help='Min fraction of sufficiently overlapping genotyped intervals that must support a call')
    parser.add_argument('--min-dangling-frac', type=float, default=0.3,
                        help='Min interval fraction of variant size for removing leftover fragments')
    parser.add_argument('--min-par-overlap', type=float, default=0.5,
                        help='Min overlap fraction for an interval to be treated as PAR')
    parser.add_argument('--min-false-negative-rescue-overlap', type=float, default=0.8,
                        help='Min overlap fraction of a variant to rescue false negatives in an existing call. '
                             'If not met, a new variant will be created')
    parser.add_argument('--min-true-positive-reciprocal-overlap', type=float, default=0.5,
                        help='Min reciprocal overlap for a variant to be considered true positive. Used for '
                             'determining whether a raw region signal needs to be rescued with a novel variant.')
    parser.add_argument('--min-region-size', type=int, default=10000, help='Min region/subdivision size for revisions')
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
    gdr_trees, gdr_del_ids, gdr_dup_ids = create_trees_from_bed_records_by_type(bed_path=args.region_bed,
                                                                                min_region_size=args.min_region_size)
    region_intervals_dict = create_region_intervals_dict(gdr_trees)
    logging.info("Reading PAR bed file...")
    par_trees = create_trees_from_bed_records(args.par_bed)
    logging.info("Reading \"median_geno\" file...")
    medians, median_vids = read_median_geno(list_path=args.median_geno_list, del_ids=gdr_del_ids, dup_ids=gdr_dup_ids,
                                            del_cutoff=args.del_median, dup_cutoff=args.dup_median,
                                            ploidy_table_dict=ploidy_table_dict, par_trees=par_trees,
                                            min_par_overlap=args.min_par_overlap)
    logging.info("Parsing VCF for variants overlapping genomic disorder regions...")
    false_positive_matches, false_negative_matches, raw_region_false_negatives, true_positives = \
        get_overlapping_samples_vcf(vcf_path=args.vcf,
                                    gdr_trees=gdr_trees,
                                    region_intervals_dict=region_intervals_dict,
                                    medians=medians,
                                    median_vids=median_vids,
                                    cutoff=args.overlap,
                                    vcf_min_size=args.vcf_min_size,
                                    min_supported_valid_overlapping_intervals_frac=args.min_frac_supporting_genotypes,
                                    min_region_overlap=args.min_region_overlap)
    revised_genotypes_tsv_path = f"{args.out}.revised_genotypes.tsv.gz"
    revised_records_before_update_path = f"{args.out}.revised_before_update.vcf.gz"
    revised_records_after_update_path = f"{args.out}.revised_after_update.vcf.gz"
    sorted_new_records_path = f"{args.out}.new_records.vcf.gz"
    revision_manifest_path = f"{args.out}.revision_manifest.tsv.gz"
    logging.info("Revising genotypes and variants...")
    with tempfile.NamedTemporaryFile(dir=args.temp, suffix=".vcf.gz") as temp_vcf:
        revise_genotypes_and_create_records(input_vcf_path=args.vcf,
                                            revised_genotypes_tsv_path=revised_genotypes_tsv_path,
                                            revised_records_before_update_path=revised_records_before_update_path,
                                            revised_records_after_update_path=revised_records_after_update_path,
                                            unsorted_new_records_path=temp_vcf.name,
                                            revision_manifest_path=revision_manifest_path,
                                            region_intervals_dict=region_intervals_dict,
                                            false_positive_matches=false_positive_matches,
                                            false_negative_matches=false_negative_matches,
                                            raw_region_false_negatives=raw_region_false_negatives,
                                            true_positives=true_positives,
                                            dangling_fraction=args.min_dangling_frac,
                                            min_false_negative_rescue_overlap=args.min_false_negative_rescue_overlap,
                                            min_true_positive_reciprocal_overlap=args.min_true_positive_reciprocal_overlap,
                                            min_supported_valid_overlapping_intervals_frac=args.min_frac_supporting_genotypes,
                                            ploidy_table_dict=ploidy_table_dict,
                                            batch=args.batch)
        pysam.tabix_index(revised_records_before_update_path, preset="vcf", force=True)
        pysam.tabix_index(revised_records_after_update_path, preset="vcf", force=True)
        with tempfile.TemporaryDirectory(dir=args.temp) as temp_dir:
            sort_vcf(vcf_path=temp_vcf.name, out_path=sorted_new_records_path, temp_dir=temp_dir)
            pysam.tabix_index(sorted_new_records_path, preset="vcf", force=True)


if __name__ == "__main__":
    main()
