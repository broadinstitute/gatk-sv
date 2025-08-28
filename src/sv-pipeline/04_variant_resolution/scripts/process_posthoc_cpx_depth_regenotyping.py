#!/usr/bin/env python

import argparse
import gzip
import logging
import sys

from collections import Counter
from statistics import median
from typing import List, Text, Dict, Optional, Set

from pysam import VariantFile
from intervaltree import IntervalTree

# Parameters
MIN_SIZE = 1000
MIN_DIFF = 0.4
MIN_SIZE_IDEL = 150


def interval_string(chrom, start, end):
    return f"{chrom}:{start}-{end}"


class VariantInfo:

    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.genotype_copy_numbers = dict()

    def __str__(self):
        return str(vars(self))

    def __repr__(self):
        return str(vars(self))


class VariantIntervalRecord:

    def __init__(self, chrom, start, end, merged_vid, carrier_del, carrier_wt, carrier_dup, control_del, control_wt,
                 control_dup, diff_case_control_del_frac, diff_case_control_dup_frac):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.merged_vid = merged_vid

        tokens = merged_vid.split(';')
        if len(tokens) < 4:
            raise ValueError(f"Encountered record with fewer than 4 semicolon-delimited tokens: {merged_vid}")
        self.vid = tokens[0]
        self.cpx_type = tokens[1]
        self.sv_type = tokens[2]
        self.intervals_str_list = tokens[3].split(',')

        self.carrier_del = carrier_del
        self.carrier_wt = carrier_wt
        self.carrier_dup = carrier_dup
        self.control_del = control_del
        self.control_wt = control_wt
        self.control_dup = control_dup
        self.diff_case_control_del_frac = diff_case_control_del_frac
        self.diff_case_control_dup_frac = diff_case_control_dup_frac

    def __str__(self):
        return str(vars(self))

    def __repr__(self):
        return str(vars(self))

    def size(self):
        return self.end - self.start


class CleanedVariantIntervalRecord:

    def __init__(self, chrom, start, end, interval_type, vid, variant_sv_type, variant_cpx_type,
                 carrier_del, carrier_wt, carrier_dup, control_del, control_wt,
                 control_dup, diff_case_control_del_frac, diff_case_control_dup_frac, cnv_assessment):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.interval_type = interval_type
        self.vid = vid

        self.variant_sv_type = variant_sv_type
        self.variant_cpx_type = variant_cpx_type

        self.carrier_del = carrier_del
        self.carrier_wt = carrier_wt
        self.carrier_dup = carrier_dup
        self.control_del = control_del
        self.control_wt = control_wt
        self.control_dup = control_dup
        self.diff_case_control_del_frac = diff_case_control_del_frac
        self.diff_case_control_dup_frac = diff_case_control_dup_frac
        self.cnv_assessment = cnv_assessment

    def __str__(self):
        return str(vars(self))

    def __repr__(self):
        return str(vars(self))

    def size(self):
        return self.end - self.start


class VcfRecord:

    def __init__(self, variant_record, end_dict):
        self.vid = variant_record.id

        # Sink is just the normal coordinates but 0-based / open
        # TODO despite implementing the same strategy as svtk vcf2bed, small differences in some CPX intervals
        #  may occur because pysam forces stop to be max(pos, end)
        self.sink_chrom = variant_record.chrom
        start, end = sorted([variant_record.pos, end_dict.get(self.vid, variant_record.stop)])
        self.sink_start = max(0, start - 1)
        if variant_record.info.get("UNRESOLVED_TYPE", "").startswith("INVERSION_SINGLE_ENDER"):
            self.sink_end = variant_record.info.get("END2", variant_record.stop)
        else:
            self.sink_end = end

        # Source interval
        if "SOURCE" in variant_record.info:
            source_tokens = variant_record.info["SOURCE"].replace("_", "\t").replace(":", "\t") \
                .replace("-", "\t").split("\t")
            if len(source_tokens) != 4:
                raise ValueError(f"Record {variant_record.vid} has SOURCE field with unexpected format: "
                                 f"{variant_record.info['SOURCE']}")
            self.source_chrom = source_tokens[1]
            self.source_start = int(source_tokens[2])
            self.source_end = int(source_tokens[3])
        else:
            self.source_chrom = None
            self.source_start = None
            self.source_end = None

    def __str__(self):
        return str(vars(self))

    def __repr__(self):
        return str(vars(self))

    def sink_size(self):
        return None if self.sink_chrom is None else self.sink_end - self.sink_start

    def source_size(self):
        return None if self.source_chrom is None else self.source_end - self.source_start

    def source_string(self):
        return interval_string(self.source_chrom, self.source_start, self.source_end)

    def sink_string(self):
        return interval_string(self.sink_chrom, self.sink_start, self.sink_end)

    def source_len(self):
        return self.source_end - self.source_start

    def sink_len(self):
        return self.sink_end - self.sink_start


class VariantAssessment:
    def __init__(self, vid, modification, reason, new_sv_type, new_cpx_type, new_cpx_intervals, new_svlen,
                 new_source, new_start, new_end):
        self.vid = vid
        self.modification = modification
        self.reason = reason
        self.new_sv_type = new_sv_type
        self.new_cpx_type = new_cpx_type
        self.new_cpx_intervals = new_cpx_intervals
        self.new_svlen = new_svlen
        self.new_source = new_source
        self.new_start = new_start
        self.new_end = new_end

    def __str__(self):
        return str(vars(self))

    def __repr__(self):
        return str(vars(self))


def has_reciprocal_overlap(chrom1, start1, end1, chrom2, start2, end2, min_reciprocal_overlap):
    overlap = get_overlap(chrom1, start1, end1, chrom2, start2, end2)
    size_1 = end1 - start1
    size_2 = end2 - start2
    return overlap / max(size_1, size_2, 1) >= min_reciprocal_overlap


def get_overlap(chrom1, start1, end1, chrom2, start2, end2):
    if chrom1 != chrom2:
        return 0
    if not (start1 <= end2 and start2 <= end1):
        return 0
    return min(end1, end2) - max(start1, start2)


def get_reciprocal_overlaps(chrom, start, end, tree_dict, min_reciprocal_overlap):
    if chrom not in tree_dict:
        return list()
    return [interval.data for interval in tree_dict[chrom].overlap(start, end)
            if has_reciprocal_overlap(chrom, start, end, chrom,
                                      interval.data.start, interval.data.end, min_reciprocal_overlap)]


def genotype_counts_per_variant(intervals_path: Text,
                                genotype_dict: Dict,
                                chr_x: Text,
                                chr_y: Text,
                                all_samples: Set[Text],
                                male_samples: Set[Text],
                                female_samples: Set[Text]):
    def _count_del_dup_wt(variants_list, samples, control_cn):
        num_del = 0
        num_wt = 0
        num_dup = 0
        for var in variants_list:
            for s in samples:
                copy_number = var.genotype_copy_numbers.get(s, None)
                if copy_number is not None:
                    if copy_number < control_cn:
                        num_del += 1
                    elif copy_number == control_cn:
                        num_wt += 1
                    else:
                        num_dup += 1
        return num_del, num_wt, num_dup

    intervals_list_dict = dict()
    record_list = list()
    with gzip.open(intervals_path, 'rt') as f:
        for line in f:
            if line.startswith("#"):
                continue
            tokens = line.strip().split('\t')
            if len(tokens) < 5:
                raise ValueError(f"Less than 5 columns in intervals record: {line.strip()}")
            chrom = tokens[0]
            start = int(tokens[1])
            end = int(tokens[2])
            vid = tokens[3]
            carrier_samples = set(tokens[4].split(','))

            # Add to intervals list for referencing later
            if vid not in intervals_list_dict:
                intervals_list_dict[vid] = list()
            intervals_list_dict[vid].append(f"{chrom}:{start}-{end}")

            # Determine carrier and control samples for this variant
            if chrom == chr_x:
                female_carriers = carrier_samples.intersection(female_samples)
                if len(female_carriers) > 0:
                    # Only use females if there are female carriers
                    carrier_samples = female_carriers
                    control_samples = female_samples - carrier_samples
                    default_median_cn = 2
                else:
                    # Otherwise use only male samples
                    carrier_samples = carrier_samples.intersection(male_samples)
                    control_samples = male_samples - carrier_samples
                    default_median_cn = 1
            elif chrom == chr_y:
                # Only use male samples
                carrier_samples = carrier_samples.intersection(male_samples)
                control_samples = male_samples - carrier_samples
                default_median_cn = 1
            else:
                # Autosomal
                control_samples = all_samples - carrier_samples
                default_median_cn = 2

            # Get valid matching variant genotype records
            variant_info_list = genotype_dict[vid]
            if end - start < 1000000:
                matching_variants = [var for var in variant_info_list
                                     if has_reciprocal_overlap(chrom, start, end, var.chrom, var.start, var.end, 0.95)]
            else:
                matching_variants = [var for var in variant_info_list
                                     if get_overlap(chrom, start, end, var.chrom, var.start, var.end) ==
                                     (var.end - var.start)]

            # Expected copy number to determine whether each genotype copy number is del/dup/wt
            control_cn = [var.genotype_copy_numbers[s] for s in control_samples for var in matching_variants
                          if var.genotype_copy_numbers.get(s, None) is not None]
            if len(control_cn) == 0:
                # If no valid controls, assume default
                median_control_cn = default_median_cn
            else:
                median_control_cn = median(control_cn)

            # Count carrier and non-carrier genotypes as del/dup/wt
            carrier_del, carrier_wt, carrier_dup = _count_del_dup_wt(variants_list=matching_variants,
                                                                     samples=carrier_samples,
                                                                     control_cn=median_control_cn)
            control_del, control_wt, control_dup = _count_del_dup_wt(variants_list=matching_variants,
                                                                     samples=control_samples,
                                                                     control_cn=median_control_cn)

            # Compute differences between carrier and control fractions for del/dup
            n_carrier = carrier_del + carrier_wt + carrier_dup
            if n_carrier == 0:
                # Edge case with no valid carriers
                diff_case_control_del_frac = 0
                diff_case_control_dup_frac = 0
            else:
                n_control = max(control_del + control_wt + control_dup, 1)
                diff_case_control_del_frac = (carrier_del / n_carrier) - (control_del / n_control)
                diff_case_control_dup_frac = (carrier_dup / n_carrier) - (control_dup / n_control)

            record_list.append(VariantIntervalRecord(chrom, start, end, vid,
                                                     carrier_del, carrier_wt, carrier_dup,
                                                     control_del, control_wt, control_dup,
                                                     diff_case_control_del_frac, diff_case_control_dup_frac))
    return intervals_list_dict, record_list


def clean_up_intervals(genotype_counts, intervals_list_dict, genotype_counts_tree_dict):
    visited_vids = set()  # Only visit each unique extended vid (w/ interval info)
    cleaned_genotype_counts = dict()
    for query_record in genotype_counts:
        if query_record.merged_vid in visited_vids:
            continue
        else:
            visited_vids.add(query_record.merged_vid)
        intervals = query_record.intervals_str_list
        if len(intervals) == 0 or (len(intervals) == 1 and intervals[0] == "NA"):
            # If not available, pull from the intervals list with matching merged vid
            if query_record.merged_vid not in intervals_list_dict:
                sys.stderr.write(f"Warning: Could not reference vid {query_record.merged_vid} in intervals\n")
            else:
                intervals = ["UNK_" + s for s in intervals_list_dict[query_record.merged_vid]]
        for interval_str in intervals:
            tokens = interval_str.replace(":", "\t").replace("_", "\t").replace("-", "\t").split("\t")
            interval_type = tokens[0]
            interval_chrom = tokens[1]
            interval_start = int(tokens[2])
            interval_end = int(tokens[3])
            if interval_chrom in genotype_counts_tree_dict:
                matches = get_reciprocal_overlaps(interval_chrom, interval_start, interval_end,
                                                  genotype_counts_tree_dict, 0.95)
                for matching_record in matches:
                    # Assign CNV type depending on frac diff support
                    if matching_record.diff_case_control_del_frac > MIN_DIFF:
                        if matching_record.diff_case_control_dup_frac > MIN_DIFF:
                            cnv_assessment = "DELDUP"
                        else:
                            cnv_assessment = "DEL"
                    elif matching_record.diff_case_control_dup_frac > MIN_DIFF:
                        cnv_assessment = "DUP"
                    elif matching_record.end - matching_record.start < MIN_SIZE:
                        cnv_assessment = "TOO_SMALL"
                    else:
                        cnv_assessment = "WT"
                    record = CleanedVariantIntervalRecord(interval_chrom, interval_start, interval_end,
                                                          interval_type, matching_record.vid, matching_record.sv_type,
                                                          matching_record.cpx_type, matching_record.carrier_del,
                                                          matching_record.carrier_wt, matching_record.carrier_dup,
                                                          matching_record.control_del, matching_record.control_wt,
                                                          matching_record.control_dup,
                                                          matching_record.diff_case_control_del_frac,
                                                          matching_record.diff_case_control_dup_frac, cnv_assessment)
                    if record.vid not in cleaned_genotype_counts:
                        cleaned_genotype_counts[record.vid] = [record]
                    else:
                        cleaned_genotype_counts[record.vid].append(record)
    return cleaned_genotype_counts


def parse_ends(vcf_path: Text) -> Dict[Text, int]:
    """
    Since pysam automatically changes invalid END fields (i.e. when less than the start position), they must
    be parsed manually.

    Parameters
    ----------
    vcf_path: Text
        input vcf path

    Returns
    -------
    header: Dict[Text, int]
        map from variant ID to END position
    """
    end_dict = dict()
    with gzip.open(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            columns = line.split('\t', 8)
            vid = columns[2]
            info = columns[7]
            info_tokens = info.split(';')
            end_field_list = [x for x in info_tokens if x.startswith("END=")]
            if len(end_field_list) > 0:
                end = int(end_field_list[0].replace("END=", ""))
            else:
                # Special case where END and POS happen to be equal
                end = int(columns[1])
            end_dict[vid] = end
    return end_dict


def read_vcf(vcf_path, cleaned_genotype_counts_vids):
    end_dict = parse_ends(vcf_path)
    with VariantFile(vcf_path) as f:
        return {r.id: VcfRecord(r, end_dict) for r in f if r.id in cleaned_genotype_counts_vids}


def final_assessment(cleaned_genotype_counts, variants_to_reclassify, min_ddup_thresh):

    def _subtract_interval(start1, end1, start2, end2):
        # Spanning case results in null (shouldn't ever do this in this context)
        if start2 <= start1 and end2 >= end1:
            return None, None
        # Non-overlapping so just return original interval
        if end2 < start1 or start2 > end1:
            return start1, end1
        if start2 <= start1:
            return end2, end1
        else:
            return start1, start2

    def _evaluate_cpx(r, records_list, new_sv_type, new_cpx_type):
        cnv_records_list = [m for m in records_list if m.interval_type == "DEL" or m.interval_type == "DUP"]
        num_cnv_intervals = 0
        num_failed_cnv_intervals = 0
        num_validated_cnv_intervals = 0
        # Check whether each CNV interval larger than min size failed to validate
        # or if at least one CNV interval validated
        # or not enough information to make a call, so keep it
        for m in cnv_records_list:
            num_cnv_intervals += 1
            if m.cnv_assessment == "WT":
                num_failed_cnv_intervals += 1
            elif m.interval_type == m.cnv_assessment or m.cnv_assessment == "DELDUP":
                num_validated_cnv_intervals += 1
        if num_cnv_intervals == 0:
            mod = "KEEP"
            reason = "NO_PREDICTED_CNV_INTERVALS"
        elif num_failed_cnv_intervals > 0:
            mod = "UNRESOLVED"
            reason = "AT_LEAST_ONE_CNV_INTERVAL_FAILED"
        elif num_validated_cnv_intervals > 0:
            mod = "KEEP"
            reason = "AT_LEAST_ONE_CNV_VALIDATED"
        else:
            mod = "KEEP"
            reason = "NO_CNVS_LARGER_THAN_MIN_SIZE"
        return VariantAssessment(
            vid=r.vid,
            modification=mod,
            reason=reason,
            new_sv_type=new_sv_type,
            new_cpx_type=new_cpx_type,
            new_cpx_intervals=None,
            new_svlen=None,
            new_source=None,
            new_start=None,
            new_end=None
        )

    # Solve translocational insertions
    # A2B: A interval inserted into B
    # B2A: B interval inserted into A
    # CTX_INS_A2B: A interval inserted into B, interchromosomal, can be inverted
    # CTX_INS_B2A: B interval inserted into A, interchromosomal, can be inverted
    # OTHERWISE: Just test for INS-iDEL
    # Note: CTX_INS_B2A variants are not resolved here, because they are indexed
    #  onto a different chromosome, and will be dealt with in that shard
    def _evaluate_ctx_ins(r, records_list, default_sv_type, default_cpx_type):
        if r.source_chrom is None:
            return VariantAssessment(
                vid=r.vid,
                modification="SKIP",
                reason="NO_SOURCE_INTERVAL_IN_CURRENT_SHARD",
                new_sv_type=default_sv_type,
                new_cpx_type=default_cpx_type,
                new_cpx_intervals=None,
                new_svlen=None,
                new_source=None,
                new_start=None,
                new_end=None
            )
        source_overlaps = [m for m in records_list
                           if has_reciprocal_overlap(r.source_chrom, r.source_start, r.source_end,
                                                     m.chrom, m.start, m.end, 0.95)]
        source_is_dup = any(m.cnv_assessment == "DUP" for m in source_overlaps)
        # TODO: could overlap against r.sink_*
        sinks = [m for m in records_list
                 if not has_reciprocal_overlap(r.source_chrom, r.source_start, r.source_end,
                                               m.chrom, m.start, m.end, 0.95)]
        if len(sinks) > 0:
            # TODO Unstable ordering
            sink_chrom = sinks[0].chrom
            sink_start = sinks[0].start
            sink_end = sinks[0].end
            sink_size = sinks[0].size()
            if sink_size >= MIN_SIZE_IDEL:
                sink_is_del = any([s.cnv_assessment == "DEL" for s in sinks])
            else:
                sink_is_del = False
        else:
            sink_chrom = None
            sink_start = None
            sink_end = None
            sink_size = None
            sink_is_del = False

        # Determine if the variant needs to be modified
        if source_is_dup and sink_is_del:
            # If source is duplicated and sink is deleted, reclassify as CPX, dDUP_iDEL
            return VariantAssessment(
                vid=r.vid,
                modification="RECLASSIFY",
                reason="DISPERSED_DUPLICATION_W_INSERT_SITE_DEL",
                new_sv_type="CPX",
                new_cpx_type="dDUP_iDEL",
                new_cpx_intervals=f"DUP_{r.source_string()},DEL_{interval_string(sink_chrom, sink_start, sink_end)}",
                new_svlen=r.source_end - r.source_start + sink_end - sink_start,
                new_source=f"DUP_{r.source_string()}",
                new_start=None,
                new_end=None
            )
        elif source_is_dup and not sink_is_del:
            # If source is duplicated and sink is not deleted, reclassify as CPX, dDUP
            return VariantAssessment(
                vid=r.vid,
                modification="RECLASSIFY",
                reason="DISPERSED_DUPLICATION",
                new_sv_type="CPX",
                new_cpx_type="dDUP",
                new_cpx_intervals=f"DUP_{r.source_string()}",
                new_svlen=None,
                new_source=f"DUP_{r.source_string()}",
                new_start=None,
                new_end=None
            )
        elif (not source_is_dup) and sink_is_del:
            # If source is not duplicated but sink is deleted, reclassify as CPX, INS_iDEL
            return VariantAssessment(
                vid=r.vid,
                modification="RECLASSIFY",
                reason="INSERT_SITE_DEL",
                new_sv_type="CPX",
                new_cpx_type="INS_iDEL",
                new_cpx_intervals=f"DEL_{interval_string(sink_chrom, sink_start, sink_end)}",
                new_svlen=r.source_end - r.source_start + sink_end - sink_start,
                new_source=f"INS_{r.source_string()}",
                new_start=None,
                new_end=None
            )
        else:
            # Otherwise, leave as canonical insertion
            return VariantAssessment(
                vid=r.vid,
                modification="KEEP",
                reason="NON_DUPLICATED_INSERTION",
                new_sv_type=default_sv_type,
                new_cpx_type=default_cpx_type,
                new_cpx_intervals=None,
                new_svlen=None,
                new_source=None,
                new_start=None,
                new_end=None
            )

    def _is_ctx_ins(t):
        return t == "INS_A2B" \
            or t == "INS_B2A" \
            or (t.startswith("CTX_") and t.endswith("INS_A2B")) \
            or (t.startswith("CTX_") and t.endswith("INS_B2A"))

    # Solve inverted dispersed duplications vs dupINV / dupINVdel or INVdup / delINVdup
    # DUP5/INS3 or dupINV / dupINVdel
    def _evaluate_dup5_ins3(r, records_list, default_sv_type, default_cpx_type, min_ddup_thresh):
        # Get duplication/insertion interval
        dup_chrom = r.source_chrom
        dup_start = r.source_start
        dup_end = r.source_end
        dup_size = dup_end - dup_start
        # Test dup interval for evidence of duplication
        # TODO unsure if the check should be == "DUP"

        dup_confirmed = any(m.cnv_assessment != "WT" and
                            has_reciprocal_overlap(dup_chrom, dup_start, dup_end, m.chrom, m.start, m.end, 0.95)
                            for m in records_list)
        # As long as dup doesn't explicitly fail depth genotyping, proceed
        if dup_confirmed:
            # If sink length >= MINSIZEiDEL, test sink interval for evidence of deletion
            sink_is_del = r.sink_size() >= MIN_SIZE_IDEL and any(
                m.cnv_assessment == "DEL" and has_reciprocal_overlap(r.sink_chrom, r.sink_start, r.sink_end, m.chrom,
                                                                     m.start, m.end, 0.95) for m in records_list)
            # Get inversion interval. Corresponds to min dup/ins coord and max sink coord
            inv_chrom = r.sink_chrom
            inv_start = dup_start
            inv_end = r.sink_end
            inv_size = inv_end - inv_start
            # If inversion length < MINdDUPTHRESH, classify as dupINV or dupINVdel
            if inv_size < min_ddup_thresh:
                # If sink is deleted, classify as dupINVdel
                if sink_is_del:
                    # Revise inv interval (subtracting del interval)
                    inv_start, inv_end = _subtract_interval(inv_start, inv_end, r.sink_start, r.sink_end)
                    if inv_start is None:
                        raise ValueError(f"Sink interval {r.sink_start}-{r.sink_end} spans inversion interval "
                                         f"{inv_start}-{inv_end} for record {r.vid}")
                    return VariantAssessment(
                        vid=r.vid,
                        modification="RECLASSIFY",
                        reason="DUP_FLANKED_INVERSION_WITH_DEL",
                        new_sv_type="CPX",
                        new_cpx_type="dupINVdel",
                        new_cpx_intervals=f"DUP_{interval_string(dup_chrom, dup_start, dup_end)},"
                                          f"INV_{interval_string(inv_chrom, inv_start, inv_end)},"
                                          f"DEL_{r.sink_string()}",
                        new_svlen=inv_size,
                        new_source=None,
                        new_start=min(dup_start, dup_end, inv_start, inv_end, r.sink_start, r.sink_end),
                        new_end=max(dup_start, dup_end, inv_start, inv_end, r.sink_start, r.sink_end)
                    )
                else:
                    # If sink is not clearly deleted (or too small), classify as dupINV (no del)
                    return VariantAssessment(
                        vid=r.vid,
                        modification="RECLASSIFY",
                        reason="DUP_FLANKED_INVERSION",
                        new_sv_type="CPX",
                        new_cpx_type="dupINV",
                        new_cpx_intervals=f"DUP_{interval_string(dup_chrom, dup_start, dup_end)},"
                                          f"INV_{interval_string(inv_chrom, inv_start, inv_end)}",
                        new_svlen=inv_size,
                        new_source=None,
                        new_start=min(dup_start, dup_end, inv_start, inv_end),
                        new_end=max(dup_start, dup_end, inv_start, inv_end)
                    )
            elif sink_is_del:
                # If sink is deleted, classify as dDUP_iDEL
                return VariantAssessment(
                    vid=r.vid,
                    modification="RECLASSIFY",
                    reason="INVERTED_DISPERSED_DUPLICATION_WITH_DELETION",
                    new_sv_type="CPX",
                    new_cpx_type="dDUP_iDEL",
                    new_cpx_intervals=f"DUP_{interval_string(dup_chrom, dup_start, dup_end)},"
                                      f"INV_{interval_string(dup_chrom, dup_start, dup_end)},"
                                      f"DEL_{r.sink_string()}",
                    new_svlen=dup_size + r.sink_size(),
                    new_source=f"DUP_{interval_string(dup_chrom, dup_start, dup_end)}",
                    new_start=None,
                    new_end=None
                )
            else:
                # If sink is not clearly deleted (or too small), classify as dDUP (no iDEL)
                return VariantAssessment(
                    vid=r.vid,
                    modification="RECLASSIFY",
                    reason="INVERTED_DISPERSED_DUPLICATION",
                    new_sv_type="CPX",
                    new_cpx_type="dDUP",
                    new_cpx_intervals=f"DUP_{interval_string(dup_chrom, dup_start, dup_end)},"
                                      f"INV_{interval_string(dup_chrom, dup_start, dup_end)}",
                    new_svlen=dup_size,
                    new_source=f"DUP_{interval_string(dup_chrom, dup_start, dup_end)}",
                    new_start=None,
                    new_end=None
                )
        else:
            # If dup explicitly fails depth genotyping, mark as unresolved
            return VariantAssessment(
                vid=r.vid,
                modification="UNRESOLVED",
                reason="PREDICTED_DUP_INTERVAL_FAILED_GT",
                new_sv_type=default_sv_type,
                new_cpx_type=default_cpx_type,
                new_cpx_intervals=None,
                new_svlen=None,
                new_source=None,
                new_start=None,
                new_end=None
            )

    # DUP3/INS5 or INVdup / delINVdup
    def _evaluate_dup3_ins5(r, records_list, default_sv_type, default_cpx_type, min_ddup_thresh):
        # Get duplication/insertion interval
        dup_chrom = r.source_chrom
        dup_start = r.source_start
        dup_end = r.source_end
        dup_size = dup_end - dup_start
        # Test dup interval for evidence of duplication
        # TODO should this check should be == "DUP" ?
        dup_confirmed = any(
            m.cnv_assessment != "WT" and has_reciprocal_overlap(dup_chrom, dup_start, dup_end, m.chrom, m.start, m.end,
                                                                0.95) for m in records_list)
        # As long as dup doesn't explicitly fail depth genotyping, proceed
        if dup_confirmed:
            # If sink length >= MINSIZEiDEL, test sink interval for evidence of deletion
            sink_is_del = r.sink_size() >= MIN_SIZE_IDEL and any(
                m.cnv_assessment == "DEL" and has_reciprocal_overlap(r.sink_chrom, r.sink_start, r.sink_end, m.chrom,
                                                                     m.start, m.end, 0.95) for m in records_list)
            # Get inversion interval. Corresponds to min sink coord and max dup/inv coord
            # Note this is different from the DUP5/INS3 case
            inv_chrom = r.sink_chrom
            inv_start = r.sink_start
            inv_end = dup_end
            inv_size = inv_end - inv_start
            # If inversion length < MINdDUPTHRESH, classify as dupINV or dupINVdel
            if inv_size < min_ddup_thresh:
                # If sink is deleted, classify as delINVdup
                if sink_is_del:
                    # Revise inv interval (subtracting del interval)
                    inv_start, inv_end = _subtract_interval(inv_start, inv_end, r.sink_start, r.sink_end)
                    if inv_start is None:
                        raise ValueError(f"Sink interval {r.sink_start}-{r.sink_end} spans inversion interval "
                                         f"{inv_start}-{inv_end} for record {r.vid}")
                    return VariantAssessment(
                        vid=r.vid,
                        modification="RECLASSIFY",
                        reason="DUP_FLANKED_INVERSION_WITH_DEL",
                        new_sv_type="CPX",
                        new_cpx_type="delINVdup",
                        new_cpx_intervals=f"DEL_{r.sink_string()},"
                                          f"INV_{interval_string(inv_chrom, inv_start, inv_end)},"
                                          f"DUP_{interval_string(dup_chrom, dup_start, dup_end)}",
                        new_svlen=inv_size,
                        new_source=None,
                        new_start=min(dup_start, dup_end, inv_start, inv_end, r.sink_start, r.sink_end),
                        new_end=max(dup_start, dup_end, inv_start, inv_end, r.sink_start, r.sink_end)
                    )
                else:
                    # If sink is not clearly deleted (or too small), classify as INVdup (no del)
                    return VariantAssessment(
                        vid=r.vid,
                        modification="RECLASSIFY",
                        reason="DUP_FLANKED_INVERSION",
                        new_sv_type="CPX",
                        new_cpx_type="INVdup",
                        new_cpx_intervals=f"INV_{interval_string(inv_chrom, inv_start, inv_end)},"
                                          f"DUP_{interval_string(dup_chrom, dup_start, dup_end)}",
                        new_svlen=inv_size,
                        new_source=None,
                        new_start=min(dup_start, dup_end, inv_start, inv_end),
                        new_end=max(dup_start, dup_end, inv_start, inv_end)
                    )
            elif sink_is_del:
                # If sink is deleted, classify as dDUP_iDEL
                return VariantAssessment(
                    vid=r.vid,
                    modification="RECLASSIFY",
                    reason="INVERTED_DISPERSED_DUPLICATION_WITH_DELETION",
                    new_sv_type="CPX",
                    new_cpx_type="dDUP_iDEL",
                    new_cpx_intervals=f"DUP_{interval_string(dup_chrom, dup_start, dup_end)},"
                                      f"INV_{interval_string(dup_chrom, dup_start, dup_end)},"
                                      f"DEL_{r.sink_string()}",
                    new_svlen=dup_size + r.sink_size(),
                    new_source=f"DUP_{interval_string(dup_chrom, dup_start, dup_end)}",
                    new_start=None,
                    new_end=None
                )
            else:
                # If sink is not clearly deleted (or too small), classify as dDUP (no iDEL)
                return VariantAssessment(
                    vid=r.vid,
                    modification="RECLASSIFY",
                    reason="INVERTED_DISPERSED_DUPLICATION",
                    new_sv_type="CPX",
                    new_cpx_type="dDUP",
                    new_cpx_intervals=f"INV_{interval_string(dup_chrom, dup_start, dup_end)},"
                                      f"DUP_{interval_string(dup_chrom, dup_start, dup_end)}",

                    new_svlen=dup_size,
                    new_source=f"DUP_{interval_string(dup_chrom, dup_start, dup_end)}",
                    new_start=None,
                    new_end=None
                )
        else:
            # If dup explicitly fails depth genotyping, mark as unresolved
            return VariantAssessment(
                vid=r.vid,
                modification="UNRESOLVED",
                reason="PREDICTED_DUP_INTERVAL_FAILED_GT",
                new_sv_type=default_sv_type,
                new_cpx_type=default_cpx_type,
                new_cpx_intervals=None,
                new_svlen=None,
                new_source=None,
                new_start=None,
                new_end=None
            )

    # Default to test for INS-iDEL
    def _evaluate_ins_idel(r, records_list, default_sv_type, default_cpx_type):
        if r.sink_chrom is not None:
            # TODO shell script uses bedtools intersect with -v flag (=>"not") that probably isn't supposed to be used
            sinks = [m for m in records_list
                     if not has_reciprocal_overlap(r.sink_chrom, r.sink_start, r.sink_end,
                                                   m.chrom, m.start, m.end, 0.95)]
            if len(sinks) > 0:
                # TODO unstable ordering
                sink_record = sinks[0]
                # TODO Both sink_is_del and sink_is_not_del could both be true. The logic here needs to be cleaned up.
                if sink_record.size() >= MIN_SIZE_IDEL:
                    # If insertion site size >= MINSIZEiDEL, check if insertion site is deleted
                    sink_is_del = any(m.cnv_assessment == "DEL" for m in sinks)
                    sink_is_not_del = any(m.cnv_assessment == "WT" for m in sinks)
                else:
                    sink_is_del = False
                    sink_is_not_del = False
                # TODO original script checks for source_is_dup but the variable is never assigned. We will assume False
                source_is_dup = False
                if (not source_is_dup) and sink_is_del:
                    # If sink is deleted, reclassify as CPX, INS_iDEL
                    return VariantAssessment(
                        vid=r.vid,
                        modification="RECLASSIFY",
                        reason="INSERT_SITE_DEL",
                        new_sv_type="CPX",
                        new_cpx_type="INS_iDEL",
                        new_cpx_intervals=f"DEL_{r.sink_string()}",
                        new_svlen=r.sink_size() + r.source_size(),
                        new_source=f"INS_{r.source_string()}",
                        new_start=None,
                        new_end=None
                    )
                elif sink_is_not_del and sink_record.size() >= MIN_SIZE:
                    # If sink is explicitly *not* deleted and is >= MINSIZE, reclassify as UNRESOLVED
                    return VariantAssessment(
                        vid=r.vid,
                        modification="UNRESOLVED",
                        reason="PREDICTED_LARGE_SINK_DEL_INTERVAL_FAILED_GT",
                        new_sv_type=default_sv_type,
                        new_cpx_type=default_cpx_type,
                        new_cpx_intervals=None,
                        new_svlen=None,
                        new_source=None,
                        new_start=None,
                        new_end=None
                    )
                else:
                    # Otherwise, leave as canonical insertion
                    return VariantAssessment(
                        vid=r.vid,
                        modification="KEEP",
                        reason="CANONICAL_INS_NO_SINK_DELETION",
                        new_sv_type=default_sv_type,
                        new_cpx_type=default_cpx_type,
                        new_cpx_intervals=None,
                        new_svlen=None,
                        new_source=None,
                        new_start=None,
                        new_end=None
                    )
            else:
                # No sinks found
                return VariantAssessment(
                    vid=r.vid,
                    modification="SKIP",
                    reason="NO_SINK_INTERVAL_IN_CURRENT_SHARD",
                    new_sv_type=default_sv_type,
                    new_cpx_type=default_cpx_type,
                    new_cpx_intervals=None,
                    new_svlen=None,
                    new_source=None,
                    new_start=None,
                    new_end=None
                )

    # Salvage inverted duplications from inversion single enders
    def _evaluate_bnd(r, records_list, default_sv_type, default_cpx_type):
        if not default_cpx_type.startswith("INVERSION_SINGLE_ENDER"):
            # Only consider inversion single enders
            # TODO original script does not set MOD/REASON for this case
            return VariantAssessment(
                vid=r.vid,
                modification="",
                reason="",
                new_sv_type=default_sv_type,
                new_cpx_type=default_cpx_type,
                new_cpx_intervals=None,
                new_svlen=None,
                new_source=None,
                new_start=None,
                new_end=None
            )
        # Get span of inversion
        dup_chrom = r.sink_chrom
        dup_start = r.sink_start
        dup_end = r.sink_end
        dup_size = dup_end - dup_start
        # Test dup interval for explicit confirmation of duplication
        dup_confirmed = any(
            m.cnv_assessment == "DUP" and has_reciprocal_overlap(dup_chrom, dup_start, dup_end,
                                                                 m.chrom, m.start, m.end, 0.95)
            for m in records_list
        )
        if dup_confirmed:
            # If dup confirms, resolve as palindromic inverted duplication
            # TODO Investigate why these are never being found
            new_cpx_type = "piDUP_FR" if default_cpx_type == "INVERSION_SINGLE_ENDER_" else "piDUP_RF"
            return VariantAssessment(
                vid=r.vid,
                modification="RECLASSIFY",
                reason="PALINDROMIC_INVERTED_DUPLICATION",
                new_sv_type="CPX",
                new_cpx_type=new_cpx_type,
                new_cpx_intervals=f"DUP_{interval_string(dup_chrom, dup_start, dup_end)},"
                                  f"INV_{interval_string(dup_chrom, dup_start, dup_end)}",
                new_svlen=dup_size,
                new_source=None,
                new_start=dup_start,
                new_end=dup_end
            )
        else:
            # Otherwise, leave as unresolved
            return VariantAssessment(
                vid=r.vid,
                modification="KEEP",
                reason="DID_NOT_RESOLVE_AS_piDUP",
                new_sv_type=default_sv_type,
                new_cpx_type=default_cpx_type,
                new_cpx_intervals=None,
                new_svlen=None,
                new_source=None,
                new_start=None,
                new_end=None
            )

    # Default case for other variant types
    def _evaluate_irrelevant(r, default_sv_type, default_cpx_type):
        return VariantAssessment(
            vid=r.vid,
            modification="KEEP",
            reason="IRRELEVANT_SV_TYPE",
            new_sv_type=default_sv_type,
            new_cpx_type=default_cpx_type,
            new_cpx_intervals=None,
            new_svlen=None,
            new_source=None,
            new_start=None,
            new_end=None
        )

    # SV and CPX class each is most common one
    for vid, records in cleaned_genotype_counts.items():
        if vid not in variants_to_reclassify:
            logging.warning(f"Variant ID {vid} was not found in vcf, ignoring...")
            continue
        vcf_record = variants_to_reclassify[vid]
        sv_type = Counter([r.variant_sv_type for r in records]).most_common(1)[0][0]
        cpx_type = Counter([r.variant_cpx_type for r in records]).most_common(1)[0][0]
        if sv_type == "CPX":
            result = _evaluate_cpx(vcf_record, records, sv_type, cpx_type)
        elif sv_type == "INS":
            if _is_ctx_ins(cpx_type):
                result = _evaluate_ctx_ins(vcf_record, records, sv_type, cpx_type)
            elif cpx_type == "DUP5/INS3":
                result = _evaluate_dup5_ins3(vcf_record, records, sv_type, cpx_type, min_ddup_thresh)
            elif cpx_type == "DUP3/INS5":
                result = _evaluate_dup3_ins5(vcf_record, records, sv_type, cpx_type, min_ddup_thresh)
            else:
                result = _evaluate_ins_idel(vcf_record, records, sv_type, cpx_type)
        elif sv_type == "BND":
            result = _evaluate_bnd(vcf_record, records, sv_type, cpx_type)
        else:
            result = _evaluate_irrelevant(vcf_record, sv_type, cpx_type)
        yield result


def write_vcf(input_vcf_path, output_path, final_assessment_list):
    def _modify_record(record, assessment):
        mod = assessment.modification
        if mod == "SKIP" or mod == "KEEP":
            return
        elif mod == "UNRESOLVED":
            # Note the filter status is now added during CleanVcf instead
            record.info["UNRESOLVED"] = True
            record.info["UNRESOLVED_TYPE"] = "POSTHOC_RD_GT_REJECTION"
        elif mod == "RECLASSIFY":
            if assessment.new_start:
                record.pos = assessment.new_start
            if assessment.new_end:
                record.stop = assessment.new_end
            if assessment.new_svlen:
                record.info["SVLEN"] = assessment.new_svlen
            if assessment.new_sv_type:
                record.info["SVTYPE"] = assessment.new_sv_type
                record.alts = ("<" + assessment.new_sv_type + ">",)
            if assessment.new_cpx_type:
                record.info["CPX_TYPE"] = assessment.new_cpx_type
            if assessment.new_source:
                record.info["SOURCE"] = assessment.new_source
            elif "SOURCE" in record.info:
                record.info.pop("SOURCE")
            if assessment.new_cpx_intervals:
                record.info["CPX_INTERVALS"] = assessment.new_cpx_intervals
            if "UNRESOLVED" in record.info:
                record.info.pop("UNRESOLVED")
            if "UNRESOLVED_TYPE" in record.info:
                record.info.pop("UNRESOLVED_TYPE")
            if "EVENT" in record.info:
                record.info.pop("EVENT")

    def _strip_cpx(record):
        svtype = record.info.get("SVTYPE", "")
        if svtype != "CPX" and svtype != "CTX":
            if "CPX_TYPE" in record.info:
                record.info.pop("CPX_TYPE")
            if "CPX_INTERVALS" in record.info:
                record.info.pop("CPX_INTERVALS")

    final_assessment_dict = {a.vid: a for a in final_assessment_list}
    with VariantFile(input_vcf_path) as fin:
        header = fin.header
        header.add_line("##CPX_TYPE_delINV=\"Complex inversion with 5' flanking deletion.\"")
        header.add_line("##CPX_TYPE_INVdel=\"Complex inversion with 3' flanking deletion.\"")
        header.add_line("##CPX_TYPE_dupINV=\"Complex inversion with 5' flanking duplication.\"")
        header.add_line("##CPX_TYPE_INVdup=\"Complex inversion with 3' flanking duplication.\"")
        header.add_line("##CPX_TYPE_delINVdel=\"Complex inversion with 5' and 3' flanking deletions.\"")
        header.add_line("##CPX_TYPE_dupINVdup=\"Complex inversion with 5' and 3' flanking duplications.\"")
        header.add_line("##CPX_TYPE_delINVdup=\"Complex inversion with 5' flanking deletion and 3' flanking "
                        "duplication.\"")
        header.add_line("##CPX_TYPE_dupINVdel=\"Complex inversion with 5' flanking duplication and 3' flanking "
                        "deletion.\"")
        header.add_line("##CPX_TYPE_piDUP_FR=\"Palindromic inverted tandem duplication, forward-reverse orientation.\"")
        header.add_line("##CPX_TYPE_piDUP_RF=\"Palindromic inverted tandem duplication, reverse-forward orientation.\"")
        header.add_line("##CPX_TYPE_dDUP=\"Dispersed duplication.\"")
        header.add_line("##CPX_TYPE_dDUP_iDEL=\"Dispersed duplication with deletion at insertion site.\"")
        header.add_line("##CPX_TYPE_INS_iDEL=\"Insertion with deletion at insertion site.\"")
        with VariantFile(output_path, 'w', header=header) as fout:
            for r in fin:
                if r.id in final_assessment_dict:
                    _modify_record(r, final_assessment_dict[r.id])
                _strip_cpx(r)
                fout.write(r)


def make_genotype_tree(genotype_counts):
    tree_dict = dict()
    for g in genotype_counts:
        if g.chrom not in tree_dict:
            tree_dict[g.chrom] = IntervalTree()
        tree_dict[g.chrom].addi(g.start, g.end, g)
    return tree_dict


def parse_genotypes(path: Text) -> Dict:
    genotype_list_dict = dict()
    with gzip.open(path, 'rt') as f:
        # Traverse depth genotypes table, containing one line per genotype per interval
        for line in f:
            if line.startswith('#') or line == "\n":
                continue
            tokens = line.strip().split('\t')
            chrom = tokens[0]
            start = int(tokens[1])
            end = int(tokens[2])
            vid = tokens[3]
            sample = tokens[4]
            copy_number = None if tokens[5] == "." else int(tokens[5])
            if vid not in genotype_list_dict:
                genotype_list_dict[vid] = list()
            matching_variant_info = None
            for variant_info in genotype_list_dict[vid]:
                # Find if there is already a variant matching on vid, chrom, start, and end
                # Note that one variant may map to multiple intervals
                if variant_info.chrom == chrom and variant_info.start == start and variant_info.end == end:
                    matching_variant_info = variant_info
                    break
            if matching_variant_info is None:
                # If match wasn't found, create new object for the site
                matching_variant_info = VariantInfo(chrom, start, end)
                genotype_list_dict[vid].append(matching_variant_info)
            # Add sample copy state
            matching_variant_info.genotype_copy_numbers[sample] = copy_number
    return genotype_list_dict


def parse_ped(path: Text) -> (Set, Set, Set):
    all_samples = set()
    male_samples = set()
    female_samples = set()
    with open(path) as f:
        for line in f:
            tokens = line.strip().split('\t')
            if len(tokens) < 5:
                raise ValueError(f"Encountered line with fewer than 5 columns: {line.strip()}")
            sample = tokens[1]
            all_samples.add(sample)
            if tokens[4] == "1":
                male_samples.add(sample)
            elif tokens[4] == "2":
                female_samples.add(sample)
    return all_samples, male_samples, female_samples


def write_reclassification_table(path, assessment):
    with open(path, 'w') as f:
        f.write("#VID\tMODIFICATION\tREASON\tNEW_SVTYPE\tNEW_CPX_TYPE\tNEW_CPX_INTERVALS\t"
                "NEW_SVLEN\tNEW_SOURCE\tNEW_START\tNEW_END\n")
        for r in assessment:
            new_cpx_intervals = r.new_cpx_intervals if r.new_cpx_intervals else "."
            new_source = r.new_source if r.new_source else "."
            new_svlen = r.new_svlen if r.new_svlen is not None else "."
            new_start = r.new_start if r.new_start is not None else "."
            new_end = r.new_end if r.new_end is not None else "."
            f.write(f"{r.vid}\t{r.modification}\t{r.reason}\t{r.new_sv_type}\t{r.new_cpx_type}\t{new_cpx_intervals}\t"
                    f"{new_svlen}\t{new_source}\t{new_start}\t{new_end}\n")


def _parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Reassign complex variant labels based on depth regenotyping",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--vcf', type=str, help='Input vcf')
    parser.add_argument('--intervals', type=str, help='BED file of genotyped intervals (gzipped)')
    parser.add_argument('--genotypes', type=str, help='Melted depth genotypes (gzipped)')
    parser.add_argument('--ped', type=str, help='PED family file')
    parser.add_argument('--out', type=str, help='Output file')
    parser.add_argument('--reclassification-table', type=str, help='Output reclassification table path', required=False)
    parser.add_argument('--min-ddup-thresh', type=int, help="Min DUP threshold", default=5000)
    parser.add_argument('--chrx', type=str, help='Chromosome X contig name', default='chrX')
    parser.add_argument('--chry', type=str, help='Chromosome Y contig name', default='chrY')
    parser.add_argument("-l", "--log-level", required=False, default="INFO",
                        help="Specify level of logging information, ie. info, warning, error (not case-sensitive). "
                             "Default: INFO")
    if len(argv) <= 1:
        parser.parse_args(["--help"])
        sys.exit(0)
    parsed_arguments = parser.parse_args(argv[1:])
    return parsed_arguments


def main(argv: Optional[List[Text]] = None):
    if argv is None:
        argv = sys.argv
    args = _parse_arguments(argv)

    # Set logging level from -l input
    log_level = args.log_level
    numeric_level = getattr(logging, log_level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % log_level)
    logging.basicConfig(level=numeric_level, format='%(asctime)s - %(levelname)s: %(message)s')

    all_samples, male_samples, female_samples = parse_ped(args.ped)
    genotype_dict = parse_genotypes(args.genotypes)
    intervals_list_dict, genotype_counts = genotype_counts_per_variant(intervals_path=args.intervals,
                                                                       genotype_dict=genotype_dict,
                                                                       chr_x=args.chrx,
                                                                       chr_y=args.chry,
                                                                       all_samples=all_samples,
                                                                       male_samples=male_samples,
                                                                       female_samples=female_samples)
    del genotype_dict
    genotype_counts_tree_dict = make_genotype_tree(genotype_counts)
    cleaned_genotype_counts = clean_up_intervals(genotype_counts, intervals_list_dict, genotype_counts_tree_dict)
    del genotype_counts, intervals_list_dict, genotype_counts_tree_dict
    variants_to_reclassify = read_vcf(args.vcf, cleaned_genotype_counts_vids=cleaned_genotype_counts.keys())
    assessment = list(final_assessment(cleaned_genotype_counts, variants_to_reclassify, args.min_ddup_thresh))
    del cleaned_genotype_counts, variants_to_reclassify
    if args.reclassification_table:
        write_reclassification_table(args.reclassification_table, assessment)
    write_vcf(args.vcf, args.out, assessment)


if __name__ == "__main__":
    main()
