#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Applies REF_PANEL_GENOTYPES filter to a single sample VCF file given the cleaned final VCF for the reference panel cohort
#
# The general objective is to remove calls in the single sample VCF which do not match a call or calls in the reference
# panel VCF AND have variant genotypes for samples from the reference panel.
#
# DELs and DUPs are filtered based on coverage by REF panel variants of the same type with matching genotypes
#
# INV and CPX are filtered based on reciprocal overlap matching with reference panel calls.
#
# Currently does not filter INS and BND SVTYPES since these require different matching rules or
# conversions to intervals.

import argparse
import sys
import pysam
import pybedtools
import subprocess
import svtk.utils as svu


def ac_filter(feature, variant_gts_allowed, sample_to_exclude):
    variant_samples = set(feature.fields[6].split(','))
    if sample_to_exclude in variant_samples:
        variant_samples.remove(sample_to_exclude)
    return len(variant_samples) > variant_gts_allowed


def filter_on_reciprocal_overlap(single_sample_vcf_file,
                                 ref_vcf_file,
                                 svtype,
                                 case_sample,
                                 overlap_frac,
                                 variant_gts_allowed):
    single_sample_vcf = single_sample_vcf_file
    ref_vcf = ref_vcf_file

    single_sample_bed = svu.vcf2bedtool(
        single_sample_vcf, annotate_ins=False, include_samples=True, svtypes=[svtype])
    ref_bed = svu.vcf2bedtool(
        ref_vcf, annotate_ins=False, include_samples=True, svtypes=[svtype])
    ref_bed.filter(ac_filter, variant_gts_allowed=variant_gts_allowed,
                   sample_to_exclude=case_sample)

    intersection = single_sample_bed.intersect(
        ref_bed, wa=True, f=overlap_frac, r=True, v=True)

    filtered_variant_ids = []

    for intx in intersection:
        filtered_variant_ids.append(intx.name)
    return filtered_variant_ids


def filter_cnv_on_coverage(single_sample_vcf_file, ref_vcf_file, svtype, case_sample, overlap_frac, variant_gts_allowed):
    single_sample_vcf = single_sample_vcf_file
    ref_vcf = ref_vcf_file

    single_sample_bed = svu.vcf2bedtool(
        single_sample_vcf, annotate_ins=False, include_samples=True, svtypes=[svtype])
    ref_bed = svu.vcf2bedtool(
        ref_vcf, annotate_ins=False, include_samples=True, svtypes=[svtype])

    # in bash bedtools this gets the results we want:
    # bedtools coverage -a single_sample_calls.bed -b ref_panel_calls.bed -d \ # compute per-base coverage of query by intervals in ref
    #  | awk '{OFS="\t"; print $1,$2,$3,$8,$9}' \ # slim down the of data by removing sample list, extra fields
    #  | bedtools groupby -g 1,2,3,5 -c 4 -o min,max \ # group together regions with the same coverage value
    #  | awk '$5 > 0 {OFS="\t"; print $1,$2+$5-1,$2+$6}' \ # make these regions into new bed intervals
    #  | bedtools intersect -a stdin -b ref_panel_calls.bed -wb \ # print out the ref intervals that overlapped these regions
    #  | bedtools groupby  -g 1,2,3 -c 10 -o distinct\ # condense the sample lists
    #  | bedtools intersect -a single_sample_calls.bed -b stdin -wao # intersect with the query, printing the amt of overlap
    #
    # pybedtools unable to handle this pipeline without blowing up disk space due to lack of working streaming support
    #
    # subprocess streaming equivalent
    single_sample_bed.saveas('single_sample_calls.bed')
    ref_bed.saveas('ref_panel_calls.bed')

    cov_hist = subprocess.Popen(['bedtools', 'coverage', '-a', 'single_sample_calls.bed', '-b', 'ref_panel_calls.bed', '-d'],
                                stdout=subprocess.PIPE)
    cov_hist_slim = subprocess.Popen(['awk', '{OFS="\t"; print $1,$2,$3,$8,$9}'],
                                     stdin=cov_hist.stdout,
                                     stdout=subprocess.PIPE)
    cov_reg_grouped = subprocess.Popen(['bedtools', 'groupby', '-g', '1,2,3,5', '-c', '4', '-o', 'min,max'],
                                       stdin=cov_hist_slim.stdout,
                                       stdout=subprocess.PIPE)
    cov_reg_grp_fix = subprocess.Popen(['awk', '$5 > 0 {OFS="\t"; print $1,$2+$5-1,$2+$6}'],
                                       stdin=cov_reg_grouped.stdout,
                                       stdout=subprocess.PIPE)
    cov_reg_ref_ovl = subprocess.Popen(['bedtools', 'intersect', '-a', 'stdin', '-b', 'ref_panel_calls.bed', '-wb'],
                                       stdin=cov_reg_grp_fix.stdout,
                                       stdout=subprocess.PIPE)
    cov_reg_ref_cds = subprocess.Popen(['bedtools', 'groupby', '-g', '1,2,3', '-c', '10', '-o', 'distinct'],
                                       stdin=cov_reg_ref_ovl.stdout,
                                       stdout=subprocess.PIPE)
    final_intersect_process = subprocess.Popen(['bedtools', 'intersect', '-a', 'single_sample_calls.bed', '-b', 'stdin', '-wao'],
                                               stdin=cov_reg_ref_cds.stdout,
                                               stdout=open('final_merged_intersection.bed', 'w'))

    final_intersect_process.communicate()[0]  # expect this to be empty
    return_code = final_intersect_process.returncode
    if return_code != 0:
        raise Exception(
            'intersection pipeline process exited with return code ' + return_code)

    intersection = pybedtools.BedTool('final_merged_intersection.bed')
    filtered_variant_ids = []

    current_case_id = ''
    has_ref_panel_gts = False
    bases_covered_by_matching_calls = 0
    current_case_length = -1

    for intx in intersection:
        new_case_id = intx.name

        if new_case_id != current_case_id:
            if current_case_id != '':
                covered_by_matching_case_calls = (
                    bases_covered_by_matching_calls / current_case_length) > overlap_frac
                if has_ref_panel_gts and not covered_by_matching_case_calls:
                    filtered_variant_ids.append(current_case_id)
            current_case_id = new_case_id
            has_ref_panel_gts = False
            bases_covered_by_matching_calls = 0
            current_case_length = intx.end - intx.start

        variant_samples = set(intx.fields[6].split(','))
        if case_sample in variant_samples:
            variant_samples.remove(case_sample)
        if len(variant_samples) > variant_gts_allowed:
            has_ref_panel_gts = True

        if intx.fields[7] != ".":
            ref_panel_gts = set(intx.fields[10].split(','))
            if len(variant_samples - ref_panel_gts) <= variant_gts_allowed:
                bases_covered_by_matching_calls += int(intx.fields[11])
    return filtered_variant_ids


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('single_sample_vcf')
    parser.add_argument('ref_vcf')
    parser.add_argument('case_sample')
    parser.add_argument('overlap_frac', type=float)
    parser.add_argument('variant_gts_allowed', type=int)
    parser.add_argument('fout')
    args = parser.parse_args()

    case_sample = args.case_sample
    variant_gts_allowed = args.variant_gts_allowed
    overlap_frac = args.overlap_frac

    VCF_FILTERS = [
        '##FILTER=<ID=REF_PANEL_GENOTYPES,Description="Variants with ref panel genotypes that conflicted with ref panel call set">'
    ]

    filtered_variant_ids = []

    filtered_variant_ids += filter_cnv_on_coverage(pysam.VariantFile(args.single_sample_vcf),
                                                   pysam.VariantFile(
                                                       args.ref_vcf),
                                                   'DEL',
                                                   case_sample,
                                                   overlap_frac,
                                                   variant_gts_allowed
                                                   )

    filtered_variant_ids += filter_cnv_on_coverage(pysam.VariantFile(args.single_sample_vcf),
                                                   pysam.VariantFile(
                                                       args.ref_vcf),
                                                   'DUP',
                                                   case_sample,
                                                   overlap_frac,
                                                   variant_gts_allowed
                                                   )

    filtered_variant_ids += filter_on_reciprocal_overlap(pysam.VariantFile(args.single_sample_vcf),
                                                         pysam.VariantFile(
                                                             args.ref_vcf),
                                                         'INV',
                                                         case_sample,
                                                         overlap_frac,
                                                         variant_gts_allowed
                                                         )

    filtered_variant_ids += filter_on_reciprocal_overlap(pysam.VariantFile(args.single_sample_vcf),
                                                         pysam.VariantFile(
                                                             args.ref_vcf),
                                                         'CPX',
                                                         case_sample,
                                                         overlap_frac,
                                                         variant_gts_allowed
                                                         )

    single_sample_vcf = pysam.VariantFile(args.single_sample_vcf)

    header = single_sample_vcf.header
    for f in VCF_FILTERS:
        header.add_line(f)

    if args.fout in '- stdout'.split():
        fout = pysam.VariantFile(sys.stdout, 'w', header=header)
    else:
        fout = pysam.VariantFile(args.fout, 'w', header=header)

    for record in single_sample_vcf:
        if record.id in filtered_variant_ids:
            record.filter.add('REF_PANEL_GENOTYPES')
        fout.write(record)


if __name__ == '__main__':
    main()
