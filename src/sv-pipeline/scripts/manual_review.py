#!/bin/env python

import argparse
import sys
import pysam
from collections import defaultdict
from typing import List, Text, Set, Dict, Optional
from itertools import chain

# Format fields to drop when converting records to CNV
CNV_FORMAT_KEYS_TO_DROP = set(['PE_GT', 'PE_GQ', 'SR_GT', 'SR_GQ'])

# For new and converted CNV records
DEFAULT_CNQ = 999
DEFAULT_GQ = 99
DEFAULT_RD_GQ = 99

# Defines how to refine genotypes spanned by erroneous deletions
DEL_SPANNED_GT_MAP = {
    (0, 0): (0, 1),
    (0, 1): (1, 1),
    (1, 1): (1, 1),
    (None, None): (None, None)
}


def _convert_to_cnv(record, fout):
    new_record = fout.new_record(id=record.id, contig=record.chrom, start=record.start, stop=record.stop,
                                 alleles=(record.alleles[0], '<CNV>'), info=record.info, filter=('MULTIALLELIC',))
    new_record.info['SVTYPE'] = 'CNV'
    if 'MINSL' in new_record.info:
        new_record.info['MINSL'] = None
    if 'NCN' in new_record.info:
        new_record.info['NCN'] = 0
    if 'NCR' in new_record.info:
        new_record.info['NCR'] = 0
    for s, gt in record.samples.items():
        new_gt = new_record.samples[s]
        new_gt['CNQ'] = DEFAULT_CNQ
        new_gt['CN'] = gt.get('RD_CN', None)
        for key, val in gt.items():
            if key == 'GT':
                new_gt['GT'] = (None,) * len(gt['GT'])
            elif key == 'EV':
                new_gt[key] = ('RD',)
            elif key == 'GQ':
                new_gt[key] = DEFAULT_GQ
            elif key == 'GT_FILTER':
                new_gt[key] = 'pass'  # currently don't filter CNVs
            elif key == 'OGQ':
                new_gt[key] = -1  # consistent with gq-recalibrator implementation
            elif key not in CNV_FORMAT_KEYS_TO_DROP:
                new_gt[key] = val
    return new_record


def _annotate_genomic_region(record, region_name):
    record.info['GD'] = region_name


def _set_filter_pass(record):
    if 'MULTIALLELIC' in record.filter:
        record.filter.clear()
        record.filter.add('MULTIALLELIC')
    else:
        record.filter.clear()
        record.filter.add('PASS')


def _annotate_manual_review(record, value):
    if record.info.get('MANUAL_REVIEW_TYPE', None) is None:
        record.info['MANUAL_REVIEW_TYPE'] = tuple()
    record.info['MANUAL_REVIEW_TYPE'] = tuple(list(record.info['MANUAL_REVIEW_TYPE']) + [value])


def _process(record, fout, remove_vids_set, spanned_del_dict, multiallelic_vids_set,
             remove_call_dict, add_call_dict, gd_dict, coords_dict, no_sex_samples, allosomes):
    if record.id in remove_vids_set:
        print(f"Deleting variant {record.id}")
        return
    if record.id in multiallelic_vids_set:
        # Convert to CNV
        print(f"Performing CNV conversion for variant {record.id}")
        record = _convert_to_cnv(record, fout)
        _annotate_manual_review(record, 'CONVERT_TO_CNV')
    if record.id in remove_call_dict:
        # Set genotypes to hom-ref
        _annotate_manual_review(record, f"DROP_{len(remove_call_dict[record.id])}_CALLS")
        _set_filter_pass(record)
        for s in remove_call_dict[record.id]:
            if s in record.samples:
                print(f"Setting genotype {s} to hom-ref for variant {record.id}")
                record.samples[s]['GT'] = (0,) * len(record.samples[s]['GT'])
    if record.id in add_call_dict:
        # Set genotypes to het
        _annotate_manual_review(record, f"ADD_{len(add_call_dict[record.id])}_CALLS")
        _set_filter_pass(record)
        for s in add_call_dict[record.id]:
            if s in record.samples:
                ploidy = len(record.samples[s]['GT'])
                print(f"Setting genotype {s} to het for variant {record.id}")
                if ploidy > 0:
                    gt = record.samples[s]['GT'] = [0] * ploidy
                    gt[-1] = 1
                    record.samples[s]['GT'] = tuple(gt)
    if record.id in gd_dict:
        print(f"Annotating genomic disorder region of variant {record.id}")
        _set_filter_pass(record)
        _annotate_genomic_region(record, gd_dict[record.id])
    if record.id in coords_dict:
        print(f"Changing coordinates of variant {record.id}")
        _set_filter_pass(record)
        _annotate_manual_review(record, 'REFINE_COORDINATES')
        value = coords_dict[record.id]
        pos = int(value[1])
        end = int(value[2])
        record.pos = pos
        record.stop = end
        record.info['SVLEN'] = end - pos + 1
    if record.id in spanned_del_dict:
        print(f"Correcting genotypes of variant {record.id} spanned an erroneous DEL")
        _set_filter_pass(record)
        _annotate_manual_review(record, f"SPANNED_DEL")
        # Flatten in case vids are repeated
        samples = [x for y in spanned_del_dict[record.id] for x in y.split(',')]
        for s in samples:
            if s in record.samples:
                gt = record.samples[s]['GT']
                record.samples[s]['GT'] = DEL_SPANNED_GT_MAP.get(gt, gt)
    if record.chrom in allosomes:
        for s in no_sex_samples:
            if s in record.samples:
                record.samples[s]['GT'] = (None, None)
    # Write variant
    fout.write(record)


def _create_new_variants(fout, new_cnv_dict):
    for key, value in new_cnv_dict.items():
        vid = key.split(';')[0]  # Sometimes semicolon delimited lists are provided
        chrom = value[0]
        pos = int(value[1])
        end = int(value[2])
        carrier_samples = set(value[3].split(','))
        is_dup = 'DUP' in vid
        is_del = 'DEL' in vid
        if is_dup:
            alt = '<DUP>'
            svtype = 'DUP'
        elif is_del:
            alt = '<DEL>'
            svtype = 'DEL'
        else:
            raise ValueError(f"Could not determine whether DEL or DUP from variant ID {vid}")
        record = fout.new_record(id=vid, contig=chrom, start=pos, stop=end,
                                 alleles=('N', alt), filter=('PASS',))
        _annotate_manual_review(record, 'NEW_VARIANT')
        record.info['SVTYPE'] = svtype
        record.info['ALGORITHMS'] = ('manual_review',)
        record.info['EVIDENCE'] = ('RD',)
        record.info['SVLEN'] = end - pos + 1
        record.info['NCN'] = 0
        record.info['NCR'] = 0
        carrier_rd_cn = 3 if is_dup else 1
        for s, gt in record.samples.items():
            if s in carrier_samples:
                gt['GT'] = (0, 1)
                gt['RD_CN'] = carrier_rd_cn
            else:
                gt['GT'] = (0, 0)
                gt['RD_CN'] = 2
            gt['EV'] = ('RD',)
            gt['GT_FILTER'] = 'pass'
            gt['GQ'] = DEFAULT_GQ
            gt['RD_GQ'] = DEFAULT_RD_GQ
        print(f"Created new variant {vid}")
        yield record


def _parse_set(path: Text) -> Set[Text]:
    if path is None:
        return set()
    with open(path, 'r') as f:
        return set([line.strip().split()[0] for line in f])


def _parse_two_column_table(path: Text) -> Dict[Text, Text]:
    if path is None:
        return dict()
    with open(path, 'r') as f:
        d = defaultdict(list)
        for line in f:
            tokens = line.strip().split()
            d[tokens[0]].append(tokens[1])
    return d


def _parse_coords_table(path: Text) -> Dict[Text, List[Text]]:
    if path is None:
        return dict()
    with open(path, 'r') as f:
        d = {}
        for line in f:
            tokens = line.strip().split()
            # { vid: [chrom, pos, end] }
            d[tokens[0]] = tokens[1:-1]
    return d


def _parse_new_cnv_table(path: Text) -> Dict[Text, List[Text]]:
    if path is None:
        return dict()
    with open(path, 'r') as f:
        d = {}
        for line in f:
            tokens = line.strip().split()
            # { vid: [chrom, pos, end, sample_list] }
            d[tokens[3]] = tokens[:3] + tokens[4:]
    return d


def _parse_no_sex_samples(path: Text) -> List[Text]:
    if path is None:
        return list()
    no_sex_samples = list()
    with open(path, 'r') as f:
        for line in f:
            tokens = line.strip().split()
            # Add sample ID if 5th column is 0
            if tokens[4] == '0':
                no_sex_samples.append(tokens[1])
    return no_sex_samples


def _parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Apply post-hoc updates to manually reviewed variants. Table columns should be delimited "
                    "with whitespace.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--vcf', type=str, help='Input vcf (defaults to stdin)')
    parser.add_argument('--out', type=str, help='Output file (defaults to stdout). May be unsorted.')

    parser.add_argument('--chr-x', type=str, default='chrX', help='Chromosome X name')
    parser.add_argument('--chr-y', type=str, default='chrY', help='Chromosome Y name')
    parser.add_argument('--ped-file', type=str, help='Ped file', required=True)

    parser.add_argument('--new-cnv-table', type=str,
                        help='Table of (1) chrom, (2) pos, (3) end, (4) unique ID possibly corresponding to IDs in '
                             '--gd-table, and (5) comma-delimited list of carrier samples, for new DEL/DUP records '
                             'to be added. The ID should contain either "DEL" or "DUP" to determine the SV type.')
    parser.add_argument('--remove-vids-list', type=str, help='List of variant IDs to remove')
    parser.add_argument('--multiallelic-vids-list', type=str,
                        help='List of variant IDs to convert to multiallelic CVNs')
    parser.add_argument('--spanned-del-table', type=str,
                        help='Table of (1) variant IDs of overlapping an erroneous spanning deletion and '
                             '(2) comma-delimited list of sample IDs whose genotypes will be corrected.')
    parser.add_argument('--remove-call-table', type=str,
                        help='Table of (1) variant ID and (2) sample ID to change to hom-ref genotypes')
    parser.add_argument('--add-call-table', type=str,
                        help='Table of (1) variant ID and (2) sample ID to change to het genotypes')
    parser.add_argument('--coords-table', type=str,
                        help='Table of (1) variant ID, (2) chrom, (3) pos, (4) end, and (5) ID of variant that led '
                             'to the new breakpoint, for changing variant coordinates.')
    parser.add_argument('--gd-table', type=str,
                        help='Table of (1) variant ID and (2) region name for annotating genomic disorder '
                             'regions (GD).')
    if len(argv) <= 1:
        parser.parse_args(["--help"])
        sys.exit(0)
    parsed_arguments = parser.parse_args(argv[1:])
    return parsed_arguments


def main(argv: Optional[List[Text]] = None):
    if argv is None:
        argv = sys.argv
    args = _parse_arguments(argv)

    if args.vcf is None:
        vcf = pysam.VariantFile(sys.stdin)
    else:
        vcf = pysam.VariantFile(args.vcf)

    header = vcf.header
    header.add_line('##INFO=<ID=GD,Number=1,Type=String,Description="Genomic disorder region">')
    header.add_line('##INFO=<ID=MANUAL_REVIEW_TYPE,Number=.,Type=String,Description="Annotation(s) for variants '
                    'modified post hoc after manual review">')
    if args.out is None:
        fout = pysam.VariantFile(sys.stdout, 'w', header=header)
    else:
        fout = pysam.VariantFile(args.out, 'w', header=header)

    remove_vids_set = _parse_set(args.remove_vids_list)
    multiallelic_vids_set = _parse_set(args.multiallelic_vids_list)
    spanned_del_dict = _parse_two_column_table(args.spanned_del_table)
    remove_call_dict = _parse_two_column_table(args.remove_call_table)
    add_call_dict = _parse_two_column_table(args.add_call_table)
    gd_dict = _parse_two_column_table(args.gd_table)
    coords_dict = _parse_coords_table(args.coords_table)
    new_cnv_dict = _parse_new_cnv_table(args.new_cnv_table)
    no_sex_samples = _parse_no_sex_samples(args.ped_file)
    print(f"Allosome genotypes of samples without an assigned sex will be set to no-call: {', '.join(no_sex_samples)}")
    for record in chain(_create_new_variants(fout=fout, new_cnv_dict=new_cnv_dict), vcf):
        _process(
            record=record,
            fout=fout,
            remove_vids_set=remove_vids_set,
            spanned_del_dict=spanned_del_dict,
            multiallelic_vids_set=multiallelic_vids_set,
            remove_call_dict=remove_call_dict,
            add_call_dict=add_call_dict,
            gd_dict=gd_dict,
            coords_dict=coords_dict,
            no_sex_samples=no_sex_samples,
            allosomes=[args.chr_x, args.chr_y]
        )
    vcf.close()
    fout.close()


if __name__ == "__main__":
    main()
