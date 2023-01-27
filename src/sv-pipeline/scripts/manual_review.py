#!/bin/env python

import argparse
import sys
import pysam
from collections import defaultdict
from typing import Any, List, Text, Set, Dict, Optional

# Format fields to drop when converting records to CNV
CNV_FORMAT_KEYS_TO_DROP = set(['PE_GT', 'PE_GQ', 'SR_GT', 'SR_GQ'])

DEFAULT_CNQ = 999
DEFAULT_GQ = 99


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


def _process(record, fout, remove_vids_set, multiallelic_vids_set, remove_call_dict, add_call_dict,
             gd_dict, coords_dict):
    if record.id in remove_vids_set:
        print(f"Deleting variant {record.id}")
        return
    elif record.id in multiallelic_vids_set:
        # Convert to CNV
        print(f"Performing CNV conversion for variant {record.id}")
        record = _convert_to_cnv(record, fout)
    elif record.id in remove_call_dict:
        # Set genotypes to hom-ref
        for s in remove_call_dict[record.id]:
            print(f"Setting genotype {s} to hom-ref for variant {record.id}")
            record.samples[s]['GT'] = (0,) * len(record.samples[s]['GT'])
    elif record.id in add_call_dict:
        # Set genotypes to het
        for s in add_call_dict[record.id]:
            ploidy = len(record.samples[s]['GT'])
            print(f"Setting genotype {s} to het for variant {record.id}")
            if ploidy > 0:
                gt = record.samples[s]['GT'] = [0] * ploidy
                gt[-1] = 1
                record.samples[s]['GT'] = tuple(gt)
    elif record.id in gd_dict:
        print(f"Annotating genomic disorder region of variant {record.id}")
        _annotate_genomic_region(record, gd_dict[record.id])
    elif record.id in coords_dict:
        print(f"Changing coordinates of variant {record.id}")
        value = coords_dict[record.id]
        pos = int(value[1])
        end = int(value[2])
        record.pos = pos
        record.stop = end
    # Write variant
    fout.write(record)


def _create_new_variants(fout, new_cnv_dict, gd_dict):
    for vid, value in new_cnv_dict.items():
        



def _parse_set(path: Text) -> Set[Text]:
    with open(path, 'r') as f:
        return set([line.strip().split('\t')[0] for line in f])


def _parse_two_column_table(path: Text) -> Dict[Text, Text]:
    with open(path, 'r') as f:
        d = defaultdict(list)
        for line in f:
            tokens = line.strip().split('\t')
            d[tokens[0]].append(tokens[1])
    return d


def _parse_coords_table(path: Text) -> Dict[Text, List[Text]]:
    with open(path, 'r') as f:
        d = {}
        for line in f:
            tokens = line.strip().split('\t')
            # { vid: [chrom, pos, end, breakpoint_vid] }
            d[tokens[0]] = tokens[1:]
    return d


def _parse_new_cnv_table(path: Text) -> Dict[Text, List[Text]]:
    with open(path, 'r') as f:
        d = {}
        for line in f:
            tokens = line.strip().split('\t')
            # { vid: [chrom, pos, end, sample_list] }
            d[tokens[3]] = tokens[:3] + tokens[4:]
    return d


def _parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Apply post-hoc updates to manually reviewed variants",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--vcf', type=str, help='Input vcf (defaults to stdin)')
    parser.add_argument('--out', type=str, help='Output file (defaults to stdout)')

    parser.add_argument('--remove-vids-list', type=str, required=True, help='List of variant IDs to remove')
    parser.add_argument('--multiallelic-vids-list', type=str, required=True,
                        help='List of variant IDs to convert to multiallelic CVNs')
    parser.add_argument('--remove-call-table', type=str, required=True,
                        help='Table of (1) variant ID and (2) sample ID to change to hom-ref genotypes')
    parser.add_argument('--add-call-table', type=str, required=True,
                        help='Table of (1) variant ID and (2) sample ID to change to het genotypes')
    parser.add_argument('--coords-table', type=str, required=True,
                        help='Table of (1) variant ID, (2) chrom, (3) pos, (4) end, and (5) ID of variant that led '
                             'to the new breakpoint, for changing variant coordinates.')
    parser.add_argument('--gd-table', type=str, required=True,
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
    if args.out is None:
        fout = pysam.VariantFile(sys.stdout, 'w', header=header)
    else:
        fout = pysam.VariantFile(args.out, 'w', header=header)

    remove_vids_set = _parse_set(args.remove_vids_list)
    multiallelic_vids_set = _parse_set(args.multiallelic_vids_list)
    remove_call_dict = _parse_two_column_table(args.remove_call_table)
    add_call_dict = _parse_two_column_table(args.add_call_table)
    gd_dict = _parse_two_column_table(args.gd_table)
    coords_dict = _parse_coords_table(args.coords_table)
    new_cnv_dict = _parse_new_cnv_table(args.new_cnv_table)
    for record in vcf:
        _process(record, fout,
                 remove_vids_set=remove_vids_set,
                 multiallelic_vids_set=multiallelic_vids_set,
                 remove_call_dict=remove_call_dict,
                 add_call_dict=add_call_dict,
                 gd_dict=gd_dict,
                 coords_dict=coords_dict)
    vcf.close()
    fout.close()


if __name__ == "__main__":
    main()
