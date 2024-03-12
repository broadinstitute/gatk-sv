#!/bin/python

import argparse
import sys
import logging
from collections import defaultdict
from typing import List, Text, Optional
import heapq
import pysam

_gt_sum_map = dict()


def _cache_gt_sum(gt):
    if gt is None:
        return 0
    s = _gt_sum_map.get(gt, None)
    if s is None:
        s = sum([1 for a in gt if a is not None and a > 0])
        _gt_sum_map[gt] = s
    return s


def _parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Deduplicates VCF records with identical coordinates and svtype",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--vcf', type=str, required=True, help='Sorted vcf path')
    parser.add_argument('--out', type=str, required=True, help='Output vcf path')
    if len(argv) <= 1:
        parser.parse_args(["--help"])
        sys.exit(0)
    parsed_arguments = parser.parse_args(argv[1:])
    return parsed_arguments


def get_record_key(record):
    return record.pos, record.chrom, record.stop, record.info.get("SVTYPE", "")


def get_non_ref_genotypes(record):
    for sample, gt in record.samples.items():
        if _cache_gt_sum(gt["GT"]) > 0:
            yield sample, gt


class RecordData:
    def __init__(self, record):
        self.record = record

    def add_record(self, record):
        for sample, gt in get_non_ref_genotypes(record):
            self_gt = self.record.samples[sample]
            for key in self_gt.keys():
                if key != "GT":
                    del self_gt[key]
            for key, val in gt.items():
                self.record.samples[sample][key] = val


def main(argv: Optional[List[Text]] = None):
    if argv is None:
        argv = sys.argv
    args = _parse_arguments(argv)
    logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

    with pysam.VariantFile(args.vcf) as fin, pysam.VariantFile(args.out, mode="w", header=fin.header) as fout:
        # VCFs are only sorted on CHROM and POS, so we must assume END and SVTYPE as not sorted
        record_data_dict = dict()
        pos_queue = list()
        current_chrom = None
        for record in fin:
            record_key = get_record_key(record)
            if record_key in record_data_dict:
                record_data_dict[record_key].add_record(record)
            else:
                if record.chrom != current_chrom:
                    for record_data in record_data_dict.values():
                        fout.write(record_data.record)
                    record_data_dict.clear()
                    pos_queue.clear()
                    current_chrom = record.chrom
                removed_keys = set()
                while len(pos_queue) > 0 and pos_queue[0][0] < record_key[0]:
                    key = heapq.heappop(pos_queue)
                    fout.write(record_data_dict[key].record)
                    removed_keys.add(key)
                for key in removed_keys:
                    del record_data_dict[key]
                record_data_dict[record_key] = RecordData(record)
                heapq.heappush(pos_queue, record_key)
        # Clean up remaining records
        for record_data in record_data_dict.values():
            fout.write(record_data.record)


if __name__ == "__main__":
    main()
