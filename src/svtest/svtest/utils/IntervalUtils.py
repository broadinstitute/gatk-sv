#!/usr/bin/env python

"""
Useful utilities for intervals and interval trees.
"""

from intervaltree import IntervalTree
import svtest.utils.VCFUtils as vu


# Creates dictionary of trees[sv_type][contig] from iterable records of VariantRecords
def create_trees_from_records(records, variant_types, contigs, padding=0):
    trees = {}
    variant_types_set = set(variant_types)
    for type in variant_types:
        trees[type] = {}
        for contig in contigs:
            trees[type][contig] = IntervalTree()
    for record in records:
        type = vu.get_sv_type(record, variant_types_set)
        contig = record.chrom
        if type == 'INS' or type == 'BND':
            length = 0
        else:
            length = vu.get_record_length(record)
        trees[type][contig].addi(
            record.start - padding, record.start + length + padding)
    return trees


# Creates dictionary of trees[sv_type][contig] from iterable records of VariantRecords
def create_trees_from_bed_records(records, variant_types, contigs, padding=0):
    trees = {}
    variant_types_set = set(variant_types)
    for type in variant_types:
        trees[type] = {}
        for contig in contigs:
            trees[type][contig] = IntervalTree()
    for record in records:
        type = record[3]
        if type not in variant_types_set:
            raise ValueError("Unexpected SVTYPE in bed file: %s" % type)
        contig = record[0]
        start = record[1]
        if type == 'INS' or type == 'BND':
            length = 0
        else:
            length = record[2] - record[1]
        trees[type][contig].addi(start - padding, start + length + padding)
    return trees


# Creates dictionary of trees[contig] from gzipped bed file
def create_trees_from_bed(f, contigs, padding):
    trees = {}
    for contig in contigs:
        trees[contig] = IntervalTree()
    contigs_set = set(contigs)
    for record in f:
        line = record.decode()
        if line.startswith('#'):
            continue
        record_tokens = line.strip().split('\t')
        contig = record_tokens[0]
        if contig not in contigs_set:
            continue
        start = int(record_tokens[1])
        end = int(record_tokens[2])
        trees[contig].addi(start - padding, end + padding, record_tokens)
    return trees


# Evaluates test tree[contig]
def evaluate_tree(test_tree, truth_tree, min_ro):
    tp = {}
    fp = {}
    fpi = {}
    for contig in test_tree:
        if contig in truth_tree:
            tp_contig, fp_contig, fpi_contig = evaluate_contig_tree(
                test_tree[contig], truth_tree[contig], min_ro)
            tp[contig] = tp_contig
            fp[contig] = fp_contig
            fpi[contig] = fpi_contig
        else:
            tp[contig] = 0
            fp[contig] = 0
            fpi[contig] = []
    return tp, fp, fpi

# Evaluates test IntervalTree


def evaluate_contig_tree(test_tree, truth_tree, min_ro):
    tp = 0
    fp_intervals = []
    for interval in test_tree:
        overlappers = truth_tree.overlap(interval[0], interval[1])
        has_overlapper = False
        for overlapper in overlappers:
            if has_reciprocal_overlap(interval, overlapper, min_ro):
                has_overlapper = True
                tp += 1
                break
        if not has_overlapper:
            fp_intervals.append(interval)
    fp = len(test_tree) - tp
    return tp, fp, fp_intervals


def has_reciprocal_overlap(interval_a, interval_b, min_ro):
    return min_ro <= reciprocal_overlap(interval_a, interval_b)


def reciprocal_overlap(interval_a, interval_b):
    return float(overlap_size(interval_a, interval_b)) / max(interval_size(interval_a), interval_size(interval_b))


def interval_size(interval):
    return interval[1] - interval[0]


def overlap_size(interval_a, interval_b):
    return max(0, min(interval_a[1], interval_b[1]) - max(interval_a[0], interval_b[0]))

# Sum tree sizes over contigs in trees[contig]


def tree_size(trees):
    return sum([len(t) for t in trees.values()])
