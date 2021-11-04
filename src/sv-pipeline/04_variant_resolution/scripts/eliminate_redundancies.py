#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Merge variants that were linked by per-sample bedtools merge
"""

import argparse
import sys
from collections import deque
import numpy as np
from scipy import sparse
from scipy.sparse import csgraph
import pysam
import svtk.utils as svu


def slink(record_links, record_map):
    """
    Single linkage cluster IDs based on whether they were merged by bedtools

    links : list of (str, str)
    """

    # Make lists of IDs and records
    # (VariantRecords aren't hashable)
    IDs = sorted(record_map.keys())
    records = np.array([record_map[ID] for ID in IDs])
    n = len(records)

    # Permit clusters of size 1
    G = sparse.eye(n, dtype=np.uint16, format='lil')

    # Add edges between linked variants
    for r1, r2 in record_links:
        idx1, idx2 = IDs.index(r1.id), IDs.index(r2.id)
        G[idx1, idx2] = 1

    # Get indices of connected components
    n_comp, comp_list = csgraph.connected_components(G, connection='weak')
    cluster_names = np.arange(n_comp)

    # Convert indices back to lists of records
    clusters = deque()
    for cname in cluster_names:
        cluster_idx = np.where(comp_list == cname)[0]
        cluster = records[cluster_idx]
        clusters.append(sorted(cluster, key=lambda r: (r.chrom, r.pos)))

    # Then sort clusters by first pair's first read's position
    return clusters


def merge_linked_depth_calls(vcf, ID_links):
    """
    vcf : pysam.VariantFile
    ID_links : list of (str, str)
    """

    # Make list of linked IDs and build map to corresponding records
    linked_IDs = sorted(set([ID for link in ID_links for ID in link]))
    record_map = {}

    # If a record wasn't linked with a bedtools merge, just return it
    for record in vcf:
        if record.id not in linked_IDs:
            yield record
        else:
            record_map[record.id] = record

    # Ignore links on other chromosomes
    linked_IDs = sorted(record_map.keys())
    ID_links = [l for l in ID_links
                if l[0] in linked_IDs and l[1] in linked_IDs]

    # Convert links from pairs of IDs to pairs of records
    record_links = np.empty([len(ID_links), 2], dtype=object)
    for i, link in enumerate(ID_links):
        record_links[i, 0] = record_map[link[0]]
        record_links[i, 1] = record_map[link[1]]

    clusters = slink(record_links, record_map)

    # Merge clusters
    for cluster in clusters:
        if len(cluster) == 1:
            yield cluster[0]
            continue

        # Take maximal region
        start = np.min([record.pos for record in cluster])
        end = np.max([record.stop for record in cluster])

        merged_record = cluster[0].copy()
        merged_record.pos = start
        merged_record.stop = end
        merged_record.info['SVLEN'] = end - start

        members = list(record.info['MEMBERS']) + [r.id for r in cluster]
        merged_record.info['MEMBERS'] = members

        # Take union of called samples
        svu.update_best_genotypes(
            merged_record, cluster, preserve_multiallelic=True)

        yield merged_record


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf')
    parser.add_argument('links', type=argparse.FileType('r'))
    parser.add_argument('fout')
    args = parser.parse_args()

    vcf = pysam.VariantFile(args.vcf)
    header = vcf.header

    if args.fout in '- stdout'.split():
        fout = pysam.VariantFile(sys.stdout, 'w', header=header)
    else:
        fout = pysam.VariantFile(args.fout, 'w', header=header)

    links = [tuple(line.strip().split()) for line in args.links.readlines()]

    for record in merge_linked_depth_calls(vcf, links):
        fout.write(record)


if __name__ == '__main__':
    main()
