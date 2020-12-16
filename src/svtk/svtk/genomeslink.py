# -*- coding: utf-8 -*-
#
"""
genomeslink.py

Implementation of single-linkage clustering for genomic coordinates.

Not actually an implementation of the SLINK algorithm.
"""

from collections import deque
from itertools import combinations
from operator import itemgetter
import numpy as np
from scipy import sparse
from scipy.sparse import csgraph

from .utils import is_smaller_chrom


class GSNode(object):
    def __init__(self, chrA, posA, chrB, posB, name='.'):
        """
        Node in graph-based single-linkage clustering of genomic coordinates.
        """
        self.chrA = str(chrA)
        self.posA = int(posA)
        self.chrB = str(chrB)
        self.posB = int(posB)
        self.name = str(name)

        self.sort_positions()

    def sort_positions(self):
        """Force chrA, posA to be upstream of chrB, posB """
        if self.chrA == self.chrB:
            self.posA, self.posB = sorted([self.posA, self.posB])

        elif not is_smaller_chrom(self.chrA, self.chrB):
            self.chrA, self.chrB = self.chrB, self.chrA
            self.posA, self.posB = self.posB, self.posA

    def is_in(self, tabixfile):
        """
        Test if breakpoints of SV fall into any region in tabix-indexed bed.

        Parameters
        ----------
        tabixfile : pysam.TabixFile

        Returns
        -------
        is_in : bool
        """

        if tabixfile is None:
            return False

        return ((self.chrA in tabixfile.contigs and
                 any(tabixfile.fetch(self.chrA, self.posA, self.posA + 1))) or
                (self.chrB in tabixfile.contigs and
                 any(tabixfile.fetch(self.chrB, self.posB, self.posB + 1))))

    @property
    def secondary(self):
        """
        Filter function
        """
        return False

    def is_allowed_chrom(self, chroms='either'):
        """
        Whitelist chromosomes.

        Parameters
        ----------
        chroms : str or list, optional
            1) List of permitted chromosomes
            2) One of 'GRCh', 'UCSC', 'either'. Allows the 22 autosomes and sex
            chromosomes in GRCh format (1, 2, ... X, Y), UCSC format
            (chr1, chr2, ... chrX, chrY), or either format.
            Defaults to 'either'.

        Returns
        -------
        is_allowed_chrom : bool
            True if both reads are on whitelisted chromosomes.
        """

        GRCh = [str(x) for x in range(1, 23)] + 'X Y'.split()
        UCSC = ['chr' + x for x in GRCh]

        if isinstance(chroms, list):
            return (self.chrA in chroms) and (self.chrB in chroms)

        if chroms == 'GRCh':
            return (self.chrA in GRCh) and (self.chrB in GRCh)
        elif chroms == 'UCSC':
            return (self.chrA in UCSC) and (self.chrB in UCSC)
        elif chroms == 'either':
            chromlist = GRCh + UCSC
            return (self.chrA in chromlist) and (self.chrB in chromlist)
        else:
            raise Exception('Invalid chromosome list: %s ' % chroms)

    def __hash__(self):
        return id(self)

    def __eq__(self, other):
        return (self.chrA == other.chrA and
                self.posA == other.posA and
                self.chrB == other.chrB and
                self.posB == other.posB and
                self.name == other.name)

    def __lt__(self, other):
        if self.chrA == other.chrA:
            return self.posA < other.posA
        else:
            return is_smaller_chrom(self.chrA, other.chrA)

    def __le__(self, other):
        if self.chrA == other.chrA:
            return self.posA <= other.posA
        else:
            return is_smaller_chrom(self.chrA, other.chrA)

    def __str__(self):
        return ('{chrA}\t{posA}\t{posB}\t{chrA}\t{name}'.format(
                **self.__dict__))


class GenomeSLINK(object):
    def __init__(self, nodes, dist, size=1, blacklist=None, single_end=False):
        """
        Graph-based single-linkage clustering of genomic coordinates.

        Parameters
        ----------
        nodes : iterable of GSNode
            GSNodes sorted by chrA, posA
        dist : int
            Maximum clustering distance.
        size : int
            Minimum cluster size. Recommended to use 1 for call/variant
            clustering, scale up for read pair clustering.
        blacklist : pysam.TabixFile, optional
            Regions to exclude from clustering. Any node with a coordinate
            inside an excluding region is omitted. (NOTE: not overlap-based.)
        single_end : bool, optional
            Require only one end to be within min dist.
        """

        self.nodes = nodes
        self.dist = dist
        self.size = size
        self.blacklist = blacklist
        self.single_end = single_end

    def is_clusterable_with(self, first, second):
        """
        Identify cluster candidates if chrA, posB meet cluster distance limit
        """
        return (first.chrA == second.chrA and
                abs(first.posA - second.posA) < self.dist)

    def clusters_with(self, first, second):
        """
        Test if candidates meet cluster distance requirement on chrB, posB
        """
        if first.chrB == second.chrB:
            if self.single_end:
                return abs(first.posA - second.posA) < self.dist or \
                        abs(first.posB - second.posB) < self.dist
            else:
                return abs(first.posA - second.posA) < self.dist and \
                        abs(first.posB - second.posB) < self.dist
        else:
            return False

    def filter_nodes(self):
        """
        Filter provided nodes. By default, remove nodes in blacklisted regions.

        Yields
        ------
        node : GSNode
        """
        for node in self.nodes:
            # Skip nodes in blacklisted regions
            if (self.blacklist and node.is_in(self.blacklist)):
                continue
            yield node

    # TODO: add parameter for filter fn to apply to each node
    # Add `filter` method to nodes for subclassing?
    def get_candidates(self):
        """
        Find batches of SVCalls eligible for clustering.

        Requires input sorted by chromosome and position of first read in each
        pair. Pairs are collected while the first read in the next pair in the
        parser is within the maximum clustering distance of the first read of
        the previous pair.

        Yields
        ------
        candidates : deque of GSNode
        """

        candidates = deque()
        prev = None

        node_count = 0
        for node in self.filter_nodes():
            node_count += 1

            if prev is None or self.is_clusterable_with(prev, node):
                candidates.append(node)

            else:
                # Permit inequality if not parallelizing by chromosome, but
                # enforce sorted order
                n, p = node, prev
                if n.chrA != p.chrA and is_smaller_chrom(n.chrA, p.chrA):
                    msg = 'Breakend with reverse CTX ordering found'
                    print(prev, list(prev.record.samples.keys())[0])
                    print(node, list(node.record.samples.keys())[0])
                    raise Exception(msg)

                yield candidates
                candidates = deque([node])

            prev = node

        yield candidates

    def cluster_candidates(self, candidates, *args, **kwargs):
        """Batch of clustering"""
        n = len(candidates)

        # Permit clusters of size 1
        # G = sparse.eye(n, dtype=np.uint8, format='lil')
        G = sparse.lil_matrix((n, n), dtype=np.uint8)

        # Add edges between nodes with overlap
        for p1, p2 in combinations(range(n), 2):
            node1, node2 = candidates[p1], candidates[p2]
            if self.clusters_with(node1, node2):
                G[p1, p2] = 1

        # Get indices of connected components
        n_comp, comp_list = csgraph.connected_components(G, connection='weak')
        cluster_names = np.arange(n_comp)

        # Remove clusters with less than minimum size
        cluster_sizes = [np.sum(comp_list == i) for i in cluster_names]
        cluster_sizes = np.array(cluster_sizes)
        cluster_names = cluster_names[np.where(cluster_sizes >= self.size)]

        # Convert indices back to lists of Nodes
        # Sort clusters internally by first read's position
        clusters = deque()
        for cname in cluster_names:
            cluster_idx = np.where(comp_list == cname)[0]
            if len(cluster_idx) == 1:
                cluster = [candidates[cluster_idx[0]]]
            else:
                cluster = itemgetter(*cluster_idx)(candidates)
            clusters.append(sorted(cluster, key=lambda v: (v.posA, v.name)))

        # Then sort clusters by first pair's first read's position
        for cluster in sorted(clusters, key=lambda c: c[0].posA):
            yield cluster

    def cluster(self, *args, **kwargs):
        """
        Perform single linkage clustering on a candidate batch of GSNodes

        Yields
        ------
        list of GSNode
            A cluster of GSNodes
        """
        for cs in self.get_candidates():
            for cluster in self.cluster_candidates(cs, *args, **kwargs):
                yield cluster
