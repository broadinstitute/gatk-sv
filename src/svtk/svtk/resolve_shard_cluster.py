# -*- coding: utf-8 -*-

"""
"""

from svtk.vcfcluster import VCFCluster

"""
Clustering class for CPX-resolve sharding

The definition of dist and frac are redefined for this case.
"""
class ResolveShardVCFCluster(VCFCluster):
    def __init__(self, vcfs,
                 dist=500, frac=0.0,
                 match_strands=True, preserve_ids=False,
                 region=None, blacklist=None,
                 preserve_genotypes=False, sample_overlap=0.0,
                 preserve_header=False):
        super().__init__(vcfs, dist, frac, match_strands=match_strands, match_svtypes=False, preserve_ids=preserve_ids,
              region=region, blacklist=blacklist, svtypes=None, preserve_genotypes=preserve_genotypes,
              sample_overlap=sample_overlap, preserve_header=preserve_header, do_cluster=True, do_merge=False)

    def is_clusterable_with(self, first, second):
        return first.chrA == second.chrA \
               and (abs(first.posA - second.posA) < self.dist or abs(first.posB - second.posB) < self.dist)

    def clusters_with(self, first, second):
        """
        Test if candidates meet cluster distance requirement on chrB, posB
        """
        return first.chrB == second.chrB and first.overlaps(second, self.frac) \
               and (abs(first.posA - second.posA) < self.dist or abs(first.posB - second.posB) < self.dist)
