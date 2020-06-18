# -*- coding: utf-8 -*-
#
# Copyright © 2015 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""
Intersect SV called by PE/SR-based algorithms.

Paired-end and split-read callers provide a reasonably precise estimation of
an SV breakpoint. This program identifies variant calls that fall within
the expected margin of error made by these programs and clusters them together.
The cluster distance defaults to 500 bp but it is recommended to use the
maximum individual clustering distance across the libraries being analyzed.
(Generally median + 7 * MAD)
"""

import heapq
import re
import pkg_resources
from pysam import VariantFile
from svtk.svfile import SVFile, SVRecordCluster, SVRecord
from svtk.genomeslink import GenomeSLINK


class VCFCluster(GenomeSLINK):
    def __init__(self, vcfs,
                 dist=500, frac=0.0,
                 match_strands=True, match_svtypes=True, preserve_ids=False,
                 region=None, blacklist=None, svtypes=None,
                 preserve_genotypes=False, sample_overlap=0.0,
                 preserve_header=False):
        """
        Clustering of VCF records.

        Records are clustered with a graph-based single linkage algorithm.
        Records are linked in the graph if their breakpoints are within a
        specified distance (default 500 bp) and if they share a minimum
        reciprocal overlap (default 0.1; records of type BND are not subjected
        to the reciprocal overlap requirement).

        VCF files must be sorted.

        Parameters
        ----------
        vcfs : list of pysam.VariantFile
            Standardized VCFs to cluster
        header : pysam.VariantHeader
            VCF header to use when creating new records
        dist : int, optional
            Clustering distance. The starts and ends of two records must both
            be within this distance in order for the records to be linked.
        frac : float, optional
            Minimum reciprocal overlap for two records to be linked.
        match_strands : bool, optional
            Two records must share strandedness in order to be linked.
        preserve_ids : bool, optional
            Keep list of constituent record IDs for each cluster.
        preserve_genotypes : bool, optional
            Report best non-reference genotype for each sample.
        region : str, optional
            Genomic region to fetch for clustering. If None, all regions
            present will be clustered.
            (chrom) or (chrom:start-end)
        blacklist : pysam.TabixFile, optional
            Blacklisted genomic regions. Records in these regions will be
            removed prior to clustering.
        svtypes : list of str, optional
            SV classes to be clustered. Records with an svtype not present in
            this list will be removed prior to clustering. If no list is
            specified, all svtypes will be clustered.
        sample_overlap : float, optional
            Minimum fraction of samples to overlap to cluster variants
        """

        # Wrap VCFs as SVFiles
        self.vcfs = vcfs
        svfiles = [SVFile(vcf) for vcf in vcfs]

        # Fetch region of interest
        if region is not None:
            chrom, start, end = parse_region(region)
            for svfile in svfiles:
                svfile.fetch(chrom, start, end)

        # Merge sorted SV files
        nodes = heapq.merge(*svfiles)

        # Make lists of unique sources and samples to construct VCF header
        sources = set()
        samples = set()

        for svfile in svfiles:
            sources = sources.union(svfile.sources)
            samples = samples.union(svfile.samples)

        # Parameterize clustering
        self.frac = frac
        self.match_strands = match_strands
        self.match_svtypes = match_svtypes
        self.svtypes = svtypes
        self.preserve_ids = preserve_ids
        self.preserve_genotypes = preserve_genotypes
        self.sample_overlap = sample_overlap
        self.preserve_header = preserve_header

        # Build VCF header for new record construction
        self.samples = sorted(samples)
        self.sources = sorted(sources)
        self.header = self.make_vcf_header()

        super().__init__(nodes, dist, 1, blacklist)

    def filter_nodes(self):
        """
        Filter records before clustering.

        In addition to default removal of records in blacklisted regions,
        filter records if:
        1) Record does not belong to one of the specified SV classes.
        2) Record is marked as SECONDARY
        3) Record is not on a whitelisted chromosome (default: 1-22,X,Y)

        Yields
        ------
        node : svfile.SVRecord
        """

        for node in super().filter_nodes():
            if self.svtypes is not None and node.svtype not in self.svtypes:
                continue
            if 'SECONDARY' in node.record.info:
                continue
            if not node.is_allowed_chrom():
                continue
            yield node

    def cluster(self, merge=True):
        """
        Yields
        ------
        record : SVRecord
        """
        clusters = super().cluster(frac=self.frac,
                                   match_strands=self.match_strands,
                                   match_svtypes=self.match_svtypes,
                                   sample_overlap=self.sample_overlap)
        
        for records in clusters:
            cluster = SVRecordCluster(records)

            if merge:
                record = self.header.new_record()
                record = cluster.merge_record_data(record)
                record = cluster.merge_record_formats(record, self.sources,
                                                      self.preserve_genotypes)
                record = cluster.merge_record_infos(record, self.header)
                if self.preserve_ids:
                    record.info['MEMBERS'] = tuple(r.record.id for r in records)

                if SVRecord(record).is_in(self.blacklist):
                    continue
                yield record
            else:
                yield cluster

    def make_vcf_header(self):
        """
        Add samples and sources to VCF template header.

        Returns
        -------
        pysam.VariantHeader
        """

        if self.preserve_header:
            header = self.vcfs[0].header
            for sample in self.samples:
                if sample not in header.samples:
                    header.add_sample(sample)

            if self.preserve_ids and 'MEMBERS' not in header.info.keys():
                info = ('##INFO=<ID=MEMBERS,Number=.,Type=String,'
                        'Description="IDs of cluster\'s constituent records.">')
                header.add_line(info)

            return header

        # Read stock template
        template = pkg_resources.resource_filename(
                        'svtk', 'data/vcfcluster_template.vcf')

        # Make header
        template = VariantFile(template)
        header = template.header

        # Add samples
        for sample in self.samples:
            header.add_sample(sample)

        # Add contigs
        contigs = []
        for vcf in self.vcfs:
            for contig in vcf.header.contigs.values():
                tup = (contig.name, contig.length)
                if tup not in contigs:
                    contigs.append(tup)

        contig_line = '##contig=<ID={0},length={1}>'
        for contig in contigs:
            header.add_line(contig_line.format(*contig))

        # Add INFO
        infos = []
        for vcf in self.vcfs:
            for tag, info in vcf.header.info.items():
                if tag in header.info.keys():
                    continue
                tup = (info.name, info.number, info.type, info.description)
                if tup not in infos:
                    infos.append(tup)
        
        info_line = '##INFO=<ID={0},Number={1},Type={2},Description="{3}">'
        for info in infos:
            header.add_line(info_line.format(*info))
        
        if self.preserve_ids and 'MEMBERS' not in header.info.keys():
            info = ('##INFO=<ID=MEMBERS,Number=.,Type=String,'
                    'Description="IDs of cluster\'s constituent records.">')
            header.add_line(info)

        # Add source
        sourcelist = sorted(set(self.sources))
        header.add_line('##source={0}'.format(','.join(sourcelist)))

        # Add source FORMAT fields
        meta = ('##FORMAT=<ID={0},Number=1,Type=Integer,'
                'Description="Called by {1}"')
        for source in self.sources:
            header.add_line(meta.format(source, source.capitalize()))

        return header


def parse_region(region):
    """
    Parameters
    ----------
    region : str
        (chrom) or (chrom:start-end)

    Returns
    -------
    chrom : str
    start : int or None
    end : int or None
    """

    # Assume only contig specified
    if ':' not in region:
        return region, None, None

    exp = re.compile(r'(.*):(\d+)-(\d+)')
    match = exp.match(region)
    chrom, start, end = match.group(1, 2, 3)

    return chrom, int(start), int(end)
