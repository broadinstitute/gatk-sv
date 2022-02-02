# -*- coding: utf-8 -*-
#

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
from svtk.utils import samples_overlap


class VCFCluster(GenomeSLINK):
    def __init__(self, vcfs,
                 dist=500, frac=0.0,
                 match_strands=True, match_svtypes=True, preserve_ids=False,
                 region=None, blacklist=None, svtypes=None,
                 preserve_genotypes=False, sample_overlap=0.0,
                 preserve_header=False,
                 do_cluster=True,
                 do_merge=True,
                 single_end=False):
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
        single_end : bool, optional
            Require only one end to be within min dist.
        """

        if (not do_cluster) and (not do_merge):
            raise ValueError('Cannot disable both clustering and merging')

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
        self.do_cluster = do_cluster
        self.do_merge = do_merge
        self.cluster_index = 0

        # Build VCF header for new record construction
        self.samples = sorted(samples)
        self.sources = sorted(sources)
        self.header = self.make_vcf_header()

        super().__init__(nodes, dist, 1, blacklist, single_end)

    def clusters_with(self, first, second):
        """
        Check if two SV cluster with each other.

        Default behavior is to check whether coordinates are within a specified
        window. Here the following additional criteria are required:
        1) SV types match
        2) SV regions share a minimum reciprocal overlap
           * Not applicable to translocations
           * Insertion "regions" are calculated as the insertion site plus the
             predicted length of the insertion.
        3) Strands of each breakpoint match (optional)
        """

        # If svtypes don't match, skip remaining calculations for efficiency
        if self.match_svtypes and first.svtype != second.svtype:
            return False

        # If both records have an INS subclass specified, require it to match
        # Otherwise, permit clustering if one or both don't have subclass
        if self.match_svtypes and first.svtype == 'INS':
            if first.record.alts[0] != second.record.alts[0]:
                if first.record.alts[0] != '<INS>' and first.record.alts[0] != '<INS>':
                    return False

        # If strands are required to match and don't, skip remaining calcs
        if self.match_svtypes and self.match_strands:
            if first.record.info['STRANDS'] != second.record.info['STRANDS']:
                return False

        clusters = (super().clusters_with(first, second) and
                    first.overlaps(second, self.frac))

        # Only compute sample overlap if a minimum sample overlap is required
        # and if records are eligible to cluster
        if clusters and self.sample_overlap > 0:
            samplesA = first.get_called_samples_set()
            samplesB = second.get_called_samples_set()
            clusters = clusters and samples_overlap(samplesA, samplesB, self.sample_overlap, self.sample_overlap)

        return clusters

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

    def cluster_by_cluster_id(self):
        """
        Identifies clusters of records with identical IDs. It is assumed node cluster records are contiguous.

        Yields: record lists, each containing a cluster
        """
        current_id = -1
        current_cluster = []
        for node in self.filter_nodes():
            node_id = node.record.info['CLUSTER']
            if node_id == current_id:
                current_cluster.append(node)
            else:
                if len(current_cluster) > 0:
                    yield current_cluster
                current_cluster = [node]
                current_id = node_id
        if len(current_cluster) > 0:
            yield current_cluster

    def cluster(self):
        """
        Yields
        ------
        record : list(SVRecord)
        """
        if self.do_cluster:
            clusters = super().cluster(frac=self.frac,
                                       match_strands=self.match_strands,
                                       match_svtypes=self.match_svtypes,
                                       sample_overlap=self.sample_overlap)
        else:
            clusters = self.cluster_by_cluster_id()

        for records in clusters:
            cluster = SVRecordCluster(records)

            if self.do_merge:
                record = self.header.new_record()
                record = cluster.merge_record_data(record)
                record = cluster.merge_record_formats(record, self.sources,
                                                      self.preserve_genotypes)
                record = cluster.merge_record_infos(record, self.header)
                if self.preserve_ids:
                    record.info['MEMBERS'] = tuple(
                        r.record.id for r in records)

                if SVRecord(record).is_in(self.blacklist):
                    continue
                yield [record]
            else:
                record_list = [r.record for r in cluster.records]
                for r in record_list:
                    r.info['CLUSTER'] = self.cluster_index
                self.cluster_index += 1
                yield record_list

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

            if (not self.do_merge) and 'CLUSTER' not in header.info.keys():
                info = ('##INFO=<ID=CLUSTER,Number=1,Type=Integer,Description="Cluster ID">')
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

        if (not self.do_merge) and 'CLUSTER' not in header.info.keys():
            info = ('##INFO=<ID=CLUSTER,Number=1,Type=Integer,Description="Cluster ID">')
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
