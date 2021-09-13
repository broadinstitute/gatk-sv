#!/usr/bin/env python

import sys
import os
import pybedtools
import pysam
import numpy
import scipy.sparse
import argparse
from typing import List, Text, Optional, Iterable, Iterator, Tuple, Set, Dict, Mapping
from types import MappingProxyType
import multiprocessing


class Keys:  # static class with re-used strings (to avoid typo errors, allow easy refactoring)
    svtype = "SVTYPE"
    ins = "INS"
    deletion = "DEL"
    dup = "DUP"
    cpx = "CPX"
    cnv = "CNV"
    unresolved = "UNRESOLVED"
    cpx_intervals = "CPX_INTERVALS"
    cpx_type = "CPX_TYPE"


class Default:  # static class with default values for kwargs
    min_cpx_reciprocal_overlap = 0.1
    cnv_cpx_reciprocal_overlap = 0.5
    cnv_cpx_sample_overlap = 0.5
    cnv_cnv_reciprocal_overlap = 0.8
    cnv_cnv_sample_overlap = 0.8
    clusterable_sv_types = frozenset({Keys.deletion, Keys.dup, Keys.cnv})
    cpx_ins_classes = frozenset({"dDUP", "dDUP_iDEL", "INS_iDEL"})
    temp_dir = "/tmp"
    num_threads = multiprocessing.cpu_count()


name_field = 3
sv_type_field = 4
is_cpx_field = 5
ref_ploidy = 2  # note, even for autosome, VCFs always have ploidy=2 calls
ref_gt = (0, 0)
non_carrier_gts = {None, (None, None), (0, 0), (0, None), (None, 0)}


def _fix_coords(start: int, end: int) -> (int, int):
    """ ensure that start preceeds end, and is >= 0 """
    start, end = (start, end) if (start <= end) else (end, start)  # ensure in sorted order
    return max(start - 1, 0), end  # convert from VCF to bed format


def _get_carrier_status(
        record: pysam.VariantRecord
) -> Tuple[numpy.ndarray, numpy.ndarray]:
    """
    Get boolean numpy arrays detailing carrier status for each sample
    Parameters
    ----------
    record: VariantRecord
        pysam record for this variant
    Returns
    -------
    is_carrier: numpy.ndarray
        boolean array that is True for samples called non-ref for this Variant, and False otherwise (including no-call)
    is_ref: numpy.ndarray
        boolean array that is True for samples called ref for this Variant, and False otherwise (including no-call)
    """
    if record.info.get(Keys.svtype, None) == Keys.cnv:  # genotype is always no-call, check info.CN
        copy_numbers = [sample_rec.get("CN") for sample_rec in record.samples.itervalues()]
        is_carrier = numpy.fromiter(
            (copy_number is not None and copy_number != ref_ploidy for copy_number in copy_numbers),
            dtype=bool, count=len(copy_numbers)
        )
        is_ref = numpy.fromiter(
            (copy_number == ref_ploidy for copy_number in copy_numbers),
            dtype=bool, count=len(copy_numbers)
        )
    else:
        genotypes = [sample_rec.get("GT") for sample_rec in record.samples.itervalues()]
        is_carrier = numpy.fromiter(
            (genotype not in non_carrier_gts for genotype in genotypes), dtype=bool, count=len(genotypes)
        )
        is_ref = numpy.fromiter(
            (genotype == ref_gt for genotype in genotypes), dtype=bool, count=len(genotypes)
        )
    return is_carrier, is_ref


def _unfiltered_vcf_records_to_bed_intervals(
        vcf_records: Iterable[pysam.VariantRecord],
        is_carrier: Dict[Text, numpy.ndarray],
        is_ref: Dict[Text, numpy.ndarray],
        clusterable_sv_types: Set[Text] = Default.clusterable_sv_types,
        cpx_ins_classes: Set[Text] = Default.cpx_ins_classes
) -> Iterator[Tuple]:
    f"""
    Iterate over input VCF, yielding records that may be redundant. Also gather is_carrier and is_ref mappings.
    Parameters
    ----------
    vcf_records: Iterable[VariantRecord]
        Iterable with pysam records from input VCF file.
    is_carrier: Dict[Text, numpy.ndarray]
        Dict from variant ID to boolean array that is True for samples called non-ref for this Variant, and False
        otherwise (including no-call). NOTE: this function *updates* is_carrier in place.
    is_ref: Dict[Text, numpy.ndarray]
        Dict from boolean array that is True for samples called ref for this Variant, and False otherwise (including
        no-call). NOTE: this function *updates* is_ref in place.
    clusterable_sv_types: Set[Text] (default={Default.clusterable_sv_types})
        SV types that may be redundant (or needed for clustering with redundant SVs).
    cpx_ins_classes: Set[Text] (default={Default.cpx_ins_classes})
        CPX SV types that should produce an INS sink (modeled as a DEL)
    Yields
    -------
    bed_tuple: Tuple
        successive records for bed object, with fields: contig, start, end, variant_id, sv_type, is_cpx
    """
    for record in vcf_records:
        sv_type = record.info[Keys.svtype]
        if sv_type == Keys.cpx:
            if Keys.unresolved in record.filter:
                continue
            # If complex, all constituent intervals are in CPX_INTERVALS
            variant_id = record.id
            is_carrier[variant_id], is_ref[variant_id] = _get_carrier_status(record)
            for cpx_interval in record.info[Keys.cpx_intervals]:
                sv_type, region = cpx_interval.split('_', 1)
                contig, coords = region.split(':', 1)
                start, end = _fix_coords(*(int(c) for c in coords.split('-', 1)))
                yield contig, start, end, variant_id, sv_type, 1
            if record.info.get(Keys.cpx_type, None) in cpx_ins_classes:
                # If complex insertion, return insertion point as 1bp DEL
                sv_type = Keys.deletion
                contig = record.contig
                end = record.pos
                start = max(0, end - 1)
                yield contig, start, end, variant_id, sv_type, 1
        elif sv_type in clusterable_sv_types:
            start, end = _fix_coords(record.pos, record.stop)
            variant_id = record.id
            is_carrier[variant_id], is_ref[variant_id] = _get_carrier_status(record)
            yield record.contig, start, end, variant_id, sv_type, 0


def _vcf_records_to_bed_intervals(
        vcf_records: Iterable[pysam.VariantRecord],
        is_carrier: Dict[Text, numpy.ndarray],
        is_ref: Dict[Text, numpy.ndarray],
        clusterable_sv_types: Set[Text] = Default.clusterable_sv_types,
        cpx_ins_classes: Set[Text] = Default.cpx_ins_classes
) -> Iterator[Tuple]:
    f"""
    Iterate over input VCF, yielding records that may be redundant. Also gather is_carrier and is_ref mappings.
    This function mainly passes results from _unfiltered_vcf_records_to_bed_intervals, but potentially filters out
    unneeded SV intervals that originated in CPX events, and duplicates SVTYPE=CNV into one DUP and one DEL.
    Parameters
    ----------
    vcf_records: Iterable[VariantRecord]
        Iterable with pysam records from input VCF file.
    is_carrier: Dict[Text, numpy.ndarray]
        Dict from variant ID to boolean array that is True for samples called non-ref for this Variant, and False
        otherwise (including no-call). NOTE: this function *updates* is_carrier in place.
    is_ref: Dict[Text, numpy.ndarray]
        Dict from boolean array that is True for samples called ref for this Variant, and False otherwise (including
        no-call). NOTE: this function *updates* is_ref in place.
    clusterable_sv_types: Set[Text] (default={Default.clusterable_sv_types})
        SV types that may be redundant (or needed for clustering with redundant SVs).
    cpx_ins_classes: Set[Text] (default={Default.cpx_ins_classes})
        CPX SV types that should produce an INS sink (modeled as a DEL)
    Yields
    -------
    bed_tuple: Tuple
        successive records for bed object, with fields: contig, start, end, variant_id, sv_type, is_cpx
    """
    for contig, start, end, variant_id, sv_type, is_cpx in _unfiltered_vcf_records_to_bed_intervals(
        vcf_records, is_carrier, is_ref, clusterable_sv_types=clusterable_sv_types, cpx_ins_classes=cpx_ins_classes
    ):
        if sv_type in clusterable_sv_types:
            # store sv_type in interval.score, is_cpx in interval.strand
            if sv_type == Keys.cnv:  # ensure CNVs cluster with both insertions and deletions
                yield contig, start, end, variant_id, Keys.deletion, is_cpx
                yield contig, start, end, variant_id, Keys.dup, is_cpx
            else:  # yield this interval normally
                yield contig, start, end, variant_id, sv_type, is_cpx


def jaccard_index(is_carrier_a: numpy.ndarray, is_carrier_b: numpy.ndarray) -> float:
    """ return Jaccard index of carrier samples based on two boolean arrays of carrier status """
    return numpy.logical_and(is_carrier_a, is_carrier_b).sum() / numpy.logical_or(is_carrier_a, is_carrier_b).sum()


def _iter_pairwise_connections(
        clusterable_bedtool: pybedtools.BedTool,
        min_reciprocal_overlap: float,
        min_sample_overlap: float = 0,
        is_carrier: Mapping[Text, numpy.ndarray] = MappingProxyType({})
) -> Iterator[Tuple[Text, Text]]:
    """
    Iterate over pairs of variant intervals that meet minimum requirement for reciprocal overlap. Exclude self-overlaps.
    Optionally impose requirement of minimum Jaccard index for carrier samples.
    Parameters
    ----------
    clusterable_bedtool: BedTool
        bed object with intervals that may overlap each other
    min_reciprocal_overlap: float
        minimum reciprocal overlap for two intervals to be connected
    min_sample_overlap: float (default=0)
        minimum Jaccard index of carrier samples for two intervals to be connected
    is_carrier: Mapping[Text, numpy.ndarray]
        map from variant ID to carrier status (array boolean True/False for each sample)
    Yields
    -------
    variant_id_1, variant_id_2: Tuple[Text, Text]
        successive pairs of variant IDs that meet the overlap requiremnts
    """
    # Cluster intervals based on reciprocal overlap
    if len(clusterable_bedtool) == 0:
        return
    overlap_bedtool = clusterable_bedtool.intersect(clusterable_bedtool, f=min_reciprocal_overlap, r=True, wa=True,
                                                    wb=True, sorted=True, nonamecheck=True)
    num_1_fields = clusterable_bedtool.field_count()
    name_1_field = name_field
    sv_type_1_field = sv_type_field
    name_2_field = num_1_fields + name_field
    sv_type_2_field = num_1_fields + sv_type_field

    if min_sample_overlap > 0:
        for overlap in overlap_bedtool:
            fields = overlap.fields
            if fields[sv_type_1_field] != fields[sv_type_2_field]:
                continue  # only cluster same sv_type
            name_1 = fields[name_1_field]
            name_2 = fields[name_2_field]
            if name_1 != name_2 and jaccard_index(is_carrier[name_1], is_carrier[name_2]) >= min_sample_overlap:
                yield name_1, name_2
    else:
        for overlap in overlap_bedtool:
            fields = overlap.fields
            if fields[sv_type_1_field] != fields[sv_type_2_field]:
                continue  # only cluster same sv_type
            name_1 = fields[name_1_field]
            name_2 = fields[name_2_field]
            if name_1 != name_2:
                yield name_1, name_2


def _get_clusters(
        clusterable_bedtool: pybedtools.BedTool,
        min_reciprocal_overlap: float,
        min_sample_overlap: float = 0,
        is_carrier: Mapping[Text, numpy.ndarray] = MappingProxyType({})
) -> List[numpy.ndarray]:
    """
    Perform single-linkage clustering of variant intervals based on reciprocal overlap. Potentially impose a clustering
    requirement of high Jaccard index for carrier samples.
    Parameters
    ----------
    clusterable_bedtool: BedTool
        bed object with intervals that may cluster with each other
    min_reciprocal_overlap: float
        minimum reciprocal overlap for two intervals to be placed in a cluster
    min_sample_overlap: float (default=0)
        minimum Jaccard index of carrier samples for two intervals to be placed in a cluster
    is_carrier: Mapping[Text, numpy.ndarray]
        map from variant ID to carrier status (array boolean True/False for each sample)
    Returns
    -------
    clusters: List[numpy.ndarray]
        each element is an object numpy array of intervals that are in a cluster
    """
    # form map from variant IDs to unique indices for this clustering
    name_to_index = {name: index for index, name in enumerate({interval.name for interval in clusterable_bedtool})}
    num_vertices = len(name_to_index)
    sparse_connections = scipy.sparse.eye(num_vertices, dtype=numpy.uint8, format="lil")
    for name_1, name_2 in _iter_pairwise_connections(
        clusterable_bedtool, min_reciprocal_overlap=min_reciprocal_overlap, min_sample_overlap=min_sample_overlap,
        is_carrier=is_carrier
    ):
        sparse_connections[(name_to_index[name_1], name_to_index[name_2])] = 1

    # Cluster graph. Use "weak" connection because bedtools will list "A overlaps B" and "B overlaps A"
    num_clusters, cluster_labels = scipy.sparse.csgraph.connected_components(sparse_connections, connection="weak")

    # Build lists of clustered Intervals
    clusters = [[] for _ in range(num_clusters)]
    for interval in clusterable_bedtool:
        cluster_label = cluster_labels[name_to_index[interval.name]]
        clusters[cluster_label].append(interval)

    # convert lists to numpy object arrays for faster indexing
    def _to_numpy_array(_cluster: List[pybedtools.Interval]) -> numpy.ndarray:
        _cluster_array = numpy.empty((len(_cluster),), dtype=numpy.object)
        _cluster_array[:] = _cluster
        return _cluster_array

    return [_to_numpy_array(cluster) for cluster in clusters]


def _is_cpx(interval: pybedtools.Interval) -> int:
    """ returns 1 if this interval originated as a CPX interval, 0 otherwise """
    return int(interval.strand)  # is_cpx is stored in strand


def _is_not_cpx(interval: pybedtools.Interval) -> bool:
    """ returns 0 if this interval originated as a CPX interval, 1 otherwise """
    return int(interval.strand) == 0  # is_cpx is stored in strand


def _get_redundant_cluster_cnv_cpx_vids(
        cluster: numpy.ndarray,
        is_carrier: Mapping[Text, numpy.ndarray],
        cnv_cpx_sample_overlap: float
) -> Iterator[Text]:
    """
    Find CNVs that are redundant with CPX events
    for each sample that participates in this interval-cluster:
        join every variant ID that the sample participates in into a sample-cluster
        if the sample-cluster contains >= 1 CPX and >= 1 non-CPX:
            find logical-or carrier status over variant IDs in sample-cluster
            for every non-CPX variant ID in sample_cluster:
                if its Jaccard index is >= cnv_cpx_sample_overlap, it's redundant
    Parameters
    ----------
    cluster: numpy.ndarray
        numpy object array of pybedtools.Interval holding intervals that cluster together
    is_carrier: Mapping[Text, numpy.ndarray]
        Map from variant ID to boolean array that is True for samples called non-ref for this Variant, and False
        otherwise (including no-call).
    cnv_cpx_sample_overlap: float
        Minimum Jaccard index for variant interval to have with sample cluster in order for it to be redundant.
    Yields
    -------
    redundant_variant_id: str
        Successive redundant variant IDs
    """
    interval_is_cpx = numpy.fromiter((_is_cpx(interval) for interval in cluster), dtype=bool, count=len(cluster))

    def _is_valid_sample_cluster(indices_in_sample_cluster: numpy.ndarray) -> bool:
        # valid sample clusters have some CPX intervals and some non-CPX intervals
        return 0 < interval_is_cpx.take(indices_in_sample_cluster).sum() < len(indices_in_sample_cluster)

    if not _is_valid_sample_cluster(numpy.arange(len(cluster))):
        # no hope of finding valid sample-clusters if the whole thing won't work
        return

    # loop over unique combinations of intervals that are all non-ref for a single sample
    is_carrier_matrix = numpy.concatenate(
        [is_carrier[interval.name].reshape(1, -1) for interval in cluster], axis=0
    )

    for interval_in_potential_cluster in numpy.unique(is_carrier_matrix, axis=1).transpose():
        indices_in_potential_cluster = numpy.nonzero(interval_in_potential_cluster)[0]
        if not _is_valid_sample_cluster(indices_in_potential_cluster):
            continue  # not a valid cluster, skip it
        # can check jaccard index a little more quickly because each interval is in the cluster, so the intersection
        # is equal to the carrier status of the interval
        num_cluster_carriers = numpy.logical_or.reduce(
            is_carrier_matrix.take(indices_in_potential_cluster, axis=0), axis=0
        ).sum()
        for index in indices_in_potential_cluster:
            if not interval_is_cpx.take(index) and \
                    is_carrier_matrix.take(index, axis=0).sum() / num_cluster_carriers >= cnv_cpx_sample_overlap:
                yield cluster.take(index).name


def _find_cnv_cpx_redundancies(
        potentially_clusterable: pybedtools.BedTool,
        is_carrier: Mapping[Text, numpy.ndarray],
        min_cpx_reciprocal_overlap: float,
        cnv_cpx_reciprocal_overlap: float,
        cnv_cpx_sample_overlap: float
) -> Set[Text]:
    """
    Subset potentially clusterable intervals to those that meet required minimum overlap with a CPX event.
    Then find clusters, and remove redundant CNVs from those clusters.
    Parameters
    ----------
    potentially_clusterable: BedTool
        bed object with intervals that could potentially be used for clustering
    is_carrier: Mapping[Text, numpy.ndarray]
        Map from variant ID to boolean array that is True for samples called non-ref for this Variant, and False
        otherwise (including no-call).
    min_cpx_reciprocal_overlap: float
        Minimum reciprocal overlap with a CPX interval for a CNV interval to be clusterable.
    cnv_cpx_reciprocal_overlap: float
        Minimum reciprocal overlap between two intervals to be part of a cluster.
    cnv_cpx_sample_overlap: float
        Minimum Jaccard index for variant interval to have with sample cluster in order for it to be redundant.
    Returns
    -------
    vids_to_remove: Set[Text]
        Set of variant IDs that are redundant and should be removed from the output VCF.
    """
    # find all potentially clusterable intervals that meet required minimum overlap with CPX
    precluster_subset = potentially_clusterable.intersect(
        potentially_clusterable.filter(_is_cpx), u=True, f=min_cpx_reciprocal_overlap, r=True, sorted=True,
        nonamecheck=True
    )

    # find clusters of intervals with high reciprocal overlap, then check each cluster for redundant variant IDs
    return {
        variant_id
        for cluster in _get_clusters(precluster_subset, min_reciprocal_overlap=cnv_cpx_reciprocal_overlap)
        for variant_id in _get_redundant_cluster_cnv_cpx_vids(cluster, is_carrier,
                                                              cnv_cpx_sample_overlap=cnv_cpx_sample_overlap)
    }


def _update_cnv_cnv_redundances(
        vids_to_remove: Set[Text],
        potentially_clusterable: pybedtools.BedTool,
        is_carrier: Mapping[Text, numpy.ndarray],
        is_ref: Mapping[Text, numpy.ndarray],
        cnv_cnv_reciprocal_overlap: float,
        cnv_cnv_sample_overlap: float
):
    """
    Update vids_to_remove by finding CNVs that are redundant with other CNVs (as opposed to CPX)
    -Find CNVs with very high reciprocal overlap, and very high carrier sample Jaccard index
    -For each CNV that is connected to any other CNVs
        Add that CNV and all its connections to vids_to_remove
        Find the "best" CNV: the maximum choosing 1st by number of carriers, 2nd by number of called refs
        Add the best CNV to set of vids that will be put back in (no matter what, even if previously or subsequently
        "removed")
    -Update vids_to_remove by removing the "best" variant IDs

    Parameters
    ----------
    vids_to_remove: Set[Text]
        set of variant IDs that are redundant and should be removed. NOTE: this function updates this set in place.
    potentially_clusterable: BedTool
        bed object with intervals that could potentially be used for clustering
    is_carrier: Mapping[Text, numpy.ndarray]
        Map from variant ID to boolean array that is True for samples called non-ref for this Variant, and False
        otherwise (including no-call).
    is_ref: Mapping[Text, numpy.ndarray]
        Map from variant ID to boolean array that is True for samples called ref for this Variant, and False otherwise
        (including no-call).
    cnv_cnv_reciprocal_overlap: float
        minimum reciprocal overlap for two CNVs to be connected
    cnv_cnv_sample_overlap: float
        minimum carrier samples Jaccard index for two CNVs to be connected
    """
    # for each non-CPX interval, find all non-CPX intervals it has sufficient reciprocal overlap and sample overlap with
    variant_pairwise_connections = {}

    non_cpx_potentially_clusterable = potentially_clusterable.filter(_is_not_cpx).saveas()
    for name_1, name_2 in _iter_pairwise_connections(
            non_cpx_potentially_clusterable, min_reciprocal_overlap=cnv_cnv_reciprocal_overlap,
            min_sample_overlap=cnv_cnv_sample_overlap, is_carrier=is_carrier
    ):
        variant_pairwise_connections[name_1] = variant_pairwise_connections.get(name_1, (name_1,)) + (name_2,)

    vids_to_remove.update(variant_pairwise_connections.keys())  # set all the clustered variants to be removed

    # for each of these variant and its direct connections
    #    - choose one "best" variant to represent it, with priority given to most carriers, followed by most ref calls
    #    - keep the "best" variant (even if it's previously or subsequently "removed") and remove all others
    num_carrier = {variant_id: variant_is_carrier.sum() for variant_id, variant_is_carrier in is_carrier.items()}
    num_ref = {variant_id: variant_is_ref.sum() for variant_id, variant_is_ref in is_ref.items()}

    def _best_variant_id(variant_id: Text) -> (int, int, str):
        return num_carrier[variant_id], num_ref[variant_id], variant_id
    # then remove the best ones
    vids_to_remove.difference_update(
        max(variant_id_cluster, key=_best_variant_id) for variant_id_cluster in variant_pairwise_connections.values()
    )


def resolve_cpx_cnv_redundancies(
        input_vcf: Text,
        output_vcf: Text,
        min_cpx_reciprocal_overlap: float = Default.min_cpx_reciprocal_overlap,
        cnv_cpx_reciprocal_overlap: float = Default.cnv_cpx_reciprocal_overlap,
        cnv_cpx_sample_overlap: float = Default.cnv_cpx_sample_overlap,
        cnv_cnv_reciprocal_overlap: float = Default.cnv_cnv_reciprocal_overlap,
        cnv_cnv_sample_overlap: float = Default.cnv_cnv_sample_overlap,
        clusterable_sv_types: Set[Text] = Default.clusterable_sv_types,
        cpx_ins_classes: Set[Text] = Default.cpx_ins_classes,
        temp_dir: str = Default.temp_dir,
        num_threads: int = Default.num_threads
):
    f"""
    From input VCF, find redundant CNVs:
        CNVs that have sufficient reciprocal overlap and carrier sample Jaccard index with a CPX
        CNVs that have sufficient reciprocal overlap and carrier sampel Jaccard index with another CNV
    Write new VCF without redundant CNVs.
    Parameters
    ----------
    input_vcf: Text
        path to input vcf
    output_vcf: Text
        path to write output vcf
    min_cpx_reciprocal_overlap: float (default={Default.min_cpx_reciprocal_overlap,})
        Minimum reciprocal overlap with a CPX interval for a CNV interval to be clusterable.
    cnv_cpx_reciprocal_overlap: float (default={Default.cnv_cpx_reciprocal_overlap})
        Minimum reciprocal overlap between two intervals to be part of a cluster.
    cnv_cpx_sample_overlap: float (default={Default.cnv_cpx_sample_overlap})
        Minimum Jaccard index for variant interval to have with sample cluster in order for it to be redundant.
    cnv_cnv_reciprocal_overlap: float (default={Default.cnv_cnv_reciprocal_overlap})
        Minimum reciprocal overlap for two CNVs to be connected
    cnv_cnv_sample_overlap: float (default={Default.cnv_cnv_sample_overlap})
        Minimum carrier samples Jaccard index for two CNVs to be connected
    clusterable_sv_types: Set[Text] (default={Default.clusterable_sv_types})
        SV types that may be redundant (or needed for clustering with redundant SVs).
    cpx_ins_classes: Set[Text] (default={Default.cpx_ins_classes})
        CPX SV types that should produce an INS sink (modeled as a DEL)
    temp_dir: str (default={Default.temp_dir})
        Base folder to create new temp folder in.
    num_threads: int (default={Default.temp_dir})
        Number of threads to use for compression/decompression of VCF files.
    """
    temp_dir = os.path.abspath(os.path.expanduser(temp_dir))
    os.makedirs(temp_dir, exist_ok=True)
    pybedtools.set_tempdir(temp_dir)
    is_carrier, is_ref = {}, {}
    with pysam.VariantFile(input_vcf, 'r', threads=num_threads) as f_in:
        header = f_in.header
        potentially_clusterable = pybedtools.BedTool(
            _vcf_records_to_bed_intervals(f_in.fetch(), is_carrier, is_ref, clusterable_sv_types=clusterable_sv_types,
                                          cpx_ins_classes=cpx_ins_classes)
        ).saveas().sort()

    # get all the potentially clusterable intervals
    vids_to_remove = _find_cnv_cpx_redundancies(
        potentially_clusterable, is_carrier, min_cpx_reciprocal_overlap=min_cpx_reciprocal_overlap,
        cnv_cpx_reciprocal_overlap=cnv_cpx_reciprocal_overlap, cnv_cpx_sample_overlap=cnv_cpx_sample_overlap
    )
    _update_cnv_cnv_redundances(
        vids_to_remove, potentially_clusterable, is_carrier, is_ref,
        cnv_cnv_reciprocal_overlap=cnv_cnv_reciprocal_overlap, cnv_cnv_sample_overlap=cnv_cnv_sample_overlap
    )

    output_folder = os.path.dirname(os.path.abspath(os.path.expanduser(output_vcf)))
    os.makedirs(output_folder, exist_ok=True)
    with pysam.VariantFile(input_vcf, 'r', threads=num_threads) as f_in, \
            pysam.VariantFile(output_vcf, 'w', header=header, threads=num_threads) as f_out:
        for record in f_in.fetch():
            if record.id not in vids_to_remove:
                f_out.write(record)


def __parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Remove CNVs that are redundant with CPX variants, or each other",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("input_vcf", type=str, help="VCF with potentially redundant CNVs")
    parser.add_argument("output_vcf", type=str, help="VCF with redundant CNVs removed")
    parser.add_argument("--min-cpx-reciprocal-overlap", type=float, default=Default.min_cpx_reciprocal_overlap,
                        help="Minimum reciprocal overlap with a CPX for an interval to be possibly redundant")
    parser.add_argument("--cnv-cpx-reciprocal-overlap", type=float, default=Default.cnv_cpx_reciprocal_overlap,
                        help="Minimum reciprocal interval overlap for clustering CNV with CPX")
    parser.add_argument("--cnv-cpx-sample-overlap", type=float, default=Default.cnv_cpx_sample_overlap,
                        help="Minimum Jaccard index (intersection/union) of samples for clustering CNV with CPX")
    parser.add_argument("--cnv-cnv-reciprocal-overlap", type=float, default=Default.cnv_cnv_reciprocal_overlap,
                        help="Minimum reciprocal interval overlap for clustering CNV with other CNV")
    parser.add_argument("--cnv-cnv-sample-overlap", type=float, default=Default.cnv_cnv_sample_overlap,
                        help="Minimum Jaccard index (intersection/union) of samples for clustering CNV with other CNV")
    parser.add_argument("--temp-dir", "-t", type=str, default=Default.temp_dir, help="directory for temp files")
    parser.add_argument("--num-threads", type=int, default=Default.num_threads,
                        help="number of threads for compressing/decompressing bgzipped files")

    if len(argv) <= 1:
        parser.parse_args(["--help"])
        sys.exit(0)
    parsed_arguments = parser.parse_args(argv[1:])
    if parsed_arguments.input_vcf is None:
        raise ValueError("Must supply input-vcf")
    if parsed_arguments.output_vcf is None:
        raise ValueError("Must supply output-vcf")

    return parsed_arguments


def main(argv: Optional[List[Text]] = None):
    if argv is None:
        argv = sys.argv
    arguments = __parse_arguments(argv)
    resolve_cpx_cnv_redundancies(**vars(arguments))


if __name__ == "__main__":
    main()
