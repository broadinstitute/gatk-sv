#!/usr/bin/env python

import sys
import os
import argparse
import json
from typing import List, Text, Mapping, Iterable, Optional, Dict, Union, Collection
from types import MappingProxyType

import numpy
import pandas

from sv_utils import genomics_io, interval_overlaps

Numeric = Union[int, float, numpy.integer, numpy.floating]


class Default:
    genome_track_split_columns = MappingProxyType({"hg38-RepeatMasker": "repClass", "hg38-gaps": "field_type"})
    genome_track_properties_save_suffix = "_track_properties.tsv.gz"
    n_jobs = interval_overlaps.Default.n_jobs
    max_overlap_task_size = interval_overlaps.Default.max_overlap_task_size
    wanted_properties = interval_overlaps.Default.wanted_properties
    # when checking overlap, pseudo-size of variants with 0 / >0 reference length = scale_factor * SVLEN + expand_bp:
    expand_point_svs_bp = interval_overlaps.Default.expand_point_svs_bp
    point_sv_scale_factor = interval_overlaps.Default.point_sv_scale_factor
    expand_non_point_svs_bp = interval_overlaps.Default.expand_non_point_svs_bp
    non_point_sv_scale_factor = interval_overlaps.Default.non_point_sv_scale_factor


class Keys:
    id = genomics_io.Keys.id
    contig = genomics_io.Keys.contig
    begin = genomics_io.Keys.begin
    end = genomics_io.Keys.end
    overlap_support = "overlap_support"
    name = genomics_io.BedKeys.name


def load_genome_tracks(
        genome_track_files: Iterable[str],
        genome_track_split_columns: Mapping[str, str] = Default.genome_track_split_columns,
) -> Dict[str, pandas.DataFrame]:
    f"""
    Load genome tracks from bed files
    Args:
        genome_track_files: Iterable[str]
            Paths to bed files containing genome tracks.
                -Must have header with columns "contig, "begin", and "end" (although they may be named "chrom",
                 "chromStart", and "chromEnd").
                -May also have columns "other_contig", "other_begin", and "other_end" (they may be named "otherChrom",
                 "otherStart", and "otherEnd"), in which case both intervals will be checked for overlap, as well as
                 spanning.
        genome_track_split_columns: Mapping[str, str] (default={Default.genome_track_split_columns})
            If a track contains multiple sub-types specified in a column, that column can be specified here and the
            track will be split into multiple tracks.
    Returns:
        genome_tracks: Dict[str, pandas.DataFrame]
            Dictionary of track_name, track table pairs
    """
    tracks = {
        os.path.basename(track_file).split('.', 1)[0]: genomics_io.bed_to_pandas(track_file)
        for track_file in genome_track_files
    }
    for track_name in list(tracks.keys()):
        split_col = genome_track_split_columns.get(track_name, None)
        if split_col is not None and split_col in tracks[track_name].columns:
            track = tracks.pop(track_name)
            for sub_track_name, sub_track in track.groupby(
                    split_col, sort=False, as_index=False, group_keys=False
            ):
                tracks[sub_track_name] = sub_track
    return tracks


def get_genome_track_overlap_properties(
        intervals_df: pandas.DataFrame,
        genome_tracks: Union[Iterable[str], Mapping[str, pandas.DataFrame]],
        genome_track_split_columns: Mapping[str, str] = Default.genome_track_split_columns,
        expand_point_svs_bp: int = Default.expand_point_svs_bp,
        point_sv_scale_factor: float = Default.point_sv_scale_factor,
        expand_non_point_svs_bp: int = Default.expand_non_point_svs_bp,
        non_point_sv_scale_factor: float = Default.non_point_sv_scale_factor,
        n_jobs: Optional[int] = Default.n_jobs,
        required_worker_memory: Optional[Numeric] = None,
        required_master_memory: Optional[Numeric] = None,
        max_task_size: Numeric = Default.max_overlap_task_size,
) -> pandas.DataFrame:
    f"""
    Calculate information about genome tracks overlapping each interval.
    Args:
        intervals_df: pandas.DataFrame
            Table of GenomeIntervals (must have columns "contig", "begin", and "end")
        genome_tracks: Union[Iterable[str], Mapping[str, pandas.DataFrame]]
            Paths to bed files containing genome tracks.
                -Must have header with columns "contig, "begin", and "end"
                -May also have columns "other_contig", "other_begin", and "other_end", in which case both intervals will
                be checked for overlap, as well as spanning.
        genome_track_split_columns: Mapping[str, str] (default={Default.genome_track_split_columns})
            If a track contains multiple sub-types specified in a column, that column can be specified here and the
            track will be split into multiple tracks.
        expand_point_svs_bp: int (Default = {Default.expand_point_svs_bp})
            when checking overlap, expand point SVs by this many BP in each direction
        point_sv_scale_factor: float (Default = {Default.point_sv_scale_factor})
            when checking overlap, expand point SVs this * SVLEN in each direction
        expand_non_point_svs_bp: int (Default = {Default.expand_non_point_svs_bp})
            when checking overlap, expand non-point SVs by this many BP in each direction
        non_point_sv_scale_factor: float (Default = {Default.non_point_sv_scale_factor})
            when checking overlap, expand non-point SVs this * SVLEN in each direction
        n_jobs: int or None (Default = {Default.n_jobs})
            Number of parallel workers to use. If None or <= 0, then use number
            of available processors.
        required_worker_memory: number or None (Default=None)
            Memory required by each worker process. Used to restrict n_jobs.
            If None, use defaults to estimate based on size of data.
        required_master_memory: number or None (Default=None)
            Memory required by master process. Used to restrict n_jobs.
            If None, use defaults to estimate based on size of data.
        max_task_size: number, (Default={Default.max_overlap_task_size})
            Chop genome up into tasks no larger than this number, with task size
            being estimated as (# evidence in one chunk) * avg overlapper_coverage(chunk)
    Returns:
        genome_annotation_info: pandas.DataFrame
            Table of information about overlap between Evidence and genome annotations with rows corresponding to
            evidence. With columns:
                -gap_overlap_fraction: proportion of each piece of Evidence
                    covered by a gap in the genome. Note, if multiple gaps
                    overlap each other, the gap_overlap_fraction can sum to > 1.
                and, for k in single_read_mapability_kmers:
                    mapability_s[k]: proportion of each piece of evidence covered
                        by single-read mappable k-mers (i.e. that location is
                        uniquely mappable for a single k-mer)
    """
    if not isinstance(genome_tracks, Mapping):
        genome_tracks = load_genome_tracks(genome_track_files=genome_tracks,
                                           genome_track_split_columns=genome_track_split_columns)

    original_index = intervals_df.index.copy()
    simple_intervals = interval_overlaps.fix_variants(
        intervals_df,
        expand_point_svs_bp=expand_point_svs_bp, point_sv_scale_factor=point_sv_scale_factor,
        expand_non_point_svs_bp=expand_non_point_svs_bp, non_point_sv_scale_factor=non_point_sv_scale_factor
    )

    def _track_interval_reducer(
            _overlap_results: pandas.DataFrame, _simple_intervals: pandas.DataFrame
    ) -> Dict[str, float]:
        return {_col: max(_overlap_results[_col]) for _col in _overlap_results.columns}

    genome_annotation_info = pandas.concat(
        tuple(
            interval_overlaps.postprocess_bnd_complex(
                interval_overlaps.apply_interval_overlap_func(
                    _overlap_support,
                    simple_intervals,
                    _genome_track,
                    n_jobs=n_jobs,
                    description=f"Getting {_track_name} overlap",
                    property_names=f"{_track_name}",
                    required_worker_memory=required_worker_memory,
                    required_master_memory=required_master_memory,
                    max_task_size=max_task_size
                ),
                simple_intervals,
                original_index,
                simple_interval_reducer=_track_interval_reducer,
                reduce_single_intervals=False
            )
            for _track_name, _genome_track in genome_tracks.items()
        ),
        axis=1
    )
    # going to save this as a TSV, so change index name to "name"
    genome_annotation_info.index.name = Keys.name

    return genome_annotation_info


def _overlap_support(
        interval: interval_overlaps.Record,
        overlappers: pandas.DataFrame,
) -> Dict[str, float]:
    """
    Worker function for interval overlaps. Find proportion of interval that is overlapped by some member of overlappers.
    Args:
        interval: interval_overlaps.Record
            namedtuple with location of test interval
        overlappers: pandas.DataFrame
            DataFrame with locations of intervals that overlap "interval"
    Returns:
        overlap_support_dict: Dict[str, float]
            Mapping with one key ("overlap_support") and its value.
    """
    if len(overlappers) == 0:
        return {Keys.overlap_support: 0.0}

    o_end = overlappers[Keys.end].values
    o_begin = overlappers[Keys.begin].values

    max_o_end = numpy.maximum.accumulate(o_end)
    break_indices = numpy.flatnonzero(o_begin[1:] > max_o_end[:-1])
    if len(break_indices) == 0:
        support = numpy.minimum(max_o_end[-1], interval.end) - numpy.maximum(o_begin[0], interval.begin)
    else:
        ends = max_o_end.take(numpy.concatenate((break_indices, [-1])))
        begins = o_begin.take(numpy.concatenate(([0], break_indices + 1)))
        support = (numpy.minimum(ends, interval.end) - numpy.maximum(begins, interval.begin)).sum()
    support_proportion = support / (interval.end - interval.begin)
    if support_proportion < 0:
        raise ValueError(f"support_proportion={support_proportion} for \n{interval}\n   and\n{overlappers}")
    return {Keys.overlap_support: support_proportion}


def _load_subset_variant_ids(subset_variant_ids_file: Optional[str] = None) -> Optional[List[str]]:
    if subset_variant_ids_file is None:
        return None
    if not os.path.isfile(subset_variant_ids_file):
        raise ValueError(f"subset_variant_ids_file {subset_variant_ids_file} is not a path to a valid file")
    with open(subset_variant_ids_file, 'r') as f_in:
        return [variant_id.rstrip() for variant_id in f_in]


def _load_genome_track_split_columns(split_columns_json: Optional[str]) -> Mapping[str, str]:
    if split_columns_json is None:
        return Default.genome_track_split_columns
    if not os.path.isfile(split_columns_json):
        raise ValueError(f"genome track split-columns json ({split_columns_json}) is not a path to a valid file")
    with open(split_columns_json, 'r') as f_in:
        return json.load(f_in)


def get_genome_track_overlaps(
        vcf: str,
        genome_track_files: Iterable[str],
        output_file: Optional[str] = None,
        genome_track_split_columns: Mapping[str, str] = Default.genome_track_split_columns,
        subset_variant_ids_file: Optional[str] = None,
        wanted_properties: Collection[str] = Default.wanted_properties,
        expand_point_svs_bp: int = Default.expand_point_svs_bp,
        point_sv_scale_factor: float = Default.point_sv_scale_factor,
        expand_non_point_svs_bp: int = Default.expand_non_point_svs_bp,
        non_point_sv_scale_factor: float = Default.non_point_sv_scale_factor
) -> pandas.DataFrame:
    # load genome tracks
    genome_tracks = load_genome_tracks(genome_track_files=genome_track_files,
                                       genome_track_split_columns=genome_track_split_columns)
    # load vcf
    available_properties = set(genomics_io.get_vcf_properties(vcf))
    wanted_properties = [prop for prop in wanted_properties if prop in available_properties]
    variant_loc = genomics_io.vcf_to_pandas(
        vcf, wanted_properties=wanted_properties, drop_trivial_multi_index=True,
        missing_properties_action=genomics_io.ErrorAction.Warn
    )
    # if only a subset of variant IDs are wanted, perform the subset
    subset_variant_ids = _load_subset_variant_ids(subset_variant_ids_file)
    if subset_variant_ids is not None:
        variant_loc = variant_loc.loc[variant_loc.index.intersection(subset_variant_ids)]
    # compute genome track overlap properties
    genome_track_overlap_properties = get_genome_track_overlap_properties(
        variant_loc, genome_tracks=genome_tracks,
        expand_point_svs_bp=expand_point_svs_bp, point_sv_scale_factor=point_sv_scale_factor,
        expand_non_point_svs_bp=expand_non_point_svs_bp, non_point_sv_scale_factor=non_point_sv_scale_factor,
    )
    # output results
    if output_file is not None:
        genomics_io.pandas_to_tsv(output_file, genome_track_overlap_properties, write_index=True)
    return genome_track_overlap_properties


def __parse_arguments(argv: List[Text]) -> argparse.Namespace:
    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Get stats of genome track overlap for each variant in a VCF",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        prog=argv[0]
    )
    parser.add_argument("--vcf", '-v', type=str, required=True,
                        help="VCF with variants to find overlaps with genome tracks")
    parser.add_argument("--track", "-t", action="extend", nargs='+',
                        help="bed file with genome track information")
    parser.add_argument("--output", "-O", type=str, required=True,
                        help="File to output results to. Will be a gzipped tab-delimited file with header.")
    parser.add_argument("--track-split-json", type=str,
                        help=f"Split one or more genome tracks into multiple tracks using the values provided in this "
                             f"JSON file: a mapping from genome track name to a column name from that track. The track"
                             f"will be split into smaller tracks that all have the same value in that column. If file "
                             f"is NOT provided, use: {Default.genome_track_split_columns}")
    parser.add_argument("--subset_variant_ids_file", "-s", type=str, required=False,
                        help="File with one variant ID per line, only find overlap properties for variants in the VCF "
                             " if they also have their variant ID in this file.")
    parser.add_argument("--expand_point_svs_bp", type=int, default=Default.expand_point_svs_bp,
                        help="When checking overlap, expand SVs with size 0 on the reference by this many BP in each"
                             " direction")
    parser.add_argument("--point-sv-scale-factor", type=float, default=Default.point_sv_scale_factor,
                        help="When checking overlap, expand SVs with size 0 on the reference by scale-factor * SVLEN in"
                             " each direction")
    parser.add_argument("--expand_non_point_svs_bp", type=int, default=Default.expand_non_point_svs_bp,
                        help="When checking overlap, expand SVs with size > 0 on the reference by this many BP in each"
                             " direction")
    parser.add_argument("--non-point-sv-scale-factor", type=float, default=Default.non_point_sv_scale_factor,
                        help="When checking overlap, expand SVs with size > 0 on the reference by scale-factor * SVLEN"
                             " in each direction")
    parsed_arguments = parser.parse_args(argv[1:] if len(argv) > 1 else ["--help"])
    return parsed_arguments


def main(argv: Optional[List[Text]] = None) -> pandas.DataFrame:
    arguments = __parse_arguments(sys.argv if argv is None else argv)
    return get_genome_track_overlaps(
        vcf=arguments.vcf, genome_track_files=arguments.track, output_file=arguments.output,
        genome_track_split_columns=_load_genome_track_split_columns(arguments.track_split_json),
        subset_variant_ids_file=arguments.subset_variant_ids_file,
        expand_point_svs_bp=arguments.expand_point_svs_bp,
        point_sv_scale_factor=arguments.point_sv_scale_factor,
        expand_non_point_svs_bp=arguments.expand_non_point_svs_bp,
        non_point_sv_scale_factor=arguments.non_point_sv_scale_factor,
    )


if __name__ == "__main__":
    main()
