#!/usr/bin/env python
import sys
import datetime
import argparse
import json
import warnings
import textwrap
import numpy
import scipy
import scipy.stats
import pandas
from matplotlib import pyplot
from matplotlib.backends.backend_pdf import PdfPages

from sv_utils import common, plotting, genomics_io, pedigree_tools, get_truth_overlap
from sv_utils.get_truth_overlap \
    import SampleConfidentVariants, ConfidentVariants, PrecisionRecallCurve, SvTypeCutoffSelector

from typing import List, Dict, Text, Iterable, Optional, Mapping, Union, Tuple, Set, Collection, Iterator, TypeVar, \
    Sequence, Container


class Keys:
    property = genomics_io.Keys.property
    sample_id = genomics_io.Keys.sample_id
    family_id = "family_id"
    contig = genomics_io.Keys.contig
    is_autosome = "is_autosome"
    svtype = genomics_io.Keys.svtype
    svlen = genomics_io.Keys.svlen
    score = "score"
    gt = genomics_io.Keys.gt
    gq = genomics_io.Keys.gq
    allele_count = genomics_io.Keys.allele_count
    rd_cn = genomics_io.Keys.rd_cn
    rd_gq = genomics_io.Keys.rd_gq
    cn = genomics_io.Keys.cn
    cnq = genomics_io.Keys.cnq

    all_variant_types = "all"
    vcf_file = "vcf"
    scores_source = "scores_source"
    scores_file = "scores_file"
    metrics_file = "metrics_file"
    passing_score = "passing_score"
    score_property = "score_property"
    proband = "proband"
    father = "father"
    mother = "mother"
    filtered = "filtered"
    unfiltered = "unfiltered"
    max_f1 = "max f1"
    num_trios = "# trios"
    mendelian = "mendelian"
    de_novo = "de novo"
    other = "other"
    bnd_types = get_truth_overlap.Keys.bnd_types
    non_bnd_types = "non-" + bnd_types
    name = genomics_io.BedKeys.name
    allele_frequency = "allele_frequency"
    het_proportion = "het_proportion"
    # keys for which figures to make
    precision_recall = "precision-recall"
    inheritance = "inheritance"
    variants_per_sample = "variants-per-sample"
    violation_curve = "violation-curve"
    scores_histogram = "scores-histogram"
    hardy_weinberg = "hardy-weinberg"


_all_make_figures = frozenset({Keys.precision_recall, Keys.inheritance, Keys.variants_per_sample, Keys.violation_curve,
                               Keys.scores_histogram, Keys.hardy_weinberg, })


class Default:
    vcf_score_property = Keys.gq
    num_precision_recall_curve_rows = 3
    num_inheritance_stats_rows = 3
    num_scores_histogram_rows = 3
    title_font_size = 6
    label_font_size = 6
    legend_font_size = 6
    tick_font_size = 4
    tick_labelrotation = -30
    plot_max_f1 = False
    replot_unfiltered = False
    use_copy_number = genomics_io.Default.use_copy_number
    use_cn = genomics_io.Default.use_cn
    # calculate and plot only these keys in Mendelian violation curve plots
    plotted_mendelian_violation_trio_types = (Keys.mendelian,)
    # only make HWE plots for these SV types:
    plotted_hardy_weinberg_categories = frozenset({Keys.non_bnd_types})
    hardy_weinberg_confidence = 0.95
    min_num_called_gt_hardy_weinberg = 5
    histogram_width = 0.9
    make_figures = _all_make_figures
    sv_selector_size_ranges = get_truth_overlap.Default.sv_selector_size_ranges


def log(text: str):
    print(f"{datetime.datetime.now().isoformat(' ')}: {text}", file=sys.stderr, flush=True)


def load_benchmark_properties_from_vcf_properties(
        vcf: str,
        wanted_properties: Collection[str],
        wanted_sample_ids: Optional[Collection[str]] = None,
        use_copy_number: bool = Default.use_copy_number,
        use_cn: bool = Default.use_cn
) -> pandas.DataFrame:
    f"""
    Load requested properties from VCF. If {Keys.allele_count} and/or {Keys.gq} is among them, load properties needed to
    compute it ({Keys.gq} requires {Keys.allele_count}, because some algorithms set a GT to no call rather than updating
    the quality score). If {Keys.is_autosome} is among them, load {Keys.contig} and compute it.
    Return table with only the requested properties
    Args:
        vcf: str
            Path to VCF
        wanted_properties: Collection[str]
            List of desired properties
        wanted_sample_ids: Optional[Collection[str]] (default=None)
            Restrict to desired samples if provided. Otherwise include data from all samples
        use_copy_number: bool (default={Default.use_copy_number})
            If True, when getting {Keys.allele_count} and/or {Keys.gq}, use copy number info if GT is no call and copy
            number adds info
        use_cn: bool (default={Default.use_cn})
            If True, when using copy number info, examine {Keys.cn} in addition to {Keys.rd_cn}. This value has no
            effect if use_copy_number is false.
    Returns:
        variants:
            Table of properties with rows corresponding to variants and (potentially multi-index) columns corresponding
            to desired properties
    """
    derived_properties = frozenset({Keys.allele_count, Keys.is_autosome})
    # figure out if we need to compute allele counts and/or quality
    if use_copy_number and Keys.gq in wanted_properties:
        need_quality = True
    else:
        need_quality = False
    if Keys.allele_count in wanted_properties:
        need_allele_count = True
    else:
        need_allele_count = False
    # determine any extra properties we need to be sure to load
    extra_load_properties = []
    if need_allele_count or need_quality:
        if use_copy_number:
            extra_load_properties = [Keys.gt, Keys.rd_cn, Keys.cn] if use_cn else [Keys.gt, Keys.rd_cn]
            if need_quality:
                extra_load_properties += ([Keys.rd_gq, Keys.cnq] if use_cn else [Keys.rd_gq])
        else:
            extra_load_properties = [Keys.gt]
    if Keys.is_autosome in wanted_properties:
        extra_load_properties.append(Keys.contig)
    load_properties = [_property for _property in wanted_properties if _property not in derived_properties] + \
                      [_property for _property in extra_load_properties if _property not in wanted_properties]
    # load the wanted properties + any extras needed for computing allele_counts or quality
    variants = genomics_io.vcf_to_pandas(vcf, samples=wanted_sample_ids, wanted_properties=load_properties,
                                         drop_trivial_multi_index=False)
    # compute any desired extra properties
    if need_quality:
        allele_counts, quality = genomics_io.get_allele_counts_and_quality(
            variants, use_copy_number=use_copy_number, use_cn=use_cn
        )
        # where the genotype is no-call, set the quality to None (perfectly bad)
        ac_is_na = allele_counts.isna()
        if ac_is_na.hasna.any(axis=None):
            quality = genomics_io.make_nullable(quality).where(~allele_counts.isna(), pandas.NA)
        variants = genomics_io.assign_dataframe_column(variants, quality, Keys.gt)
        if need_allele_count:
            variants = genomics_io.assign_dataframe_column(variants, allele_counts, Keys.allele_count)
    elif need_allele_count:
        allele_counts = genomics_io.get_allele_counts(variants, use_copy_number=use_copy_number, use_cn=use_cn)
        variants = genomics_io.assign_dataframe_column(variants, allele_counts, Keys.allele_count)
    if Keys.is_autosome in wanted_properties:
        variants = genomics_io.assign_dataframe_column(
            variants, genomics_io.get_is_autosome(variants), Keys.is_autosome
        )
    # keep only the wanted properties, drop multi-index if it's not needed
    return genomics_io.drop_trivial_columns_multi_index(variants.loc[:, (slice(None), wanted_properties)])


class TruthData:
    """
    Class for holding truth data
        confident_variants: Mapping[str, SampleConfidentVariants]
            Mapping from sample ID to SampleConfidentVariants (class holding good_variant_ids and bad_variant_ids) that
            we are confident are good (definitely present) or bad (definitely not present) for that particular sample
        pedigree_file_info: PedigreeFileInfo
            Class with info on relatedness of samples
        overall_confident_variants: SampleConfidentVariants:
            SampleConfidentVariants object that holds variant IDs for variants that overall we are sure are good or bad.
    """
    __slots__ = ("confident_variants", "pedigree_file_info", "_overall_confident_variants")

    def __init__(self, overlap_results_file: str, pedigree_files: Optional[Collection[str]] = None):
        self.confident_variants = TruthData.load_confident_variants(overlap_results_file)
        self.pedigree_file_info = None if (pedigree_files is None or not pedigree_files) \
            else pedigree_tools.PedigreeFileInfo.load(pedigree_files)
        self._overall_confident_variants = None

    @staticmethod
    def load_confident_variants(overlap_results_file: str) -> ConfidentVariants:
        with open(overlap_results_file, 'r') as f_in:
            confident_variants = {
                sample_id: SampleConfidentVariants(**confident_variants)
                for sample_id, confident_variants in json.load(f_in).items()
            }

        return confident_variants

    @property
    def truth_sample_ids(self) -> Set[str]:
        trio_truth_sample_ids = {
            sample_id for pedigree_line in self.pedigree_file_info.pedigree_lines
            for sample_id in pedigree_line.trio_ids
            if not pedigree_line.any_unknown
        }
        return set(self.confident_variants.keys()).union(trio_truth_sample_ids)

    @property
    def all_truth_variants(self) -> Set[str]:
        return {
            v_id
            for confident_variants in self.confident_variants.values()
            for v_id in confident_variants.all_confident_variant_ids
        }

    @property
    def overall_confident_variants(self) -> SampleConfidentVariants:
        """
        Get list of variants that are confidently good (real variants in at least some sample) and confidently bad
        (never real). Good variants are confidently good in at least one sample. Bad variants are not confidently good
        in any sample, and are confidently bad in at least one sample. NOTE: it's okay for good variant IDs to be bad in
        some sample: we don't expect every sample to have a variant, even if it's real.
        Returns:
            overall_confident_variants: SampleConfidentVariants
                Class holding variant IDs that we are confident are good/bad overall (i.e. real vs not real variants).
        """
        if self._overall_confident_variants is None:
            good_variant_ids = set()
            bad_variant_ids = set()
            # form unions to find set of variant IDs that are good or bad in ANY sample
            for sample_id, sample_test_truth_overlaps in self.confident_variants.items():
                good_variant_ids.update(sample_test_truth_overlaps.good_variant_ids)
                bad_variant_ids.update(sample_test_truth_overlaps.bad_variant_ids)
            # use setdiff to restrict bad variant IDs to those that are bad in some sample, but never good
            bad_variant_ids = bad_variant_ids.difference(good_variant_ids)

            self._overall_confident_variants = get_truth_overlap.SampleConfidentVariants(
                good_variant_ids=good_variant_ids, bad_variant_ids=bad_variant_ids
            )
        return self._overall_confident_variants


def _ordered_intersection(
        keep_order_members: Iterable[str],
        check_members: Iterable[str]
):
    check_members = set(check_members)
    return [member for member in keep_order_members if member in check_members]


def load_variants_from_pickle(
        pickle_file: str,
        wanted_variant_ids: Optional[Iterable[str]] = None,
        wanted_sample_ids: Optional[Iterable[str]] = None,
        wanted_properties: Optional[Iterable[str]] = None
) -> pandas.DataFrame:
    variants = pandas.read_pickle(pickle_file)
    if wanted_variant_ids is None and wanted_sample_ids is None and wanted_properties is None:
        return variants
    variants_select = slice(None) if wanted_variant_ids is None \
        else _ordered_intersection(wanted_variant_ids, variants.index)
    # noinspection PyTypeChecker
    samples_select = slice(None) if wanted_sample_ids is None or Keys.sample_id not in variants.columns.names \
        else [None] + _ordered_intersection(wanted_sample_ids, variants.columns.get_level_values(level=Keys.sample_id))
    properties_select = slice(None) if wanted_properties is None or Keys.property not in variants.columns.names \
        else _ordered_intersection(wanted_properties, variants.columns.get_level_values(level=Keys.property))
    columns_select = samples_select if Keys.property not in variants.columns.names \
        else properties_select if Keys.sample_id not in variants.columns.names \
        else (samples_select, properties_select)
    return variants.loc[variants_select, columns_select]


class ScoresSource:
    __slots__ = ("score_files", "score_property", "passing_score")

    score_files: Optional[Collection[str]]
    score_property: Optional[str]
    passing_score: float

    def __init__(
            self,
            passing_score: float,
            score_property: Optional[str] = None,
            scores_file: Union[None, str, Iterable[str]] = None
    ):
        self.score_files = (scores_file,) if isinstance(scores_file, str) else scores_file
        self.score_property = score_property
        self.passing_score = passing_score

    @property
    def is_vcf_source(self) -> bool:
        _ends_vcf = [True] if self.score_files is None \
            else [score_file.endswith(".vcf") or score_file.endswith(".vcf.gz") for score_file in self.score_files]
        if any(_ends_vcf):
            if all(_ends_vcf):
                return True
            else:
                raise ValueError("Currently only ScoresSources that are all VCFs or all not VCFs are supported.")
        else:
            return False

    @property
    def is_pickle_source(self) -> bool:
        _ends_pickle = [False] if self.score_files is None \
            else [score_file.endswith(".pickle") for score_file in self.score_files]
        if any(_ends_pickle):
            if all(_ends_pickle):
                return True
            else:
                raise ValueError("Currently only ScoresSources that are all pickle files or all not pickle files are "
                                 "supported.")
        else:
            return False

    @property
    def has_empty_files(self) -> bool:
        return True if self.score_files is None or len(self.score_files) == 0 else False

    @property
    def failing_score(self) -> float:
        # If we know the passing score, the failing score is one less. Otherwise return -1, which is a failing score in
        # any scoring scheme I'm aware of.
        return -1 if self.passing_score is None else self.passing_score - 1

    def _load_scores_file(
            self,
            scores_file: str,
            wanted_variant_ids: Optional[Collection[str]] = None,
            wanted_sample_ids: Optional[Iterable[str]] = None,
            use_copy_number: bool = Default.use_copy_number,
            use_cn: bool = Default.use_cn
    ) -> pandas.DataFrame:
        try:
            if scores_file.endswith(".vcf") or scores_file.endswith(".vcf.gz"):
                # it's a VCF
                scores = load_benchmark_properties_from_vcf_properties(
                    vcf=scores_file, wanted_properties=(self.score_property,),
                    use_copy_number=use_copy_number, use_cn=use_cn
                )
                return scores if wanted_variant_ids is None \
                    else scores.loc[_ordered_intersection(wanted_variant_ids, scores.index)]
            elif ".pickle" in scores_file:
                # it's a pre-processed pickled table
                scores = load_variants_from_pickle(scores_file, wanted_variant_ids=wanted_variant_ids,
                                                   wanted_sample_ids=wanted_sample_ids,
                                                   wanted_properties=(self.score_property,))
                if scores.columns.nlevels > 1:
                    scores.columns = scores.columns.droplevel(Keys.property)
                return scores
            else:
                # it's a tsv
                if wanted_sample_ids is not None:
                    raise ValueError("Non-VCF scores files cannot have per-genotype scores")
                _score_property = Keys.score if self.score_property is None else self.score_property
                scores = genomics_io.tsv_to_pandas(
                    scores_file, columns=(_score_property,)
                ).set_index(Keys.name)
                return scores if wanted_variant_ids is None \
                    else scores.loc[_ordered_intersection(wanted_variant_ids, scores.index)]
        except Exception as exception:
            common.add_exception_context(exception, f"loading scores file {scores_file}")
            raise

    def load_scores(
            self,
            vcf_data: pandas.DataFrame,
            wanted_variant_ids: Optional[Collection[str]] = None,
            wanted_sample_ids: Optional[Iterable[str]] = None,
    ) -> Union[pandas.Series, pandas.DataFrame]:
        if self.has_empty_files:
            return genomics_io.get_format_property(vcf_data, self.score_property)
        else:
            scores = pandas.concat(
                (
                    self._load_scores_file(scores_file=scores_file, wanted_variant_ids=wanted_variant_ids,
                                           wanted_sample_ids=wanted_sample_ids)
                    for scores_file in self.score_files
                ),
                axis=0
            )
            return scores if wanted_variant_ids is None \
                else scores.loc[_ordered_intersection(wanted_variant_ids, scores.index)]


def get_dtype_info(dt: type) -> Union[numpy.iinfo, numpy.finfo]:
    numpy_dt = dt.numpy_dtype if hasattr(dt, "numpy_dtype") else dt
    return numpy.iinfo(numpy_dt) if numpy.issubdtype(numpy_dt, numpy.integer) else numpy.finfo(numpy_dt)


class ScoresDataSet:
    """
    Class to hold data set with scores, metrics, allele counts, and inheritance stats for validation/benchmarking
    Can be initialized in unloaded state (only files, label, and default passing score specified) or loaded state
    (all data necessary for benchmarking included). Typically initialized in unloaded state and then loaded, using
    static from_json() function.
    """
    __slots__ = ("vcfs", "scores_sources", "label", "scores", "allele_counts", "inheritance_stats", "metrics",
                 "category", "_num_called_genotypes")

    def __init__(
            self,
            vcfs: Union[str, Iterable[str]],
            scores_sources: Mapping[str, ScoresSource],
            label: str,
            scores: Optional[pandas.DataFrame] = None,
            allele_counts: Optional[pandas.DataFrame] = None,
            inheritance_stats: Optional[pandas.DataFrame] = None,
            metrics: Optional[pandas.DataFrame] = None,
            category: str = Keys.all_variant_types,
            num_called_genotypes: Optional[pandas.Series] = None
    ):
        """

        Args:
            vcfs: Unions[str, Iterable[str]]
                One or more vcfs that are interpreted as all being part of the same data set.
            scores_sources: Mapping[str, ScoresSource]
            label:
            scores:
            allele_counts:
            inheritance_stats:
            metrics:
            category:
        """
        self.vcfs = (vcfs,) if isinstance(vcfs, str) else tuple(vcfs)
        self.scores_sources = scores_sources
        self.label = label
        self.scores = scores
        self.allele_counts = allele_counts
        self.inheritance_stats = inheritance_stats
        self.metrics = metrics
        self.category = category
        self._num_called_genotypes = num_called_genotypes

    @staticmethod
    def from_json(
            json_file: str,
            truth_data: Optional[TruthData] = None,
    ) -> List["ScoresDataSet"]:
        with open(json_file, 'rb') as f_in:
            all_data_sets_info = json.load(f_in)
        unloaded_data_sets = [
            ScoresDataSet(
                vcfs=data_set_info[Keys.vcf_file],
                scores_sources={source_label: ScoresSource(**source_info)
                                for source_label, source_info in data_set_info[Keys.scores_source].items()},
                label=label
            )
            for label, data_set_info in all_data_sets_info.items()
        ]
        if truth_data is None:
            return unloaded_data_sets
        else:
            return [
                loaded_data_set
                for unloaded_data_set in unloaded_data_sets
                for loaded_data_set in unloaded_data_set.load_data()
            ]

    @staticmethod
    def _load_vcf_data(
            vcf: str,
            wanted_properties: Tuple[str, ...],
            wanted_variant_ids: Optional[Iterable[str]] = None,
            wanted_sample_ids: Optional[Iterable[str]] = None
    ) -> pandas.DataFrame:
        try:
            if vcf.endswith(".vcf") or vcf.endswith(".vcf.gz"):
                vcf_data = load_benchmark_properties_from_vcf_properties(vcf, wanted_properties=wanted_properties,
                                                                         wanted_sample_ids=wanted_sample_ids)
                return vcf_data if wanted_variant_ids is None \
                    else vcf_data.loc[_ordered_intersection(wanted_variant_ids, vcf_data.index)]
            elif ".pickle" in vcf:
                return load_variants_from_pickle(vcf, wanted_variant_ids=wanted_variant_ids,
                                                 wanted_sample_ids=wanted_sample_ids,
                                                 wanted_properties=wanted_properties)
            else:
                raise ValueError("vcf data must be a VCF file ('.vcf' or '.vcf.gz') or a .pickle file")
        except Exception as exception:
            common.add_exception_context(exception, f"loading vcf data {vcf}")
            raise

    def load_data(
            self,
            wanted_variant_ids: Optional[Iterable[str]] = None,
            wanted_sample_ids: Optional[Iterable[str]] = None
    ) -> Iterator["ScoresDataSet"]:
        if self.scores is None or self.allele_counts is None:
            # not already loaded, or need to load genotypes
            wanted_properties = (Keys.svtype, Keys.svlen, Keys.allele_count, Keys.is_autosome)
            for scores_source in self.scores_sources.values():
                if scores_source.has_empty_files:
                    wanted_properties += (scores_source.score_property,)
            vcf_data = pandas.concat(
                (
                    ScoresDataSet._load_vcf_data(
                        vcf, wanted_properties=wanted_properties, wanted_variant_ids=wanted_variant_ids,
                        wanted_sample_ids=wanted_sample_ids
                    )
                    for vcf in self.vcfs
                ),
                axis=0
            )
            # extract allele_counts and metrics, and drop now-unneeded columns from vcf_data
            vcf_data.drop(Keys.gt, axis=1, level=Keys.property, inplace=True)
            allele_counts = genomics_io.get_format_property(vcf_data, Keys.allele_count)
            vcf_data.drop(Keys.allele_count, axis=1, level=Keys.property, inplace=True)
            metrics = genomics_io.get_info_property(vcf_data, [Keys.svtype, Keys.svlen, Keys.is_autosome])
            vcf_data.drop([Keys.svtype, Keys.svlen, Keys.is_autosome], axis=1, level=Keys.property, inplace=True)

            # yield one ScoresDataSet for each ScoresSource
            for source_label, scores_source in self.scores_sources.items():
                # noinspection PyTypeChecker
                try:
                    scores = scores_source.load_scores(vcf_data, wanted_variant_ids=vcf_data.index,
                                                       wanted_sample_ids=wanted_sample_ids)
                except Exception as exception:
                    common.add_exception_context(exception, f"Loading {source_label}")
                    raise

                if not scores.index.equals(allele_counts.index):
                    # this scores_source has different variants from the vcfs
                    # 1) could be a superset, caused by passing in a larger VCF for comparison
                    #    - in which case remove "extra" variants
                    # 2) most likely be a subset, caused by rejecting this missing variants and filtering them out
                    #    - in which case, add failing scores for these variants
                    extra_variants = scores.index.difference(allele_counts.index)
                    if extra_variants.any():
                        scores.drop(labels=extra_variants, inplace=True)
                    missing_variants = allele_counts.index.difference(scores.index)
                    if missing_variants.any():
                        scores = pandas.concat(
                            (
                                genomics_io.make_nullable(scores),
                                pandas.Series(index=missing_variants) if isinstance(scores, pandas.Series) else
                                pandas.DataFrame(index=missing_variants)
                            )
                        )
                    # now scores can be arranged with the same index as allele_counts
                    scores = scores.loc[allele_counts.index]

                yield ScoresDataSet(
                    vcfs=self.vcfs,
                    scores_sources={source_label: scores_source},
                    label=f"{self.label} {source_label}" if source_label else self.label,
                    scores=scores,
                    allele_counts=allele_counts,
                    inheritance_stats=None,
                    metrics=metrics,
                    category=self.category
                )
        else:
            yield self

    def filter_variants(self, variant_is_wanted: pandas.Series, new_category: Optional[str] = None) -> "ScoresDataSet":
        return ScoresDataSet(
            vcfs=self.vcfs, scores_sources=self.scores_sources, label=self.label,
            scores=self.scores.loc[variant_is_wanted],
            allele_counts=None if self.allele_counts is None else self.allele_counts[variant_is_wanted],
            inheritance_stats=None if self.inheritance_stats is None else self.inheritance_stats[variant_is_wanted],
            metrics=self.metrics.loc[variant_is_wanted],
            category=self.category if new_category is None else new_category,
            num_called_genotypes=None if self._num_called_genotypes is None
            else self._num_called_genotypes.lloc[variant_is_wanted]
        )

    def select_category(
            self,
            category: str,
            sv_type_cutoff_selector: SvTypeCutoffSelector
    ) -> "ScoresDataSet":
        return self.filter_variants(
            variant_is_wanted=sv_type_cutoff_selector(self.metrics), new_category=category
        )

    def has_category_data(self, sv_type_cutoff_selector: SvTypeCutoffSelector) -> bool:
        return sv_type_cutoff_selector(self.metrics).any()

    def to_autosome(
            self,
            min_num_called_gt_hardy_weinberg: int = Default.min_num_called_gt_hardy_weinberg
    ) -> "ScoresDataSet":
        wanted = self.metrics[Keys.is_autosome] & \
            (self.num_ref_genotypes >= min_num_called_gt_hardy_weinberg) & \
            (self.num_non_ref_genotypes >= min_num_called_gt_hardy_weinberg) \
            if min_num_called_gt_hardy_weinberg > 0 else self.metrics[Keys.is_autosome]
        return self.filter_variants(variant_is_wanted=wanted)

    @property
    def num_samples(self) -> int:
        return self.allele_counts.shape[1]

    @property
    def is_variant_filter(self) -> bool:
        return isinstance(self.scores, pandas.Series)

    @property
    def is_genotype_filter(self) -> bool:
        return isinstance(self.scores, pandas.DataFrame)

    @property
    def passing_score(self) -> float:
        # assuming data is loaded, there should be only one ScoresSource. Get its passing score
        return next(iter(self.scores_sources.values())).passing_score

    @property
    def variants_per_sample(self) -> pandas.Series:
        # noinspection PyTypeChecker,PyUnresolvedReferences
        return (self.allele_counts > 0).sum(axis=0, skipna=True)

    @property
    def num_ref_genotypes(self) -> pandas.Series:
        return (self.allele_counts == 0).sum(axis=1, skipna=True)

    @property
    def num_het_genotypes(self) -> pandas.Series:
        return (self.allele_counts == 1).sum(axis=1, skipna=True)

    @property
    def num_homvar_genotypes(self) -> pandas.Series:
        return (self.allele_counts == 2).sum(axis=1, skipna=True)

    @property
    def num_non_ref_genotypes(self) -> pandas.Series:
        return (self.allele_counts > 0).sum(axis=1, skipna=True)

    @property
    def num_called_genotypes(self) -> pandas.Series:
        if self._num_called_genotypes is None:
            self._num_called_genotypes = self.allele_counts.shape[1] - (self.allele_counts.isna()).sum(axis=1)
        return self._num_called_genotypes

    @property
    def num_variant_alleles(self) -> pandas.Series:
        return self.allele_counts.sum(axis=1, skipna=True)

    @property
    def allele_frequency(self) -> pandas.Series:
        return self.num_variant_alleles / (2 * self.num_called_genotypes)

    @property
    def het_proportion(self) -> pandas.Series:
        return self.num_het_genotypes / self.num_called_genotypes

    def apply_hard_filter(self, passing_score: Optional[float] = None) -> "ScoresDataSet":
        if passing_score is None:
            # get default passing score from ScoresSource, should be only one since data is loaded
            passing_score = self.passing_score
        # noinspection PyTypeChecker
        if self.is_variant_filter:
            # noinspection PyTypeChecker
            return self.filter_variants(variant_is_wanted=(self.scores >= passing_score).fillna(False).astype(bool))
        else:  # this is a genotype filter
            # noinspection PyTypeChecker
            # keep allele counts where score >= passing score, otherwise set to NA
            # NOTE: it is possible for some VCFs to have missing (i.e. NULL) scores, in which case we will filter
            # those alleles as well.
            allele_counts = self.allele_counts.where(
                (self.scores >= passing_score).fillna(False).astype(bool),
                pandas.NA
            )
            return ScoresDataSet(
                vcfs=self.vcfs, scores_sources=self.scores_sources, label=self.label, scores=self.scores,
                allele_counts=allele_counts, inheritance_stats=None, metrics=self.metrics, category=self.category,
                num_called_genotypes=None
            )

    @staticmethod
    def _dataframe_to_trio_3d_numpy(
            truth_data: "TruthData",
            samples_df: pandas.DataFrame,
            numpy_dtype: Union[numpy.dtype, type]
    ) -> (numpy.ndarray, List[str]):
        """
        briefly use 3D numpy array instead of pandas DataFrame, because we want to use vectorized reduction along the
        3rd dimension (individual trios) and pandas provides no way to do this. Switching to numpy gives a big
        speed-up, and if needed we can move the final results into a pandas DataFrame basically for free
        Args:
            truth_data: TruthData

            samples_df:
            numpy_dtype:

        Returns:

        """
        sample_ids = set(samples_df.columns)
        num_variants = samples_df.shape[0]
        trio_probands = [
            pedigree_line.proband_id for pedigree_line in truth_data.pedigree_file_info.pedigree_lines
            if not pedigree_line.any_unknown and all(sample_id in sample_ids for sample_id in pedigree_line.trio_ids)
        ]
        num_trios = len(trio_probands)
        numpy_trios_tensor = numpy.empty((num_variants, num_trios, 3), dtype=numpy_dtype)
        trio_ind = 0
        null_value = get_dtype_info(numpy_dtype).min  # set null value to the min allowed value for this dtype
        for pedigree_line in truth_data.pedigree_file_info.pedigree_lines:
            if pedigree_line.any_unknown or not all(sample_id in sample_ids for sample_id in pedigree_line.trio_ids):
                continue
            numpy_trios_tensor[:, trio_ind, 0] = numpy.where(
                samples_df[pedigree_line.proband_id].isna(), null_value, samples_df[pedigree_line.proband_id]
            )
            numpy_trios_tensor[:, trio_ind, 1] = numpy.where(
                samples_df[pedigree_line.father_id].isna(), null_value, samples_df[pedigree_line.father_id]
            )
            numpy_trios_tensor[:, trio_ind, 2] = numpy.where(
                samples_df[pedigree_line.mother_id].isna(), null_value, samples_df[pedigree_line.mother_id]
            )
            trio_ind += 1
        return numpy_trios_tensor, trio_probands

    def _get_trio_categories(
            self,
            truth_data: "TruthData"
    ) -> Dict[str, pandas.DataFrame]:
        f"""
        return boolean pandas.DataFrame (with variant_ids for index and trio_proband_ids as columns) that indicates
        whether each trio is
            {Keys.num_trios}: "testable", aka every member is called
            {Keys.de_novo}: testable and inheritance pattern is de-novo
            {Keys.mendelian}: testable and inheritance pattern is mendelian
            {Keys.other}: testable and inheritance pattern is non-mendelian and not de-novo
        """
        # briefly use 3D numpy array instead of pandas DataFrame, because we want to use vectorized reduction along the
        # 3rd dimension (individual trios) and pandas provides no way to do this. Switching to numpy gives a big
        # speed-up, and we can move the final results into a pandas DataFrame basically for free
        # form num_variants x num_trios x 3 tensor of allele counts, reserving -1 for no-calls
        trio_acs, trio_probands = ScoresDataSet._dataframe_to_trio_3d_numpy(
            truth_data=truth_data, samples_df=self.allele_counts, numpy_dtype=numpy.int8
        )

        # testable if all members of a trio are defined, and at least one is called
        testable = numpy.logical_and(
            numpy.logical_and.reduce(trio_acs >= 0, axis=2),
            numpy.logical_or.reduce(trio_acs > 0, axis=2)
        )
        # denovo if proband is called non-ref and parents are both called ref
        de_novo = numpy.logical_and(
            trio_acs[:, :, 0] > 0,
            numpy.logical_and(trio_acs[:, :, 1] == 0, trio_acs[:, :, 2] == 0)
        )
        # mendelian if proband has neither too low nor too high allele count based on parents
        mendelian = numpy.logical_and(
            testable,
            numpy.logical_and(
                trio_acs[:, :, 1] // 2 + trio_acs[:, :, 2] // 2 <= trio_acs[:, :, 0],
                trio_acs[:, :, 0] <= (trio_acs[:, :, 1] > 0) + (trio_acs[:, :, 2] > 0)
            )
        )
        other_non_mendelian = numpy.logical_and(testable, numpy.logical_not(numpy.logical_or(de_novo, mendelian)))
        return {
            Keys.num_trios: pandas.DataFrame(testable, columns=trio_probands, index=self.allele_counts.index),
            Keys.de_novo: pandas.DataFrame(de_novo, columns=trio_probands, index=self.allele_counts.index),
            Keys.mendelian: pandas.DataFrame(mendelian, columns=trio_probands, index=self.allele_counts.index),
            Keys.other: pandas.DataFrame(other_non_mendelian, columns=trio_probands, index=self.allele_counts.index)
        }

    def get_inheritance_stats(
        self,
        truth_data: "TruthData"
    ) -> "ScoresDataSet":
        if truth_data.pedigree_file_info is None:
            return self  # can't calculate anything

        # sum across trios and put into DataFrame
        self.inheritance_stats = pandas.DataFrame(
            {trio_category: trio_is_category.sum(axis=1)
             for trio_category, trio_is_category in self._get_trio_categories(truth_data).items()}
        )
        return self

    def get_mendelian_violation_curve(
            self,
            truth_data: "TruthData",
            plotted_mendelian_violation_keys: Optional[Iterable[str]] = Default.plotted_mendelian_violation_trio_types
    ) -> Optional["MendelianViolationCurve"]:
        f"""
        Given supplied truth data compute trio_category_scores and from them a MendelianViolationCurve.
        trio_category_scores are a dictionary with
            keys: str
                trio category ({Keys.num_trios}, {Keys.de_novo}, {Keys.mendelian}, {Keys.other}) as described in
                    ScoresDataSet._get_trio_categories
            values: pandas.Serise
                Series with index = minimum score in the trio (hence the trio will be eliminated if the threshold goes
                above that score), and values = counts of trios in this trio_category with this minimum score
        Args:
            truth_data: TruthData
                Class containing truth data. Relevant info is the PedigreeFileInfo.
            plotted_mendelian_violation_keys:  Optional[Iterable[str]]
                (default={Default.plotted_mendelian_violation_trio_types})
                If supplied, restrict plotting (and computing) curves to only the relevant keys.
        Returns:
            mendelian_violation_curve: MendelianViolationCurve
                Object with enough information to plot the wanted mendelian violation rates vs quality score threshold
        """
        if truth_data.pedigree_file_info is None:
            return None
        column_infos = [get_dtype_info(dt) for dt in set(self.scores.dtypes)]
        min_value = min(column_info.min for column_info in column_infos)
        # for speed purposes, we're not using masked array, just setting nulled values to a low value and using that as
        # null. So set null to be the min allowed value, and insist that it be at least as small as -1
        min_type = genomics_io.IntPropertyCollator.range_to_min_int_dtype(
            min_value=min(min_value, -1), max_value=max(column_info.max for column_info in column_infos), nullable=False
        )
        # noinspection PyArgumentList
        min_scores = ScoresDataSet._dataframe_to_trio_3d_numpy(
            truth_data=truth_data, samples_df=self.scores, numpy_dtype=min_type
        )[0].min(axis=2)

        def _numpy_values_to_value_counts(_values: numpy.ndarray) -> pandas.Series:
            """ Get unique set of all non-null score values, and count of number of times each occurs """
            # first get rid of null scores (set to minimum allowed value by _dataframe_to_trio_3d_numpy)
            _values = _values.compress(_values != numpy.iinfo(_values.dtype).min)
            if min_value >= 0:
                # if original scores data was strictly non-negative, we can use this faster routine
                _value_counts = pandas.Series(numpy.bincount(_values))
                return _value_counts.loc[_value_counts > 0]
            else:
                # use general purpose routine
                return pandas.Series(_values).value_counts().sort_index()

        return MendelianViolationCurve.from_trio_category_scores(
            {
                category: _numpy_values_to_value_counts(
                    min_scores.ravel().compress(trios_are_category.values.ravel())
                )
                for category, trios_are_category in self._get_trio_categories(truth_data).items()
            },
            plotted_mendelian_violation_keys=plotted_mendelian_violation_keys
        )

    @staticmethod
    def get_good_bad_thresholds(
            confident_variants: SampleConfidentVariants, quality_scores: pandas.DataFrame
    ) -> (numpy.ndarray, numpy.ndarray):
        """ Return good quality_scores and bad_quality scores, by getting quality scores of confident variants """
        # noinspection PyTypeChecker
        return (
            quality_scores.loc[
                _ordered_intersection(quality_scores.index, confident_variants.good_variant_ids)
            ].astype("float").values,
            quality_scores.loc[
                _ordered_intersection(quality_scores.index, confident_variants.bad_variant_ids)
            ].astype("float").values
        )

    def get_precision_recall_curve(self, truth_data: TruthData) -> PrecisionRecallCurve:
        if self.is_variant_filter:
            return PrecisionRecallCurve.from_good_bad_thresholds(
                *ScoresDataSet.get_good_bad_thresholds(truth_data.overall_confident_variants, self.scores)
            )
        else:  # genotype filter
            # noinspection PyTypeChecker
            return PrecisionRecallCurve.from_good_bad_thresholds(
                *(
                    numpy.concatenate(sample_thresholds)
                    for sample_thresholds in zip(
                        *(
                            ScoresDataSet.get_good_bad_thresholds(truth_data.confident_variants[sample_id],
                                                                  sample_scores)
                            for sample_id, sample_scores in self.scores.items()
                            if sample_id in truth_data.confident_variants
                        )
                    )
                )
            )

    @property
    def all_sv_types(self) -> Set[str]:
        return set(self.metrics[Keys.svtype].unique())


class MendelianViolationCurve:
    __slots__ = ("dataframe",)
    num_trios_key = Keys.num_trios
    mendelian_key = Keys.mendelian
    de_novo_key = Keys.de_novo
    other_bad_key = Keys.other

    def __init__(self, dataframe: pandas.DataFrame):
        self.dataframe = dataframe

    @property
    def thresholds(self) -> numpy.ndarray:
        """ return quality scores where the number of passing trios change """
        return self.dataframe.index.to_numpy()

    @property
    def proportion_testable(self) -> numpy.ndarray:
        """ return number of trios that are testable (indexed by quality threshold) as a proportion of max number """
        testable = self.dataframe[MendelianViolationCurve.num_trios_key].values
        return testable / testable.max()

    @property
    def proportion_mendelian(self) -> numpy.ndarray:
        """ return proportion of trios that are mendelian (indexed by quality threshold) """
        return self.dataframe[MendelianViolationCurve.mendelian_key].values \
            / self.dataframe[MendelianViolationCurve.num_trios_key].values

    @property
    def proportion_de_novo(self) -> numpy.ndarray:
        """ return proportion of trios that are de-novo (indexed by quality threshold) """
        return self.dataframe[MendelianViolationCurve.de_novo_key].values \
            / self.dataframe[MendelianViolationCurve.num_trios_key].values

    @property
    def proportion_other_bad(self) -> numpy.ndarray:
        """ return proportion of trios that are non-mendelian but not de-novo (indexed by quality threshold) """
        return self.dataframe[MendelianViolationCurve.other_bad_key].values \
            / self.dataframe[MendelianViolationCurve.num_trios_key].values

    @staticmethod
    def _num_passing_lookup(score_counts: pandas.Series) -> pandas.Series:
        """
        Given number of trios indexed by lowest quality score, compute cumulative: number of trios with score >= a given
        quality score.
        Args:
            score_counts: pandas.Series
                Number of trios that are in some category (mendelian, de-novo, etc) indexed by the lowest quality score
                in that trio
        Returns:
            cumulative_score_counts: pandas.Series
                number of trios with score >= a given quality score.
        """
        # get reversed cumulative sum of counts, i.e. the number of trios that will have every member pass at each
        # threshold
        score_counts = score_counts[::-1].sort_index(ascending=False)
        try:
            if score_counts.index[0] != numpy.inf:
                # add a point at infinity with no counts
                score_counts = score_counts.reindex(score_counts.index.insert(0, numpy.inf), fill_value=0)
        except IndexError:  # score_counts was empty
            score_counts[numpy.inf] = 0
        # do cumsum then flip to ascending order of scores
        return score_counts.cumsum().sort_index(ascending=True)

    @staticmethod
    def _num_passing_at_thresholds(score_counts: pandas.Series, thresholds: numpy.ndarray) -> numpy.ndarray:
        """
        Compute the number of trios that pass at each provided quality threshold
        Args:
            score_counts: pandas.Series
                Number of trios that are in some category (mendelian, de-novo, etc) indexed by the lowest quality score
                in that trio
            thresholds: numpy.ndarray
                Quality scores to use a thresholds for a potential filter
        Returns:
            num_passing: numpy.ndarray
                Number of trios that pass at each point of the supplied trios
        """
        num_passing = MendelianViolationCurve._num_passing_lookup(score_counts)
        return num_passing.values.take(numpy.searchsorted(num_passing.index.values, thresholds))

    @staticmethod
    def from_trio_category_scores(
            trio_category_score_counts: Mapping[str, pandas.Series],
            plotted_mendelian_violation_keys: Optional[Iterable[str]] = Default.plotted_mendelian_violation_trio_types
    ) -> Optional["MendelianViolationCurve"]:
        f"""
        {Keys.num_trios}: "testable", aka every member is called
            {Keys.de_novo}: testable and inheritance pattern is de-novo
            {Keys.mendelian}: testable and inheritance pattern is mendelian
            {Keys.other}: testable and inheritance pattern is non-mendelian and not de-novo
        Form MendelianViolationCurve from counts of trio categories at different quality scores.
        Args:
            trio_category_score_counts: Mapping[str, pandas.Series]
                Map from trio category ({Keys.num_trios}, {Keys.de_novo}, {Keys.mendelian}, {Keys.other}) to pandas
                Series with count of trios in that category, indexed by lowest quality score in that trio
            plotted_mendelian_violation_keys: Optional[Iterable[str]]
                Compute MendelianViolationCurve for each category passed in. If None, compute all curves.
        Returns:
            mendelian_violation_curve: Optional[MendelianViolationCurve]
                Computed MendelianViolationCurve. If there are no scores available, return None
        """
        # get thresholds where Mendelian violation rate may change. Note these indices are guaranteed to be unique
        # because they are formed from pandas.Series.value_counts(), which returns the count of items with a unique
        # value
        thresholds = trio_category_score_counts[MendelianViolationCurve.num_trios_key].index.values
        if not numpy.isfinite(thresholds).any():
            return None  # the curve is empty

        if plotted_mendelian_violation_keys is not None and Keys.num_trios not in plotted_mendelian_violation_keys:
            # always need num_trios, because that's needed for normalization
            plotted_mendelian_violation_keys += (Keys.num_trios,)
        return MendelianViolationCurve(
            dataframe=pandas.DataFrame(
                {
                    trio_category: MendelianViolationCurve._num_passing_at_thresholds(score_counts, thresholds)
                    for trio_category, score_counts in trio_category_score_counts.items()
                    if plotted_mendelian_violation_keys is None or trio_category in plotted_mendelian_violation_keys
                },
                index=pandas.Index(thresholds, name="thresholds")
            )
        )

    def get_point_at_threshold(self, threshold: float, get_proportions: bool) -> pandas.Series:
        """
        Get mendelian/violation rates at point on curve corresponding to given threshold. This will be the first row
        with threshold >= requested_threshold
        Args:
            threshold: float
                Requested threshold on curve
            get_proportions: bool
                If true, scale rates by number of trios
        Returns:
            mendelian_violation_point: pandas.Series
                Series with threshold, and mendelian/violation number. Due to how pandas.Series work, threshold will be
                accessed as mendelian_violation_point.name
        """
        row_index = self.dataframe.index.searchsorted(threshold, side="left")

        if row_index < 0 or row_index >= len(self.dataframe):
            raise ValueError(f"No data at requested threshold ({threshold})")
        point = self.dataframe.iloc[row_index].rename(threshold)
        if get_proportions:
            num_trios = point[MendelianViolationCurve.num_trios_key]
            return pandas.Series(
                {key: (value if key == MendelianViolationCurve.num_trios_key else value / num_trios)
                 for key, value in point.items()},
                name=point.name
            )
        else:
            return point


def plot_precision_recall_curves(
        ax: pyplot.Axes, precision_recall_curves: Mapping[str, PrecisionRecallCurve],
        best_thresholds: Mapping[str, float], default_thresholds: Mapping[str, float],
        title_str: str = "",
        title_font_size: int = Default.title_font_size,
        label_font_size: int = Default.label_font_size,
        tick_font_size: int = Default.tick_font_size,
        plot_max_f1: bool = Default.plot_max_f1
):
    # scale curves that may be drawn from different numbers of variants
    num_good = 0
    num_bad = 0
    for label, precision_recall_curve in precision_recall_curves.items():
        if precision_recall_curve.num_good > num_good:
            num_good = precision_recall_curve.num_good
            num_bad = precision_recall_curve.num_bad
    ax.plot([], [])  # make empty line with no legend label, to sync colors across plot types
    for label, precision_recall_curve in precision_recall_curves.items():
        recall_scale = precision_recall_curve.num_good / num_good
        curve_line, = ax.plot(recall_scale * precision_recall_curve.recall, precision_recall_curve.precision,
                              label=label)
        if precision_recall_curve.is_empty:
            continue

        best_threshold = best_thresholds[label]
        default_threshold = default_thresholds[label]
        # noinspection PyTypeChecker
        p_default = precision_recall_curve.get_point_at_threshold(default_threshold)
        ax.plot(recall_scale * p_default[PrecisionRecallCurve.recall_key],
                p_default[PrecisionRecallCurve.precision_key], 's', color=curve_line.get_color(),
                label=f"thresh={p_default.name: .2f}, f={p_default[PrecisionRecallCurve.f_score_key]: .2f}")
        if plot_max_f1:
            p_best = precision_recall_curve.get_point_at_threshold(best_threshold)
            ax.plot(recall_scale * p_best[PrecisionRecallCurve.recall_key], p_best[PrecisionRecallCurve.precision_key],
                    'o', color=curve_line.get_color(),
                    label=f"thresh={p_best.name: .2f}, f={p_best[PrecisionRecallCurve.f_score_key]: .2f}")

    ax.tick_params(axis='both', which='major', labelsize=tick_font_size)
    ax.set_xlabel("Recall", fontsize=label_font_size)
    ax.set_ylabel("Precision", fontsize=label_font_size)
    # noinspection PyTypeChecker
    ax.set_xlim([0.0, 1.0])
    # noinspection PyTypeChecker
    ax.set_ylim([0.0, 1.0])
    if title_str:
        title_str = f"{title_str} ({num_good} good, {num_bad} bad)"
        ax.set_title(title_str, fontsize=title_font_size, verticalalignment="top")


def _draw_extra_legend_axis(
        fig: pyplot.Figure,
        gridspec: pyplot.GridSpec,
        row: int,
        column: int,
        line_labels: Iterable[str],
        markers: Iterable[Tuple[str, str]],
        legend_font_size: int = Default.legend_font_size,
        advance_color_cycle: int = 0
):
    # make extra axis with legend
    ax = fig.add_subplot(gridspec[row, column])
    ax.set_facecolor('w')
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    bbox = ax.get_tightbbox(renderer=fig.canvas.get_renderer())
    approximate_ax_width_chars = bbox.width / legend_font_size
    wrap_width = int(approximate_ax_width_chars)
    for __ in range(advance_color_cycle):
        # unlabeled lines to advance color cycle
        ax.plot([numpy.nan, numpy.nan], [numpy.nan, numpy.nan])
    for line_label in line_labels:
        ax.plot([numpy.nan, numpy.nan], [numpy.nan, numpy.nan], label='\n'.join(textwrap.wrap(line_label, wrap_width)))
    for marker, marker_label in markers:
        ax.plot(numpy.nan, numpy.nan, marker, label='\n'.join(textwrap.wrap(marker_label, wrap_width)))
    ax.legend(fontsize=legend_font_size)


def make_fig_and_gridspec(num_axes: int, num_rows: int) -> (pyplot.Figure, pyplot.GridSpec, int):
    fig = pyplot.figure(num=plotting.next_fig(), constrained_layout=True)
    num_columns = int(numpy.ceil(num_axes / num_rows))
    gridspec = fig.add_gridspec(num_rows, num_columns)
    return fig, gridspec, num_columns


def get_precision_recall_curves_figure(
        truth_data: TruthData,
        scores_data_sets: Collection[ScoresDataSet],
        sv_selectors: Mapping[str, SvTypeCutoffSelector],
        num_precision_recall_curve_rows: int = Default.num_precision_recall_curve_rows,
        title_font_size: int = Default.title_font_size,
        label_font_size: int = Default.label_font_size,
        legend_font_size: int = Default.legend_font_size,
        tick_font_size: int = Default.tick_font_size,
        plot_max_f1: bool = Default.plot_max_f1,
        pdf_writer: Optional[PdfPages] = None
) -> (pyplot.Figure, Dict[str, float]):
    # get the PrecisionRecallCurves broken down by category, keeping only the categories that have at least one
    # non-empty curve
    precision_recall_curves = {
        category: category_pr_curves
        for category, category_pr_curves in get_precision_recall_curves(
            truth_data=truth_data, scores_data_sets=scores_data_sets, sv_selectors=sv_selectors
        ).items()
        if not all(precision_recall_curve.is_empty for precision_recall_curve in category_pr_curves.values())
    }
    # get the overall best thresholds
    best_thresholds = {
        data_set_label: precision_recall_curve.get_point_at_max_f_score().name
        for data_set_label, precision_recall_curve in precision_recall_curves[Keys.non_bnd_types].items()
    }
    # get the default thresholds
    default_thresholds = {
        data_set.label: data_set.passing_score
        for data_set in scores_data_sets
    }

    # construct figure add axis for legend
    fig, gridspec, num_columns = make_fig_and_gridspec(num_axes=len(precision_recall_curves) + 1,
                                                       num_rows=num_precision_recall_curve_rows)
    row, column = 0, 0
    for sv_category, category_pr_curves in precision_recall_curves.items():
        ax = fig.add_subplot(gridspec[row, column])
        plot_precision_recall_curves(ax, category_pr_curves, best_thresholds, default_thresholds, title_str=sv_category,
                                     title_font_size=title_font_size, label_font_size=label_font_size,
                                     tick_font_size=tick_font_size, plot_max_f1=plot_max_f1)
        column += 1
        if column == num_columns:
            row, column = row + 1, 0

    # make extra axis with legend
    _draw_extra_legend_axis(
        fig, gridspec, row, column,
        line_labels=precision_recall_curves[Keys.all_variant_types].keys(),
        markers=[("ks", "default"), ("ko", "f1-max")] if plot_max_f1 else [("ks", "default")],
        legend_font_size=legend_font_size, advance_color_cycle=1
    )

    if pdf_writer is not None:
        pdf_writer.savefig(fig, bbox_inches="tight")
        return None, best_thresholds
    else:
        return fig, best_thresholds


def iter_categories(sv_selectors: Mapping[str, SvTypeCutoffSelector]) -> Iterator[Tuple[str, SvTypeCutoffSelector]]:
    yield Keys.all_variant_types, lambda metrics: pandas.Series(common.true(metrics.shape[0]),
                                                                index=metrics.index)
    if Keys.bnd_types in sv_selectors:
        bnd_selector = sv_selectors[Keys.bnd_types]
        yield Keys.non_bnd_types, lambda info: ~bnd_selector(info)

    for category in sorted(sv_selectors.keys()):
        yield category, sv_selectors[category]


def _filter_null_curves_categories(
        precision_recall_curves: Dict[str, Dict[str, PrecisionRecallCurve]]
) -> Dict[str, Dict[str, PrecisionRecallCurve]]:
    return {
        category: category_pr_curves
        for category, category_pr_curves in precision_recall_curves.items()
        if not all(precision_recall_curve.is_empty for precision_recall_curve in category_pr_curves.values())
    }


def get_precision_recall_curves(
        truth_data: TruthData,
        scores_data_sets: Collection[ScoresDataSet],
        sv_selectors: Mapping[str, SvTypeCutoffSelector]
) -> Dict[str, Dict[str, PrecisionRecallCurve]]:
    """
    return dict from SV category to another mapping, that being a dict from data set label to precision-recall curve
    """
    log("get_precision_recall_curves() ...")
    precision_recall_curves = _filter_null_curves_categories(
        {
            category: {
                data_set.label: data_set.select_category(
                    category=category, sv_type_cutoff_selector=selector
                ).get_precision_recall_curve(truth_data)
                for data_set in scores_data_sets
            }
            for category, selector in iter_categories(sv_selectors)
        }
    )
    log("COMPLETE: get_precision_recall_curves()")
    return precision_recall_curves


def _filter_null_stats_categories(
        stats_by_category: Optional[Mapping[str, pandas.DataFrame]]
) -> Optional[Dict[str, pandas.DataFrame]]:
    if stats_by_category is None:
        return None
    # noinspection PyUnresolvedReferences
    return {
        category: stats for category, stats in stats_by_category.items()
        if not (stats == 0).all().all()
    }


def get_mendelian_inheritance_stats(
        truth_data: TruthData,
        scores_data_sets: Collection[ScoresDataSet],
        sv_selectors: Mapping[str, SvTypeCutoffSelector],
        max_f1_thresholds: Mapping[str, float],
        plot_max_f1: bool,
        replot_unfiltered: bool,
) -> Optional[Dict[str, pandas.DataFrame]]:
    log("get_mendelian_inheritance_stats() ...")
    if truth_data.pedigree_file_info is None:
        log("SKIPPED: get_mendelian_inheritance_stats(), no truth data")
        return None

    inheritance_stats_dict = {}
    calculate_unfiltered = True
    for data_set in scores_data_sets:
        for filter_label, filter_score in {
            Keys.unfiltered: None,
            Keys.filtered: data_set.passing_score,
            Keys.max_f1: max_f1_thresholds.get(data_set.label, None)
        }.items():
            if filter_label == Keys.unfiltered and not calculate_unfiltered:
                continue
            elif filter_label == Keys.max_f1 and not plot_max_f1:
                continue  # don't need/want max-f1 so skip calculation
            log(f"get_inheritance_stats({data_set.label}:{filter_label}) ...")
            filtered_data_set = data_set.get_inheritance_stats(truth_data=truth_data) if filter_score is None \
                else data_set.apply_hard_filter(filter_score).get_inheritance_stats(truth_data=truth_data)

            for category, selector in iter_categories(sv_selectors):
                log(f"\tselect_variants: {category} ...")
                selected_filtered_data_set = filtered_data_set.select_category(category=category,
                                                                               sv_type_cutoff_selector=selector)
                log("\tsum inheritance stats ...")
                inheritance_stats_dict[(data_set.label, filter_label, category)] = \
                    selected_filtered_data_set.inheritance_stats.astype(pandas.UInt64Dtype()).sum(axis=0)
            log(f"get_inheritance_stats({data_set.label}:{filter_label}) COMPLETE")
            if filter_label == Keys.unfiltered and not replot_unfiltered:
                # don't re-compute unfiltered stats
                if calculate_unfiltered:
                    calculate_unfiltered = False

    filter_keys = (Keys.unfiltered, Keys.filtered, Keys.max_f1) if plot_max_f1 \
        else (Keys.unfiltered, Keys.filtered)
    mendelian_inheritance_stats = _filter_null_stats_categories(
        {
            category: pandas.DataFrame.from_dict(
                {
                    (data_set.label, filter_label): inheritance_stats_dict[(data_set.label, filter_label, category)]
                    for data_set in scores_data_sets
                    for filter_label in filter_keys
                    if (data_set.label, filter_label, category) in inheritance_stats_dict
                },
                orient="index"
            )
            for category, __ in iter_categories(sv_selectors)
        }
    )
    log("COMPLETE: get_mendelian_inheritance_stats()")
    return mendelian_inheritance_stats


def get_variants_per_sample(
        scores_data_sets: Collection[ScoresDataSet],
        sv_selectors: Mapping[str, SvTypeCutoffSelector],
        max_f1_thresholds: Mapping[str, float],
        plot_max_f1: bool,
        replot_unfiltered: bool,
) -> Dict[str, pandas.DataFrame]:
    log("get_variants_per_sample() ...")

    variants_per_sample_dict = {}
    calculate_unfiltered = True
    for data_set in scores_data_sets:
        for filter_label, filter_score in {
            Keys.unfiltered: None,
            Keys.filtered: data_set.passing_score,
            Keys.max_f1: max_f1_thresholds.get(data_set.label, None)
        }.items():
            if filter_label == Keys.unfiltered and not calculate_unfiltered:
                continue
            if filter_label == Keys.max_f1 and not plot_max_f1:
                continue  # don't need/want max-f1 so skip calculation
            log(f"get variants per sample({data_set.label}:{filter_label}) ...")
            filtered_data_set = data_set if filter_score is None else data_set.apply_hard_filter(filter_score)

            for category, selector in iter_categories(sv_selectors):
                log(f"\tselect_variants: {category} ...")
                selected_filtered_data_set = filtered_data_set.select_category(category=category,
                                                                               sv_type_cutoff_selector=selector)
                log("\tget variants per sample ...")
                variants_per_sample_dict[(data_set.label, filter_label, category)] = \
                    selected_filtered_data_set.variants_per_sample
            if filter_label == Keys.unfiltered and not replot_unfiltered:
                # don't re-compute unfiltered stats
                if calculate_unfiltered:
                    calculate_unfiltered = False

    filter_keys = (Keys.unfiltered, Keys.filtered, Keys.max_f1) if plot_max_f1 \
        else (Keys.unfiltered, Keys.filtered)

    variants_per_sample = _filter_null_stats_categories(
        {
            category: pandas.DataFrame.from_dict(
                {
                    (data_set.label, filter_label): variants_per_sample_dict[(data_set.label, filter_label, category)]
                    for data_set in scores_data_sets
                    for filter_label in filter_keys
                    if (data_set.label, filter_label, category) in variants_per_sample_dict
                },
                orient="index"
            )
            for category, __ in iter_categories(sv_selectors)
        }
    )
    log("COMPLETE: get_variants_per_sample()")
    return variants_per_sample


def get_mendelian_violation_curves(
        truth_data: TruthData,
        scores_data_sets: Collection[ScoresDataSet],
        sv_selectors: Mapping[str, SvTypeCutoffSelector],
        plotted_mendelian_violation_keys: Optional[Iterable[str]] = Default.plotted_mendelian_violation_trio_types
) -> Optional[Dict[str, Dict[str, MendelianViolationCurve]]]:
    if truth_data.pedigree_file_info is None:
        return None
    log("get_mendelian_violation_curves() ...")

    mendelian_violation_curves_dict = {category: {} for category, __ in iter_categories(sv_selectors)}
    for data_set in scores_data_sets:
        log(f"get_inheritance_stats({data_set.label}) ...")
        filtered_data_set = data_set.get_inheritance_stats(truth_data=truth_data)
        for category, selector in iter_categories(sv_selectors):
            log(f"\tselect_variants: {category} ...")
            selected_filtered_data_set = filtered_data_set.select_category(category=category,
                                                                           sv_type_cutoff_selector=selector)
            log("\tget violation curve ...")
            violation_curve = selected_filtered_data_set.get_mendelian_violation_curve(
                truth_data=truth_data, plotted_mendelian_violation_keys=plotted_mendelian_violation_keys
            )
            if violation_curve is not None:
                mendelian_violation_curves_dict[category][data_set.label] = violation_curve

    log("COMPLETE: get_mendelian_violation_curves()")
    return {
        sv_category: category_violation_curves
        for sv_category, category_violation_curves in mendelian_violation_curves_dict.items()
        if len(category_violation_curves) > 0
    }


def scale_inheritance_stats(
        mendelian_inheritance_stats: pandas.DataFrame
) -> (pandas.DataFrame, pandas.DataFrame, pandas.Series):
    num_trios = mendelian_inheritance_stats[Keys.num_trios].copy()
    # scale all stats by num_trios_stat
    scaled_stats = mendelian_inheritance_stats.div(num_trios, axis=0)
    # replace num_trios stat by num_trios scaled to unfiltered num_trios for each data set
    unfiltered_num_trios = num_trios.xs(Keys.unfiltered, level=1, drop_level=True)
    scaled_stats[Keys.num_trios] = num_trios / unfiltered_num_trios.max()
    return scaled_stats, num_trios.astype("int"), unfiltered_num_trios.astype("int")


def plot_mendelian_inheritance_stats(
        ax: pyplot.Axes,
        mendelian_inheritance_stats: pandas.DataFrame,
        title: str = "",
        title_font_size: int = Default.title_font_size,
        label_font_size: int = Default.label_font_size,
        tick_font_size: int = Default.tick_font_size,
        tick_labelrotation: float = Default.tick_labelrotation,
        plot_max_f1: bool = Default.plot_max_f1
):
    # scale inheritance stats for plotting
    scaled_inheritance_stats, num_trios, unfiltered_num_trios = scale_inheritance_stats(mendelian_inheritance_stats)

    # Draw an invisible line on the plot for each data set label, and collect the colors. This also preps the legend
    data_colors = {
        data_label: ax.plot([numpy.nan, numpy.nan], [numpy.nan, numpy.nan],
                            label=f"{data_label}: {num_trios} trios")[0].get_color()
        for data_label, num_trios in num_trios.items()
        if plot_max_f1 or not "".join(data_label).endswith(Keys.max_f1)
    }

    x = 0
    stat_tick_marks = {}
    for stat_name, stat_values in scaled_inheritance_stats.items():
        x += 1  # add padding between different stats
        stat_x_start = x
        for data_label, scaled_stat_value in stat_values.items():
            if "".join(data_label).endswith(Keys.max_f1) and not plot_max_f1:
                continue
            if not pandas.isnull(scaled_stat_value):
                x_range = [x, x + 1]
                y_range = [0, scaled_stat_value]
                plotting.rectangle(ax, x_range, y_range, facecolor=data_colors[data_label],
                                   edgecolor=data_colors[data_label], label="_nolegend_")
            x += 1
        stat_tick_marks[stat_name] = (stat_x_start + x) / 2

    ax.tick_params(axis='both', which='major', labelsize=tick_font_size)
    ax.set_xticks(list(stat_tick_marks.values()))
    ax.set_xticklabels(list(stat_tick_marks.keys()), horizontalalignment="left", fontsize=tick_font_size)
    ax.tick_params(axis='x', labelrotation=tick_labelrotation)
    ax.autoscale(enable=True, tight=True)
    ax.set_ylabel("Proportion", fontsize=label_font_size)
    # noinspection PyTypeChecker
    ax.set_xlim([0.5, x + 0.5])
    # noinspection PyTypeChecker
    ax.set_ylim([0.0, 1.0])
    ax.set_yscale("symlog", linthresh=0.01)
    ax.set_yticks([0, 0.01, 0.1, 1.0])
    ax.set_yticklabels(["0%", "1%", "10%", "100%"], fontsize=tick_font_size)
    unfiltered_num_trios = unfiltered_num_trios.values.max()
    title_str = f"{title}: {unfiltered_num_trios:d} trios" if title else f"{unfiltered_num_trios:d} trios"
    ax.set_title(title_str, fontsize=title_font_size, verticalalignment="top")


def plot_mendelian_violation_curves(
        ax: pyplot.Axes,
        mendelian_violation_curves: Mapping[str, MendelianViolationCurve],
        best_thresholds: Mapping[str, float],
        default_thresholds: Mapping[str, float],
        title_str: str = "",
        title_font_size: int = Default.title_font_size,
        label_font_size: int = Default.label_font_size,
        tick_font_size: int = Default.tick_font_size,
        plot_max_f1: bool = Default.plot_max_f1,
        plotted_mendelian_violation_keys: Optional[Iterable[str]] = Default.plotted_mendelian_violation_trio_types
):
    min_x = None
    max_x = None
    max_y = 0
    ax.plot([], [])  # make empty line with no legend label, to sync colors across plot types

    # noinspection PyPep8Naming
    FloatType = TypeVar("FloatType", float, numpy.ndarray)

    def _get_bad_percent_and_label(curve_y_value: FloatType, violation_curve_key: str) -> (FloatType, str):
        if violation_curve_key == Keys.mendelian:
            return 100 * (1.0 - curve_y_value), "non-mendelian-rate"
        elif violation_curve_key == Keys.num_trios:
            return 100 * (1.0 - curve_y_value), "non-testable-rate"
        elif violation_curve_key == Keys.de_novo:
            return 100 * curve_y_value, "de-novo-rate"
        else:
            return 100 * curve_y_value, "other-bad-rate"

    for label, mendelian_violation_curve in mendelian_violation_curves.items():
        def _get_curve_data(violation_curve_key: str) -> FloatType:
            if violation_curve_key == Keys.mendelian:
                return mendelian_violation_curve.proportion_mendelian
            elif violation_curve_key == Keys.num_trios:
                return mendelian_violation_curve.proportion_testable
            elif violation_curve_key == Keys.de_novo:
                return mendelian_violation_curve.proportion_de_novo
            else:
                return mendelian_violation_curve.proportion_other_bad

        if min_x is None:
            dt_info = get_dtype_info(mendelian_violation_curve.thresholds.dtype)
            min_x = dt_info.max  # set these values to totally wrong side of things, so final result is tight
            max_x = dt_info.min  # on the data
        min_x = mendelian_violation_curve.thresholds.min(initial=min_x)
        max_x = mendelian_violation_curve.thresholds.max(initial=max_x)
        best_threshold = best_thresholds[label] if plot_max_f1 else None
        default_threshold = default_thresholds[label]
        # noinspection PyTypeChecker
        p_default = mendelian_violation_curve.get_point_at_threshold(default_threshold, get_proportions=True)
        p_best = mendelian_violation_curve.get_point_at_threshold(best_threshold, get_proportions=True) if plot_max_f1 \
            else None

        for curve_key in p_default.index if plotted_mendelian_violation_keys is None \
                else plotted_mendelian_violation_keys:
            curve_y_values, curve_label = _get_bad_percent_and_label(_get_curve_data(curve_key), curve_key)
            max_y = curve_y_values.max(initial=max_y)
            curve_line, = ax.plot(mendelian_violation_curve.thresholds, curve_y_values, label=curve_label)

            y_default, __ = _get_bad_percent_and_label(p_default[curve_key], curve_key)
            ax.plot(
                p_default.name, y_default, 's', color=curve_line.get_color(),
                label=f"default thresh={p_default.name: .2f}, {curve_label}={y_default: .2f}"
            )
            if plot_max_f1:
                y_best, __ = _get_bad_percent_and_label(p_best[curve_key], curve_key)
                ax.plot(p_best.name, y_best, 'o', color=curve_line.get_color(),
                        label=f"best thresh={p_best.name: .2f}, {curve_label}={y_best: .2f}")

    ax.tick_params(axis='both', which='major', labelsize=tick_font_size)
    ax.set_xlabel("Threshold", fontsize=label_font_size)
    ax.set_ylabel("% non-mendelian", fontsize=label_font_size)
    # noinspection PyTypeChecker
    ax.set_xlim([min_x, max_x])
    ax.set_xscale("symlog")
    # noinspection PyTypeChecker
    ax.set_ylim([0.0, max_y])
    ax.set_yscale("symlog", linthresh=5)  # linear to 5 %, then log scale
    if title_str:
        ax.set_title(title_str, fontsize=title_font_size, verticalalignment="top")


def _set_component_color(_plot_component, _color, _alpha=0.7):
    if hasattr(_plot_component, "set_color"):
        _plot_component.set_color(_color)
    else:
        if hasattr(_plot_component, "set_facecolor"):
            _plot_component.set_facecolor(_color)
        if hasattr(_plot_component, "set_edgecolor"):
            _plot_component.set_edgecolor(_color)
    if hasattr(_plot_component, "set_alpha"):
        _plot_component.set_alpha(_alpha)


def _set_violin_obj_color(violin_obj, color):
    for pc in violin_obj["bodies"]:
        _set_component_color(pc, color)
    for component_name, component in violin_obj.items():
        if isinstance(component, Iterable):
            for pc in component:
                _set_component_color(pc, color)
        else:
            _set_component_color(component, color)


def plot_variants_per_sample(
        ax: pyplot.Axes,
        variants_per_sample: pandas.DataFrame,
        title: str = "",
        title_font_size: int = Default.title_font_size,
        label_font_size: int = Default.label_font_size,
        tick_font_size: int = Default.tick_font_size,
        plot_max_f1: bool = Default.plot_max_f1,
        histogram_width: float = Default.histogram_width
):
    median_variants_per_sample = variants_per_sample.median(axis=1)
    # noinspection PyTypeChecker
    # Draw an invisible line on the plot for each data set label, and collect the colors.
    data_colors = {
        data_label: ax.plot([], [], label=f"{data_label}: {median_value} variants_per_sample")[0].get_color()
        for data_label, median_value in median_variants_per_sample.items()
        if plot_max_f1 or not "".join(data_label).endswith(Keys.max_f1)
    }

    x = 0
    max_y = 0
    for data_label, num_variants in variants_per_sample.iterrows():
        # noinspection PyUnresolvedReferences,PyTypeChecker
        if "".join(data_label).endswith(Keys.max_f1) and not plot_max_f1:
            continue
        x += 1
        max_y = max(max_y, num_variants.max())
        violin_obj = ax.violinplot(num_variants, positions=[x], widths=[histogram_width],
                                   showmeans=False, showmedians=True, showextrema=True)
        _set_violin_obj_color(violin_obj, data_colors[data_label])

    ax.tick_params(axis='both', which='major', labelsize=tick_font_size)
    ax.set_xticks([])
    ax.autoscale(enable=True, tight=True)
    ax.set_ylabel("variants/sample", fontsize=label_font_size)
    # noinspection PyTypeChecker
    ax.set_xlim([0.5, x + 0.5])
    # noinspection PyTypeChecker
    ax.set_ylim([0, max_y])
    # noinspection PyUnresolvedReferences
    if title:
        ax.set_title(title, fontsize=title_font_size, verticalalignment="top")


def get_mendelian_inheritance_figure(
        truth_data: TruthData,
        scores_data_sets: Collection[ScoresDataSet],
        sv_selectors: Mapping[str, SvTypeCutoffSelector],
        best_thresholds: Optional[Mapping[str, float]] = None,
        num_inheritance_stats_rows: int = Default.num_inheritance_stats_rows,
        title_font_size: int = Default.title_font_size,
        label_font_size: int = Default.label_font_size,
        legend_font_size: int = Default.legend_font_size,
        tick_font_size: int = Default.tick_font_size,
        tick_labelrotation: float = Default.tick_labelrotation,
        plot_max_f1: bool = Default.plot_max_f1,
        replot_unfiltered: bool = Default.replot_unfiltered,
        pdf_writer: Optional[PdfPages] = None
) -> Optional[pyplot.Figure]:
    mendelian_inheritance_stats = get_mendelian_inheritance_stats(
        truth_data=truth_data, scores_data_sets=scores_data_sets, sv_selectors=sv_selectors,
        max_f1_thresholds=best_thresholds, plot_max_f1=plot_max_f1, replot_unfiltered=replot_unfiltered,
    )
    log("get_mendelian_inheritance_figure() ...")
    if mendelian_inheritance_stats is None:
        log("Skipped: get_mendelian_inheritance_figure(), no inheritance data")
        return None

    # construct mendelian inheritance figure and variants per sample figure in the same loop,
    # add extra axis for legend
    mendelian_inheritance_fig, mendelian_inheritance_gridspec, num_columns = make_fig_and_gridspec(
        num_axes=len(mendelian_inheritance_stats) + 1, num_rows=num_inheritance_stats_rows
    )

    row, column = 0, 0
    for sv_category, category_inheritance_stats in mendelian_inheritance_stats.items():
        inheritance_ax = mendelian_inheritance_fig.add_subplot(mendelian_inheritance_gridspec[row, column])
        plot_mendelian_inheritance_stats(
            inheritance_ax, category_inheritance_stats, title=sv_category, title_font_size=title_font_size,
            label_font_size=label_font_size, tick_labelrotation=tick_labelrotation, tick_font_size=tick_font_size,
            plot_max_f1=plot_max_f1
        )
        column += 1
        if column == num_columns:
            row, column = row + 1, 0

    filter_line_labels = [
        ' '.join(data_label)
        for data_label in mendelian_inheritance_stats[Keys.all_variant_types].index
        if plot_max_f1 or not "".join(data_label).endswith(Keys.max_f1)
    ]
    _draw_extra_legend_axis(
        mendelian_inheritance_fig, mendelian_inheritance_gridspec, row, column,
        line_labels=filter_line_labels, markers=[],
        legend_font_size=legend_font_size
    )

    log("COMPLETE: get_mendelian_inheritance_figure()")
    if pdf_writer is not None:
        pdf_writer.savefig(mendelian_inheritance_fig, bbox_inches="tight")
        return None
    else:
        return mendelian_inheritance_fig


def get_variants_per_sample_figure(
        scores_data_sets: Collection[ScoresDataSet],
        sv_selectors: Mapping[str, SvTypeCutoffSelector],
        best_thresholds: Optional[Mapping[str, float]] = None,
        num_inheritance_stats_rows: int = Default.num_inheritance_stats_rows,
        title_font_size: int = Default.title_font_size,
        label_font_size: int = Default.label_font_size,
        legend_font_size: int = Default.legend_font_size,
        tick_font_size: int = Default.tick_font_size,
        plot_max_f1: bool = Default.plot_max_f1,
        replot_unfiltered: bool = Default.replot_unfiltered,
        pdf_writer: Optional[PdfPages] = None
) -> Optional[pyplot.Figure]:
    variants_per_sample = get_variants_per_sample(
        scores_data_sets=scores_data_sets, sv_selectors=sv_selectors,
        max_f1_thresholds=best_thresholds, plot_max_f1=plot_max_f1, replot_unfiltered=replot_unfiltered,
    )

    log("get_variants_per_sample_figure() ...")

    variants_per_sample_fig, variants_per_sample_gridspec, num_columns = make_fig_and_gridspec(
        num_axes=len(variants_per_sample) + 1, num_rows=num_inheritance_stats_rows
    )

    row, column = 0, 0
    for sv_category, category_variants_per_sample in variants_per_sample.items():
        variants_per_sample_ax = variants_per_sample_fig.add_subplot(variants_per_sample_gridspec[row, column])
        plot_variants_per_sample(
            variants_per_sample_ax, category_variants_per_sample, title=sv_category, title_font_size=title_font_size,
            label_font_size=label_font_size, tick_font_size=tick_font_size, plot_max_f1=plot_max_f1
        )
        column += 1
        if column == num_columns:
            row, column = row + 1, 0

    filter_line_labels = [
        ' '.join(data_label)
        for data_label in variants_per_sample[Keys.all_variant_types].index
        if plot_max_f1 or not "".join(data_label).endswith(Keys.max_f1)
    ]
    _draw_extra_legend_axis(
        variants_per_sample_fig, variants_per_sample_gridspec, row, column,
        line_labels=filter_line_labels, markers=[],
        legend_font_size=legend_font_size
    )

    log("COMPLETE: variants_per_sample_figure()")
    if pdf_writer is not None:
        pdf_writer.savefig(variants_per_sample_fig, bbox_inches="tight")
        return None
    else:
        return variants_per_sample_fig


def get_mendelian_violation_curve_figure(
        truth_data: TruthData,
        scores_data_sets: Collection[ScoresDataSet],
        sv_selectors: Mapping[str, SvTypeCutoffSelector],
        best_thresholds: Optional[Mapping[str, float]] = None,
        num_inheritance_stats_rows: int = Default.num_inheritance_stats_rows,
        title_font_size: int = Default.title_font_size,
        label_font_size: int = Default.label_font_size,
        legend_font_size: int = Default.legend_font_size,
        tick_font_size: int = Default.tick_font_size,
        plot_max_f1: bool = Default.plot_max_f1,
        plotted_mendelian_violation_keys: Optional[Iterable[str]] = Default.plotted_mendelian_violation_trio_types,
        pdf_writer: Optional[PdfPages] = None
) -> Optional[pyplot.Figure]:
    mendelian_violation_curves = get_mendelian_violation_curves(
        truth_data=truth_data, scores_data_sets=scores_data_sets, sv_selectors=sv_selectors,
        plotted_mendelian_violation_keys=plotted_mendelian_violation_keys
    )

    log("get_mendelian_violation_curve_figure() ...")
    if mendelian_violation_curves is None:
        # no data to plot
        log("SKIPPED: get_mendelian_violation_curve_figure(), no inheritance data")
        return None

    mendelian_violation_fig, mendelian_violation_gridspec, num_columns = make_fig_and_gridspec(
        num_axes=len(mendelian_violation_curves) + 1, num_rows=num_inheritance_stats_rows
    )

    # get the default thresholds
    default_thresholds = {
        data_set.label: data_set.passing_score
        for data_set in scores_data_sets
    }

    row, column = 0, 0
    violation_curve_axes = []
    for sv_category, category_violation_curves in mendelian_violation_curves.items():
        violation_curve_ax = mendelian_violation_fig.add_subplot(mendelian_violation_gridspec[row, column])
        violation_curve_axes.append(violation_curve_ax)
        category_violation_curves = mendelian_violation_curves[sv_category]
        plot_mendelian_violation_curves(
            violation_curve_ax, category_violation_curves, title_str=sv_category, title_font_size=title_font_size,
            label_font_size=label_font_size, tick_font_size=tick_font_size, plot_max_f1=plot_max_f1,
            default_thresholds=default_thresholds, best_thresholds=best_thresholds,
            plotted_mendelian_violation_keys=plotted_mendelian_violation_keys
        )
        column += 1
        if column == num_columns:
            row, column = row + 1, 0

    plotting.sync_axes(violation_curve_axes, sync_x=True, sync_y=True, zoom_to_fit_all=True)
    _draw_extra_legend_axis(
        mendelian_violation_fig, mendelian_violation_gridspec, row, column,
        line_labels=mendelian_violation_curves[Keys.all_variant_types].keys(),
        markers=[("ks", "default"), ("ko", "f1-max")] if plot_max_f1 else [("ks", "default")],
        legend_font_size=legend_font_size, advance_color_cycle=1
    )

    log("COMPLETE: get_mendelian_violation_curve_figure()")
    if pdf_writer is not None:
        pdf_writer.savefig(mendelian_violation_fig, bbox_inches="tight")
        return None
    else:
        return mendelian_violation_fig


def get_scores_histogram_figure(
        scores_data_sets: Collection[ScoresDataSet],
        sv_selectors: Mapping[str, SvTypeCutoffSelector],
        num_scores_histogram_rows: int = Default.num_scores_histogram_rows,
        histogram_width: float = Default.histogram_width,
        title_font_size: int = Default.title_font_size,
        label_font_size: int = Default.label_font_size,
        legend_font_size: int = Default.legend_font_size,
        tick_font_size: int = Default.tick_font_size,
        tick_labelrotation: float = Default.tick_labelrotation,
        pdf_writer: Optional[PdfPages] = None
) -> Optional[pyplot.Figure]:
    log("get_scores_histogram_figure() ...")

    # Need to compute which sv_categories have scores
    log("\tFind which SV categories have scores ....")
    category_has_data = {
        sv_category: any(scores_data_set.has_category_data(sv_selector) for scores_data_set in scores_data_sets)
        for sv_category, sv_selector in iter_categories(sv_selectors)
    }
    # allocate one more axes than there are sv categories with data (for legend)
    histogram_fig, histogram_gridspec, num_columns = make_fig_and_gridspec(
        num_axes=sum(category_has_data.values(), start=1), num_rows=num_scores_histogram_rows
    )

    # get the colors from the color cycle, starting on index 1 (since there's no "unfiltered" scores set)
    colors = pyplot.rcParams["axes.prop_cycle"].by_key()["color"]
    data_colors = {scores_data_set.label: colors[ind] for ind, scores_data_set in enumerate(scores_data_sets, start=1)}
    row, column = 0, 0
    for sv_category, sv_selector in iter_categories(sv_selectors):
        if not category_has_data[sv_category]:
            continue
        ax = histogram_fig.add_subplot(histogram_gridspec[row, column])
        # noinspection PyUnboundLocalVariable
        plot_scores_histogram(
            ax=ax, sv_category=sv_category, sv_selector=sv_selector, scores_data_sets=scores_data_sets,
            data_colors=data_colors, histogram_width=histogram_width, title_font_size=title_font_size,
            label_font_size=label_font_size, tick_font_size=tick_font_size, tick_labelrotation=tick_labelrotation
        )
        column += 1
        if column == num_columns:
            row, column = row + 1, 0

    _draw_extra_legend_axis(
        histogram_fig, histogram_gridspec, row, column,
        line_labels=[scores_data_set.label for scores_data_set in scores_data_sets], markers=[],
        legend_font_size=legend_font_size, advance_color_cycle=1
    )

    log("COMPLETE: get_scores_histogram_figure()")
    if pdf_writer is not None:
        pdf_writer.savefig(histogram_fig, bbox_inches="tight")
        return None
    else:
        return histogram_fig


def _get_scores_hisogram(
        scores: pandas.DataFrame,
        scores_label: str,
        null_sentinal: float,
        num_bins: int = 100
) -> (numpy.ndarray, numpy.ndarray):
    # return counts, edges
    if scores.empty:
        log(f"\t\t{scores_label}: no data")
        # no data, so return empty histogram
        return numpy.empty(0), numpy.empty(0)
    # basically want to:
    # 1) convert scores to int64 numpy array and call ravel to get as a 1D array
    # 2) call histogram to get the histogram
    # however there are some fiddly bits with point 1, because:
    #    -scores may or may not be a nullable dtype, and
    #    -scores may be signed or unsigned
    # These problems interact if we need to change a nulled value to a negative number, but the values are unsigned

    if any(pandas.api.types.is_float_dtype(_dt) for _dt in scores.dtypes):
        # At least some scores are float scores that may contain NaN. Convert all to float and then use histogram
        scores = scores.astype(numpy.float64).fillna(null_sentinal)
    elif any(isinstance(_dt, pandas.api.extensions.ExtensionDtype) for _dt in scores.dtypes):
        # scores are nullable integers. If null_sentinal is non-negative, there's no problem, otherwise we'll need to
        # convert all unsigned nullable types to signed
        if null_sentinal < 0:
            for _column, _dt in list(scores.dtypes.items()):
                if isinstance(_dt, pandas.api.extensions.ExtensionDtype) and \
                        pandas.api.types.is_unsigned_integer_dtype(_dt):
                    # convert nullable unsigned int to nullable signed int
                    # noinspection PyUnresolvedReferences
                    _signed_dt = pandas.core.dtypes.cast.pandas_dtype(_dt.name.replace('UInt', 'Int'))
                    scores.loc[:, _column] = scores.loc[:, _column].astype(_signed_dt)

        scores = scores.fillna(int(numpy.floor(null_sentinal))).astype(numpy.int64)

    log(f"\t\t{scores_label}: {scores.size} scores")
    return numpy.histogram(scores.values.ravel(), bins=num_bins, density=False)


def plot_scores_histogram(
        ax: pyplot.Axes,
        sv_category: str,
        sv_selector: SvTypeCutoffSelector,
        scores_data_sets: Collection[ScoresDataSet],
        data_colors: Mapping[str, str],
        histogram_width: float = Default.histogram_width,
        title_font_size: int = Default.title_font_size,
        label_font_size: int = Default.label_font_size,
        tick_font_size: int = Default.tick_font_size,
        tick_labelrotation: float = Default.tick_labelrotation,
) -> pyplot.Axes:
    log(f"\tPlot scores for {sv_category} ...")
    min_y, max_y = 0, 0
    pad = (1.0 - histogram_width) / 2
    x_ticks = []
    x_tick_labels = []
    for data_set_index, scores_data_set in enumerate(scores_data_sets):
        counts, edges = _get_scores_hisogram(
            scores=scores_data_set.select_category(sv_category, sv_selector).scores,
            scores_label=scores_data_set.label,
            null_sentinal=scores_data_set.passing_score - 1, num_bins=100
        )
        if counts.size > 0:
            log_max_counts = numpy.log10(counts.max(initial=0))
            if log_max_counts >= 2:  # can make three ticks at nice even numbers
                ticks = [1, 10 ** round(log_max_counts / 2), 10 ** int(numpy.floor(log_max_counts))]
                ticks_str = ["1", f"10^{round(log_max_counts / 2)}", f"10^{int(numpy.floor(log_max_counts))}"]
            elif log_max_counts >= 1:
                ticks = [1, 10 ** int(numpy.floor(log_max_counts))]
                ticks_str = ["1", f"10^{int(numpy.floor(log_max_counts))}"]
            elif log_max_counts >= 0:
                ticks = [1, round(10 ** log_max_counts)]
                ticks_str = ["1", str(ticks[-1])]
            else:
                ticks = [1]
                ticks_str = ["1"]
            x_tick_labels.extend(ticks_str)

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                values = numpy.log(counts)
            scale = (histogram_width / values.max(initial=0))
            values = scale * values
            x_baseline = data_set_index - 0.5 + pad
            x_ticks.extend(x_baseline + scale * numpy.log(tick) for tick in ticks)

            plotting.barh(ax, values=values, edges=edges, x_baseline=x_baseline,
                          color=data_colors[scores_data_set.label], fill=True)
            min_y = min(edges[0], min_y)
            max_y = max(edges[-1], max_y)

    ax.tick_params(axis='both', which='major', labelsize=tick_font_size)
    ax.set_title(sv_category, fontsize=title_font_size)
    ax.set_ylabel("score", fontsize=label_font_size)
    ax.set_xlabel("variant density", fontsize=label_font_size)
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_tick_labels)
    ax.tick_params(axis='x', labelrotation=tick_labelrotation)
    ax.set_xlim(-0.5, len(scores_data_sets) - 0.5)
    if max_y > min_y:
        ax.set_ylim(min_y, max_y)
    else:
        ax.set_ylim(-0.5, 0.5)
    return ax


def get_hardy_weinberg_figures(
        scores_data_sets: Sequence[ScoresDataSet],
        sv_selectors: Mapping[str, SvTypeCutoffSelector],
        title_font_size: int = Default.title_font_size,
        label_font_size: int = Default.label_font_size,
        tick_font_size: int = Default.tick_font_size,
        plotted_hardy_weinberg_categories: Optional[Container[str]] = Default.plotted_hardy_weinberg_categories,
        min_num_called_gt_hardy_weinberg: int = Default.min_num_called_gt_hardy_weinberg,
        hardy_weinberg_confidence: float = Default.hardy_weinberg_confidence,
        pdf_writer: Optional[PdfPages] = None
) -> Tuple[pyplot.Figure, ...]:
    return tuple(
        get_hardy_weinberg_figure(
            figure_title=f"Hardy-Weinberg: {category}",
            scores_data_sets=[scores_data_set.select_category(category, selector)
                              for scores_data_set in scores_data_sets],
            title_font_size=title_font_size,
            label_font_size=label_font_size,
            tick_font_size=tick_font_size,
            min_num_called_gt_hardy_weinberg=min_num_called_gt_hardy_weinberg,
            hardy_weinberg_confidence=hardy_weinberg_confidence,
            pdf_writer=pdf_writer
        )
        for category, selector in iter_categories(sv_selectors)
        if plotted_hardy_weinberg_categories is None or category in plotted_hardy_weinberg_categories
    )


def get_hardy_weinberg_figure(
        figure_title: str,
        scores_data_sets: Sequence[ScoresDataSet],
        title_font_size: int = Default.title_font_size,
        label_font_size: int = Default.label_font_size,
        tick_font_size: int = Default.tick_font_size,
        min_num_called_gt_hardy_weinberg: int = Default.min_num_called_gt_hardy_weinberg,
        hardy_weinberg_confidence: float = Default.hardy_weinberg_confidence,
        pdf_writer: Optional[PdfPages] = None
) -> Optional[pyplot.Figure]:
    log(f"Getting figure for {figure_title}")

    # construct figure add axes for legend. We want these to be close to square, so dynamically scale num_rows based
    # on the number of data sets
    num_axes = len(scores_data_sets) + 2  # one for each ScoresDataSet, one for unfiltered, one for legend
    num_rows = int(numpy.ceil(numpy.sqrt(num_axes)))
    fig, gridspec, num_columns = make_fig_and_gridspec(num_axes=num_axes, num_rows=num_rows)
    fig.suptitle(figure_title, fontsize=title_font_size)
    # create axes for colorbar
    subgrid_spec = gridspec[num_rows - 1, num_columns - 1].subgridspec(1, max(2, 24 // num_columns))
    colorbar_axes = fig.add_subplot(subgrid_spec[0, 0])
    colorbar_axes.set_facecolor('w')
    colorbar_axes.grid(False)
    colorbar_axes.set_xticks([])
    colorbar_axes.set_yticks([])
    legend_axes = fig.add_subplot(subgrid_spec[0, 1:])
    legend_axes.set_facecolor('w')
    legend_axes.grid(False)
    legend_axes.set_xticks([])
    legend_axes.set_yticks([])

    # plot each HWE axis
    row, column = 0, 0
    ax = fig.add_subplot(gridspec[row, column])
    ax.grid(visible=False, which="both")
    plot_hardy_weinberg_proportions(
        ax,
        scores_data_sets[0].to_autosome(min_num_called_gt_hardy_weinberg=min_num_called_gt_hardy_weinberg),
        title=f"{scores_data_sets[0].label} {Keys.unfiltered}", title_font_size=title_font_size,
        label_font_size=label_font_size, tick_font_size=tick_font_size, colorbar_axes=colorbar_axes,
        legend_axes=legend_axes, min_num_called_gt_hardy_weinberg=min_num_called_gt_hardy_weinberg,
        hardy_weinberg_confidence=hardy_weinberg_confidence
    )
    column += 1
    if column == num_columns:
        row, column = row + 1, 0
    for scores_data_set in scores_data_sets:
        ax = fig.add_subplot(gridspec[row, column])
        plot_hardy_weinberg_proportions(
            ax,
            scores_data_set.to_autosome(min_num_called_gt_hardy_weinberg=min_num_called_gt_hardy_weinberg)
            .apply_hard_filter(),
            title=scores_data_set.label, title_font_size=title_font_size, label_font_size=label_font_size,
            tick_font_size=tick_font_size, min_num_called_gt_hardy_weinberg=min_num_called_gt_hardy_weinberg,
            hardy_weinberg_confidence=hardy_weinberg_confidence
        )
        column += 1
        if column == num_columns:
            row, column = row + 1, 0

    log(f"COMPLETE: Getting figure for {figure_title}")
    if pdf_writer is not None:
        pdf_writer.savefig(fig, bbox_inches="tight")
        return None
    else:
        return fig


def _hardy_weinberg_chi_squared_to_num_hets_boundary(
        hardy_weinberg_confidence: float,
        alt_allele_frequency: numpy.ndarray,
        num_ref_genotypes: numpy.ndarray,
        num_called_genotypes: Union[numpy.ndarray, int]
) -> (numpy.ndarray, numpy.ndarray):
    """
    Given a likelihood threshold and existing observed allele / genotype data, compute confidence interval on the number
    of HET genotypes compatible with Hardy-Weinberg Equilibrium (HWE). Uses chi2 distribution as model of variability.
    Args:
        hardy_weinberg_confidence: float
            Proportion of calls that should be in confidence interval for number of HET calls, assuming HWE
        alt_allele_frequency: numpy.ndarray
            num_variants length array of allele frequency for ALT allele
        num_ref_genotypes: numpy.ndarray
            number of genotypes called REF for each variant
        num_called_genotypes: numpy.ndarray or int
            number of genotypes called for each variant
    Returns:
        num_expected_het: numpy.ndarray
            number of expected HET calls, if each variant was in HWE
        het_delta: numpy.ndarray
            +-delta to high/low boundary of confidence interval for HET calls in HWE
    """
    # get chi-squared threshold at desired rareness, with one degree of freedom (only #HETs varying)
    chi_squared = scipy.stats.chi2.ppf(hardy_weinberg_confidence, df=1)
    ref_allele_frequency = 1.0 - alt_allele_frequency
    num_expected_ref = num_called_genotypes * ref_allele_frequency ** 2
    num_expected_het = num_called_genotypes * 2 * ref_allele_frequency * alt_allele_frequency
    # MATH and logic below
    # num_expected_homvar = num_called_genotypes * alt_allele_frequency ** 2
    # chi_squared = (
    #     (num_ref_genotypes - num_expected_ref) ** 2 +
    #     (num_homvar_genotypes - num_expected_homvar) ** 2 +
    #     (num_het_genotypes - num_expected_het) ** 2
    # ) / num_called_genotypes
    # num_homvar_genotypes = (num_alt_alleles - num_het_genotypes) / 2
    # num_homvar_genotypes = num_called_genotypes * alt_allele_frequency - num_het_genotypes / 2
    # num_homvar_genotypes = num_expected_homvar / alt_allele_frequency - num_het_genotypes / 2
    #      => num_homvar_genotypes - num_expected_homvar =
    #               = num_expected_homvar * (1 - alt_allele_frequency) / alt_allele_frequency - num_het_genotypes / 2
    #               = num_expected_homvar * ref_allele_frequency / alt_allele_frequency - num_het_genotypes / 2
    #               = num_expected_het / 2 - num_het_genotypes / 2
    # chi_squared = (
    #     (num_ref_genotypes - num_expected_ref) ** 2 +
    #     1/4 * (num_het_genotypes - num_expected_het) ** 2 +
    #     (num_het_genotypes - num_expected_het) ** 2
    # ) / num_called_genotypes
    # het_delta = +/- abs(num_het_genotypes - num_expected_het)
    #           = sqrt(4/5 * (chi_squared * num_called_genotypes - (num_ref_genotypes - num_expected_ref) ** 2)
    het_delta = numpy.sqrt(
        numpy.maximum(
            0.8 * (chi_squared * num_called_genotypes - (num_ref_genotypes - num_expected_ref) ** 2),
            0
        )
    )
    return num_expected_het, het_delta


def _compute_hardy_weinberg_proportion_accepted(
        num_het_genotypes: numpy.ndarray,
        num_ref_genotypes: numpy.ndarray,
        alt_allele_frequency: numpy.ndarray,
        num_called_genotypes: numpy.ndarray,
        hardy_weinberg_confidence: float
) -> float:
    """ compute proportion of calls that are within """
    # get confidence interval for number of HET calls for each variant
    num_expected_het, het_delta = _hardy_weinberg_chi_squared_to_num_hets_boundary(
        hardy_weinberg_confidence=hardy_weinberg_confidence, alt_allele_frequency=alt_allele_frequency,
        num_ref_genotypes=num_ref_genotypes, num_called_genotypes=num_called_genotypes
    )
    # compute proportion of calls inside their confidence interval
    # noinspection PyUnresolvedReferences
    return (numpy.abs(num_het_genotypes - num_expected_het) <= het_delta).mean()


def _compute_hardy_weinberg_plotted_boundary(
        proportion_het_grid: numpy.ndarray,
        alt_allele_frequency_grid: numpy.ndarray,
        num_samples: int,
        hardy_weinberg_confidence: float
) -> (numpy.ndarray, numpy.ndarray):
    """
    Compute confidence boundary for proportion of genotypes called HET on the Hardy-Weinberge Equilibrium (HWE) heat-map
    plot.
    Args:
        proportion_het_grid: numpy.ndarray
            Ed
        alt_allele_frequency_grid:
        num_samples:
        hardy_weinberg_confidence:

    Returns:
        proportion_het_low: numpy.ndarray
        proportion_het_high: numpy.ndarray
    """
    # We're working on a grid, where the number of REF calls is not defined. To visualize approximate confidence
    # boundary, treat all no-calls as REF, in which case:
    #   num_called_genotypes = num_samples
    #   num_ref_genotypes = num_samples - num_het_genotypes - num_homvar_genotypes,
    #           num_het_genotypes = num_samples * proportion_het
    #           num_homvar_genotypes = (num_alt_alleles - num_het_genotypes) / 2
    #           num_homvar_genotypes = (2 * num_samples * alt_allele_frequency - num_het_genotypes) / 2
    #   => num_ref_genotypes = num_samples  - num_het_genotypes / 2 - alt_allele_frequency * num_samples
    #                        = num_samples * (1 - proportion_het / 2 - alt_allele_frequency)

    def _num_ref_genotypes(_proportion_het: numpy.ndarray) -> numpy.ndarray:
        return num_samples * (1.0 - _proportion_het / 2 - alt_allele_frequency_grid)

    log("\t_compute_hardy_weinberg_plotted_boundary() ...")

    # don't know proportion_het at the boundary (that's what we're computing) so iterate
    # start with low boundary, guess midway between the lowest possible and expected values
    proportion_het_low = alt_allele_frequency_grid * (1.0 - alt_allele_frequency_grid)
    previous = None
    tol = 0.1 * numpy.diff(proportion_het_grid).mean()  # set tol to about 10% of a grid width
    while previous is None or numpy.abs(proportion_het_low - previous).max(initial=0) >= tol:
        # to aid convergence, only take a half step each iteration:
        # noinspection PyUnresolvedReferences
        num_expected_het, het_delta = _hardy_weinberg_chi_squared_to_num_hets_boundary(
            hardy_weinberg_confidence=hardy_weinberg_confidence,
            alt_allele_frequency=alt_allele_frequency_grid,
            num_ref_genotypes=_num_ref_genotypes(
                proportion_het_low if previous is None else (previous + proportion_het_low) / 2
            ),
            num_called_genotypes=num_samples
        )
        previous = proportion_het_low
        proportion_het_low = numpy.minimum(numpy.maximum(num_expected_het - het_delta, 0) / num_samples, 1)

    # now do proportion_het_high, start midway between the highest possible and expected values
    proportion_het_high = (
        2 * alt_allele_frequency_grid * (1.0 - alt_allele_frequency_grid) +
        1.0 - numpy.abs(2 * alt_allele_frequency_grid - 1)
    ) / 2.0
    previous = None
    while previous is None or numpy.abs(proportion_het_high - previous).max() >= tol:
        # to aid convergence, only take a half step each iteration:
        # noinspection PyTypeChecker
        num_expected_het, het_delta = _hardy_weinberg_chi_squared_to_num_hets_boundary(
            hardy_weinberg_confidence=hardy_weinberg_confidence,
            alt_allele_frequency=alt_allele_frequency_grid,
            num_ref_genotypes=_num_ref_genotypes(
                proportion_het_high if previous is None else (previous + proportion_het_high) / 2
            ),
            num_called_genotypes=num_samples
        )
        previous = proportion_het_high
        proportion_het_high = numpy.minimum(numpy.maximum(num_expected_het + het_delta, 0) / num_samples, 1)

    log("\tCOMPLETE: _compute_hardy_weinberg_plotted_boundary()")
    return proportion_het_low, proportion_het_high


def plot_hardy_weinberg_proportions(
        ax: pyplot.Axes,
        scores_data_set: ScoresDataSet,
        title: str,
        title_font_size: int = Default.title_font_size,
        label_font_size: int = Default.label_font_size,
        tick_font_size: int = Default.tick_font_size,
        colorbar_axes: Optional[pyplot.Axes] = None,
        legend_axes: Optional[pyplot.Axes] = None,
        hardy_weinberg_confidence: float = Default.hardy_weinberg_confidence,
        min_num_called_gt_hardy_weinberg: int = Default.min_num_called_gt_hardy_weinberg
):
    log(f"\tplotting {title}")

    num_het_genotypes = scores_data_set.num_het_genotypes
    num_called_genotypes = scores_data_set.num_called_genotypes
    alt_allele_frequency = scores_data_set.num_variant_alleles / (2 * num_called_genotypes)
    het_proportion = num_het_genotypes / num_called_genotypes
    __, histogram_mat, x_edges, y_edges = plotting.hist2d(
        ax=ax, x=alt_allele_frequency, y=het_proportion,
        xlabel="alt allele frequency", ylabel="proportion het", title=title, cmap="Reds", xlim=(0, 1), ylim=(0, 1),
        mask_condition=lambda hist_mat: hist_mat == 0, mask_color='w', colorbar=colorbar_axes is not None,
        colorbar_axes=colorbar_axes, zscale="log", zlabel="variant counts",
        title_fontsize=title_font_size, label_fontsize=label_font_size, tick_fontsize=tick_font_size,
        x_bins=50, y_bins=50
    )

    proportion_accepted = _compute_hardy_weinberg_proportion_accepted(
        num_het_genotypes=num_het_genotypes.values, num_ref_genotypes=scores_data_set.num_ref_genotypes.values,
        alt_allele_frequency=alt_allele_frequency, num_called_genotypes=num_called_genotypes.values,
        hardy_weinberg_confidence=hardy_weinberg_confidence
    )

    proportion_het_low, proportion_het_high = _compute_hardy_weinberg_plotted_boundary(
        proportion_het_grid=y_edges, alt_allele_frequency_grid=x_edges,
        num_samples=scores_data_set.num_samples, hardy_weinberg_confidence=hardy_weinberg_confidence
    )

    ax.plot(x_edges, proportion_het_low, ":b", label="acceptance region")
    ax.plot(x_edges, proportion_het_high, ":b", label="acceptance region")
    ax.text(0.99, 0.99, f"{100 * proportion_accepted:.1f}% accepted",
            horizontalalignment="right", verticalalignment="top", fontsize=label_font_size)

    if legend_axes is not None:
        legend_axes.plot([], [], "b:", label="acceptance region")
        legend_axes.legend(fontsize=label_font_size, loc="upper left")


def get_filter_quality_figures(
        truth_data: TruthData,
        scores_data_sets: Sequence[ScoresDataSet],
        sv_selectors: Mapping[str, SvTypeCutoffSelector],
        num_precision_recall_curve_rows: int = Default.num_precision_recall_curve_rows,
        num_inheritance_stats_rows: int = Default.num_inheritance_stats_rows,
        title_font_size: int = Default.title_font_size,
        label_font_size: int = Default.label_font_size,
        legend_font_size: int = Default.legend_font_size,
        tick_font_size: int = Default.tick_font_size,
        tick_labelrotation: float = Default.tick_labelrotation,
        plot_max_f1: bool = Default.plot_max_f1,
        plotted_mendelian_violation_keys: Optional[Iterable[str]] = Default.plotted_mendelian_violation_trio_types,
        plotted_hardy_weinberg_categories: Optional[Container[str]] = Default.plotted_hardy_weinberg_categories,
        min_num_called_gt_hardy_weinberg: int = Default.min_num_called_gt_hardy_weinberg,
        hardy_weinberg_confidence: float = Default.hardy_weinberg_confidence,
        make_figures: Set[str] = Default.make_figures,
        pdf_writer: Optional[PdfPages] = None
) -> Tuple[Optional[pyplot.Figure], ...]:
    if Keys.precision_recall in make_figures:
        precision_recall_fig, best_thresholds = get_precision_recall_curves_figure(
            truth_data=truth_data,
            scores_data_sets=scores_data_sets,
            sv_selectors=sv_selectors,
            num_precision_recall_curve_rows=num_precision_recall_curve_rows,
            title_font_size=title_font_size,
            label_font_size=label_font_size,
            legend_font_size=legend_font_size,
            tick_font_size=tick_font_size,
            plot_max_f1=plot_max_f1,
            pdf_writer=pdf_writer
        )
    else:
        precision_recall_fig = None
        if plot_max_f1:
            log("get max-f1 thresholds ...")
            # need to calculate the max_f1_thresholds
            # in this case, only need the non-bnd PRC, because we just want a point off of it
            precision_recall_curves = get_precision_recall_curves(
                truth_data=truth_data, scores_data_sets=scores_data_sets,
                sv_selectors={Keys.all_variant_types: sv_selectors[Keys.non_bnd_types]}
            )
            # get the overall best thresholds
            best_thresholds = {
                data_set_label: precision_recall_curve.get_point_at_max_f_score().name
                for data_set_label, precision_recall_curve in
                precision_recall_curves.pop(Keys.all_variant_types).items()
            }
            log("COMPLETE: get max-f1 thresholds")
        else:
            best_thresholds = {}

    if Keys.inheritance in make_figures:
        mendelian_inheritance_fig = \
            get_mendelian_inheritance_figure(
                truth_data=truth_data,
                scores_data_sets=scores_data_sets,
                sv_selectors=sv_selectors,
                best_thresholds=best_thresholds,
                num_inheritance_stats_rows=num_inheritance_stats_rows,
                title_font_size=title_font_size,
                label_font_size=label_font_size,
                legend_font_size=legend_font_size,
                tick_font_size=tick_font_size,
                tick_labelrotation=tick_labelrotation,
                plot_max_f1=plot_max_f1,
                pdf_writer=pdf_writer
            )
    else:
        mendelian_inheritance_fig = None

    if Keys.variants_per_sample in make_figures:
        variants_per_sample_fig = \
            get_variants_per_sample_figure(
                scores_data_sets=scores_data_sets,
                sv_selectors=sv_selectors,
                best_thresholds=best_thresholds,
                num_inheritance_stats_rows=num_inheritance_stats_rows,
                title_font_size=title_font_size,
                label_font_size=label_font_size,
                legend_font_size=legend_font_size,
                tick_font_size=tick_font_size,
                plot_max_f1=plot_max_f1,
                pdf_writer=pdf_writer
            )
    else:
        variants_per_sample_fig = None

    if Keys.violation_curve in make_figures:
        mendelian_violation_curve_fig = get_mendelian_violation_curve_figure(
            truth_data=truth_data,
            scores_data_sets=scores_data_sets,
            sv_selectors=sv_selectors,
            best_thresholds=best_thresholds,
            num_inheritance_stats_rows=num_inheritance_stats_rows,

            title_font_size=title_font_size,
            label_font_size=label_font_size,
            legend_font_size=legend_font_size,
            tick_font_size=tick_font_size,
            plot_max_f1=plot_max_f1,
            plotted_mendelian_violation_keys=plotted_mendelian_violation_keys,
            pdf_writer=pdf_writer
        )
    else:
        mendelian_violation_curve_fig = None

    if Keys.scores_histogram in make_figures:
        scores_histogram_fig = get_scores_histogram_figure(
            scores_data_sets=scores_data_sets,
            sv_selectors=sv_selectors,
            title_font_size=title_font_size,
            label_font_size=label_font_size,
            tick_font_size=tick_font_size,
            tick_labelrotation=tick_labelrotation,
            pdf_writer=pdf_writer
        )
    else:
        scores_histogram_fig = None

    if Keys.hardy_weinberg in make_figures:
        hardy_weinberg_figs = get_hardy_weinberg_figures(
            scores_data_sets=scores_data_sets,
            sv_selectors=sv_selectors,
            title_font_size=title_font_size,
            label_font_size=label_font_size,
            tick_font_size=tick_font_size,
            plotted_hardy_weinberg_categories=plotted_hardy_weinberg_categories,
            min_num_called_gt_hardy_weinberg=min_num_called_gt_hardy_weinberg,
            hardy_weinberg_confidence=hardy_weinberg_confidence,
            pdf_writer=pdf_writer
        )
    else:
        hardy_weinberg_figs = (None,)

    return precision_recall_fig, mendelian_inheritance_fig, mendelian_violation_curve_fig, variants_per_sample_fig, \
        scores_histogram_fig, *hardy_weinberg_figs


def load_quality_data(
        overlap_results_file: str,
        pedigree_files: Optional[Collection[str]],
        optimal_overlap_cutoffs_file: str,
        scores_data_json: Optional[str] = None,
        size_ranges: Mapping[str, Tuple[float, float]] = Default.sv_selector_size_ranges
) -> (TruthData, Dict[str, SvTypeCutoffSelector], List[ScoresDataSet]):
    truth_data = TruthData(overlap_results_file=overlap_results_file, pedigree_files=pedigree_files)
    scores_data_sets = ScoresDataSet.from_json(
        scores_data_json, truth_data=truth_data
    )
    all_sv_types = set.union(*(data_set.all_sv_types for data_set in scores_data_sets))
    sv_selectors = get_truth_overlap.get_sv_selectors(all_sv_types=all_sv_types, size_ranges=size_ranges) \
        if optimal_overlap_cutoffs_file is None \
        else {
            sv_category: cutoff.selector
            for sv_category, cutoff in get_truth_overlap.load_optimal_overlap_cutoffs(
                optimal_overlap_cutoffs_file
            ).items()
        }  # noqa: E123

    return truth_data, sv_selectors, scores_data_sets


class SplitArgs(argparse.Action):
    """ Helper class to split (potentially) comma-separated input str values and store as a set """
    def __call__(self, parser, namespace, values, option_string=None):
        new_values = values.split(',')
        values = getattr(namespace, self.dest, set())
        if isinstance(values, frozenset):
            # this is the default, which we're overridding
            setattr(namespace, self.dest, set(new_values))
        else:
            # updating because user supplied the argument more than once
            setattr(namespace, self.dest, values.union(new_values))


def __parse_arguments(argv: List[Text]) -> argparse.Namespace:

    # noinspection PyTypeChecker
    parser = argparse.ArgumentParser(
        description="Plot precision/recall curve and inheritance stats for trained variant filters",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        prog=argv[0]
    )
    parser.add_argument("--overlap-results", '-o', type=str, required=True,
                        help="file with saved mapping from sample ID to lists of good and bad variant IDs, as "
                             "determined by get_truth_overlap and overlaps between test and truth VCFs or VaPoR")
    parser.add_argument("--optimal-overlap-cutoffs-file", "-c", type=str, required=False,
                        help="file with pickled optimal overlap cutoffs, produced by get_truth_overlap")
    parser.add_argument("--ped-file", "-p", type=str, action="extend", nargs="+", required=False,
                        help="pedigree file for samples")
    parser.add_argument("--scores-data-json", "-s", type=str, nargs='?', required=True,
                        help="path to json file with object mapping from data set label to object with keys:\n"
                             f"{Keys.scores_file}: array of paths to files with variant scores\n"
                             f"{Keys.metrics_file}: array of paths to files with variant metrics\n"
                             f"{Keys.passing_score}: minimum score considered passing by default")
    parser.add_argument("--figure-save-file", "-f", type=str, required=True, help="path to save output .pdf")
    parser.add_argument("--plot-max-f1", action="store_true", help="Add point at maximum f1 score to plots")
    # noinspection PyTypeChecker
    parser.add_argument("--make-figures", choices=_all_make_figures, default=Default.make_figures, type=str,
                        action=SplitArgs, help="Comma-separated list of which figures to make.")
    parser.add_argument("--size-ranges", type=str, default=_size_ranges_to_arg(Default.sv_selector_size_ranges),
                        help="comma-separated list of size ranges to break down summary of SVs. Each size range is of "
                             "form description:low_size:high_size")

    parsed_arguments = parser.parse_args(argv[1:] if len(argv) > 1 else ["--help"])

    if not parsed_arguments.figure_save_file.lower().endswith(".pdf"):
        raise ValueError("--figure-save-file must be a .pdf file")
    return parsed_arguments


def _parse_size_ranges(size_ranges_arg: str) -> dict[str, tuple[float, float]]:
    def _parse_arg(_arg: str) -> Iterator[tuple[str, tuple[float, float]]]:
        for _word in _arg.split(','):
            _split_word = _word.split(':', 1)
            yield _split_word[0], tuple(float(_x) for _x in _split_word[1].split(':', 1))

    return {description: size_range for description, size_range in _parse_arg(size_ranges_arg)}


def _size_ranges_to_arg(size_ranges: dict[str, tuple[float, float]]) -> str:
    return ','.join(f"{key}:{val[0]}-{val[1]}" for key, val in size_ranges.items())


def main(argv: Optional[List[Text]] = None) -> (pyplot.Figure, pyplot.Figure):
    arguments = __parse_arguments(sys.argv if argv is None else argv)
    truth_data, sv_selectors, scores_data_sets = load_quality_data(
        overlap_results_file=arguments.overlap_results,
        pedigree_files=arguments.ped_file,
        optimal_overlap_cutoffs_file=arguments.optimal_overlap_cutoffs_file,
        scores_data_json=arguments.scores_data_json,
        size_ranges=_parse_size_ranges(arguments.size_ranges)
    )

    with PdfPages(arguments.figure_save_file) as pdf_writer:
        get_filter_quality_figures(
            truth_data=truth_data, sv_selectors=sv_selectors, scores_data_sets=scores_data_sets,
            plot_max_f1=arguments.plot_max_f1, make_figures=set(arguments.make_figures), pdf_writer=pdf_writer
        )


if __name__ == "__main__":
    main()
