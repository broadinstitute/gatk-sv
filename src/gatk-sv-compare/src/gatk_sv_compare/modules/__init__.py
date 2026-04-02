"""Analysis module interfaces and registry helpers."""

from .allele_freq import AlleleFreqModule
from .base import AnalysisModule
from .binned_counts import BinnedCountsModule
from .counts_per_genome import CountsPerGenomeModule
from .family_analysis import FamilyAnalysisModule
from .genotype_dist import GenotypeDistModule
from .genotype_concordance import GenotypeConcordanceModule
from .genotype_exact_match import GenotypeExactMatchModule
from .genotype_quality import GenotypeQualityModule
from .overall_counts import OverallCountsModule
from .site_overlap import SiteOverlapModule
from .size_signatures import SizeSignaturesModule

SITE_LEVEL_MODULES = [
    BinnedCountsModule,
    OverallCountsModule,
    SiteOverlapModule,
    AlleleFreqModule,
    GenotypeDistModule,
    GenotypeQualityModule,
    CountsPerGenomeModule,
    SizeSignaturesModule,
]

GENOTYPE_LEVEL_MODULES = [
    GenotypeExactMatchModule,
    GenotypeConcordanceModule,
    FamilyAnalysisModule,
]

ALL_MODULES = SITE_LEVEL_MODULES + GENOTYPE_LEVEL_MODULES

__all__ = [
    "ALL_MODULES",
    "SITE_LEVEL_MODULES",
    "AnalysisModule",
    "AlleleFreqModule",
    "BinnedCountsModule",
    "CountsPerGenomeModule",
    "FamilyAnalysisModule",
    "GenotypeDistModule",
    "GenotypeConcordanceModule",
    "GenotypeExactMatchModule",
    "GenotypeQualityModule",
    "GENOTYPE_LEVEL_MODULES",
    "OverallCountsModule",
    "SiteOverlapModule",
    "SizeSignaturesModule",
]

