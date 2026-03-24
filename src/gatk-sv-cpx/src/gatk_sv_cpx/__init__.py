"""
gatk-sv-cpx: Consolidated complex SV resolution, evidence evaluation, and
genotype refinement for the GATK-SV pipeline.

Subcommands
-----------
resolve              Identify and annotate complex (CPX) SV linkages from
                     clustered breakpoints, replacing RCV.
evaluate-evidence    Gather PE and RD evidence for putative CPX/CTX variants,
                     replacing the fragmented evidence collection in GCV and
                     RefCV.
genotype-and-refine  Apply depth-based genotyping and PE/RD refinement in a
                     single pass, replacing GCV genotyping + RefCV revision.
"""

__version__ = "0.1.0"
