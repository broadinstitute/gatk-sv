"""
gatk-sv-gd: Genomic Disorder CNV detection from binned read counts.

Modules:
    preprocess - Preprocess read-depth data for GD CNV inference
    infer      - Bayesian CNV inference at GD loci (Pyro model)
    call       - Call GD CNVs from model posterior probabilities
    plot       - Generate visualisation plots for GD CNV calls
    eval       - Evaluate GD CNV calls against a truth table
    models     - Shared data structures (GDLocus, GDTable)
"""

__version__ = "0.1.0"

from gatk_sv_gd.models import GDLocus, GDTable

__all__ = ["GDLocus", "GDTable", "__version__"]
