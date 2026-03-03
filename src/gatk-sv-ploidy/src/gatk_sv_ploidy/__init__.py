"""
gatk-sv-ploidy: Whole-genome aneuploidy detection from binned read counts.

Modules
-------
preprocess  Normalise and filter depth data
infer       Bayesian CN inference (Pyro model)
call        Assign sex karyotype and aneuploidy type
plot        Generate diagnostic and summary plots
eval        Evaluate predictions against a truth set
"""

__version__ = "0.1.0"

from gatk_sv_ploidy.data import DepthData

__all__ = ["DepthData", "__version__"]
