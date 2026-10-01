# TODO: this import is the reason svtk needs an old setuptools. Recent setuptools
# releases dropped pkg_resources, so every svtk command dies on this line unless
# something holds setuptools back. On this branch the only guard is the cap in
# dockerfiles/sv-pipeline/Dockerfile; PR #959 adds 'setuptools<81' to
# install_requires on main. Replacing this import removes the constraint instead
# of renewing it: importlib.metadata for the version below, importlib.resources
# for the resource_filename() call sites in cli/rdtest2vcf.py,
# cli/standardize_vcf.py, utils/rdtest.py and vcfcluster.py.
from pkg_resources import get_distribution

__all__ = []

__version__ = get_distribution('svtk').version
