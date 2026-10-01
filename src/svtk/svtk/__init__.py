# TODO: this import is the reason svtk needs an old setuptools. Recent setuptools
# releases dropped pkg_resources, so every svtk command dies on this line unless
# something holds setuptools back. That is 'setuptools<81' in setup.py
# install_requires, plus the cap in dockerfiles/sv-pipeline/Dockerfile. Replacing
# it removes the constraint instead of renewing it: importlib.metadata for the
# version below, importlib.resources for the resource_filename() call sites in
# cli/rdtest2vcf.py, cli/standardize_vcf.py, utils/rdtest.py and vcfcluster.py.
# See PR #959, which put the cap in install_requires.
from pkg_resources import get_distribution

__all__ = []

__version__ = get_distribution('svtk').version
