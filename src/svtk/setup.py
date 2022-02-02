#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""

"""

from setuptools import setup, find_packages
from Cython.Build import cythonize
import pysam

setup(
    name='svtk',
    version='0.1',
    description='Structural variation toolkit',
    author='Matthew Stone',
    author_email='mstone5@mgh.harvard.edu',
    packages=find_packages(),
    package_data={'svtk': ['data/*_template.vcf']},
    scripts=['scripts/svtk'],
    ext_modules=cythonize('svtk/utils/helpers.pyx'),
    include_dirs=pysam.get_include(),
    install_requires=[
        'numpy',
        'scipy',
        'pysam>=0.11.2.2',
        'pybedtools',
        'cython',
        'natsort',
        'pandas',
    ]
)
