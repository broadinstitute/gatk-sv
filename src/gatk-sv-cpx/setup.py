#!/usr/bin/env python

from setuptools import setup, find_packages
from glob import glob
import os

setup(
    name="gatk-sv-cpx",
    version="0.1.0",
    description="Consolidated complex SV resolution, evidence evaluation, and genotype refinement",
    author="GATK-SV Team",
    package_dir={"": "src"},
    packages=find_packages("src"),
    py_modules=[
        os.path.splitext(os.path.basename(path))[0]
        for path in glob("src/*.py")
    ],
    entry_points={
        "console_scripts": [
            "gatk-sv-cpx=gatk_sv_cpx.cli:main",
        ]
    },
    python_requires=">=3.9",
    install_requires=[
        "pysam",
        "numpy",
        "intervaltree",
    ],
    extras_require={
        "dev": ["pytest", "pytest-cov"],
    },
    include_package_data=True,
    zip_safe=False,
)
