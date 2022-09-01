#!/usr/bin/env python


from distutils.core import setup
from setuptools import find_packages
from glob import glob
import os


setup(
    name="sv_utils",
    version="1.0.0",
    description="Utility and machine-learning scripts for processing structural variant data",
    author="Ted Brookings",
    author_email="tbrookin@broadinstitute.org",
    package_dir={"": "src"},
    packages=find_packages("src"),
    py_modules=[os.path.splitext(os.path.basename(path))[0]
                for path in glob("src/*.py")],
    entry_points={
        "console_scripts": [
            "sv-utils=sv_utils.command_line:main",
            "sv_utils=sv_utils.command_line:main",
            "gq-recalibrator=sv_utils.command_line:main",
            "gq_recalibrator=sv_utils.command_line:main"
        ]
    },
    python_requires=">=3.10.4",
    install_requires=[
        "attrs>=22.1.0",
        "dask>=2023.7.0",
        "dill>=0.3.5.1",
        "matplotlib>=3.5.1",
        "numpy>=1.22.4",
        "pandas>=1.5.3",
        "psutil>=5.8.0",
        "pyarrow>=12.0.0",
        "pympler>=0.9",
        "pysam",
        "scipy>=1.10.1",
        "seaborn>=0.11.2",
        "tqdm>=4.62.3",
    ],
    extras_require={
        "tests": ["pytest>=7.1.2", "pytest-cov>=3.0.0"],
        "gq-recalibrator": ["pytorch>=1.13.1"],
        "gq-recalibrator-gpu":  []
    },
    include_package_data=True,
    zip_safe=False
)
