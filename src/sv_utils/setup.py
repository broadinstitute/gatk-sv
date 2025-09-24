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
            "sv_utils=sv_utils.command_line:main"
        ]
    },
    python_requires=">3.8",
    install_requires=[
        "numpy",
        "scipy",
        "pandas",
        "tqdm",
        "psutil",
        "pysam>=0.23.3",
        "dill",
        "pympler",
        "matplotlib",
        "seaborn"
    ],
    extras_require={
        "tests": ["pytest", "pytest-cov"]
    },
    include_package_data=True,
    zip_safe=False
)
