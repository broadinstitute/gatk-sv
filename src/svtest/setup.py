#!/usr/bin/env python

from setuptools import setup

setup(name='svtest',
      version='0.1',
      description='Test package for the GATK SV pipeline',
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Programming Language :: Python :: 3.8',
      ],
      url='https://github.com/talkowski-lab/gatk-sv-v1',
      author='Mark Walker',
      author_email='markw@broadinsitute.org',
      packages=['svtest'],
      include_package_data=True,
      zip_safe=False,
      scripts=['scripts/svtest'],
      install_requires=['numpy', 'matplotlib', 'pandas', 'intervaltree', 'pysam>=0.23.3'])
