#!/usr/bin/env python

from setuptools import setup

setup(name='svqc',
      version='0.1',
      description='QC package for the GATK SV pipeline',
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Programming Language :: Python :: 2.7',
      ],
      url='https://github.com/talkowski-lab/gatk-sv-v1',
      author='Mark Walker',
      author_email='markw@broadinsitute.org',
      packages=['svqc'],
      include_package_data=True,
      zip_safe=False,
      entry_points={
          'console_scripts': ['svqc=svqc.command_line:main'],
      },
      test_suite='nose.collector',
      tests_require=['nose'],
      install_requires=['pandas']
      )
