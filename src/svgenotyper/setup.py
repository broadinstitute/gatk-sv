#!/usr/bin/env python

from setuptools import setup

setup(name='svgenotyper',
      version='0.1',
      description='GATK-SV genotyping package',
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Programming Language :: Python :: 3.8',
      ],
      url='https://github.com/talkowski-lab/gatk-sv-v1',
      author='Mark Walker',
      author_email='markw@broadinsitute.org',
      packages=['svgenotyper'],
      include_package_data=True,
      zip_safe=False,
      scripts=['scripts/svgenotyper'],
      install_requires=['numpy', 'pyro-ppl', 'torch'])
