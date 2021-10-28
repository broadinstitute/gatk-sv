#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Simple classes for parsing and interpreting .fam files
"""

from collections import defaultdict


class Fam:
    def __init__(self, samples, families):
        self.samples = samples
        self.families = families


class FamSample:
    def __init__(self, family, ID, father, mother, sex, phenotype):
        self.family = family
        self.ID = ID
        self.father = father if father != '0' else None
        self.mother = mother if mother != '0' else None
        self.sex = sex
        self.phenotype = phenotype
        self.is_parent = False
        self.children = []

    @property
    def is_male(self):
        return self.sex == 1

    @property
    def is_female(self):
        return self.sex == 2

    @property
    def has_parents(self):
        return self.father is not None and self.mother is not None


def parse_famfile(famfile):
    """
    Parse a fam file.

    Returns
    -------
    fam : dict of {str: FamSample}
        Mapping of sample IDs to FamSample entries
    """

    samples = {}
    families = defaultdict(list)

    for line in famfile:
        data = line.strip().split()

        sample = FamSample(*data)
        samples[data[1]] = sample
        families[data[0]].append(sample)

    # Label parents
    for ID, sample in samples.items():
        if sample.father is not None:
            samples[sample.father].is_parent = True
            samples[sample.father].children.append(ID)
        if sample.mother is not None:
            samples[sample.mother].is_parent = True
            samples[sample.mother].children.append(ID)

    return Fam(samples, families)
