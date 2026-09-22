#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""

"""

import numpy as np


class TrainingLabeler:
    def __init__(self):
        pass

    def label(self, metrics):
        return metrics.apply(self.label_row, axis=1)

    def label_row(self, row):
        return 'Unlabeled'


class BAF1TrainingLabeler(TrainingLabeler):
    def label_row(self, row):
        if 0 <= row.RD_MEDIAN_SEPARATION < 0.1:
            return 'Fail'
        elif 0.4 <= row.RD_MEDIAN_SEPARATION < 1.0:
            return 'Pass'
        else:
            return 'Unlabeled'


class SR1TrainingLabeler(TrainingLabeler):
    def label_row(self, row):
        if (row.RD_MEDIAN_SEPARATION < 0.15 and
                row.BAF1_prob < 0.4 and
                row.PEQ < -10 * np.log10(0.05)):
            return 'Fail'
        elif row.RD_MEDIAN_SEPARATION >= 0.4 and row.BAF1_prob >= 0.9:
            return 'Pass'
        else:
            return 'Unlabeled'


class RDTrainingLabeler(TrainingLabeler):
    def label_row(self, row):
        if 'depth' not in row['name'] and row.svsize >= 1000:
            if row.BAF1_prob < 0.4 and row.SR1_prob < 0.4:
                return 'Fail'
            elif row.BAF1_prob >= 0.9 and row.SR1_prob >= 0.9:
                return 'Pass'
            else:
                return 'Unlabeled'
        elif 'depth' not in row['name'] and row.svsize < 1000:
            if row.SR1_prob < 0.4:
                return 'Fail'
            elif row.SR1_prob >= 0.9:
                return 'Pass'
            else:
                return 'Unlabeled'
        else:
            if row.BAF1_prob < 0.4:
                return 'Fail'
            elif row.BAF1_prob >= 0.9:
                return 'Pass'
            else:
                return 'Unlabeled'


class PETrainingLabeler(TrainingLabeler):
    def label_row(self, row):
        if row.SR1_prob < 0.4 and row.RD_prob < 0.4 and row.svsize >= 1000:
            return 'Fail'
        elif row.SR1_prob >= 0.9 and row.RD_prob >= 0.9 and row.svsize >= 1000:
            return 'Pass'
        else:
            return 'Unlabeled'


class BAF2TrainingLabeler(TrainingLabeler):
    def label_row(self, row):
        if (row.RD_prob < 0.4 and row.PE_prob < 0.4 and
                row.SR1_prob < 0.4 and row.BAF1_prob < 0.4):
            return 'Fail'
        elif (row.RD_prob >= 0.9 and row.PE_prob >= 0.9 and
                row.SR1_prob >= 0.9 and row.BAF1_prob >= 0.4):
            return 'Pass'
        else:
            return 'Unlabeled'


class SR2TrainingLabeler(TrainingLabeler):
    def label_row(self, row):
        if (row.RD_prob < 0.4 or row.PE_prob < 0.4):
            return 'Fail'
        elif (row.RD_prob >= 0.9 and row.PE_prob >= 0.9 and
                row.SR1_prob >= 0.4):
            return 'Pass'
        else:
            return 'Unlabeled'


class PESRTrainingLabeler(TrainingLabeler):
    def label_row(self, row):
        if (row.RD_prob < 0.4 and row.PE_prob < 0.4 and
                row.SR1_prob < 0.4):
            return 'Fail'
        elif (row.RD_prob >= 0.9 and row.PE_prob >= 0.9 and
                row.SR1_prob >= 0.9):
            return 'Pass'
        else:
            return 'Unlabeled'
