#!/usr/bin/env python

######################################################
#
# Sample QC evaluation for the GATK SV pipeline
# written by Mark Walker (markw@broadinstitute.org)
#
######################################################

import pandas as pd
import svqc.constants as const

# Tests whether numbers are on a specified interval (inclusive)


class Criterion:
    def __init__(self, name, min=None, max=None):
        self.name = name
        self.min = min
        self.max = max

    def test(self, value):
        if not self.is_testable():
            raise ValueError(
                "Criterion %s could not be tested because it does not have assigned bounds" % self.name)
        return self.min <= value <= self.max

    def __repr__(self):
        return "[" + self.name + ":" + self.min_string() + "," + self.max_string() + "]"

    def min_string(self):
        if self.min is None:
            return const.NA_STRING
        return str(self.min)

    def max_string(self):
        if self.max is None:
            return const.NA_STRING
        return str(self.max)

    def is_testable(self):
        return self.min is not None and self.max is not None


# Determines the result of testing a given value against a given Criterion
# Defines behavior for missing values or missing criterion
class Test:
    def __init__(self, criterion, value=None):
        self.value = value
        self.criterion = criterion
        if not self.is_testable():
            self.result = const.NA_STRING
        else:
            if self.criterion.test(self.value):
                self.result = const.PASS_STRING
            else:
                self.result = const.FAIL_STRING

    def __repr__(self):
        return "(" + self.get_value_string() + "," + str(self.criterion) + "," + self.result + ")"

    def is_testable(self):
        return self.value is not None and self.criterion.is_testable()

    def get_value_string(self):
        if self.value is None:
            return const.NA_STRING
        return str(self.value)

    def pretty_string(self):
        return "Metric \"%s\" = %s evaluated against range [%s, %s] and resulted in %s\n" \
               % (self.criterion.name,
                  self.get_value_string(),
                  self.criterion.min_string(),
                  self.criterion.max_string(),
                  self.result)


# Runs tests on given metrics against a given set of criteria
class QCEvaluator:
    def __init__(self, metrics_file_path, criteria_file_path):
        self._load_metrics(metrics_file_path)
        self._load_criteria(criteria_file_path)

    def _load_metrics(self, filename):
        cols = [const.METRIC_FILE_NAME_COLUMN, const.METRIC_FILE_VALUE_COLUMN]
        self.data = pd.read_csv(
            filename, sep=const.TSV_DELIMITER, names=cols, header=None)

    def _load_criteria(self, filename):
        cols = [const.CRITERIA_METRIC_COLUMN,
                const.CRITERIA_MIN_COLUMN, const.CRITERIA_MAX_COLUMN]
        raw_data = pd.read_csv(
            filename, sep=const.TSV_DELIMITER, names=cols, header=0)
        self.criteria = {}
        for i in range(raw_data.shape[0]):
            name = raw_data.at[i, const.CRITERIA_METRIC_COLUMN]
            min = raw_data.at[i, const.CRITERIA_MIN_COLUMN]
            max = raw_data.at[i, const.CRITERIA_MAX_COLUMN]
            self.criteria[name] = Criterion(name, min=min, max=max)

    def get_tests(self):
        tests = {}
        name_idx = self.data.columns.get_loc(const.METRIC_FILE_NAME_COLUMN)
        val_idx = self.data.columns.get_loc(const.METRIC_FILE_VALUE_COLUMN)
        for row in self.data.iterrows():
            name = row[1][name_idx]
            value = row[1][val_idx]
            if name not in self.criteria:
                tests[name] = Test(Criterion(name), value=value)
            else:
                criterion = self.criteria[name]
                tests[name] = Test(criterion, value=value)
        # Criteria that were defined but there were no metrics for
        missed_tests = [Test(self.criteria[x])
                        for x in self.criteria if x not in tests]
        for test in missed_tests:
            name = test.criterion.name
            tests[name] = test
        return tests.values()

    def write_tests(self, filename, tests):
        with open(filename, 'w') as f:
            f.write("#METRIC\tMIN\tMAX\tVALUE\tRESULT\n")
            f.writelines(["\t".join([x.criterion.name, x.criterion.min_string(
            ), x.criterion.max_string(), x.get_value_string(), x.result]) + "\n" for x in tests])
