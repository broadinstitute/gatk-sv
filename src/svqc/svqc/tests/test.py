#!/usr/bin/env python

from unittest import TestCase
from svqc.evaluator import QCEvaluator
import svqc.tests.constants as const
import os


class TestSVQC(TestCase):

    def tearDown(self):
        os.remove(const.OUT_FILE_PATH)

    def integration_test(self):
        eval = QCEvaluator(const.METRICS_FILE_PATH, const.CRITERIA_FILE_PATH)
        tests = eval.get_tests()
        self.assertIsNot(tests, None)
        self.assertEqual(len(tests), const.EXPECTED_NUM_TESTS)
        eval.write_tests(const.OUT_FILE_PATH, tests)
