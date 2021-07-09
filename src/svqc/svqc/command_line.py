#!/usr/bin/env python

import argparse
import sys
from svqc.evaluator import QCEvaluator


def main():
    usage = """Detects QC failures"""
    parser = argparse.ArgumentParser(usage=usage)
    parser.add_argument("metrics_file", help="Table of metrics")
    parser.add_argument("criteria_file", help="Table of QC criteria")
    parser.add_argument("output", help="Output file")
    parser.add_argument(
        "--verbose", help="Enable logging to stderr", action='store_true')
    args = parser.parse_args()

    eval = QCEvaluator(args.metrics_file, args.criteria_file)
    tests = eval.get_tests()
    eval.write_tests(args.output, tests)

    if args.verbose:
        for test in tests:
            sys.stderr.write(test.pretty_string())


if __name__ == "__main__":
    main()
