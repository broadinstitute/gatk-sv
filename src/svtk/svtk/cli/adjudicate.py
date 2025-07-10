#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Random-forest adjudication of metrics supporting putative SV calls.
"""

import argparse
import sys
import pandas as pd
from svtk.adjudicate import adjudicate_SV


def main(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        prog='svtk adjudicate',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('metrics')
    parser.add_argument('scores', type=argparse.FileType('w'),
                        help='Final variant scores.')
    parser.add_argument('cutoffs', type=argparse.FileType('w'),
                        help='Learned cutoffs.')
    parser.add_argument('--remove-wham', action='store_true',
                        help='Remove WHAM deletions from metrics table.')

    # Print help if no arguments specified
    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    metrics = pd.read_table(args.metrics)
    scores, cutoffs = adjudicate_SV(metrics, remove_wham=args.remove_wham)
    scores.to_csv(args.scores, index=False, sep='\t')
    cutoffs.to_csv(args.cutoffs, index=False, sep='\t')


if __name__ == '__main__':
    main()
