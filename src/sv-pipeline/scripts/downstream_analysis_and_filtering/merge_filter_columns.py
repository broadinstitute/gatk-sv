#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Sanitize two filter columns stripped from paired VCFs
"""

import argparse


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('file1', help='Input FILTER values from VCF 1.')
    parser.add_argument('file2', help='Input FILTER values from VCF 2.')
    parser.add_argument('fout', help='Output FILTER values.')

    args = parser.parse_args()

    fout = open(args.fout, 'w')

    with open(args.file1) as f1, open(args.file2) as f2:
        for x, y in zip(f1, f2):
            x = x.strip().split(';')
            y = y.strip().split(';')
            # Only return PASS if both are PASS with no other filters
            if x == ['PASS'] and y == ['PASS']:
                newfilt = 'PASS'
            else:
                x = [f for f in x if f != 'PASS']
                y = [f for f in y if f != 'PASS']
                newfilt = ';'.join(sorted(list(set(x + y))))

            fout.write(newfilt + '\n')

    fout.close()


if __name__ == '__main__':
    main()
