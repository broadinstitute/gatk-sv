#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""

"""

import sys
import itertools


def make_link_combos():
    pass


def main():
    # parser = argparse.ArgumentParser(
    #     description=__doc__,
    #     formatter_class=argparse.RawDescriptionHelpFormatter)
    # args = parser.parse_args()

    for line in sys.stdin:
        links = line.strip().split(',')
        for pair in itertools.combinations(links, 2):
            sys.stdout.write('{0}\t{1}\n'.format(*pair))


if __name__ == '__main__':
    main()
