#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""

"""

import argparse
import numpy as np
import scipy.stats as ss


def convert_poisson_p(qual):
    """
    Convert poisson phred-scaled quality score to count, assuming 0 background

    Can't use ss.poisson.ppf because k/mu are flipped in pe/sr test calculation
    """
    count = 0
    while True:
        pval_cmp = -10. * np.log10(ss.poisson.cdf(0, count))

        if pval_cmp > qual:
            return max(1, count - 1)

        count += 1


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('qual', type=float)
    args = parser.parse_args()

    print(convert_poisson_p(args.qual))


if __name__ == '__main__':
    main()
