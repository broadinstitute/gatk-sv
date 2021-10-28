#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""

"""

import argparse
import numpy as np
import scipy.stats as ss


def convert_poisson_p(log_pval):
    """
    Convert poisson p-value to count, assuming 0 background

    Can't use ss.poisson.ppf because k/mu are flipped in pe/sr test calculation
    """
    count = 0
    while True:
        pval_cmp = -np.log10(ss.poisson.cdf(0, count))

        if pval_cmp > log_pval:
            return count - 1

        count += 1


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('pval', type=float)
    args = parser.parse_args()

    print(convert_poisson_p(args.pval))


if __name__ == '__main__':
    main()
