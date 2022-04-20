#!/bin/python

# Synopsis:
#  Creates a new QC definitions file (valid ranges) from a given Batch run's QC metrics output
#

import argparse
import math


def get_number(val):
    try:
        return int(val)
    except:
        try:
            return float(val)
        except:
            return val


def get_metric_range(val, args):
    val = get_number(val)
    if isinstance(val, float):
        if math.isnan(val):
            return None
        else:
            return val * (1.0 - args.range), val * (1.0 + args.range)
    elif isinstance(val, int):
        if val < args.integer_zero_lower:
            return 0, args.integer_zero_upper
        else:
            return round(val * (1.0 - args.range)), round(val * (1.0 + args.range))
    else:
        return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("metrics_file_batch", type=str,
                        help="GATKSVPipelineBatch metrics output file")
    parser.add_argument("--integer-zero-lower", type=int, default=10,
                        help="Integer metrics below this value will have their minimum "
                             "value set to 0 and maximum value set to --integer-zero-upper")
    parser.add_argument("--integer-zero-upper", type=int, default=11,
                        help="See --integer-zero-lower")
    parser.add_argument("--range", type=float, default=0.1,
                        help="Fraction of metric value to define valid range (+/-)")
    args = parser.parse_args()

    with open(args.metrics_file_batch, 'r') as f:
        for line in f:
            key, val = line.rstrip().split('\t')
            range = get_metric_range(val, args)
            if range is not None:
                print(f"{key}\t{range[0]}\t{range[1]}")


if __name__ == "__main__":
    main()
