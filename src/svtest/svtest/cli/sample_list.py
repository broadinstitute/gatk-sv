#!/usr/bin/env python

"""
Collect sample list metrics. Writes stats to stdout.

Metrics:
  sample_list_count   : Number of samples

"""

import argparse
import sys
import svtest.utils.IOUtils as iou

KEY_PREFIX = "sample_list_"


def main(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        prog='svtest sample-list',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('test_sample_list', type=str)
    parser.add_argument('--valid-sample-list', type=str, default=None,
                        help='Sample ids not found in this list will cause an error')
    parser.add_argument('--prefix', type=str, default=None,
                        help='Prefix to add to metric names')

    # Print help if no arguments specified
    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    samples = iou.read_samples_list(args.test_sample_list, fail_on_empty=False)
    if args.valid_sample_list is not None:
        valid_samples = iou.read_samples_list(args.valid_sample_list)
    else:
        valid_samples = None

    # Get metrics
    metrics = get_metrics(samples, valid_samples, args.prefix)

    # Write metrics
    write_metrics(metrics)


def write_metrics(metrics):
    for key in metrics:
        sys.stdout.write("%s\t%s\n" % (key, str(metrics[key])))


def get_metrics(samples, valid_samples, metric_prefix):
    if valid_samples is not None:
        unexpected_samples = set(samples) - set(valid_samples)
        if len(unexpected_samples) > 0:
            raise ValueError('Unexpected samples: %s' % unexpected_samples)

    if metric_prefix is None:
        pfx = KEY_PREFIX
    else:
        pfx = metric_prefix + "_" + KEY_PREFIX
    metrics = {
        pfx + "count": len(samples)
    }
    return metrics


if __name__ == '__main__':
    main()
