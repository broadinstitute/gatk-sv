#!/usr/bin/env python

"""
Collect ped file metrics. Writes stats to stdout.

Metrics:
  ped_file_count : Number of records
  ped_file_<family> : Number of families of each type

"""

import argparse
import sys
import pandas as pd
import svtest.utils.IOUtils as iou

KEY_PREFIX = "ped_file_"

SINGLETON_STR = "singletons"
DUO_STR = "duos"
TRIO_STR = "family_size_3"
QUAD_STR = "family_size_4"
QUINTET_PLUS_STR = "family_size_5_or_larger"

MALE = "1"
FEMALE = "2"
MALE_METRIC = "male"
FEMALE_METRIC = "female"
OTHER_METRIC = "other"


def main(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        prog='svtest ped-file',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('test_ped_file', type=str)
    parser.add_argument('--sample-list', type=str, default=None,
                        help='Sample ids not found in this list will cause an error')
    parser.add_argument('--prefix', type=str, default=None,
                        help='Prefix to add to metric names')

    # Print help if no arguments specified
    if len(argv) == 0:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args(argv)

    if args.sample_list is not None:
        samples = iou.read_samples_list(args.sample_list)
    else:
        samples = None

    # Get metrics
    df = pd.read_csv(args.test_ped_file, sep='\t', names=range(6))
    metrics = get_metrics(df, valid_samples=samples, metric_prefix=args.prefix)

    # Write metrics
    write_metrics(metrics)


def write_metrics(metrics):
    for key in metrics:
        sys.stdout.write("%s\t%s\n" % (key, str(metrics[key])))


def get_metrics(df, valid_samples=None, metric_prefix=None):
    check_samples(df, valid_samples)
    if metric_prefix is None:
        pfx = KEY_PREFIX
    else:
        pfx = metric_prefix + "_" + KEY_PREFIX
    metrics = {
        pfx + "count": df.shape[0]
    }
    metrics = add_family_count_metrics(metrics, df, pfx)
    metrics = add_sex_metrics(metrics, df, pfx)
    return metrics


def add_sex_metrics(metrics, df, prefix):
    num_male = count_sex(MALE, df)
    num_female = count_sex(FEMALE, df)
    num_other = df.shape[0] - num_male - num_female
    metrics[prefix + MALE_METRIC] = num_male
    metrics[prefix + FEMALE_METRIC] = num_female
    metrics[prefix + OTHER_METRIC] = num_other
    return metrics


def count_sex(type, df):
    sex_col = df[4].astype('str')
    return sex_col[sex_col == type].size


def add_family_count_metrics(metrics, df, prefix):
    counts_by_size = get_family_counts(df)
    for key in counts_by_size:
        metrics[prefix + key] = counts_by_size[key]
    return metrics


def get_family_counts(df):
    family_id_col = df[0]
    family_ids = list(set(family_id_col))
    counts_by_id = {}
    for id in family_ids:
        counts_by_id[id] = family_id_col[family_id_col == id].size
    counts_by_size = {SINGLETON_STR: 0, DUO_STR: 0,
                      TRIO_STR: 0, QUAD_STR: 0, QUINTET_PLUS_STR: 0}
    for id in counts_by_id:
        if counts_by_id[id] == 1:
            counts_by_size[SINGLETON_STR] += 1
        elif counts_by_id[id] == 2:
            counts_by_size[DUO_STR] += 1
        elif counts_by_id[id] == 3:
            counts_by_size[TRIO_STR] += 1
        elif counts_by_id[id] == 4:
            counts_by_size[QUAD_STR] += 1
        else:
            counts_by_size[QUINTET_PLUS_STR] += 1
    return counts_by_size


def check_samples(df, valid_samples):
    samples = set(df[1])
    if len(samples) < df.shape[0]:
        raise ValueError('There are duplicate sample ids in the ped file')
    if valid_samples is not None:
        unexpected_samples = samples - set(valid_samples)
        if len(unexpected_samples) > 0:
            raise ValueError('Unexpected samples: %s' % unexpected_samples)


if __name__ == '__main__':
    main()
