#!/bin/bash
#
# split.sh
#

set -e

batch=$1
source=$2
chrom=$3

bed="../algorithm_integration/rdtest_beds/${batch}.${source}.${chrom}.bed"

sed '/^#/d' $bed \
  | split -a 5 -d -l 100 - "split_beds/${batch}.${source}.${chrom}."
