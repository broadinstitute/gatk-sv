#!/bin/bash
#
# split.sh
#
# 
#
# Copyright (C) 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

set -e

batch=$1
source=$2
chrom=$3

bed="../algorithm_integration/rdtest_beds/${batch}.${source}.${chrom}.bed"

sed '/^#/d' $bed \
  | split -a 5 -d -l 100 - "split_beds/${batch}.${source}.${chrom}."
