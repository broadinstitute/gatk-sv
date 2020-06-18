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


cat split_rdtest/${batch}.${source}.${chrom}.*.metrics \
  | sed -r -e '/^chr\s/d' \
  | sort -k1,1V -k2,2n \
  | cat <(head -n1 split_rdtest/Phase1.delly.1.00000.metrics) - \
  > rdtest/${batch}.${source}.${chrom}.metrics
