#!/bin/bash
#
# make_bkpt_rmsk.sh
#
# 
#
# Copyright (C) 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

set -e

cat \
    /data/talkowski/rlc47/src/GRCh37.*.RMSK.merged.bed \
    /data/talkowski/rlc47/src/GRch37.segdups_gaps_abParts_heterochrom.lumpy.exclude.bed \
  | sort -k1,1V -k2,2n \
  > breakpoint.rmsk.bed
