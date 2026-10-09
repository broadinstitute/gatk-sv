#!/bin/bash
#
# reformat_metrics.sh
#
# TODO: do this as part of metric aggregation
#

set -e

blacklisted=$1
metricfile=$2

  ##remove duplicate lines##	\
  ##create metric file remove variants where all samples are called to have a CNV## \
  ##add poor region coverage and size as metrics at the end## \
  ##add size and coverage NA for none CNV sv types## \
  ##remove chr X and Y \
  ##get rid of straggler header line## \

  # BAF p-value/log-scale and inf handling moved to aggregate/preprocess. Do not re-add a sign flip
  # for BAF_DEL_LOGLIK here: the producer (GATK BafHetRatioTester) already negates the log-likelihood.

cat ${@:3} \
  > $metricfile

