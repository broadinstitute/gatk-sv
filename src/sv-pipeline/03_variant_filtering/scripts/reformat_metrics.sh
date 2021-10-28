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
  ####fix BAF p-value to be on log scale to match others### \
  ####flip sign for BAF DEL log likiehood so postive is more significant ### \
  ##replace inf from -log10 p-value with max -log10(P) of 300## \
  ##add poor region coverage and size as metrics at the end## \
  ##add size and coverage NA for none CNV sv types## \
  ##remove chr X and Y \
  ##get rid of straggler header line## \

  # moved to aggregate/preprocess
  # | awk '{if ($34!="NA" && $34!="BAF_KS_pval") {$34=-log($34)/log(10)} print}' \
  # | awk '{if ($31!="NA" && $31!="BAF_del_loglik") {$31=-$31} print}' \
  # | sed 's/inf/300/g' \

cat ${@:3} \
  > $metricfile

