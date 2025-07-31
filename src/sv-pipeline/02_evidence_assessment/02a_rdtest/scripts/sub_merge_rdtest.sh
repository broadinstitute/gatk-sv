#!/bin/bash
#
# sub_rdtest_splits.sh
#

set -e

for batch in Phase1 Pilot; do
  for source in delly dragen lumpy manta wham depth; do
    for chrom in $(seq 1 22) X Y; do
      bsub -q normal -o merge_logs/${batch}.${source}.${chrom}.out -sla miket_sc -J "merge_${batch}_${source}_${chrom}" "
        ./merge.sh $batch $source $chrom" > /dev/null
    done
  done
done

