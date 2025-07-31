#!/bin/bash
#
# sub_rdtest_splits.sh
#

set -e

coveragefile=/data/talkowski/Samples/common-mind/matrices/CMC.all.binCov.bed.gz
medianfile=/data/talkowski/Samples/common-mind/matrices/CMC.all.binCov.median
famfile=/data/talkowski/Samples/common-mind/ref/CMC.fam

for batch in CMC; do
  for source in delly dragen lumpy manta wham depth; do
    for chrom in 1; do
    # for chrom in $(seq 1 22) X Y; do
      for bed in split_beds/${batch}.${source}.${chrom}.*; do
        split=$(basename $bed | cut -d"." -f4)
        # echo $bed
        # echo $split
        bsub -q short -o split_logs/${batch}.${source}.${chrom}.out -sla miket_sc -J "rdtest_${bed}" "
          Rscript scripts/RdTest.R \
            -b $bed \
            -o split_rdtest \
            -n ${batch}.${source}.${chrom}.${split} \
            -c $coveragefile \
            -m $medianfile \
            -f $famfile " > /dev/null
      done
    done
  done
done

