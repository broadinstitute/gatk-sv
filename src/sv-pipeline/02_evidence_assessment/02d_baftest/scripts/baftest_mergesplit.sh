#!bin/bash
# This script will merge split rdtest statistics by split number
# The output will be there be 1 file for each chrom+batch+source for allosomes
# For sex chromosome the result will be one file per sex chromosome for each sex
# The usage format is ```bash rdtest_mergesplit.sh {batch} {source} {chrom}

set -e

batch=$1
source=$2
chrom=$3
inputfolder="baftest_split"
outputfolder="baftest"
input=()
for item in  $inputfolder/$batch.$source.$chrom.*.metrics; do
        input+=($item)
done

    cat $inputfolder/$batch.$source.$chrom.*.metrics \
          | sed -r -e '/^chr/d' \
          | sort -k1,1V -k2,2n \
          | cat <(head -n1 $input) - \
          > $outputfolder/$batch.$source.$chrom.metrics
fi

