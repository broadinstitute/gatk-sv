#!/bin/bash
#
# label.sh
#

set -e

bed=$1
name=$2

cut -f -3 $bed \
  | sort -k1,1V -k2,2n \
  | bedtools merge -i stdin \
  | awk -v OFS="\t" -v name=$name '{print $0, name}' \
  > cleaned/${name}.bed
