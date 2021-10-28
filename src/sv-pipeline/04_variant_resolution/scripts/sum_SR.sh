#!/bin/bash
#
# sum_SR.sh
#

set -e

cat $1 \
  | fgrep -v -e "name" \
  | awk '{print $1"@"$3 "\t" $2 "\t" $4}' \
  | awk '{a[$1]+=$3;}END{for(i in a)print i"\t"a[i];}' \
  | tr '@' '\t' \
  | sort -k1,1 \
  | gzip \
  > $2
