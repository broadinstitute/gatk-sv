#!/bin/bash
#
# clean_VCF_part3.sh
#
# 
#
# Copyright (C) 2018 Harrison Brand<hbrand1@mgh.harvard.edu>
# Distributed under terms of the MIT license.

set -euo pipefail

combined_file=$1

awk '{print $1}' $combined_file \
  | sort \
  | uniq -c \
  | sort -nrk1,1  \
  > variant.count.txt

final=0
prev=0
var=0

while read line
do
  i=$(echo $line|awk -v prev=$prev '{print $1+prev}' )
  let "var=$var+1"
  if [ $i -gt 5000  ] || [ $var -gt 100  ]
  then
    final=$(echo $final|awk '{print $1+1}')
    prev=0
    var=0
  else
    prev=$i
  fi
done < variant.count.txt

j=0
prev=0
mkdir shards

while read line
do
  i=$(echo $line|awk -v prev=$prev '{print $1+prev}' )
  let "var=$var+1"
  if [ $i -gt 5000  ] || [ $var -gt 100  ]
  then
    j=$(echo $j|awk '{print $1+1}')
    prev=0
    var=0
  else
    prev=$i
  fi
  out=$(echo $j"_"$final)
  echo $line|awk '{print $2}'|fgrep -wf - $combined_file >>shards/out.$out.txt || true
done < variant.count.txt
