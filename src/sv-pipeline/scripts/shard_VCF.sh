#!/bin/bash

# Shard VCF into sequential, evenly-sized splits


###USAGE
usage(){
cat <<EOF

usage: shard_VCF.sh [-h] VCF LINES PREFIX

Shard a VCF into sequential, evenly-sized splits

Positional arguments:
  VCF             VCF from sv-pipeline
  LINES           Number of records per shard
  PREFIX          Output prefix for all shards

Optional arguments:
  -h  HELP        Show this help message and exit

EOF
}


###PARSE ARGS
while getopts "h" opt; do
	case "$opt" in
		h)
			usage
			exit 0
			;;
	esac
done
shift $(( ${OPTIND} - 1))
VCF=$1
LINES=$2
PREFIX=$3


###PROCESS ARGS & OPTIONS
#Check for required input
if [ -z ${VCF} ]; then
  echo -e "\nERROR: input VCF not specified\n"
  usage
  exit 0
fi
if ! [ -s ${VCF} ]; then
  echo -e "\nERROR: input VCF either empty or not found\n"
  usage
  exit 0
fi


###SHARD VCF
zcat ${VCF} | sed -n '1,2000p' | fgrep "#" > header.tmp
zcat ${VCF} | fgrep -v "#" | \
split -a 6 -d -l ${LINES} - ${PREFIX}
nsplits=$( find `pwd` -name "$( echo ${PREFIX} | sed 's/\//\t/g' | \
           awk -v FS="\t" '{ print $NF }' )*" | wc -l )
for i in $( seq -w 000000 $(( ${nsplits} - 1 )) ); do
  cat header.tmp ${PREFIX}${i} | vcf-sort > ${PREFIX}${i}.vcf
  bgzip -f ${PREFIX}${i}.vcf
  tabix -f ${PREFIX}${i}.vcf.gz
done
