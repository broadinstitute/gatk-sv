#!/bin/bash

# Copyright (c) 2018 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# Benchmarks the sensitivity of one SV callset against another

set -e

###USAGE
usage(){
cat <<EOF

usage: compare_callsets.sh [options] SET1 SET2 CONTIGS

Helper tool to compare the sensitivity of one SV callset vs another

Positional arguments:
  SET1            SV callset to be benchmarked
  SET2            Ground truth SV callset to use for benchmarking
  CONTIGS         List of contigs to evaluate

Optional arguments:
  -h  HELP        Show this help message and exit
  -O  OUTFILE     Output file (default: stdout)
  -d  DISTANCE    Maximum DISTance between breakpoints during comparisons (default: 250bp)
  -p  PREFIX      Prefix for benchmarking variant IDs (default: Benchmarking_SV)
  -C  CARRIER     Report carrier frequencies (default: report allele frequencies)

Notes:
  1) SET1/SET2 are expected to be BED3+ formatted files
  2) SET1 must have at least the following five columns, in order:
     chr, start, end, SV ID, SV type
  3) The last column of SET1 must contain allele frequency
  4) SET2 must have exactly six columns as follows, in order:
     chr, start, end, SV type, SV size, allele frequency

EOF
}


###PARSE ARGS
OUTFILE=/dev/stdout
PREFIX="Benchmarking_SV"
DIST=250
CARRIER=0
while getopts ":O:d:p:Ch" opt; do
	case "$opt" in
		h)
			usage
			exit 0
			;;
    O)
      OUTFILE=${OPTARG}
      ;;
    d)
      DIST=${OPTARG}
      ;;
    p)
      PREFIX=${OPTARG}
      ;;
    C)
      CARRIER=1
      ;;
	esac
done
shift $(( ${OPTIND} - 1))
SET1=$1
SET2=$2
CONTIGS=$3


###SET BIN
BIN=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )


###PROCESS ARGS & OPTIONS
#Check for required input SET1
if [ -z ${SET1} ]; then
  echo -e "\nERROR: input SET1 not specified\n"
  usage
  exit 0
elif ! [ -s ${SET1} ]; then
  echo -e "\nERROR: input SET1 either empty or not found\n"
  usage
  exit 0
fi
#Check for required input SET2
if [ -z ${SET2} ]; then
  echo -e "\nERROR: input SET2 not specified\n"
  usage
  exit 0
elif ! [ -s ${SET2} ]; then
  echo -e "\nERROR: input SET2 either empty or not found\n"
  usage
  exit 0
fi
#Check for required list of contigs
if [ -z ${CONTIGS} ]; then
  echo -e "\ncompare_callsets_perSample.sh ERROR: input contig list not specified\n"
  usage
  exit 1
elif ! [ -s ${CONTIGS} ]; then
  echo -e "\ncompare_callsets_perSample.sh ERROR: input contig list either empty or not found\n"
  usage
  exit 1
fi
#Checks output file
if [ -z ${OUTFILE} ]; then
  echo -e "\nERROR: output file not specified\n"
  usage
  exit 0
fi
#Checks prefix
if [ -z ${PREFIX} ]; then
  PREFIX="Benchmarking_SV"
fi


###PREP INPUT FILES
OVRTMP=`mktemp -d`
#Unzip SET1, if gzipped, restrict to contigs in $CONTIGS, and automatically set SVs with size <1 to 1
if [ $( file ${SET1} | fgrep " gzip " | wc -l ) -gt 0 ] || \
   [ $( echo ${SET1} | awk -v FS="." '{ if ($NF ~ /gz|bgz/) print "TRUE" }' ) ]; then
  zcat ${SET1} | fgrep "#" > ${OVRTMP}/set1.bed
  zcat ${SET1} | fgrep -v "#" | awk -v OFS="\t" '{ if ($3<$2) $3=$2; print }' | \
  grep -f <( awk '{ print "^"$1"\t" }' ${CONTIGS} ) | \
  sort -Vk1,1 -k2,2n -k3,3n | uniq >> ${OVRTMP}/set1.bed
else
  cat ${SET1} | fgrep "#" > ${OVRTMP}/set1.bed
  cat ${SET1} | fgrep -v "#" | awk -v OFS="\t" '{ if ($3<$2) $3=$2; print }' | \
  grep -f <( awk '{ print "^"$1"\t" }' ${CONTIGS} ) | \
  sort -Vk1,1 -k2,2n -k3,3n | uniq >> ${OVRTMP}/set1.bed
fi
#Set carrierFrequency as final column, if optioned
if [ ${CARRIER} == 1 ]; then
  idx=$( head -n1 ${OVRTMP}/set1.bed | sed 's/\t/\n/g' | \
         awk '{ if ($1=="carrierFreq") print NR }' )
  awk -v FS="\t" -v OFS="\t" -v idx=${idx} \
  '{ print $0, $(idx) }' ${OVRTMP}/set1.bed > \
  ${OVRTMP}/set1.bed2
  mv ${OVRTMP}/set1.bed2 ${OVRTMP}/set1.bed
fi
#Unzip & format SET2, if gzipped
if [ $( file ${SET2} | fgrep " gzip " | wc -l ) -gt 0 ] || \
   [ $( echo ${SET2} | awk -v FS="." '{ if ($NF ~ /gz|bgz/) print "TRUE" }' ) ]; then
  zcat ${SET2} | fgrep "#" | awk -v OFS="\t" \
  '{ print $1, $2, $3, "VID", $4, $5, $6 }' > ${OVRTMP}/set2.bed
  zcat ${SET2} | fgrep -v "#" | \
  awk -v OFS="\t" -v PREFIX=${PREFIX} \
  '{ if ($3<$2) $3=$2; print $1, $2, $3, PREFIX"_"NR, $4, $5, $6 }' | \
  grep -f <( awk '{ print "^"$1"\t" }' ${CONTIGS} ) | \
  sort -Vk1,1 -k2,2n -k3,3n | uniq >> ${OVRTMP}/set2.bed
else
  cat ${SET2} | fgrep "#" | awk -v OFS="\t" \
  '{ print $1, $2, $3, "VID", $4, $5, $6 }' > ${OVRTMP}/set2.bed
  cat ${SET2} | fgrep -v "#" | \
  awk -v OFS="\t" -v PREFIX=${PREFIX} \
  '{ if ($3<$2) $3=$2; print $1, $2, $3, PREFIX"_"NR, $4, $5, $6 }' | \
  grep -f <( awk '{ print "^"$1"\t" }' ${CONTIGS} ) | \
  sort -Vk1,1 -k2,2n -k3,3n | uniq >> ${OVRTMP}/set2.bed
fi


###RUN INTERSECTIONS
# Check if any variants are left in set 2 (benchmarking set) after subsetting to contigs of interest
if [ $( cat ${OVRTMP}/set2.bed | fgrep -v "#" | wc -l ) -gt 0 ]; then

  #Intersect method 1 data
  bedtools intersect -loj -r -f 0.5 \
    -a <( awk -v small_cutoff=5000 -v OFS="\t" '{ if ($6>=small_cutoff) print $0 }' ${OVRTMP}/set2.bed ) \
    -b ${OVRTMP}/set1.bed > \
    ${OVRTMP}/OVR1.raw.bed
  bedtools intersect -loj -r -f 0.1 \
    -a <( awk -v small_cutoff=5000 -v OFS="\t" '{ if ($6<small_cutoff) print $0 }' ${OVRTMP}/set2.bed ) \
    -b ${OVRTMP}/set1.bed >> \
    ${OVRTMP}/OVR1.raw.bed
  #Intersect method 1a: 50% reciprocal overlap (10% for small SV), matching SV types
  awk -v FS="\t" -v OFS="\t" \
  '{ if ($5==$12 || $5=="DUP" && $12 ~ /"CNV"/ || $12=="DUP" && $5=="MCNV" || $5=="DEL" && $12=="MCNV" || $5=="MCNV" && $12=="DEL") print $4, $NF; else if ($12==".") print $4, "NO_OVR" }' \
  ${OVRTMP}/OVR1.raw.bed | sort -Vk1,1 -k2,2n | uniq | \
  awk -v OFS="\t" '{ if ($2=="NA") $2="1"; print $1, $2 }' > \
  ${OVRTMP}/OVR1a.raw.txt
  cut -f1 ${OVRTMP}/OVR1a.raw.txt | fgrep -wvf - ${OVRTMP}/OVR1.raw.bed | \
  awk -v OFS="\t" '{ print $4, "NO_OVR" }' | sort -Vk1,1 -k2,2n | uniq >> \
  ${OVRTMP}/OVR1a.raw.txt
  sort -Vk1,1 -k2,2n ${OVRTMP}/OVR1a.raw.txt | uniq > ${OVRTMP}/OVR1a.raw.txt2
  mv ${OVRTMP}/OVR1a.raw.txt2 ${OVRTMP}/OVR1a.raw.txt
  #Intersect method 1b: 50% reciprocal overlap (10% for small SV), any SV types
  awk -v FS="\t" -v OFS="\t" \
  '{ if ($12==".") print $4, "NO_OVR"; else print $4, $NF }' \
  ${OVRTMP}/OVR1.raw.bed | sort -Vk1,1 -k2,2n | uniq | \
  awk -v OFS="\t" '{ if ($2=="NA") $2="1"; print $1, $2 }' > \
  ${OVRTMP}/OVR1b.raw.txt
  cut -f1 ${OVRTMP}/OVR1b.raw.txt | fgrep -wvf - ${OVRTMP}/OVR1.raw.bed | \
  awk -v OFS="\t" '{ print $4, "NO_OVR" }' | sort -Vk1,1 -k2,2n | uniq >> \
  ${OVRTMP}/OVR1b.raw.txt
  sort -Vk1,1 -k2,2n ${OVRTMP}/OVR1b.raw.txt | uniq > ${OVRTMP}/OVR1b.raw.txt2
  mv ${OVRTMP}/OVR1b.raw.txt2 ${OVRTMP}/OVR1b.raw.txt
  #Intersect method 2 data
  bedtools intersect -loj -a ${OVRTMP}/set2.bed -b ${OVRTMP}/set1.bed > \
  ${OVRTMP}/OVR2.raw.bed
  #Intersect method 2a: any overlap, breakpoints within $DIST, matching SV types
  awk -v FS="\t" -v OFS="\t" -v DIST=${DIST} \
  '{ if ($12!="." && ($2-$9<=DIST && $2-$9>=-DIST) && ($3-$10<=DIST && $3-$10>=-DIST) && ($5==$12 || $5=="DUP" && $12=="MCNV" || $12=="DUP" && $5=="MCNV" || $5=="DEL" && $12=="MCNV" || $5=="MCNV" && $12=="DEL")) print $4, $NF }' \
  ${OVRTMP}/OVR2.raw.bed | sort -Vk1,1 -k2,2n | uniq | \
  awk -v OFS="\t" '{ if ($2=="NA") $2="1"; print $1, $2 }' > \
  ${OVRTMP}/OVR2a.raw.txt
  cut -f1 ${OVRTMP}/OVR2a.raw.txt | sort | uniq | fgrep -wvf - ${OVRTMP}/OVR2.raw.bed | \
  awk -v OFS="\t" '{ print $4, "NO_OVR" }' | sort -Vk1,1 -k2,2n | uniq >> \
  ${OVRTMP}/OVR2a.raw.txt
  sort -Vk1,1 -k2,2n ${OVRTMP}/OVR2a.raw.txt | uniq > ${OVRTMP}/OVR2a.raw.txt2
  mv ${OVRTMP}/OVR2a.raw.txt2 ${OVRTMP}/OVR2a.raw.txt
  #Intersect method 2b: any overlap, breakpoints within $DIST, any SV types
  awk -v FS="\t" -v OFS="\t" -v DIST=${DIST} \
  '{ if ($12!="." && ($2-$9<=DIST && $2-$9>=-DIST) && ($3-$10<=DIST && $3-$10>=-DIST)) print $4, $NF }' \
  ${OVRTMP}/OVR2.raw.bed | sort -Vk1,1 -k2,2n | uniq | \
  awk -v OFS="\t" '{ if ($2=="NA") $2="1"; print $1, $2 }' > \
  ${OVRTMP}/OVR2b.raw.txt
  cut -f1 ${OVRTMP}/OVR2b.raw.txt | sort | uniq | fgrep -wvf - ${OVRTMP}/OVR2.raw.bed | \
  awk -v OFS="\t" '{ print $4, "NO_OVR" }' | sort -Vk1,1 -k2,2n | uniq >> \
  ${OVRTMP}/OVR2b.raw.txt
  sort -Vk1,1 -k2,2n ${OVRTMP}/OVR2b.raw.txt | uniq > ${OVRTMP}/OVR2b.raw.txt2
  mv ${OVRTMP}/OVR2b.raw.txt2 ${OVRTMP}/OVR2b.raw.txt
  #Intersect method 3: any overlap, buffer ± $DIST, any svtype
  bedtools intersect -loj -a ${OVRTMP}/set2.bed \
  -b <( awk -v OFS="\t" -v DIST=${DIST} '{ $2=$2-DIST; $3=$3+DIST; print }' \
        ${OVRTMP}/set1.bed | awk -v OFS="\t" '{ if ($2<0) $2=0; print }' ) | \
  sort -Vk1,1 -k2,2n | uniq | awk -v FS="\t" -v OFS="\t" \
  '{ if ($12==".") print $4, "NO_OVR"; else print $4, $NF }' \
  | sort -Vk1,1 -k2,2n | uniq | \
  awk -v OFS="\t" '{ if ($2=="NA") $2="1"; print $1, $2 }' > \
  ${OVRTMP}/OVR3.raw.txt

  ###CONVERT INTERSECTIONS TO FINAL TABLE
  ${BIN}/compare_callsets_helper.R \
    ${OVRTMP}/set2.bed \
    ${OVRTMP}/OVR1a.raw.txt \
    ${OVRTMP}/OVR1b.raw.txt \
    ${OVRTMP}/OVR2a.raw.txt \
    ${OVRTMP}/OVR2b.raw.txt \
    ${OVRTMP}/OVR3.raw.txt \
    ${OUTFILE}

# If no variants remain in set 2 after subsetting to contigs of interest, output empty OUTFILE
else
  echo -e "#chr\tstart\tend\tVID\tsvtype\tlength\tAF\tovr1a\tovr1b\tovr2a\tovr2b\tovr3" > ${OUTFILE}
fi

###CLEAN UP
rm -rf ${QCTMP}

