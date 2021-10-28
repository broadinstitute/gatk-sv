#!/bin/bash

# Contact: Ryan Collins <rlcollins@g.harvard.edu>

# Intelligently shards a VCF prior to complex resolution (for parallelization)

set -euo pipefail

###USAGE
usage(){
cat <<EOF
usage: shardVCF_preResolveCPX.sh [-h] [-L MIN_LINES_PER_SHARD] [-S MAX_SHARDS]
                                 [-N NONCLUSTER_SHARDS] [-P PREFIX] [-T VCF_INDEX] VCF
Shard VCF prior to complex resolution (for parallelization)
Positional arguments:
  VCF                      VCF from sv-pipeline prior to svtk resolve
  OUTDIR                   Output directory for all shards
Optional arguments:
  -h  HELP                 Show this help message and exit
  -L  MIN_LINES_PER_SHARD  Minimum variants per shard [default: 10]
  -S  MAX_SHARDS           Maximum number of shards [default: 100]
  -N  NONCLUSTER_SHARDS    Add an additional N shards to evenly split the 
                           remaining records in the VCF that do not cluster with
                           any other variants [default: 30]
  -P  PREFIX               String appended to the filename of each output VCF 
                           shard [default: "vcf_shard"]
  -O  OUTDIR               Output directory [default: pwd]
  -T  VCF_index            Tabix index corresponding to input VCF 
                           [default: VCF path + .tbi]
Notes:
  1. The total number of shards created will be no more 
     than MAX_SHARDS + NONCLUSTER_SHARDS.
EOF
}


###PARSE ARGS
MIN_LINES_PER_SHARD=10
MAX_SHARDS=100
NONCLUSTER_SHARDS=30
PREFIX="vcf_shard"
OUTDIR=`pwd`
while getopts ":L:S:N:P:O:T:h" opt; do
  case "$opt" in
    h)
      usage
      exit 0
      ;;
    L)
      MIN_LINES_PER_SHARD=${OPTARG}
      ;;
    S)
      MAX_SHARDS=${OPTARG}
      ;;
    N)
      NONCLUSTER_SHARDS=${OPTARG}
      ;;
    P)
      PREFIX=${OPTARG}
      ;;
    O)
      OUTDIR=${OPTARG}
      ;;
    T)
      VCF_INDEX=${OPTARG}
      ;;
  esac
done
shift $(( ${OPTIND} - 1))
VCF=$1
if [ -z ${VCF_INDEX} ]; then
  VCF_INDEX=${VCF}.tbi
fi


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
if [ $( file ${VCF} | fgrep "gzip" | wc -l ) -lt 1 ]; then
  echo -e "\nERROR: input VCF must be bgzipped\n"
  usage
  exit 0
fi
if [ -z ${VCF_INDEX} ]; then
  echo -e "\nERROR: input VCF index not specified\n"
  usage
  exit 0
fi
if ! [ -s ${VCF_INDEX} ]; then
  echo -e "\nERROR: input VCF index either empty or not found\n"
  usage
  exit 0
fi
if [ -z ${OUTDIR} ]; then
  echo -e "\nERROR: output directory not specified\n"
  usage
  exit 0
fi
#Creates $OUTDIR if necessary
if ! [ -e ${OUTDIR} ]; then
  mkdir ${OUTDIR}
fi
#Checks for numeric L and S options
if [ ${MIN_LINES_PER_SHARD} -lt 1 ] || [ ${MIN_LINES_PER_SHARD} -gt 1000000 ]; then
  echo -e "\nERROR: MIN_LINES_PER_SHARD must be an integer between 1 and 1000000\n"
  usage
  exit 0
fi
if [ ${MAX_SHARDS} -lt 1 ] || [ ${MAX_SHARDS} -gt 10000 ]; then
  echo -e "\nERROR: MAX_SHARDS must be an integer between 1 and 10000\n"
  usage
  exit 0
fi
#Set path to execution directory
BIN=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )


###IDENTIFY CANDIDATE COMPLEX CLUSTERS
#Make temporary directory
SHARD_VCF_TMP=`mktemp -d`
#Cut vcf to single-sample for improved clustering speed
zcat ${VCF} | cut -f1-10 | bgzip -c > ${SHARD_VCF_TMP}/single_sample_input.vcf.gz
#Identify all candidate complex variant clusters (generous 1kb clustering)
svtk vcfcluster \
  -d 1000 \
  -f 0 \
  -p candidate_complex_clusters \
  --ignore-svtypes \
  -o 0 \
  --preserve-ids \
  <( echo "${SHARD_VCF_TMP}/single_sample_input.vcf.gz" ) \
  ${SHARD_VCF_TMP}/input_vcf.clustered.vcf
#Convert clustered variants to bed
svtk vcf2bed \
  --no-samples \
  --info ALL \
  ${SHARD_VCF_TMP}/input_vcf.clustered.vcf \
  ${SHARD_VCF_TMP}/input_vcf.clustered.bed
#Write list of clusters with >1 constituent variant
mem_idx=$( head -n1 ${SHARD_VCF_TMP}/input_vcf.clustered.bed | \
           sed 's/\t/\n/g' | awk '{ if ($1=="MEMBERS") print NR }' )
awk -v idx=${mem_idx} -v OFS="\t" \
'$idx ~ /,/ { print $1, $2, $3, $idx }' \
${SHARD_VCF_TMP}/input_vcf.clustered.bed | fgrep -v "#" > \
${SHARD_VCF_TMP}/candidate_complex_clusters.bed
#Add all non-CNV single-record variants
class_idx=$( head -n1 ${SHARD_VCF_TMP}/input_vcf.clustered.bed | \
           sed 's/\t/\n/g' | awk '{ if ($1=="SVTYPE") print NR }' )
awk -v idx=${class_idx} -v mem_idx=${mem_idx} -v OFS="\t" \
'$mem_idx !~ /,/ { print $1, $2, $3, $idx, $mem_idx }' \
${SHARD_VCF_TMP}/input_vcf.clustered.bed | \
awk -v OFS="\t" '$4 !~ /DEL|DUP|CNV|MCNV|mCNV/ { print $1, $2, $3, $5 }' | fgrep -v "#" >> \
${SHARD_VCF_TMP}/candidate_complex_clusters.bed
#Get min/max coordinates of all variants in list of VIDs
cat <( zcat ${VCF} | fgrep "#" | cut -f1-10 ) \
<( cut -f4 ${SHARD_VCF_TMP}/candidate_complex_clusters.bed | \
   sed 's/\,/\n/g' | sort -Vk1,1 | uniq | { fgrep -wf - <( zcat ${VCF} ) || true; } \
   | cut -f1-10 ) | \
svtk vcf2bed --no-samples /dev/stdin /dev/stdout > \
${SHARD_VCF_TMP}/candidate_complex_clusters.variant_coordinates.bed


#Split into breakpoints and pad all breakpoints by ±0kb 
#DEV NOTE: padding breakpoints for large chromosomes & many samples was causing
# issues where tens of thousands of breakpoints would end up in the same shard
# and take >36h to resolve, defeating the purpose of sharding
fgrep -v "#" ${SHARD_VCF_TMP}/candidate_complex_clusters.variant_coordinates.bed | \
awk -v OFS="\t" -v buffer=0 \
'{ print $1, $2-buffer, $2+buffer, $4"\n"$1, $3-buffer, $3+buffer, $4 }' | \
awk -v OFS="\t" '{ if ($2<1) $2=1; print $1, $2, $3, $4 }' | \
sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - -c 4 -o distinct > \
${SHARD_VCF_TMP}/breakpoint_intervals.bed
#Iterate over breakpoint intervals and write list of maximum nonredundant intervals
in_cluster=`mktemp`
remaining=`mktemp`
cp ${SHARD_VCF_TMP}/breakpoint_intervals.bed ${remaining}
while read chr start end VIDs; do
  #Get all lines associated with current VIDs
  echo -e "${VIDs}" | sed 's/,/\n/g' | { fgrep -wf - ${remaining} || true; } > ${in_cluster}
  #Only run if at least one line added to ${in_cluster}
  if [ $( cat ${in_cluster} | wc -l ) -gt 0 ]; then
    #Exclude all lines in ${in_cluster} from ${remaining}
    bedtools intersect -v -a ${remaining} -b ${in_cluster} > ${remaining}2
    mv ${remaining}2 ${remaining}
    #Iterate until no more related VIDs are present in ${remaining}
    until [ $( cut -f4 ${in_cluster} | sed 's/\,/\n/g' | { fgrep -wf - ${remaining} || true; } | wc -l ) -eq 0 ]; do
      #Add new lines to ${in_cluster}
      cut -f4 ${in_cluster} | sed 's/\,/\n/g' | { fgrep -wf - ${remaining} || true; } >> ${in_cluster}
      #Exclude all lines in ${in_cluster} from ${remaining}
      bedtools intersect -v -a ${remaining} -b ${in_cluster} > ${remaining}2
      mv ${remaining}2 ${remaining}
    done
    #Write out final interval
    for wrapper in 1; do
      #Print list of coordinates
      cut -f1-3 ${in_cluster} | sort -Vk1,1 -k2,2 -k3,3 | bedtools merge -i - | \
      awk '{ print $1":"$2"-"$3 }' | paste -s -d\;
      #Print list of involved VIDs
      cut -f4 ${in_cluster} | sed 's/,/\n/g' | sort | uniq | paste -s -d,
    done | paste -s
  fi
done < ${SHARD_VCF_TMP}/breakpoint_intervals.bed > \
${SHARD_VCF_TMP}/complex_intervals_to_test.final.txt
#Pull out exceptionally large clusters to the side to be placed in their own shards
while read ints VIDs; do
  if [ $( echo ${VIDs} | sed 's/,/\n/g' | wc -l ) -ge ${MIN_LINES_PER_SHARD} ]; then
    echo -e "${ints}\t${VIDs}"
  fi
done < ${SHARD_VCF_TMP}/complex_intervals_to_test.final.txt \
  > ${SHARD_VCF_TMP}/large_intervals_to_test.final.txt
if [ $( cat ${SHARD_VCF_TMP}/large_intervals_to_test.final.txt | wc -l ) -gt 0 ]; then
  cut -f2 ${SHARD_VCF_TMP}/large_intervals_to_test.final.txt \
  | sed 's/,/\n/g' \
  | { fgrep -wvf - ${SHARD_VCF_TMP}/complex_intervals_to_test.final.txt || true; } \
  > ${SHARD_VCF_TMP}/complex_intervals_to_test.final.txt2
  mv ${SHARD_VCF_TMP}/complex_intervals_to_test.final.txt2 \
  ${SHARD_VCF_TMP}/complex_intervals_to_test.final.txt
fi


###DETERMINE COORDINATES FOR EACH SHARD
#Split variants into shards based on number of variants
#If total number of intervals/MAX_SHARDS < MIN_LINES_PER_SHARD, evenly split into MIN_LINES_PER_SHARD sites per shard
if [ $(( $( cat ${SHARD_VCF_TMP}/complex_intervals_to_test.final.txt | wc -l ) / ${MAX_SHARDS} )) -lt ${MIN_LINES_PER_SHARD} ]; then
  ${BIN}/evenSplitter.R \
  -L ${MIN_LINES_PER_SHARD} \
  ${SHARD_VCF_TMP}/complex_intervals_to_test.final.txt \
  ${SHARD_VCF_TMP}/${PREFIX}.shard_intervals_
#Otherwise, split into MAX_SHARDS evenly-sized shards
else
  ${BIN}/evenSplitter.R \
  -S ${MAX_SHARDS} \
  ${SHARD_VCF_TMP}/complex_intervals_to_test.final.txt \
  ${SHARD_VCF_TMP}/${PREFIX}.shard_intervals_
fi
#Writes exceptionally large clusters to their own shards
n_shards=$( find ${SHARD_VCF_TMP} -name "${PREFIX}*" | wc -l )
if [ $( cat ${SHARD_VCF_TMP}/large_intervals_to_test.final.txt | wc -l ) -gt 0 ]; then
  while read ints VIDs; do
    n_shards=$(( ${n_shards} + 1 ))
    echo -e "${ints}\t${VIDs}" > ${SHARD_VCF_TMP}/${PREFIX}.shard_intervals_${n_shards}
  done < ${SHARD_VCF_TMP}/large_intervals_to_test.final.txt
fi
#Reformat interval shards
for i in $( seq 1 ${n_shards} ); do
  cut -f1 ${SHARD_VCF_TMP}/${PREFIX}.shard_intervals_${i} | \
  sed -e 's/\;/\n/g' -e 's/\:/\t/g' -e 's/\-/\t/g' | \
  sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - > \
  ${SHARD_VCF_TMP}/${PREFIX}.shard_intervals_${i}.bed
  rm ${SHARD_VCF_TMP}/${PREFIX}.shard_intervals_${i}
done


###SHARD VCF
#Convert full, original VCF to BED
svtk vcf2bed --no-samples \
  ${VCF} int.bed
#Harrison's patch for sharding
awk '{if ($1!~"#") print $1,$2,$2+1,$4,$5 \
  "\n" $1,$3-1,$3,$4,$5;else print}' OFS='\t' int.bed \
  > ${SHARD_VCF_TMP}/input_vcf.vcf2bed.bed
rm int.bed
#Iterate over all sharded intervals
for i in $( seq 1 $(( ${n_shards} )) ); do
  #Write exclusion list of VIDs already used in earlier shards
  touch ${SHARD_VCF_TMP}/used_VIDs.tmp
  if [ ${i} -gt 1 ]; then
    for j in $( seq 1 $(( ${i} - 1 )) ); do
      cat ${SHARD_VCF_TMP}/${PREFIX}.shard_${j}.VIDs.list
    done | sort | uniq > ${SHARD_VCF_TMP}/used_VIDs.tmp
  fi
  #Get list of IDs to be used in shard
  bedtools intersect \
  -a ${SHARD_VCF_TMP}/input_vcf.vcf2bed.bed \
  -b ${SHARD_VCF_TMP}/${PREFIX}.shard_intervals_${i}.bed | \
  cut -f4 | { fgrep -wvf ${SHARD_VCF_TMP}/used_VIDs.tmp || true; } > \
  ${SHARD_VCF_TMP}/${PREFIX}.shard_${i}.VIDs.list
  #Print header
  zcat ${VCF} | head -n1000 | fgrep "#" > \
  ${OUTDIR}/${PREFIX}.shard_${i}.vcf
  #Shard based on VIDs (slower than tabix, but avoids omitting variants)
  zcat ${VCF} | { fgrep -wf ${SHARD_VCF_TMP}/${PREFIX}.shard_${i}.VIDs.list || true; } >> \
  ${OUTDIR}/${PREFIX}.shard_${i}.vcf
  #Sanity check shard
  if [ $( { fgrep -v "#" ${OUTDIR}/${PREFIX}.shard_${i}.vcf || true; } | wc -l ) -gt 0 ]; then
    #Bgzip & tabix shard
    bgzip -f ${OUTDIR}/${PREFIX}.shard_${i}.vcf
    tabix -f ${OUTDIR}/${PREFIX}.shard_${i}.vcf.gz
  else
    rm ${OUTDIR}/${PREFIX}.shard_${i}.vcf
  fi
  #Clean up used VID list
  rm ${SHARD_VCF_TMP}/used_VIDs.tmp
done
#Write list of all VIDs used in cluster shards
zcat ${OUTDIR}/${PREFIX}.shard_*.vcf.gz \
  | cut -f1-3 | fgrep -v "#" | cut -f3 \
  > ${SHARD_VCF_TMP}/used_VIDs.cluster_shards.list
#Write list of eligible VIDs
fgrep -v "#" ${SHARD_VCF_TMP}/input_vcf.vcf2bed.bed \
  | cut -f4 \
  | { fgrep -wvf ${SHARD_VCF_TMP}/used_VIDs.cluster_shards.list || true; } \
  > ${SHARD_VCF_TMP}/remaining_VIDs.list
#Shard remaining records into no more than $NONCLUSTER_SHARDS shards
#If total number of records/NONCLUSTER_SHARDS < MIN_LINES_PER_SHARD, evenly split into MIN_LINES_PER_SHARD sites per shard
if [ $(( $( cat ${SHARD_VCF_TMP}/remaining_VIDs.list | wc -l ) / ${NONCLUSTER_SHARDS} )) -lt ${MIN_LINES_PER_SHARD} ]; then
  ${BIN}/evenSplitter.R \
  -L ${MIN_LINES_PER_SHARD} \
  ${SHARD_VCF_TMP}/remaining_VIDs.list \
  ${SHARD_VCF_TMP}/${PREFIX}.remainder_VIDs_
#Otherwise, split into MAX_SHARDS evenly-sized shards
else
  ${BIN}/evenSplitter.R \
  -S ${NONCLUSTER_SHARDS} \
  ${SHARD_VCF_TMP}/remaining_VIDs.list \
  ${SHARD_VCF_TMP}/${PREFIX}.remainder_VIDs_
fi
#Iterate over all non-cluster shards and generate VCF shards
n_noncluster_shards=$( find ${SHARD_VCF_TMP} -name "${PREFIX}.remainder_VIDs_*" | wc -l )
for i in $( seq 1 ${n_noncluster_shards} ); do
  idx=$(( ${n_shards} + ${i} ))
  #Print header
  zcat ${VCF} | head -n1000 | fgrep "#" > \
  ${OUTDIR}/${PREFIX}.shard_${idx}.vcf
  #Shard based on VIDs (slower than tabix, but avoids omitting variants)
  zcat ${VCF} | { fgrep -wf ${SHARD_VCF_TMP}/${PREFIX}.remainder_VIDs_${i} || true; } >> \
  ${OUTDIR}/${PREFIX}.shard_${idx}.vcf
  #Sanity check shard
  if [ $( { fgrep -v "#" ${OUTDIR}/${PREFIX}.shard_${idx}.vcf || true; } | wc -l ) -gt 0 ]; then
    #Bgzip & tabix shard
    bgzip -f ${OUTDIR}/${PREFIX}.shard_${idx}.vcf
    tabix -f ${OUTDIR}/${PREFIX}.shard_${idx}.vcf.gz
  else
    rm ${OUTDIR}/${PREFIX}.shard_${idx}.vcf
  fi
done


###SANITY CHECK SHARDS
while read shard; do 
  zcat ${shard} | cut -f1 | { fgrep -v "#" || true; } | wc -l
done < <( find ${OUTDIR} -name "${PREFIX}.shard_*.vcf.gz" ) \
  | sort -nrk1,1 \
  > ${SHARD_VCF_TMP}/vars_per_shard.txt
echo -e "FINISHED SHARDING VCF. RESULTING RECORDS PER SHARD FOR LARGEST 100 SHARDS:"
head -n100 ${SHARD_VCF_TMP}/vars_per_shard.txt | paste -s -d','
#If shard with most variants is >10-fold more than next-largest shard, exit with code 1
if [ $( find ${OUTDIR} -name "${PREFIX}.shard_*.vcf.gz" | wc -l ) -gt 1 ]; then
  first=$( sed -n '1p' ${SHARD_VCF_TMP}/vars_per_shard.txt )
  second=$( sed -n '2p' ${SHARD_VCF_TMP}/vars_per_shard.txt )
  if [ ! -z ${second} ] && [ ${second} -gt 0 ]; then
    if [ $(( ${first} / ${second} )) -ge 10 ]; then
      echo -e "CRITICAL WARNING: LARGEST SHARD IS AT LEAST 10 TIMES LARGER THAN SECOND-LARGEST SHARD"
      exit 1
    fi
  fi
fi


###CLEAN UP
rm -rf ${SHARD_VCF_TMP}
