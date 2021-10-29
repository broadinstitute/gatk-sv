#!/bin/bash

# Original author: Ryan Collins

# Intelligently shards a VCF prior to complex resolution (for parallelization)

# Subsetted to first half, to just output lists of which variants should go in each shard

set -Eeu -o pipefail

# ARGS defaults and hard-coded values
DIST=${DEFAULT_DIST:-1000}
RECIP=${DEFAULT_RECIP:-0.1}
MIN_LINES_PER_SHARD=${DEFAULT_MIN_LINES_PER_SHARD:-10}
MAX_SHARDS=${DEFAULT_MAX_SHARDS:-100}
NONCLUSTER_SHARDS=${DEFAULT_NONCLUSTER_SHARDS:-30}
PREFIX=${DEFAULT_PREFIX:-"vcf_shard"}
BREAKPOINT_PADDING=${DEFAULT_BREAKPOINT_PADDING:-5000}
IGNORE_SV_TYPES=${DEFAULT_IGNORE_SV_TYPES:-false}
ADD_SINGLE_REC=${DEFAULT_ADD_SINGLE_REC:-false}
SHARD_LARGE_CLUSTERS=${DEFAULT_SHARD_LARGE_CLUSTERS:-true}
SCRIPT_NAME=${SCRIPT_NAME:-$(basename "${BASH_SOURCE[0]}")}
#Set path to execution directory
BIN=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
#Set default output folder
OUTDIR=`pwd`

###USAGE
usage(){
cat <<EOF
usage: $SCRIPT_NAME [-h] [-D DIST] [-R RECIP] [-L MIN_LINES_PER_SHARD]
                                       [-S MAX_SHARDS] [-N NONCLUSTER_SHARDS] [-P PREFIX]
                                       [-T VCF_INDEX] VCF
Shard VCF prior to complex resolution (for parallelization)
Positional arguments:
  VCF                      VCF from sv-pipeline prior to svtk resolve
  OUTDIR                   Output directory for all shards
Optional arguments:
  -h  HELP                 Show this help message and exit
  -D  DIST                 Breakpoint distance used in clustering [default: $DIST]
  -R  RECIP                Reciprocal overlap used in clustering [default: $((100 * RECIP))%]
  -L  MIN_LINES_PER_SHARD  Minimum variants per shard [default: $MIN_LINES_PER_SHARD]
  -S  MAX_SHARDS           Maximum number of shards [default: $MAX_SHARDS]
  -N  NONCLUSTER_SHARDS    Add an additional N shards to evenly split the 
                           remaining records in the VCF that do not cluster with
                           any other variants [default: $NONCLUSTER_SHARDS]
  -P  PREFIX               String appended to the filename of each output VCF 
                           shard [default: $PREFIX]
  -B  BREAKPOINT_PADDING   Distance breakpoints are padded for clustering
                           [default: $BREAKPOINT_PADDING]
  -C  SHARD_LARGE_CLUSTERS Put exceptionally large clusters in their own shard
                           [default: $SHARD_LARGE_CLUSTERS]
  -I  IGNORE_SV_TYPES      If true, ignore sv types when clustering [default: $IGNORE_SV_TYPES]
  -A  ADD_SINGLE_REC       If true, add single record non-CNVs to candidate
                           complex clusters [default: $ADD_SINGLE_REC]
  -O  OUTDIR               Output directory [default: pwd]
  -T  VCF_index            Tabix index corresponding to input VCF 
                           [default: VCF path + .tbi]
Notes:
  1. The total number of shards created will be no more 
     than MAX_SHARDS + NONCLUSTER_SHARDS.
EOF
}

###PARSE ARGS
while getopts ":D:R:L:S:N:P:B:C:I:A:O:T:h" opt; do
  case "$opt" in
    h)
      usage
      exit 0
      ;;
    D)
      DIST=${OPTARG}
      ;;
    R)
      RECIP=${OPTARG}
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
    B)
      BREAKPOINT_PADDING=${OPTARG}
      ;;
    C)
      SHARD_LARGE_CLUSTERS=${OPTARG}
      ;;
    I)
      IGNORE_SV_TYPES=${OPTARG}
      ;;
    A)
      ADD_SINGLE_REC=${OPTARG}
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
VCF_INDEX=${VCF_INDEX:-"${VCF}.tbi"}


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
if ! file ${VCF} | fgrep -q "gzip"; then
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
#Prep string to pass to vcfcluster
if IGNORE_SV_TYPES; then
    IGNORE_SV_TYPES_ARG="--ignore-svtypes"
else
    IGNORE_SV_TYPES_ARG=""
fi

###IDENTIFY CANDIDATE CLUSTERS
#Make temporary directory
SHARD_VCF_TMP=`mktemp -d`
#Cut vcf to single-sample for improved clustering speed
zcat ${VCF} | cut -f1-10 | bgzip -c > ${SHARD_VCF_TMP}/single_sample_input.vcf.gz
#Identify all candidate complex variant clusters (generous 1kb clustering)
svtk vcfcluster \
  -d ${DIST} \
  -f ${RECIP} \
  -p candidate_complex_clusters \
  ${IGNORE_SV_TYPES_ARG} \
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
mem_idx=$( head -n1 ${SHARD_VCF_TMP}/input_vcf.clustered.bed \
           | sed 's/\t/\n/g' | awk '{ if ($1=="MEMBERS") print NR }' )
awk -v idx=${mem_idx} -v OFS="\t" '$idx ~ /,/ { print $1, $2, $3, $idx }' \
    ${SHARD_VCF_TMP}/input_vcf.clustered.bed \
  | (grep -v "^#" || printf "") \
  > ${SHARD_VCF_TMP}/candidate_complex_clusters.bed
if ADD_SINGLE_REC; then
  #Add all non-CNV single-record variants
  class_idx=$( head -n1 ${SHARD_VCF_TMP}/input_vcf.clustered.bed | \
              sed 's/\t/\n/g' | awk '{ if ($1=="SVTYPE") print NR }' )
  awk -v idx=${class_idx} -v mem_idx=${mem_idx} -v OFS="\t" \
    '$mem_idx !~ /,/ { print $1, $2, $3, $idx, $mem_idx }' \
    ${SHARD_VCF_TMP}/input_vcf.clustered.bed \
    | awk -v OFS="\t" '$4 !~ /DEL|DUP|CNV|MCNV|mCNV/ { print $1, $2, $3, $5 }' \
    | (grep -v "^#" || printf "") \
    >> ${SHARD_VCF_TMP}/candidate_complex_clusters.bed
fi

#Get min/max coordinates of all variants in list of VIDs
{
  zcat ${VCF} \
    | (grep "^#" || printf "") \
    | cut -f1-10;
  cut -f4 ${SHARD_VCF_TMP}/candidate_complex_clusters.bed \
    | sed 's/\,/\n/g' \
    | sort -Vk1,1 \
    | uniq \
    | (fgrep -wf - <(zcat ${VCF}) || printf "") \
    | cut -f1-10;
} \
  | svtk vcf2bed --no-samples /dev/stdin /dev/stdout \
  > ${SHARD_VCF_TMP}/candidate_complex_clusters.variant_coordinates.bed

###ONLY PERFORM CLUSTER-BASED SHARDING IF ANY VARIANTS PREDICTED TO CLUSTER
if grep -vq "^#" ${SHARD_VCF_TMP}/candidate_complex_clusters.variant_coordinates.bed; then
  #Split into breakpoints and pad all breakpoints by ±BREAKPOINT_PADDING
  ###DETERMINE SET OF NONREDUNDANT INTERVALS FOR ALL CLUSTERS
  grep -v "^#" ${SHARD_VCF_TMP}/candidate_complex_clusters.variant_coordinates.bed \
    | awk -v OFS="\t" -v buffer=$BREAKPOINT_PADDING \
        '{ print $1, $2-buffer, $2+buffer, $4;
           print $1, $3-buffer, $3+buffer, $4 }' \
    | awk -v OFS="\t" '{ if ($2<1) $2=1; print $1, $2, $3, $4 }' \
    | sort -Vk1,1 -k2,2n -k3,3n \
    | bedtools merge -i - -c 4 -o distinct \
    > ${SHARD_VCF_TMP}/breakpoint_intervals.bed
  #Iterate over breakpoint intervals and write list of maximum nonredundant intervals
  in_cluster=`mktemp`
  remaining=`mktemp`
  cp ${SHARD_VCF_TMP}/breakpoint_intervals.bed ${remaining}
  while read chr start end VIDs; do
    #Get all lines associated with current VIDs
    echo -e "${VIDs}" \
      | sed 's/,/\n/g' \
      | (fgrep -wf - ${remaining} || printf "") \
      > ${in_cluster}
    #Only run if at least one line added to ${in_cluster}
    if [ -s "${in_cluster}" ]; then
      #Exclude all lines in ${in_cluster} from ${remaining}
      bedtools intersect -v -a ${remaining} -b ${in_cluster} > ${remaining}2
      mv ${remaining}2 ${remaining}
      #Iterate until no more related VIDs are present in ${remaining}
      until ! cut -f4 ${in_cluster} | sed 's/\,/\n/g' | fgrep -qwf - ${remaining}; do
        #Add new lines to ${in_cluster}
        cut -f4 ${in_cluster} \
          | sed 's/\,/\n/g' \
          | (fgrep -wf - ${remaining} || printf "") \
          >> ${in_cluster}
        #Exclude all lines in ${in_cluster} from ${remaining}
        bedtools intersect -v -a ${remaining} -b ${in_cluster} > ${remaining}2
        mv ${remaining}2 ${remaining}
      done
      #Write out final interval
      for wrapper in 1; do
        #Print list of coordinates
        cut -f1-3 ${in_cluster} \
          | sort -Vk1,1 -k2,2 -k3,3 \
          | bedtools merge -i - \
          | awk '{ print $1":"$2"-"$3 }' \
          | paste -s -d\;
        #Print list of involved VIDs
        cut -f4 ${in_cluster} \
          | sed 's/,/\n/g' \
          | sort \
          | uniq \
          | paste -s -d,
      done | paste -s
    fi
  done < ${SHARD_VCF_TMP}/breakpoint_intervals.bed \
    > ${SHARD_VCF_TMP}/complex_intervals_to_test.final.txt

  if SHARD_LARGE_CLUSTERS; then
      #Pull out exceptionally large clusters to the side to be placed in their own shards
      while read ints VIDs; do
        if [ $( echo ${VIDs} | sed 's/,/\n/g' | wc -l ) -ge ${MIN_LINES_PER_SHARD} ]; then
          echo -e "${ints}\t${VIDs}"
        fi
      done < ${SHARD_VCF_TMP}/complex_intervals_to_test.final.txt \
        > ${SHARD_VCF_TMP}/large_intervals_to_test.final.txt
      if [ -s "${SHARD_VCF_TMP}/large_intervals_to_test.final.txt" ]; then
        cut -f2 ${SHARD_VCF_TMP}/large_intervals_to_test.final.txt \
          | sed 's/,/\n/g' \
          | (fgrep -wvf - ${SHARD_VCF_TMP}/complex_intervals_to_test.final.txt || printf "") \
          > ${SHARD_VCF_TMP}/complex_intervals_to_test.final.txt2
        mv ${SHARD_VCF_TMP}/complex_intervals_to_test.final.txt2 \
          ${SHARD_VCF_TMP}/complex_intervals_to_test.final.txt
      fi
  fi

  ###DETERMINE COORDINATES FOR EACH SHARD
  #Split variants into shards based on number of variants
  #If total number of intervals/MAX_SHARDS < MIN_LINES_PER_SHARD, evenly split into MIN_LINES_PER_SHARD sites per shard
  if [ $(( $( wc -l < "${SHARD_VCF_TMP}/complex_intervals_to_test.final.txt" ) / ${MAX_SHARDS} )) -lt ${MIN_LINES_PER_SHARD} ]; then
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
  #Determine number of shards
  n_shards=$( find ${SHARD_VCF_TMP} -name "${PREFIX}.shard_intervals_*" | wc -l )
  #Writes exceptionally large clusters to their own shards
  if [ -s "${SHARD_VCF_TMP}/large_intervals_to_test.final.txt" ]; then
    while read ints VIDs; do
      ((++n_shards))
      echo -e "${ints}\t${VIDs}" > ${SHARD_VCF_TMP}/${PREFIX}.shard_intervals_${n_shards}
    done < ${SHARD_VCF_TMP}/large_intervals_to_test.final.txt
  fi
  #Reformat interval shards
  for i in $( seq 1 ${n_shards} ); do
    cut -f1 ${SHARD_VCF_TMP}/${PREFIX}.shard_intervals_${i} \
      | sed -e 's/\;/\n/g' -e 's/\:/\t/g' -e 's/\-/\t/g' \
      | sort -Vk1,1 -k2,2n -k3,3n \
      | bedtools merge -i - \
      > ${SHARD_VCF_TMP}/${PREFIX}.shard_intervals_${i}.bed
    rm ${SHARD_VCF_TMP}/${PREFIX}.shard_intervals_${i}
  done


  ###SHARD CLUSTERABLE VCF
  #Convert full, original VCF to BED
  svtk vcf2bed --no-samples ${VCF} int.bed
  #Harrison's patch for sharding
  awk '{ if ($1!~"#") { print $1,$2,$2+1,$4,$5;
                        print $1,$3-1,$3,$4,$5 }
         else print }' OFS='\t' int.bed \
    > ${SHARD_VCF_TMP}/input_vcf.vcf2bed.bed
  rm int.bed

  #Create exclusion list of VIDs already used in earlier shards
  touch ${SHARD_VCF_TMP}/used_VIDs.tmp
  #Iterate over all sharded intervals
  for i in $( seq 1 ${n_shards} ); do
    if [ ${i} -gt 1 ]; then
      cat ${OUTDIR}/${PREFIX}.shard_*.VIDs.list \
        | sort \
        | uniq \
        > ${SHARD_VCF_TMP}/used_VIDs.tmp
    else
      touch ${SHARD_VCF_TMP}/used_VIDs.tmp
    fi
    #Get list of IDs to be used in shard
    bedtools intersect \
        -a ${SHARD_VCF_TMP}/input_vcf.vcf2bed.bed \
        -b ${SHARD_VCF_TMP}/${PREFIX}.shard_intervals_${i}.bed \
        | cut -f4 \
        | (fgrep -wvf ${SHARD_VCF_TMP}/used_VIDs.tmp || printf "")\
        | sort \
        | uniq \
        > ${OUTDIR}/${PREFIX}.shard_${i}.VIDs.list

    #Update exclusion list of VIDs already used in earlier shards
    cat ${OUTDIR}/${PREFIX}.shard_*.VIDs.list \
      | sort \
      | uniq \
      > ${SHARD_VCF_TMP}/used_VIDs.tmp
  done

  #Write list of eligible VIDs
  (grep -v "^#" ${SHARD_VCF_TMP}/input_vcf.vcf2bed.bed || printf "") \
    | cut -f4 \
    | (fgrep -wvf ${SHARD_VCF_TMP}/used_VIDs.tmp || printf "") \
    | sort \
    | uniq \
    > ${SHARD_VCF_TMP}/remaining_VIDs.list
else
  n_shards=0
  zcat ${VCF} \
    | (grep -v "^#" || printf "") \
    | cut -f3 \
    | sort \
    | uniq \
    > ${SHARD_VCF_TMP}/remaining_VIDs.list
fi


###SHARD NONCLUSTERABLE VCF
#Shard remainder intervals into no more than $NONCLUSTER_SHARDS shards
#If total number of variants/NONCLUSTER_SHARDS < MIN_LINES_PER_SHARD, evenly split into MIN_LINES_PER_SHARD sites per shard
if [ -s "${SHARD_VCF_TMP}/remaining_VIDs.list" ]; then
  if [ $(( $(wc -l < ${SHARD_VCF_TMP}/remaining_VIDs.list) / ${NONCLUSTER_SHARDS} )) -lt ${MIN_LINES_PER_SHARD} ]; then
    ${BIN}/evenSplitter.R \
      -L ${MIN_LINES_PER_SHARD} \
      ${SHARD_VCF_TMP}/remaining_VIDs.list \
      ${SHARD_VCF_TMP}/${PREFIX}.remaining_variants_
  #Otherwise, split into MAX_SHARDS evenly-sized shards
  else
    ${BIN}/evenSplitter.R \
      -S ${NONCLUSTER_SHARDS} \
      ${SHARD_VCF_TMP}/remaining_VIDs.list \
      ${SHARD_VCF_TMP}/${PREFIX}.remaining_variants_
  fi
  n_nonclusterable_shards=$( find ${SHARD_VCF_TMP} -name "${PREFIX}.remaining_variants_*" | wc -l )
  #Iterate over all sharded variant lists
  for i in $( seq 1 ${n_nonclusterable_shards} ); do
    idx=$(( ${n_shards} + ${i} ))
    mv ${SHARD_VCF_TMP}/${PREFIX}.remaining_variants_${i} \
      ${OUTDIR}/${PREFIX}.shard_${idx}.VIDs.list
  done
fi

###CLEAN UP
rm -rf ${SHARD_VCF_TMP}

1>&2 echo "$(basename $0): Success"
