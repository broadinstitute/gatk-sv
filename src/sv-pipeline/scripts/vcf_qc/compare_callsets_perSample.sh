#!/bin/bash

# Compares per-sample concordance between two SV callsets

set -e

###USAGE
usage(){
cat <<EOF

usage: compare_callsets_perSample.sh [options] SET1.tar.gz SET2.tar.gz SAMPLES CONTIGS OUTDIR

Helper tool to compare the sensitivity of one SV callset vs another

Positional arguments:
  SET1.tar.gz            SV callset to be benchmarked
  SET2.tar.gz            Ground truth SV callset to use for benchmarking
  SAMPLES                List of samples to be considered
  CONTIGS                List of contigs to evaluate
  OUTDIR                 Output directory for per-sample comparison results

Optional arguments:
  -h  HELP        Show this help message and exit
  -d  DISTANCE    Maximum distance between breakpoints during comparisons (default: 250bp)
  -p  PREFIX      Prefix for benchmarking variant IDs (default: Benchmarking_SV)
  -q  QUIET       Silence all status updates

Notes:
  1) Both SET1 and SET2 are expected to be provided as tarred, gzipped directories
     containing one identically-formatted file per sample.
  2) SET1/SET2 are expected to be BED3+ formatted files
  3) SET1 must have at least the following five columns, in order:
     chr, start, end, SV ID, SV type
  4) The last column of SET1 must contain allele frequency
  5) SET2 must have exactly six columns as follows, in order:
     chr, start, end, SV type, SV size, allele frequency
  6) Assumes all per-sample files in SET1 and SET2 follow the same file naming convention.
     Filename scheme: [sample_ID].[callset_name].[other_suffixes].SV_calls.bed.gz

EOF
}


###PARSE ARGS
PREFIX="Benchmarking_SV"
DIST=250
QUIET=0
while getopts ":p:d:qh" opt; do
	case "$opt" in
		h)
			usage
			exit 1
			;;
    p)
      PREFIX=${OPTARG}
      ;;
    d)
      DIST=${OPTARG}
      ;;
    q)
      QUIET=1
      ;;
	esac
done
shift $(( ${OPTIND} - 1))
SET1=$1
SET2=$2
SAMPLES=$3
CONTIGS=$4
OUTDIR=$5


###SET BIN
BIN=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )


###PROCESS ARGS & OPTIONS
#Check for required input SET1
if [ -z ${SET1} ]; then
  echo -e "\ncompare_callsets_perSample.sh ERROR: input SET1.tar.gz not specified\n"
  usage
  exit 1
elif ! [ -s ${SET1} ]; then
  echo -e "\ncompare_callsets_perSample.sh ERROR: input SET1.tar.gz either empty or not found\n"
  usage
  exit 1
fi
#Check for required input SET2
if [ -z ${SET2} ]; then
  echo -e "\ncompare_callsets_perSample.sh ERROR: input SET2.tar.gz not specified\n"
  usage
  exit 1
elif ! [ -s ${SET2} ]; then
  echo -e "\ncompare_callsets_perSample.sh ERROR: input SET2.tar.gz either empty or not found\n"
  usage
  exit 1
fi
#Check for required list of samples
if [ -z ${SAMPLES} ]; then
  echo -e "\ncompare_callsets_perSample.sh ERROR: input sample list not specified\n"
  usage
  exit 1
elif ! [ -s ${SAMPLES} ]; then
  echo -e "\ncompare_callsets_perSample.sh ERROR: input sample list either empty or not found\n"
  usage
  exit 1
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
#Checks output directory
if [ -z ${OUTDIR} ]; then
  echo -e "\ncompare_callsets_perSample.sh ERROR: output directory not specified\n"
  usage
  exit 1
fi
if ! [ -e ${OUTDIR} ]; then
  mkdir ${OUTDIR}
fi
#Checks prefix
if [ -z ${PREFIX} ]; then
  PREFIX="Benchmarking_SV"
fi


###PREP INPUT FILES
OVRTMP=`mktemp -d`
#Print status
if [ ${QUIET} == 0 ]; then
  echo -e "$( date ) - PER-SAMPLE COMPARISON STATUS: Preparing input files."
fi
#Untar SET1
mkdir ${OVRTMP}/SET1/
mkdir ${OVRTMP}/SET1/original/
tar -xzvf ${SET1} --directory ${OVRTMP}/SET1/original/
#Get list of all samples in SET1
find ${OVRTMP}/SET1/original/ -name "*.SV_calls.bed.gz" > \
${OVRTMP}/SET1/original/SET1_paths.list
sed 's/\//\t/g' ${OVRTMP}/SET1/original/SET1_paths.list | \
awk '{ print $NF }' | sed 's/\./\t/g' | awk '{ print $1 }' > \
${OVRTMP}/SET1/original/SET1_samples.list
paste ${OVRTMP}/SET1/original/SET1_samples.list \
  ${OVRTMP}/SET1/original/SET1_paths.list > \
  ${OVRTMP}/SET1/original/SET1_samples_and_paths.list
#Untar SET2
mkdir ${OVRTMP}/SET2/
mkdir ${OVRTMP}/SET2/original/
tar -xzvf ${SET2} --directory ${OVRTMP}/SET2/original/
#Get list of all samples in SET2
#Print status
if [ ${QUIET} == 0 ]; then
  echo -e "$( date ) - PER-SAMPLE COMPARISON STATUS: Determining samples present in SET1, SET2, and sample list."
fi
find ${OVRTMP}/SET2/original/ -name "*.SV_calls.bed.gz" > \
${OVRTMP}/SET2/original/SET2_paths.list
sed 's/\//\t/g' ${OVRTMP}/SET2/original/SET2_paths.list | \
awk '{ print $NF }' | sed 's/\./\t/g' | awk '{ print $1 }' > \
${OVRTMP}/SET2/original/SET2_samples.list
paste ${OVRTMP}/SET2/original/SET2_samples.list \
  ${OVRTMP}/SET2/original/SET2_paths.list > \
  ${OVRTMP}/SET2/original/SET2_samples_and_paths.list
#Find list of samples in SET1, SET2, and SAMPLES
fgrep -wf ${OVRTMP}/SET1/original/SET1_samples.list \
          ${OVRTMP}/SET2/original/SET2_samples.list | \
  fgrep -wf ${SAMPLES} > \
${OVRTMP}/int_samples.list || true


###ONLY RUN COMPARISONS IF ANY SAMPLES FOUND IN SET1, SET2, AND SAMPLE LIST
nsamps=$( cat ${OVRTMP}/int_samples.list | wc -l )
if [ ${nsamps} -eq 0 ]; then
  #Print status
  if [ ${QUIET} == 0 ]; then
    echo -e "$( date ) - PER-SAMPLE COMPARISON STATUS: No matching samples found in SET1, SET2, and sample list."
    echo -e "$( date ) - PER-SAMPLE COMPARISON STATUS: Exiting..."
  fi
  exit 0
else
  #Print status
  if [ ${QUIET} == 0 ]; then
    echo -e "$( date ) - PER-SAMPLE COMPARISON STATUS: Found ${nsamps} samples matching between SET1, SET2, and sample list."
    echo -e "$( date ) - PER-SAMPLE COMPARISON STATUS: Starting sample comparisons..."
  fi

  #Create list of samples with SET1 and SET2 paths
  while read ID; do
    echo -e "${ID}"
    awk -v ID=${ID} '{ if ($1==ID) print $2 }' \
    ${OVRTMP}/SET1/original/SET1_samples_and_paths.list
    awk -v ID=${ID} '{ if ($1==ID) print $2 }' \
    ${OVRTMP}/SET2/original/SET2_samples_and_paths.list
  done < ${OVRTMP}/int_samples.list | paste - - - > \
  ${OVRTMP}/int_samples_and_paths.list

  #Iterate over samples and run compare_callsets.sh for sensitivity & specificity
  i=0
  j=0
  while read ID SET1s SET2s; do
    #Clean callsets for sensitivity analysis
    zcat ${SET1s} | awk -v OFS="\t" '{ if ($3!="end" && $3<$2) $3=$2; print }' > \
    ${OVRTMP}/${ID}.SET1.sens.cleaned.bed
    bgzip -f ${OVRTMP}/${ID}.SET1.sens.cleaned.bed
    zcat ${SET2s} | awk -v OFS="\t" '{ if ($3!="end" && $3<$2) $3=$2; print }' > \
    ${OVRTMP}/${ID}.SET2.sens.cleaned.bed
    bgzip -f ${OVRTMP}/${ID}.SET2.sens.cleaned.bed

    #Run sensitivity analysis
    ${BIN}/compare_callsets.sh \
      -d ${DIST} \
      -p ${ID}_${PREFIX}_sensitivity \
      -O ${OUTDIR}/${ID}.sensitivity.bed \
      ${OVRTMP}/${ID}.SET1.sens.cleaned.bed.gz \
      ${OVRTMP}/${ID}.SET2.sens.cleaned.bed.gz \
      ${CONTIGS}
    bgzip -f ${OUTDIR}/${ID}.sensitivity.bed
    rm ${OVRTMP}/${ID}.SET1.sens.cleaned.bed.gz \
       ${OVRTMP}/${ID}.SET2.sens.cleaned.bed.gz

    #Clean callsets for specificity analysis
    echo -e "#chr\tstart\tend\tsvtype\tlength\tAF" > \
    ${OVRTMP}/${ID}.SET1.spec.cleaned.bed
    zcat ${SET1s} | fgrep -v "#" | awk -v OFS="\t" \
    '{ if ($3!="end" && $3<$2) $3=$2; print $1, $2, $3, $5, $6, $NF }' >> \
    ${OVRTMP}/${ID}.SET1.spec.cleaned.bed
    bgzip -f ${OVRTMP}/${ID}.SET1.spec.cleaned.bed
    echo -e "#chr\tstart\tend\tVID\tsvtype\tlength\tAF" > \
    ${OVRTMP}/${ID}.SET2.spec.cleaned.bed
    zcat ${SET2s} | fgrep -v "#" | awk -v OFS="\t" -v ID=${ID} -v PREFIX=${PREFIX} \
    '{ if ($3!="end" && $3<$2) $3=$2; print $1, $2, $3, ID"_"PREFIX"_"NR, $4, $5, $6 }' >> \
    ${OVRTMP}/${ID}.SET2.spec.cleaned.bed
    bgzip -f ${OVRTMP}/${ID}.SET2.spec.cleaned.bed

    #Run specificity analysis
    ${BIN}/compare_callsets.sh \
      -d ${DIST} \
      -p ${ID}_${PREFIX}_specificity \
      -O ${OUTDIR}/${ID}.specificity.bed \
      ${OVRTMP}/${ID}.SET2.spec.cleaned.bed.gz \
      ${OVRTMP}/${ID}.SET1.spec.cleaned.bed.gz \
      ${CONTIGS}
    bgzip -f ${OUTDIR}/${ID}.specificity.bed
    rm ${OVRTMP}/${ID}.SET1.spec.cleaned.bed.gz \
       ${OVRTMP}/${ID}.SET2.spec.cleaned.bed.gz

    #Report counter, if relevant
    i=$(( ${i} + 1 ))
    j=$(( ${j} + 1 ))
    if [ ${i} -eq 10 ]; then
      i=0
      if [ ${QUIET} == 0 ]; then
        echo -e "$( date ) - PER-SAMPLE COMPARISON STATUS: Finished comparisons for ${j}/${nsamps} samples..."
      fi
    fi
  done < ${OVRTMP}/int_samples_and_paths.list
fi


###CLEAN UP
rm -rf ${OVRTMP}

