#!/bin/bash

# Collects QC data for SV VCF output by SV pipeline

set -e

###USAGE
usage(){
cat <<EOF

usage: collectQC.sh [-h] [-F FAM] [-R] VCF SVTYPES OUTDIR

Helper tool to collect QC data for a VCF output by sv-pipeline

Positional arguments:
  VCF             VCF from sv-pipeline
  SVTYPES         List of SV types to evaluate. Two-column, tab-delimited file.
                  First column: sv type. Second column: HEX color for sv type.
  OUTDIR          Output directory for all QC data

Optional arguments:
  -h  HELP        Show this help message and exit
  -F  FAM         FAM file (for Mendelian violation analysis)
  -R  RESTRICT    Restrict analysis to samples present in FAM file

EOF
}


###PARSE ARGS
unset FAM
RESTRICT=0
while getopts ":F:Rh" opt; do
	case "$opt" in
		h)
			usage
			exit 0
			;;
    F)
      FAM=${OPTART}
      ;;
    R)
      RESTRICT=1
      ;;
	esac
done
shift $(( ${OPTIND} - 1))
VCF=$1
SVTYPES=$2
OUTDIR=$2


###PROCESS ARGS & OPTIONS
#Check for required input
if [ -z ${VCF} ]; then
  echo -e "\nERROR: input VCF not specified\n"
  usage
  exit 0
fi
if [ -s ${VCF} ]; then
  echo -e "\nERROR: input VCF either empty or not found\n"
  usage
  exit 0
fi
if [ -z ${OUTDIR} ]; then
  echo -e "\nERROR: output directory not specified\n"
  usage
  exit 0
fi
#Checks for FAM file if RESTRICT is optioned
if [ ${RESTRICT} != "0" ]; then
  if [ -z ${FAM} ] || [ ! -s ${FAM} ]; then
    echo -e "\nERROR: FAM file must be provided if RESTRICT option is selected\n"
    usage
    exit 0
  fi
fi


###SET BIN
BIN=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )


###PREP INPUT FILES
#Prep directories
QCTMP=`mktemp -d` && cd ${QCTMP}
if ! [ -e ${OUTDIR} ]; then
  mkdir ${OUTDIR}
fi
cd ${QCTMP}
mkdir ${QCTMP}/perSample
mkdir ${QCTMP}/perFamily
#Unzip VCF, if gzipped
if [ $( file ${VCF} | fgrep " gzip " | wc -l ) -gt 0 ]; then
  zcat ${VCF} > ${QCTMP}/input.vcf
else
  cp ${VCF} ${QCTMP}/input.vcf
fi
#Gather SV types to process
cut -f1 ${SVTYPES} | sort | uniq > ${QCTMP}/svtypes.txt


###DETERMINE SAMPLES TO PROCESS
#If RESTRICT option is not set, write list of all samples in VCF header
if [ ${RESTRICT} == 0 ]; then
  fgrep "#" ${QCTMP}/input.vcf | fgrep -v "##" | cut -f10- | \
  sed 's/\t/\n/g' | sort | uniq > ${QCTMP}/vcf_samples.list
fi
#If FAM file not provided, use all samples in VCF header
if [ -z ${FAM} ] || [ ! -s ${FAM} ]; then
  cp ${QCTMP}/vcf_samples.list \
  ${QCTMP}/analysis_samples.list
#Otherwise, action depends on RESTRICT
else
  #If RESTRICT is not selected, take the union of samples in the FAM file & VCF
  if [ ${RESTRICT} == 0 ]; then
    cat <( fgrep -v "#" ${FAM} | cut -f2 ) ${QCTMP}/vcf_samples.list | \
    sort | uniq > ${QCTMP}/analysis_samples.list
  #If RESTRICT is selected, use only samples in FAM file
  else
    fgrep -v "#" ${FAM} | cut -f2 | sort | uniq > \
    ${QCTMP}/analysis_samples.list
  fi
fi
#Get number of samples
nsamp=$( cat ${QCTMP}/analysis_samples.list | wc -l )


###PROCESS VCF
#Run svtk vcf2bed with all info fields
svtk vcf2bed --split-bnd -i ALL \
${QCTMP}/input.vcf ${QCTMP}/vcf2bed_unsorted_unfiltered.bed
#Get genotype counts per variant
echo -e "VID\tnsamp_gt\thomref\thet\thomalt\tother\tunknown\tAC\tAN\tAF" > \
${QCTMP}/genotype_counts_per_SV.txt
grep -v ^# ${QCTMP}/input.vcf | \
awk -v FS="\t" -v OFS="\t" -v nsamp=${nsamp} '{
 unknown=0; homref=0; het=0; homalt=0; other=0; for(i=10;i<=NF;i++) {
  split($i,a,":"); split(a[1],GT,"[/|]");
  if (GT[1]=="." && GT[2]==".") {unknown++}
  else if (GT[1]==0 && GT[2]==0) {homref++}
  else if (GT[1]==1 && GT[2]==1) {homalt++}
  else if ((GT[1]==0 && GT[2]==1) || (GT[1]==1 && GT[2]==0)) {het++}
  else {other++} };
  if (other > 0 || (nsamp==other+unknown)) {AC="NA"; AN="NA"; AF="NA"}
  else {AC=(2*homalt)+het; AN=2*(nsamp-(unknown+other)); AF=AC/AN};
  print $3, nsamp-unknown, homref, het, homalt, other, unknown, AC, AN, AF }' >> \
${QCTMP}/genotype_counts_per_SV.txt
#Write header
head -n1 ${QCTMP}/vcf2bed_unsorted_unfiltered.bed > \
${QCTMP}/vcf2bed_cleaned.bed
#Sort & filter 
fgrep -wf ${QCTMP}/analysis_samples.list \
${QCTMP}/vcf2bed_unsorted_unfiltered.bed | \
fgrep -wf ${QCTMP}/svtypes.txt | \
sort -Vk1,1 -k2,2n -k3,3n >> \
${QCTMP}/vcf2bed_cleaned.bed
#Run Rscript to clean VCF stats
${BIN}/clean_vcf2bed_output.R \
-N ${nsamp} \
${QCTMP}/vcf2bed_cleaned.bed \
${QCTMP}/genotype_counts_per_SV.txt \
${QCTMP}/vcf2bed_cleaned.simple.bed


###PLOT VCF SUMMARY DISTRIBUTIONS
${BIN}/plot_sv_vcf_distribs.R \
-N ${nsamp} -S ${SVTYPES} \
${QCTMP}/vcf2bed_cleaned.simple.bed \
${OUTDIR}/


###COLLECT LISTS OF VARIANTS PER SAMPLE
#Cut list of samples & variant ID for faster grepping
fgrep -v "#" ${QCTMP}/vcf2bed_cleaned.bed | cut -f4,6 | sed 's/\,/\t/g' > \
${QCTMP}/VIDs_samples.list
#Iterate over samples & write list of VIDs to file
while read ID; do
  fgrep -w ${ID} ${QCTMP}/VIDs_samples.list | cut -f1 | sort -Vk1,1 | uniq > \
  ${QCTMP}/perSample/${ID}.VIDs.list
done < ${QCTMP}/analysis_samples.list


###CLEAN UP
rm -rf ${QCTMP}

