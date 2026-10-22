#!/bin/bash

# Collects QC data for SV VCF output by SV pipeline
# Summary, VCF-wide statistics only

set -e

###USAGE
usage(){
cat <<EOF

usage: collectQC.vcf_wide.sh [-h] [-F FAM] [-R] [-r REF] VCF SVTYPES OUTDIR

Helper tool to collect vcf-wide QC data for a VCF output by sv-pipeline

Positional arguments:
  VCF             VCF from sv-pipeline
  SVTYPES         List of SV types to evaluate. Two-column, tab-delimited file.
                  First column: sv type. Second column: HEX color for sv type.
  OUTDIR          Output directory for all QC data

Optional arguments:
  -h  HELP        Show this help message and exit
  -F  FAM         FAM file (for Mendelian violation analysis)
  -R  RESTRICT    Restrict analysis to samples present in FAM file
  -r  REFERENCE   Reference build (default: GRCh37/hg19)
  -q  QUIET       Silence all status updates

EOF
}


###PARSE ARGS
unset FAM
RESTRICT=0
REF="hg19"
QUIET=0
while getopts ":F:Rr:qh" opt; do
	case "$opt" in
		h)
			usage
			exit 0
			;;
    F)
      FAM=${OPTARG}
      ;;
    R)
      RESTRICT=1
      ;;
    r)
      REF=${OPTARG}
      ;;
    q)
      QUIET=1
      ;;
	esac
done
shift $(( ${OPTIND} - 1))
VCF=$1
SVTYPES=$2
OUTDIR=$3


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
#Checks assembly
if [ ${REF} != "hg19" ] && [ ${REF} != "GRCh37" ] && [ ${REF} != "h37" ]; then
  echo -e "\nALERT: VCF QC workflow currently only supports hg19/GRCh37. Exiting.\n"
  usage
  exit 0
fi


###SET BIN
BIN=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )


###PREP INPUT FILES
#Print status
if [ ${QUIET} == 0 ]; then
  echo -e "$( date ) - VCF QC STATUS: Preparing input files for QC metric generation"
fi
#Prep directories
QCTMP=`mktemp -d`
if ! [ -e ${OUTDIR} ]; then
  mkdir ${OUTDIR}
fi
if ! [ -e ${OUTDIR}/data ]; then
  mkdir ${OUTDIR}/data
fi
mkdir ${QCTMP}/perSample
#Unzip VCF
zcat ${VCF} > ${QCTMP}/input.vcf
#Gather SV types to process
cut -f1 ${SVTYPES} | sort | uniq > ${QCTMP}/svtypes.txt


###DETERMINE SAMPLES TO PROCESS
#Write list of all samples in VCF header
fgrep "#" ${QCTMP}/input.vcf | fgrep -v "##" | cut -f10- | \
sed 's/\t/\n/g' | sort | uniq > ${QCTMP}/vcf_samples.list
#If FAM file not provided, use all samples in VCF header
if [ -z ${FAM} ] || [ ! -s ${FAM} ]; then
  cp ${QCTMP}/vcf_samples.list \
  ${QCTMP}/analysis_samples.list
#Otherwise, action depends on RESTRICT
else
  #If RESTRICT is not selected, use all samples in VCF header
  if [ ${RESTRICT} == 0 ]; then
    cp ${QCTMP}/vcf_samples.list \
    ${QCTMP}/analysis_samples.list
  #If RESTRICT is selected, use only samples in FAM file
  else
    fgrep -v "#" ${FAM} | cut -f2 | sort | uniq | \
    fgrep -wf - ${QCTMP}/vcf_samples.list > \
    ${QCTMP}/analysis_samples.list
  fi
fi
#Get number of samples
nsamp=$( cat ${QCTMP}/analysis_samples.list | wc -l )
#Exit if no samples left for analysis
if [ -z ${nsamp} ] | [ ${nsamp} -eq 0 ]; then
  echo -e "$( date ) - VCF QC ERROR: No samples included for analysis. Check input VCF, fam file, and -R flag and rerun."
  exit 1
fi
#Move list of analysis samples to $OUTDIR
cp ${QCTMP}/analysis_samples.list ${OUTDIR}/


###PROCESS VCF
#Print status
if [ ${QUIET} == 0 ]; then
  echo -e "$( date ) - VCF QC STATUS: Processing VCF"
fi
#Run svtk vcf2bed with all info fields
svtk vcf2bed -i ALL \
${QCTMP}/input.vcf ${QCTMP}/vcf2bed_unsorted_unfiltered.bed

#Get genotype counts per variant for all types of SVs except for multi-alleleic CNVs
echo -e "VID\tnsamp_gt\thomref\thet\thomalt\tother\tunknown\tAC\tAN\tAF" > \
${QCTMP}/genotype_counts_per_SV.txt

grep -v ^# ${QCTMP}/input.vcf | grep -v "<CNV>" | \
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

#Get genotype counts per variant for multi-alleleic CNVs
grep -v ^# ${QCTMP}/input.vcf | grep "<CNV>" | \
awk -v FS="\t" -v OFS="\t" -v nsamp=${nsamp} '{
 unknown=0; homref=0; het=0; homalt=0; other=0;
 split($9,format,":"); for (f=1;f<=length(format);f++) { if (format[f] == "CN") CN=f };
 for(i=10;i<=NF;i++) {
  split($i,a,":");
    if (a[CN]==2) {homref++}
    else if (a[CN]==0 || a[CN]>=4) {homalt++}
    else {het++} };
  AC=(2*homalt)+het; AN=2*(nsamp-(unknown+other)); AF=AC/AN;
  print $3, nsamp-unknown, homref, het, homalt, other, unknown, AC, AN, AF }' >> \
${QCTMP}/genotype_counts_per_SV.txt


#Write header
head -n1 ${QCTMP}/vcf2bed_unsorted_unfiltered.bed > \
${QCTMP}/vcf2bed_cleaned.bed
#Sort & filter 
fgrep -wf ${QCTMP}/analysis_samples.list \
${QCTMP}/vcf2bed_unsorted_unfiltered.bed | \
sort -Vk1,1 -k2,2n -k3,3n >> \
${QCTMP}/vcf2bed_cleaned.bed
#Run Rscript to clean VCF stats
${BIN}/clean_vcf2bed_output.R \
  -N ${nsamp} \
  ${QCTMP}/vcf2bed_cleaned.bed \
  ${QCTMP}/genotype_counts_per_SV.txt \
  ${SVTYPES} \
  ${QCTMP}/vcf2bed_cleaned.simple.bed
#Reformat cleaned VCF stats
head -n1 ${QCTMP}/vcf2bed_cleaned.simple.bed | awk '{ print "#"$0 }' > \
${QCTMP}/vcf2bed_cleaned.simple.bed2
sed '1d' ${QCTMP}/vcf2bed_cleaned.simple.bed | sed 's/mCNV/MCNV/g' | \
sort -Vk1,1 -k2,2n -k3,3n | uniq >> \
${QCTMP}/vcf2bed_cleaned.simple.bed2
mv ${QCTMP}/vcf2bed_cleaned.simple.bed2 ${QCTMP}/vcf2bed_cleaned.simple.bed
bgzip -f ${QCTMP}/vcf2bed_cleaned.simple.bed
tabix -f ${QCTMP}/vcf2bed_cleaned.simple.bed.gz
#Moves cleaned VCF stats to $OUTDIR
cp ${QCTMP}/vcf2bed_cleaned.simple.bed.gz \
   ${OUTDIR}/data/VCF_sites.stats.bed.gz
cp ${QCTMP}/vcf2bed_cleaned.simple.bed.gz.tbi \
   ${OUTDIR}/data/VCF_sites.stats.bed.gz.tbi


###CLEAN UP
rm -rf ${QCTMP}

