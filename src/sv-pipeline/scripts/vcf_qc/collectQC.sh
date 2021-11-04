#!/bin/bash

# Collects QC data for SV VCF output by SV pipeline

set -e

###USAGE
usage(){
cat <<EOF

usage: collectQC.sh [-h] [-F FAM] [-R] VCF SVTYPES BENCHDIR OUTDIR

Helper tool to collect QC data for a VCF output by sv-pipeline

Positional arguments:
  VCF             VCF from sv-pipeline
  SVTYPES         List of SV types to evaluate. Two-column, tab-delimited file.
                  First column: sv type. Second column: HEX color for sv type.
  BENCHDIR        Directory containing benchmark archives
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
BENCHDIR=$3
OUTDIR=$4


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
if [ ${REF} != "hg19" ] || [ ${REF} != "GRCh37" ] || [ ${REF} != "h37" ]; then
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
QCTMP=`mktemp -d` && cd ${QCTMP}
if ! [ -e ${OUTDIR} ]; then
  mkdir ${OUTDIR}
fi
if ! [ -e ${OUTDIR}/data ]; then
  mkdir ${OUTDIR}/data
fi
cd ${QCTMP}
mkdir ${QCTMP}/perSample
#Unzip VCF, if gzipped
if [ $( file ${VCF} | fgrep " gzip " | wc -l ) -gt 0 ]; then
  zcat ${VCF} > ${QCTMP}/input.vcf
else
  cp ${VCF} ${QCTMP}/input.vcf
fi
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
${SVTYPES} \
${QCTMP}/vcf2bed_cleaned.simple.bed
#Reformat cleaned VCF stats
head -n1 ${QCTMP}/vcf2bed_cleaned.simple.bed | awk '{ print "#"$0 }' > \
${QCTMP}/vcf2bed_cleaned.simple.bed2
sed '1d' ${QCTMP}/vcf2bed_cleaned.simple.bed | sort -Vk1,1 -k2,2n -k3,3n | uniq >> \
${QCTMP}/vcf2bed_cleaned.simple.bed2
mv ${QCTMP}/vcf2bed_cleaned.simple.bed2 ${QCTMP}/vcf2bed_cleaned.simple.bed
bgzip -f ${QCTMP}/vcf2bed_cleaned.simple.bed
tabix -f ${QCTMP}/vcf2bed_cleaned.simple.bed.gz
#Moves cleaned VCF stats to $OUTDIR
cp ${QCTMP}/vcf2bed_cleaned.simple.bed.gz \
   ${OUTDIR}/data/VCF_sites.stats.bed.gz
cp ${QCTMP}/vcf2bed_cleaned.simple.bed.gz.tbi \
   ${OUTDIR}/data/VCF_sites.stats.bed.gz.tbi


###GATHER EXTERNAL BENCHMARKING
#Print status
if [ ${QUIET} == 0 ]; then
  echo -e "$( date ) - VCF QC STATUS: Starting external benchmarking"
fi
#1000G (Sudmant) with allele frequencies
for pop in ALL AFR AMR EAS EUR SAS; do
  ${BIN}/compare_callsets.sh \
    -O ${QCTMP}/1000G_Sudmant.SV.${pop}.overlaps.bed \
    -p 1000G_Sudmant_${pop}_Benchmarking_SV \
    ${QCTMP}/vcf2bed_cleaned.simple.bed.gz \
    ${BENCHDIR}/1000G_Sudmant.SV.${pop}.bed.gz
done
cp ${QCTMP}/1000G_Sudmant.SV.${pop}.overlaps.bed \
   ${OUTDIR}/data/1000G_Sudmant.SV.overlaps.bed
bgzip -f ${OUTDIR}/data/1000G_Sudmant.SV.overlaps.bed
tabix -f ${OUTDIR}/data/1000G_Sudmant.SV.overlaps.bed.gz
#ASC (Werling) with carrier frequencies
for pop in ALL EUR OTH; do
  ${BIN}/compare_callsets.sh -C \
    -O ${QCTMP}/ASC_Werling.SV.${pop}.overlaps.bed \
    -p ASC_Werling_${pop}_Benchmarking_SV \
    ${QCTMP}/vcf2bed_cleaned.simple.bed.gz \
    ${BENCHDIR}/ASC_Werling.SV.${pop}.bed.gz
done
cp ${QCTMP}/ASC_Werling.SV.ALL.overlaps.bed \
   ${OUTDIR}/data/ASC_Werling.SV.overlaps.bed
bgzip -f ${OUTDIR}/data/ASC_Werling.SV.overlaps.bed
tabix -f ${OUTDIR}/data/ASC_Werling.SV.overlaps.bed.gz
#HGSV (Chaisson) with carrier frequencies
for pop in ALL AFR AMR EAS; do
  ${BIN}/compare_callsets.sh -C \
    -O ${QCTMP}/HGSV_Chaisson.SV.${pop}.overlaps.bed \
    -p HGSV_Chaisson_${pop}_Benchmarking_SV \
    ${QCTMP}/vcf2bed_cleaned.simple.bed.gz \
    ${BENCHDIR}/HGSV_Chaisson.SV.hg19_liftover.${pop}.bed.gz
done
cp ${QCTMP}/HGSV_Chaisson.SV.ALL.overlaps.bed \
   ${OUTDIR}/data/HGSV_Chaisson.SV.overlaps.bed
bgzip -f ${OUTDIR}/data/HGSV_Chaisson.SV.overlaps.bed
tabix -f ${OUTDIR}/data/HGSV_Chaisson.SV.overlaps.bed.gz


###COLLECT LISTS OF VARIANTS PER SAMPLE
#Print status
if [ ${QUIET} == 0 ]; then
  echo -e "$( date ) - VCF QC STATUS: Starting per-sample QC metric generation"
fi
#Gather list of VIDs and genotypes per sample
fgrep -v "##" ${QCTMP}/input.vcf | \
${BIN}/perSample_vcf_parsing_helper.R \
/dev/stdin \
${QCTMP}/analysis_samples.list \
${QCTMP}/perSample/
#Gzip all output lists
for file in ${QCTMP}/perSample/*.VIDs_genotypes.txt; do
  gzip -f ${file}
done
#Create tarball of all sample-variant genotype pairs and move to $OUTDIR
tar -czvf ${OUTDIR}/data/variants_and_genotypes_per_sample.tar.gz \
${QCTMP}/perSample


###GATHER FILES PER FAMILY
#Only evaluate if fam file is supplied
if ! [ -z ${FAM} ]; then
  #Break down families by complete duos and trios
  #Write header
  echo -e "#FAM_ID\tPROBAND\tFATHER\tMOTHER\tPROBAND_SEX\tPHENOTYPE" > \
  ${QCTMP}/input_fam.cleaned.fam
  #Write revised list of duos and trios after masking on analysis_samples.list
  while read fam pro fa mo sex pheno; do
    #Only consider families with child in analysis set
    if [ $( fgrep -w ${pro} ${QCTMP}/analysis_samples.list | wc -l ) -gt 0 ]; then
      #Check for father
      if [ $( fgrep -w ${fa} ${QCTMP}/analysis_samples.list | wc -l ) -eq 0 ]; then
        fa="."
      fi
      
      #Check for mother
      if [ $( fgrep -w ${mo} ${QCTMP}/analysis_samples.list | wc -l ) -eq 0 ]; then
        mo="."
      fi

      #Only consider families with at least one parent
      if [ ${fa} != "." ] || [ ${mo} != "." ]; then
        #Write as duo if one parent is missing
        if [ ${fa} == "." ] || [ ${mo} == "." ]; then
          echo -e "DUO_\t${pro}\t${fa}\t${mo}\t${sex}\t${pheno}"
        #Otherwise, write as trio
        else
          echo -e "TRIO_\t${pro}\t${fa}\t${mo}\t${sex}\t${pheno}"
        fi
      fi
    fi
  done < ${FAM} | sort -Vk1,1 -k2,2V > ${QCTMP}/input_fam.cleaned.fam.tmp

  #Clean output fam file
  for f in DUO TRIO; do
    grep -e "^${f}" ${QCTMP}/input_fam.cleaned.fam.tmp | \
    sort -Vk2,2 | awk -v OFS="\t" '{ print $1NR, $2, $3, $4, $5, $6 }'
  done >> ${QCTMP}/input_fam.cleaned.fam
  rm ${QCTMP}/input_fam.cleaned.fam.tmp
fi


###CLEAN UP
rm -rf ${QCTMP}

