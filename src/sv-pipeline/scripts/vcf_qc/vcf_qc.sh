#!/bin/bash

# Master wrapper to collect QC data & generate plots for SV VCF output by SV pipeline

set -e

###USAGE
usage(){
cat <<EOF

usage: vcf_qc.sh [-h] [-F FAM] [-R] VCF SVTYPES BENCHDIR OUTDIR

Perform several quality control (QC) analyses on an SV VCF from sv-pipeline

Positional arguments:
  VCF             VCF from sv-pipeline
  SVTYPES         List of SV types to evaluate. Two-column, tab-delimited file.
                  First column: sv type. Second column: HEX color for sv type.
  BENCHDIR        Directory containing benchmark archives
  OUTDIR          Output directory for all QC data & plots

Optional arguments:
  -h  HELP        Show this help message and exit
  -F  FAM         FAM file (for Mendelian violation analysis)
  -R  RESTRICT    Restrict analysis to samples present in FAM file
  -q  QUIET       Silence all status updates

EOF
}


###PARSE ARGS
unset FAM
RESTRICT=0
QUIET=0
while getopts ":F:Rqh" opt; do
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


###SET BIN
BIN=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )


###COLLECT QC METRICS
#Print status
if [ ${QUIET} == 0 ]; then
  echo -e "$( date ) - VCF QC STATUS: Beginning QC workflow"
fi
#Set command line options for collectQC.sh
opts=""
if ! [ -z ${FAM} ]; then
  opts="${opts} -F ${FAM}"
fi
if [ ${RESTRICT} -eq 1 ]; then
  opts="${opts} -R "
fi
if [ ${QUIET} == 1 ]; then
  opts="${opts} -q "
fi
#Run collectQC.sh
${BIN}/collectQC.sh ${opts} ${VCF} ${SVTYPES} ${BENCHDIR} ${OUTDIR}


###PLOT VCF SUMMARY DISTRIBUTIONS
#Print status
if [ ${QUIET} == 0 ]; then
  echo -e "$( date ) - VCF QC STATUS: Generating VCF-wide QC plots"
fi
#Generate plots
${BIN}/plot_sv_vcf_distribs.R \
-N ${nsamp} -S ${SVTYPES} \
${QCTMP}/vcf2bed_cleaned.simple.bed.gz \
${OUTDIR}/


###PLOT EXTERNAL BENCHMARKING DATA FOR THREE STANDARD DATASETS
#Print status
if [ ${QUIET} == 0 ]; then
  echo -e "$( date ) - VCF QC STATUS: Plotting external benchmarking"
fi
#1000 Genomes Project phase 3 samples (Sudmant et al.) with allele frequencies
${BIN}/plot_callset_comparison.R \
-p 1000G_Sudmant \
${OUTDIR}/data/1000G_Sudmant.SV.overlaps.bed.gz \
${OUTDIR}/supporting_plots/1000G_Sudmant_comparisons/
cp ${OUTDIR}/supporting_plots/1000G_Sudmant_comparisons/main_plots/VCF_QC.1000G_Sudmant.callset_benchmarking.png \
${OUTDIR}/main_plots/
for pop in ALL AFR AMR EAS EUR SAS; do
  ${BIN}/plot_callset_comparison.R \
  -p 1000G_${pop} \
  ${QCTMP}/1000G_Sudmant.SV.${pop}.overlaps.bed \
  ${OUTDIR}/supporting_plots/1000G_Sudmant_comparisons/${pop}_samples/
done
#Autism Sequencing Consortium pilot & phase 1 families (Werling et al.) with carrier frequencies
${BIN}/plot_callset_comparison.R -C \
-p ASC_Werling \
${OUTDIR}/data/ASC_Werling.SV.overlaps.bed.gz \
${OUTDIR}/supporting_plots/ASC_Werling_comparisons/
cp ${OUTDIR}/supporting_plots/ASC_Werling_comparisons/main_plots/VCF_QC.ASC_Werling.callset_benchmarking.png \
${OUTDIR}/main_plots/
for pop in ALL EUR OTH; do
  ${BIN}/plot_callset_comparison.R -C \
  -p ASC_${pop} \
  ${QCTMP}/ASC_Werling.SV.${pop}.overlaps.bed \
  ${OUTDIR}/supporting_plots/ASC_Werling_comparisons/${pop}_samples/
done
#Human Genome Structural Variation Consortium trios (Chaisson et al.) with carrier frequencies
${BIN}/plot_callset_comparison.R -C \
-p HGSV_Chaisson \
${OUTDIR}/data/HGSV_Chaisson.SV.overlaps.bed.gz \
${OUTDIR}/supporting_plots/HGSV_Chaisson_comparisons/
cp ${OUTDIR}/supporting_plots/HGSV_Chaisson_comparisons/main_plots/VCF_QC.HGSV_Chaisson.callset_benchmarking.png \
${OUTDIR}/main_plots/
for pop in ALL AFR AMR EAS; do
  ${BIN}/plot_callset_comparison.R -C \
  -p HGSV_${pop} \
  ${QCTMP}/HGSV_Chaisson.SV.${pop}.overlaps.bed \
  ${OUTDIR}/supporting_plots/HGSV_Chaisson_comparisons/${pop}_samples/
done


###PLOT PER-SAMPLE DISTRIBUTIONS
${BIN}/plot_sv_perSample_distribs.R \
-S ${SVTYPES} \
${QCTMP}/vcf2bed_cleaned.simple.bed.gz \
${QCTMP}/analysis_samples.list \
${QCTMP}/perSample/ \
${OUTDIR}/


###PLOT INHERITANCE ANALYSES
#Only run if any complete duos or trios exist in the dataset
if [ $( fgrep -v "#" ${QCTMP}/input_fam.cleaned.fam | wc -l ) -gt 0 ]; then
  ${BIN}/analyze_fams.R \
  -S ${SVTYPES} \
  ${QCTMP}/vcf2bed_cleaned.simple.bed.gz \
  ${QCTMP}/input_fam.cleaned.fam \
  ${QCTMP}/perSample/ \
  ${OUTDIR}/
fi








