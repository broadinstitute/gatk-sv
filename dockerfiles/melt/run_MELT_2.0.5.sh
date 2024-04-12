#!/usr/bin/env bash

#############################
#    gnomAD SV Discovery    #
#############################

# gnomAD credits: http://gnomad.broadinstitute.org/

# Wrapper to run MELT (Gardner et al., Unpublished)

# Updated by ShuangBroad
# Several differences between this and Ryan L. Collins' version
#   1): read_len AND mean_is are now required, positional arguments
#   2): unified support for HG19 and HG38
#   3): allows specifying JVM max mem

set -euo pipefail

##### Usage statement
usage(){
cat <<EOF
  usage: runMELT.sh bam ref cov read_len mean_is MELT_DIR REF_VER
  Runs MELT for mobile element detection. Requires paired-end deep (>10X) WGS data mapped with BWA.
  Positional arguments (all required):
    bam        Full path to mapped bam file (bam.bai assumed accompanying the BAM)
    ref        Full path to reference FASTA
    cov        Approximate nucleotide coverage of bam file (+/- 2x is ok)
    read_len   Mean read length of library
    mean_is    Mean insert size of library
    MELT_DIR   Full path to MELT install directory
    REF_VER    Reference version (19|38)
EOF
}

##### Args parsing and validation
if [[ "$#" -eq 0 ]]; then
  usage && exit 0;
elif [[ "$#" -lt 7 ]]; then
  echo "At least one of the required parameters is not properly set by the given command:"
  temp_args="$@" && echo "$0 ${temp_args}" && exit 1;
fi

JVM_MAX_MEM=${JVM_MAX_MEM:-12G}
# see http://melt.igs.umaryland.edu/manual.php, for option "-d"
MIN_CHR_LENGTH=${MIN_CHR_LENGTH:-40000000}

bam=$1
ref=$2
cov=$3
read_len=$4
mean_is=$5
MELT_DIR=$6
REF_VER=$7

##### Check for required input (unset or empty)
if [ -z "${bam}" ] || [ -z "${ref}" ] || [ -z "${cov}" ] || [ -z "${read_len}" ] || [ -z "${mean_is}" ] || [ -z "${MELT_DIR}" ] || [ -z "${REF_VER}" ]; then
  echo "At least one of the required parameters is not properly set by the given command:"
  temp_args="$@" && echo "$0 ${temp_args}" && exit 1; # non-zero exit because it indicates user errror
fi

if [[ ! -f "${bam}" ]]; then
  echo "Provided bam file doesn't exist: ${bam}" && exit 1
elif [[ ! -f "${bam}.bai" ]]; then
  echo "Provided bam file doesn't have accompanying index file" && exit 1
elif [[ ! -f "${ref}" ]]; then
  echo "Provided reference file doesn't exist: ${ref}" && exit 1
elif [[ ! -d "${MELT_DIR}" ]]; then
  echo "Provided MELT directory doesn't exist: ${MELT_DIR}" && exit 1
fi

##### Floor coverage value, read length, and insert size
cov=$( echo "${cov}" | cut -f1 -d\. )
read_len=$( echo "${read_len}" | cut -f1 -d\. )
mean_is=$( echo "${mean_is}" | cut -f1 -d\. )

##### remove trailing slash just to make sure
MELT_DIR="${MELT_DIR%/}"

##### Create transposons reference list
if [[ "$REF_VER" == "38" ]]; then
  ls "${MELT_DIR}"/me_refs/Hg38/*zip | sed 's/\*//g' > "transposon_reference.list"
  GENE_BED_FILE="${MELT_DIR}/add_bed_files/Hg38/Hg38.genes.bed"
elif [[ "$REF_VER" == "19" ]]; then
  ls "${MELT_DIR}"/me_refs/1KGP_Hg19/*zip | sed 's/\*//g' > "transposon_reference.list"
  GENE_BED_FILE="${MELT_DIR}/add_bed_files/1KGP_Hg19/hg19.genes.bed"
fi

##### Run MELT Single
java -Xmx"${JVM_MAX_MEM}" \
  -jar "${MELT_DIR}/MELT.jar" \
  Single \
  -bamfile "${bam}" \
  -h "${ref}" \
  -c "${cov}" \
  -r "${read_len}" \
  -e "${mean_is}" \
  -d "${MIN_CHR_LENGTH}" \
  -t "transposon_reference.list" \
  -n "${GENE_BED_FILE}" \
  -w "."
