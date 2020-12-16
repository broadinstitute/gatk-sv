#!/bin/bash

# Resolve redundancies between simple CNVs and unbalanced complex SV in mod04b

set -eu -o pipefail

###USAGE
usage(){
cat <<EOF

usage: resolve_CPX_CNV_redundancies.sh [-h] INVCF OUTVCF

Resolve redundancies between simple CNVs and unbalanced complex SV in mod04b

Positional arguments:
  INVCF                    Original input VCF prior to regenotyping
  OUTVCF                   Full path to output VCF after relabeling

Optional arguments:
  -h  HELP                 Show this help message and exit
  -t  TEMPBASE              Use TEMPBASE as the base for temporary directory, instead of /tmp
Notes:
  1. All input files must be compressed with bgzip.

EOF
}


###PARSE ARGS
TEMPBASE=${TEMPBASE:-""}
while getopts ":ht:" opt; do
  case "$opt" in
    h)
      usage
      exit 0
      ;;
    t)
      TEMPBASE=${OPTARG}
      echo "TEMPBASE=$TEMPBASE"
      mkdir -p $TEMPBASE
      ;;
    *)
      echo "opt=$opt"
  esac
done
shift $(( ${OPTIND} - 1))
INVCF=$1
OUTVCF=$2

###PROCESS ARGS & OPTIONS
#Check for required input
if [ -z ${INVCF} ]; then
  echo -e "\nERROR: input VCF not specified\n"
  usage
  exit 0
fi
if ! [ -s ${INVCF} ]; then
  echo -e "\nERROR: input VCF either empty or not found\n"
  usage
  exit 0
fi
if [ $( file ${INVCF} | fgrep "gzip" | wc -l ) -lt 1 ]; then
  echo -e "\nERROR: input VCF must be bgzipped\n"
  usage
  exit 0
fi
if [ -z ${OUTVCF} ]; then
  echo -e "\nERROR: path to output VCF not specified\n"
  usage
  exit 0
fi
#Set path to execution directory
BIN=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
#Prepares temporary directory
PROCDIR=`mktemp ${TEMPBASE:+"-p" "$TEMPBASE"} -d`

###PREP FILES
#Convert full VCF to BED intervals
#Ignore CPX events with UNRESOLVED filter status
svtk vcf2bed --split-cpx --info SVTYPE \
  <(bcftools view -e 'INFO/SVTYPE == "CPX" && FILTER == "UNRESOLVED"' ${INVCF}) - \
  | grep -e '^#\|DEL\|DUP\|CNV\|CPX' \
  | awk -v OFS="\t" '{ if ($5=="CNV") print $1, $2, $3, $4, $6, "DEL", $7"\n"$1, $2, $3, $4, $6, "DUP", $7; \
                       else if ($5=="DEL" || $5=="DUP") print $1, $2, $3, $4, $6, $5, $7 }' \
  | sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
  > ${PROCDIR}/intervals.preclustered.bed

awk '$7 == "CPX" {print $4}' ${PROCDIR}/intervals.preclustered.bed | sort -u > ${PROCDIR}/cpx_vids.list
awk '$7 != "CPX" {print $4}' ${PROCDIR}/intervals.preclustered.bed | sort -u > ${PROCDIR}/non_cpx_vids.list
# get rid of extra column and zip
cut -f1-6 ${PROCDIR}/intervals.preclustered.bed | bgzip -c > ${PROCDIR}/intervals.preclustered.bed.gz
# rm ${PROCDIR}/intervals.preclustered.bed

###REMOVE CNVS REDUNDANT WITH COMPLEX EVENTS
#Subset to only variants that share some overlap (at least 10% recip) with at least one CPX variant
bedtools intersect -wa -r -f 0.1 \
  -a ${PROCDIR}/intervals.preclustered.bed.gz \
  -b <( zgrep -wf ${PROCDIR}/cpx_vids.list ${PROCDIR}/intervals.preclustered.bed.gz ) \
  | sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
  | uniq \
  | bgzip -c \
  > ${PROCDIR}/intervals.preclustered.subset.bed.gz
#Melt subsetted variants
while read chr start end VID samples CNV; do
  echo -e "${samples}" \
    | sed 's/,/\n/g' \
    | awk -v OFS="\t" -v chr=${chr} -v start=${start} -v end=${end} -v VID=${VID} -v CNV=${CNV} \
      '{ print chr, start, end, VID, $1, CNV}'
done < <( zcat ${PROCDIR}/intervals.preclustered.subset.bed.gz ) \
  | bgzip -c \
  > ${PROCDIR}/intervals.preclustered.subset.melted.bed.gz
#Cluster BED intervals (50% RO)
svtk bedcluster -f 0.5 \
  ${PROCDIR}/intervals.preclustered.subset.melted.bed.gz - \
  | bgzip -c > \
  ${PROCDIR}/intervals.clustered.bed.gz

#Get list of all variants that cluster with a complex variant,
# evaluate sample overlap from original intervals file, 
# and, if overlap >50%, write that ID to be stripped from the output VCF
while read VIDs; do
  #Get nonredundant list of sample IDs involved in any clustered variant
  num_nonredundant_samples=$(echo -e "${VIDs}" | sed 's/,/\n/g' \
    | fgrep -wf - <( zcat ${PROCDIR}/intervals.preclustered.bed.gz ) \
    | cut -f5 | tr , '\n' | sort | uniq | wc -l)
  #Iterate over VIDs and print non-CPX VID if sample overlap >50%
  while read VID samples; do
    #Get list of samples in variant
    nsamp=$(echo -e "${samples}" | sed 's/,/\n/g' | sort | uniq | wc -l)

    #Compare
    frac=$((nsamp * 100 / num_nonredundant_samples))
    if [ ${frac} -ge 50 ]; then
      echo "${VID}"
    fi

  done < <( echo -e "${VIDs}" | sed 's/,/\n/g'  | fgrep -wf ${PROCDIR}/non_cpx_vids.list \
              | fgrep -wf - <( zcat ${PROCDIR}/intervals.preclustered.bed.gz ) \
              | cut -f4,5 | sort | uniq )

done < <( zcat ${PROCDIR}/intervals.clustered.bed.gz \
            | cut -f7 | fgrep -wf ${PROCDIR}/cpx_vids.list | fgrep -wf ${PROCDIR}/non_cpx_vids.list ) \
  | sort -V | uniq \
  > ${PROCDIR}/VIDs_to_remove.list

###FIND REMAINING REDUNDANT CNVS WITH STRONG (80%) OVERLAP IN SAMPLES AND SIZE
#Find CNV intervals that have 80% reciprocal overlap
zgrep -wFf ${PROCDIR}/non_cpx_vids.list ${PROCDIR}/intervals.preclustered.bed.gz \
  > ${PROCDIR}/intervals.non_cpx.preclustered.bed
bedtools intersect -wa -wb -r -f 0.8 \
  -a ${PROCDIR}/intervals.non_cpx.preclustered.bed \
  -b ${PROCDIR}/intervals.non_cpx.preclustered.bed \
  | awk -v FS="\t" '{ if ($4!=$10 && $6==$12) print $0 }' \
  | bgzip -c \
  > ${PROCDIR}/step2.intervals.preclustered.subset.bed.gz
#Determine which events share 80% sample overlap
while read VIDa sa VIDb sb; do
  na=$( echo -e "${sa}" | sed 's/,/\n/g' | sort | uniq | wc -l )
  nb=$( echo -e "${sb}" | sed 's/,/\n/g' | sort | uniq | wc -l )
  denom=$( echo -e "${sa},${sb}" | sed 's/,/\n/g' | sort | uniq | wc -l )
  numer=$( echo -e "${sa}" | sed 's/,/\n/g' | (fgrep -wf - \
            <( echo -e "${sb}" | sed 's/,/\n/g' ) || true) \
            | sort | uniq | wc -l )
  if [ ${denom} -gt 0 ]; then
    ovr=$(( 100 * ${numer} / ${denom} ))
  fi
  if [ -z ${ovr} ]; then
    ovr=0
  fi
  if [ ${ovr} -ge 80 ]; then
    echo -e "${VIDa}\n${VIDb}" \
      | sort | uniq | paste -s -d,
  fi
done < <( zcat ${PROCDIR}/step2.intervals.preclustered.subset.bed.gz \
            | cut -f4,5,10,11 ) \
  | sort | uniq \
  > ${PROCDIR}/step2.variants_to_resolve.list
#Iterate over variants, pick info & coords from variant with largest N, 
# and consolidate genotypes
sed 's/,/\n/g' ${PROCDIR}/step2.variants_to_resolve.list \
  | sort | uniq \
  > ${PROCDIR}/step2.variants_to_resolve.melted.list
if [ -e ${PROCDIR}/records_to_add.vcf ]; then
  rm ${PROCDIR}/records_to_add.vcf
fi
touch ${PROCDIR}/records_to_add.vcf

function get_best_vcf_line() {
  python <(cat <<'END'
import sys
from typing import List
ref_gts = {"0", "0/0", "0|0"}
ref_or_nocall_gts = {"0", "0/0", "0|0", '.', "./.", "0/.", "./0", ".|.", "0|.", ".|0"}

def _get_cn_score(split_line: List[str], cn_index: int) -> (int, int, str):
  cns = [column.split(':', cn_index + 1)[cn_index] for column in split_line[9:]]
  ref_cn = "1" if split_line[0] in {"chrY"} else "2"
  num_carrier = sum(cn not in {ref_cn, '.'} for cn in cns)
  num_ref = sum(cn == ref_cn for cn in cns)
  return num_carrier, num_ref, split_line[2]

def _get_gt_score(split_line: List[str], gt_index: int) -> (int, int, str):
  gts = [column.split(':', gt_index + 1)[gt_index] for column in split_line[9:]]
  num_carrier = sum(gt not in ref_or_nocall_gts for gt in gts)
  num_ref = sum(gt in ref_gts for gt in gts)
  return num_carrier, num_ref, split_line[2]

def _get_line_score(vcf_line: str) -> (int, int, str):
  split_line = vcf_line.split()
  format_split = split_line[8].split(':')
  if "CN" in format_split:
    return _get_cn_score(split_line, format_split.index("CN"))
  else:
    return _get_gt_score(split_line, format_split.index("GT"))

print(max((vcf_line.strip() for vcf_line in sys.stdin), key=_get_line_score))
END
)
}

until [ $( cat ${PROCDIR}/step2.variants_to_resolve.melted.list | wc -l ) -eq 0 ]; do
  #get next variant
  VID=$( head -n1 ${PROCDIR}/step2.variants_to_resolve.melted.list )
  #get all other variants from clusters containing this variant
  fgrep -w ${VID} ${PROCDIR}/step2.variants_to_resolve.list \
    | sed 's/,/\n/g' | sort | uniq \
    > ${PROCDIR}/step2.partners.tmp
  #Add new record with greatest total number of samples to final append tmp file
  zgrep $INVCF -wFf ${PROCDIR}/step2.partners.tmp | get_best_vcf_line >> ${PROCDIR}/records_to_add.vcf
  #Write list of variants to exclude from original VCF
  cat ${PROCDIR}/step2.partners.tmp >> ${PROCDIR}/VIDs_to_remove.list
  #Exclude variants from list of VIDs to resolve
  fgrep -wvf ${PROCDIR}/step2.partners.tmp \
    ${PROCDIR}/step2.variants_to_resolve.melted.list \
    > ${PROCDIR}/step2.variants_to_resolve.melted.list2 \
    || true
  mv ${PROCDIR}/step2.variants_to_resolve.melted.list2 \
    ${PROCDIR}/step2.variants_to_resolve.melted.list
done


###CLEAN UP FINAL OUTPUT
zcat ${INVCF} \
  | fgrep -wvf ${PROCDIR}/VIDs_to_remove.list \
  | cat - ${PROCDIR}/records_to_add.vcf \
  | vcf-sort \
  | bgzip -c \
  > ${OUTVCF}


###CLEAN UP
rm -rf ${PROCDIR}
