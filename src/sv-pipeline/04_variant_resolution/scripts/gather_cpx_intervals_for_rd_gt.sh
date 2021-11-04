#!/bin/bash

# Gather relevant complex SV intervals for testing

set -e

###USAGE
usage(){
cat <<EOF

usage: gather_cpx_intervals_for_rd_gt.sh [-h] [-D <integer>] [-I integer] INVCF OUT

Reassign variant labels based on depth regenotyping in mod04b

Positional arguments:
  INVCF                    Input VCF prior to regenotyping
  OUT                      Output BED file of intervals to be genotyped

Optional arguments:
  -h  HELP                 Show this help message and exit
  -D  <integer>            Minimum insertion site size (in bp) to be considered for 
                           distinguishing insertion site deletions [default: 150 bp]
  -I  <integer>            Minimum size for considering palindromic inverted 
                           duplication [default: 1000 bp]

Notes:
  1. All input files must be compressed with bgzip.

EOF
}


###PARSE ARGS
MINSINK=150
MININV=1000
while getopts ":D:I:h" opt; do
  case "$opt" in
    h)
      usage
      exit 0
      ;;
    D)
      MINSINK=${OPTARG}
      ;;
    I)
      MININV=${OPTARG}
      ;;
  esac
done
shift $(( ${OPTIND} - 1))
INVCF=$1
OUT=$2


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
if [ -z ${OUT} ]; then
  echo -e "\nERROR: path to output BED not specified\n"
  usage
  exit 0
fi
#Prepares temporary directory
SCDIR=`mktemp -d`

# #Old implementation
# /opt/sv-pipeline/04_variant_resolution/scripts/CPX_CNV_Parser.py \
# ${INVCF} ${OUT}

#New implementation
#Convert to bed
zcat ${INVCF} | grep -e '^#\|CPX' | fgrep -v UNRESOLVED | bgzip -c > ${SCDIR}/cleaned.in.vcf.gz
svtk vcf2bed --split-cpx \
-i SVTYPE -i CPX_TYPE -i CPX_INTERVALS -i SOURCE \
${SCDIR}/cleaned.in.vcf.gz ${SCDIR}/intervals.bed

#Convert unresolved inversion single enders to bed
zcat ${INVCF} | grep -e '^#\|INVERSION_SINGLE_ENDER' \
  | grep -e '^#\|UNRESOLVED' | bgzip -c > ${SCDIR}/cleaned.inv_se.in.vcf.gz
svtk vcf2bed --split-cpx \
  -i SVTYPE -i UNRESOLVED_TYPE -i CPX_INTERVALS -i SOURCE -i END \
  ${SCDIR}/cleaned.inv_se.in.vcf.gz ${SCDIR}/cleand.inv_se.in.bed
awk -F "\t" -v OFS="\t" '{ print $1, $2, $NF, $4, $5, $6, $7, $8, $9, $10 }' \
  ${SCDIR}/cleand.inv_se.in.bed \
  > ${SCDIR}/inversion_singleender_intervals.bed

#Write header
echo -e "#chr\tstart\tend\tVID\tsamples\tCNV" \
  > ${SCDIR}/output.bed

#Wrap all lines of output in same output for sorting & uniq
for wrapper in 1; do
  #Gather basic intervals to test
  fgrep -v "#" ${SCDIR}/intervals.bed \
    | awk -F "\t" -v OFS="\t" '$8 !~ /INS:ME:|MEI_INS|CTX_P|CTX_Q/ { if ($8!="INV" && $8!="DEL" && $8!="DUP") \
    print $1, $2, $3, $4";"$8";"$7";"$9, $6, "DUP" }'

  #Supplement with explicit CPX_INTERVALS & SOURCE intervals
  while read chr start end VID alt samples svtype cpxtype cpxint sourceint; do
    #Wrap CPX_INTERVALS and SOURCE intervals in same stream
    for wrapperB in 1; do
      #Split complex intervals
      if [ ${cpxint} != "NA" ]; then
        echo "${cpxint}" \
          | sed -e 's/\_/\t/g' -e 's/\:/\t/g' -e 's/\-/\t/g' -e 's/,/\n/g' \
          | awk -F "\t" -v OFS="\t" '{ print $2, $3-1, $4, $1 }'
      fi
      #Split insertion intervals
      if [ ${sourceint} != "NA" ]; then
        echo "${sourceint}" \
          | sed -e 's/\_/\t/g' -e 's/\:/\t/g' -e 's/\-/\t/g' -e 's/,/\n/g' \
          | awk -F "\t" -v OFS="\t" '{ print $2, $3-1, $4, "INS" }'
      fi
    done \
      | sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
      | awk -F "\t" -v OFS="\t" -v VID=${VID} -v svtype=${svtype} -v cpxtype=${cpxtype} \
        -v cpxint=${cpxint} \
        '{ print $1, $2, $3, VID";"cpxtype";"svtype";"cpxint }' \
      | sed "s/$/\t$samples\tDUP/g" \
      | uniq
  done < <( fgrep -v "#" ${SCDIR}/intervals.bed \
            | awk -F "\t" -v FS="\t" '{ if ($7!="INV" && ($9!="NA" || $10!="NA")) print $0 }' )

  #Add inversion single ender intervals if interval >= $MININV 
  fgrep -v "#" ${SCDIR}/inversion_singleender_intervals.bed \
    | awk -F "\t" -v OFS="\t" -v MININV=${MININV} \
      '{ if ($3-$2>=MININV) print $1, $2, $3, $4";"$8";"$7";"$9, $6, "DUP" }'

  #Add insertion sinks if >= $MINSINK
  awk -F "\t" -v FS="\t" -v OFS="\t" -v MINSINK=${MINSINK} \
    '{ if ($7=="INS" && $3-$2>=MINSINK) print $1, $2, $3, $4";"$8";"$7";"$9, $6, "DUP" }' \
    ${SCDIR}/intervals.bed

done \
  | sort -Vk1,1 -k2,2n -k3,3n -k4,4V -k5,5V \
  | uniq \
  >> ${SCDIR}/output.bed

#Compress output & localize
bgzip -f ${SCDIR}/output.bed
mv ${SCDIR}/output.bed.gz ${OUT}

#Clean up
rm -rf ${SCDIR}

