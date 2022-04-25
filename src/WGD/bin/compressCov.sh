#!/bin/bash

# Wrapper for bedtools map to compress raw binCov output into
# larger bin sizes

set -e 

#Usage statement
usage(){
cat <<EOF

usage: compressCov.sh [-h] [-N] [-s] [-z] [-o OUTFILE] INPUT RATIO

Helper tool to automate compression of raw binCov.py output bed files or
bed-style coverage matrices into larger bin sizes

Positional arguments:
  INPUT     path to binCov.py bed file or bed-stype matrix
  RATIO     compression ratio

Optional arguments:
  -h  HELP        Show this help message and exit
  -N  NO HEADER   Input BED file does not have a header (default: has header)
  -s  SUM         Report sum (default: report median)
  -z  BGZIP       Bgzip and tabix output file
  -o  OUTFILE     Output file (default: stdout)

EOF
}

#Parse arguments
OUTFILE=/dev/stdout
HEADER=1
BEDOP="median"
GZ=0
while getopts ":hNszo:" opt; do
  case "$opt" in
    h)
      usage
      exit 0
      ;;
    N)
      HEADER=0
      ;;
    s)
      BEDOP="sum"
      ;;
    z)
      GZ=1
      ;;
    o)
      OUTFILE=${OPTARG}
      ;;
  esac
done
shift $(( ${OPTIND} - 1))
INPUT=$1
RATIO=$2

#Check for required input
if [ -z ${INPUT} ] || [ -z ${RATIO} ]; then
  usage
  exit 0
fi

#Check for gzip optioned only if output file specified
if [ ${OUTFILE} == "/dev/stdout" ] && [ ${GZ} == 1 ]; then
  echo -e "\nOUTFILE required for zip-compressed output"
  usage
  exit 0
fi

#Scrub ".gz" from output filename if provided by user
if [ ${GZ} == 1 ] && [ ${OUTFILE: -3} == ".gz" ]; then
  OUTFILE=$( echo "${OUTFILE}" | sed 's/\.gz//g' )
fi

#Unzip input file if gzipped
GZI=0
if file ${INPUT} | grep -qE "gzip (compressed|compatible)"; then
  GZI=1
  TMPI=`mktemp`; mv ${TMPI} ${TMPI}.gz; TMPI=${TMPI}.gz
  cp ${INPUT} ${TMPI}
  gunzip ${TMPI}
  INPUT=$( echo "${TMPI}" | sed "s/\.gz/\t/g" | cut -f1 )
fi

#Gather number of samples present in file
NSAMP=$( head -n1 ${INPUT} | awk '{ print NF-3 }' )
COLS=$( seq 4 $((${NSAMP}+3)) | paste -s -d, )

#Calculate new bin size
OBIN=$( fgrep -v "#" ${INPUT} | head -n2 | tail -n1 | awk '{ print $3-$2 }' )
NBIN=$(( ${RATIO}*${OBIN} ))

#Print header from original file to OUTFILE
if [ ${HEADER} -eq 1 ]; then
  head -n1 ${INPUT} > ${OUTFILE}
fi

#Iterate over contigs present in input file
while read CONTIG; do
  #Get min and max coordinate present in input file for contig
  MIN=$( awk -v CONTIG=${CONTIG} '{ if ($1==CONTIG) print $2 }' ${INPUT} | head -n1 )
  MAX=$( awk -v CONTIG=${CONTIG} '{ if ($1==CONTIG) print $3 }' ${INPUT} | tail -n1 )
  #Perform bedtools map function (operation dependent on norm/non-norm)
  bedtools map -c ${COLS} -o ${BEDOP} \
  -a <( paste <( seq ${MIN} ${NBIN} ${MAX} ) \
              <( seq $((${MIN}+${NBIN})) ${NBIN} $((${MAX}+${NBIN})) ) | \
        awk -v CONTIG=${CONTIG} -v OFS="\t" '{ print CONTIG, $1, $2 }' ) \
 -b <( awk -v CONTIG=${CONTIG} -v OFS="\t" '{ if ($1==CONTIG) print $0 }' \
  ${INPUT} ) | awk '{ if ($4!=".") print $0 }'
done < <( fgrep -v "#" ${INPUT} | cut -f1 | sort -Vk1,1 | uniq ) >> ${OUTFILE}

#Gzip output if optioned
if [ ${GZ} == 1 ]; then
  bgzip -f ${OUTFILE}
  tabix -f ${OUTFILE}.gz
fi

#Clean up
# if [ GZI == 1 ]; then
#   rm -rf ${TMPI}
#   rm -rf ${INPUT}
# fi
