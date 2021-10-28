#!/bin/bash

# Batch script to iterate binCov over a supplied list of contigs

set -e

#Usage statement
usage(){
cat <<EOF

usage: binCov_batch.sh [-h] [-b BINSIZE] [-m MODE] [-n] [-z] [-C] [-I INDEX]
                       [-L CONTIGS] [-x BLACKLIST] [-v OVERLAP] 
                       BAM ID OUTDIR

Wrapper for serialized execution of binCov.py across multiple chromosomes

Positional arguments:
  BAM     Input bam
  ID      Sample ID
  OUTDIR  Output directory

Optional arguments:
  -h  HELP         Show this help message and exit
  -b  BINSIZE      Bin size in bp (default: 1000)
  -m  MODE         Evaluate physical or nucleotide coverage (default: nucleotide)
  -n  NORMALIZED   Also generate normalized coverage values
  -z  GZIP         Attempt to tar & gzip the output directory
  -C  CRAM         Input file is in CRAM format
  -I  INDEX        Path to BAM/CRAM index
  -L  CONTIGS      List of contigs to evaluate (default: all contigs in bam header)
  -x  BLACKLIST    BED file of regions to ignore
  -v  OVERLAP      Maximum tolerated blacklist overlap before excluding bin

EOF
}

#Parse arguments
binsize=1000
mode=nucleotide
contigs=DEFAULT
blist=NONE
v=0.05
norm=0
tgz=0
CRAM=0
index_path=0
while getopts ":b:m:nzCI:L:x:v:h" opt; do
  case "$opt" in
    b)
      binsize=${OPTARG}
      ;;
    m)
      mode=${OPTARG}
      ;;
    n)
      norm=1
      ;;
    z)
      tgz=1
      ;;
    C)
      CRAM=1
      ;;
    I)
      index_path=${OPTARG}
      ;;
    L)
      contigs=${OPTARG}
      ;;
    x)
      blist=${OPTARG}
      ;;
    v)
      v=${OPTARG}
      ;;
    h)
      usage
      exit 0
      ;;
  esac
done
shift $(( ${OPTIND} - 1))
bam=$1
ID=$2
OUTDIR=$3

#Check positional arguments
if [ -z ${bam} ] || [ -z ${ID} ] || [ -z ${OUTDIR} ]; then
  usage
  exit 0
fi

#Create output directory if it doesn't exist
if ! [ -e ${OUTDIR} ]; then
  mkdir ${OUTDIR}
fi

#Determine list of contigs to use (note: requires samtools)
if [ ${contigs} == "DEFAULT" ]; then
  contigs_list=mktemp
  samtools view -H ${bam} | fgrep -w "@SQ" | \
  awk '{ print $2 }' | cut -d\: -f2 > \
  ${contigs_list}
else
  contigs_list=${contigs}
fi

#Run binCov.py on all contigs
spath=$( dirname $( readlink -f $0 ) )
while read contig; do
  #Concatenate command line options
  binCovOptions=$( echo -e "-z -b ${binsize} -m ${mode} -v ${v} ${bam} ${contig} \
                            ${OUTDIR}/${ID}.${contig}.${binsize}bpBins.${mode}.rawCov.bed" )
  if [ ${norm} -eq 1 ]; then
    binCovOptions=$( echo -e "-n ${OUTDIR}/${ID}.${contig}.${binsize}bpBins.${mode}.normCov.bed ${binCovOptions}" )
  fi
  if ! [ ${blist} == "NONE" ]; then
    binCovOptions=$( echo -e "-x ${blist} ${binCovOptions}" )
  fi
  if [ ${CRAM} == 1 ]; then
    binCovOptions=$( echo -e "-C ${binCovOptions}" )
  fi
  if [ ${index_path} != 0 ]; then
    binCovOptions=$( echo -e "-I ${index_path} ${binCovOptions}" )
  fi
  #Run binCov
  ${spath}/binCov.py ${binCovOptions}
done < ${contigs_list}

#Tar & gzip output if optioned
if [ ${tgz} -eq 1 ]; then
  OUTDIR_name=$( echo ${OUTDIR} | sed 's/\/$//g' | sed "s/\//\t/g" | awk '{ print $NF }' )
  cd ${OUTDIR}/../; tar -czvf ${OUTDIR_name}.tar.gz ${OUTDIR_name}
fi
