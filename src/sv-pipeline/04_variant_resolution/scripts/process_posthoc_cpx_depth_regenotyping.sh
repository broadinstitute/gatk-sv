#!/bin/bash

# Reassign variant labels based on depth regenotyping in mod04b

set -eo pipefail
###USAGE
usage(){
cat <<EOF

usage: process_posthoc_cpx_depth_regenotyping.sh [-h] [-s MINSIZE] [-d MINDIFF] [-D integer] [-T integer]
                                                 INVCF INTERVALS GENOTYPES OUTVCF

Reassign variant labels based on depth regenotyping in mod04b

Positional arguments:
  INVCF                    Original input VCF prior to regenotyping
  INTERVALS                BED file of genotyped intervals
  GENOTYPES                Melted depth genotypes
  FAMFILE                  .fam file for all samples to be considered
  OUTVCF                   Full path to output VCF after relabeling

Optional arguments:
  -h  HELP                 Show this help message and exit
  -s  MINSIZE              Minimum size (in bp) of CNV interval to be considered
                           during variant reclassification [default: 1000]
  -d  MINDIFF              Minimum difference in non-ref CNV genotype frequency
                           between predicted carriers and noncarriers to consider
                           a CNV interval to have adequate depth support
                           [default: 0.4]
  -D  <integer>            Minimum insertion site size (in bp) to be considered for
                           distinguishing insertion site deletions [default: 150 bp]
  -T  <integer>            Minimum size (in bp) at which to prioritize an inverted
                           dDUP classification over a dupINV or INVdup classification
                           [default: 1000000 bp]
  -R  <path>               Path to table containing the final reclassification
                           decision made per variant. [default: no table output]
  -G  <path>               Path to table containing the raw genotype counts table
                           per interval per variant. [default: no table output]


Notes:
  1. All input files other than FAMFILE must be compressed with bgzip.
  2. OUTVCF will be automatically bgzipped.
  3. Rationale for minimum insertion site size: reference alignment bias
     can make insertion sites appear deleted, so -D should be set to a size
     at which a confident depth assessment can be made on insertion site. Should
     be at least as large as the library fragment size.

EOF
}


###PARSE ARGS
MINSIZE=1000
MINDIFF=0.4
MINSIZEiDEL=150
MINdDUPTHRESH=1000000
unset RTABLE
unset GTCOUNTSTABLE
while getopts ":s:d:D:T:R:G:h" opt; do
  case "$opt" in
    h)
      usage
      exit 0
      ;;
    s)
      MINSIZE=${OPTARG}
      ;;
    d)
      MINDIFF=${OPTARG}
      ;;
    D)
      MINSIZEiDEL=${OPTARG}
      ;;
    T)
      MINdDUPTHRESH=${OPTARG}
      ;;
    R)
      RTABLE=${OPTARG}
      ;;
    G)
      GTCOUNTSTABLE=${OPTARG}
      ;;
  esac
done
shift $(( ${OPTIND} - 1))
INVCF=$1
INTERVALS=$2
GENOTYPES=$3
FAMFILE=$4
OUTVCF=$5


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
if [ -z ${INTERVALS} ]; then
  echo -e "\nERROR: intervals file not specified\n"
  usage
  exit 0
fi
if ! [ -s ${INTERVALS} ]; then
  echo -e "\nERROR: intervals file either empty or not found\n"
  usage
  exit 0
fi
if [ $( file ${INTERVALS} | fgrep "gzip" | wc -l ) -lt 1 ]; then
  echo -e "\nERROR: intervals file must be bgzipped\n"
  usage
  exit 0
fi
if [ -z ${GENOTYPES} ]; then
  echo -e "\nERROR: melted genotypes file not specified\n"
  usage
  exit 0
fi
if ! [ -s ${GENOTYPES} ]; then
  echo -e "\nERROR: melted genotypes file either empty or not found\n"
  usage
  exit 0
fi
if [ $( file ${GENOTYPES} | fgrep "gzip" | wc -l ) -lt 1 ]; then
  echo -e "\nERROR: melted genotypes file must be bgzipped\n"
  usage
  exit 0
fi
if [ -z ${FAMFILE} ]; then
  echo -e "\nERROR: input .fam file not specified\n"
  usage
  exit 0
fi
if ! [ -s ${FAMFILE} ]; then
  echo -e "\nERROR: input .fam file either empty or not found\n"
  usage
  exit 0
fi
if [ -z ${OUTVCF} ]; then
  echo -e "\nERROR: path to output VCF not specified\n"
  usage
  exit 0
fi
if [ -z ${MINSIZE} ]; then
  echo -e "\nERROR: minimum CNV size not set\n"
  usage
  exit 0
fi
if [ -z ${MINDIFF} ]; then
  echo -e "\nERROR: minimum CNV carrier freq difference not set\n"
  usage
  exit 0
fi
#Set path to execution directory
BIN=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
#Prepares temporary directory
GTDIR=`mktemp -d`


###PREPARE SAMPLE lists
#Master lists of samples extracted from famfile
if [ $( file ${FAMFILE} | fgrep "gzip" | wc -l ) -gt 0 ]; then
  #All samples
  zcat ${FAMFILE} | fgrep -v "#" | cut -f2 | sort | uniq \
    > ${GTDIR}/all.samples.list
  #Males
  zcat ${FAMFILE} | fgrep -v "#" \
    | awk '{ if ($5=="1") print $2 }' | sort | uniq \
    > ${GTDIR}/male.samples.list
  #Females
  zcat ${FAMFILE} | fgrep -v "#" \
    | awk '{ if ($5=="2") print $2 }' | sort | uniq \
    > ${GTDIR}/female.samples.list
else
  #All samples
  fgrep -v "#" ${FAMFILE} | cut -f2 | sort | uniq \
    > ${GTDIR}/all.samples.list
  #Males
  fgrep -v "#" ${FAMFILE} \
    | awk '{ if ($5=="1") print $2 }' | sort | uniq \
    > ${GTDIR}/male.samples.list
  #Females
  fgrep -v "#" ${FAMFILE} \
    | awk '{ if ($5=="2") print $2 }' | sort | uniq \
    > ${GTDIR}/female.samples.list
fi


###GATHER DATA PER VARIANT
#Header
echo -e "#chr\tstart\tend\tVID\tcarrier_del\tcarrier_wt\tcarrier_dup\
\tcontrol_del\tcontrol_wt\tcontrol_dup\
\tdiff_case_control_del_frac\tdiff_case_control_dup_frac" \
  > ${GTDIR}/genotype_counts_per_variant.bed
#Count copy states of carriers & noncarriers per variant
while read chr start end VID samps trash; do
  #Wrap one line per variant
  for wrapper in 1; do
    #Print variant info
    echo -e "${chr}\t${start}\t${end}\t${VID}"

    #Clear existing files
    if [ -e ${GTDIR}/carrier_samples.tmp ]; then
      rm ${GTDIR}/carrier_samples.tmp
    fi
    if [ -e ${GTDIR}/control_samples.tmp ]; then
      rm ${GTDIR}/control_samples.tmp
    fi
    unset medCN

    ###Get list of samples & reference CN to consider, dependent on chr of call
    #ChrX: use only diploid females if possible, otherwise use haploid males
    if [ ${chr} == "X" ] || [ ${chr} == "chrX" ]; then

      #Try to get female carrier samples
      echo -e "${samps}" | sed 's/,/\n/g' \
        | fgrep -wf - ${GTDIR}/female.samples.list \
        > ${GTDIR}/carrier_samples.tmp ||true

      #Use only females if there are any female carriers
      if [ $( cat ${GTDIR}/carrier_samples.tmp | wc -l ) -gt 0 ]; then

        #Get non-carrier female controls
        fgrep -wvf ${GTDIR}/carrier_samples.tmp \
          ${GTDIR}/female.samples.list \
          > ${GTDIR}/control_samples.tmp ||true

        #Set default medCN
        medCNdefault=2

      #Otherwise, use males
      else

        #Get male carrier samples
        echo -e "${samps}" | sed 's/,/\n/g' \
          | fgrep -wf - ${GTDIR}/male.samples.list \
          > ${GTDIR}/carrier_samples.tmp ||true

        #Get non-carrier male controls
        fgrep -wvf ${GTDIR}/carrier_samples.tmp \
          ${GTDIR}/male.samples.list \
          > ${GTDIR}/control_samples.tmp ||true

        #Set default medCN
        medCNdefault=1
      fi

    #ChrY: use only haploid males
    else
      if [ ${chr} == "Y" ] || [ ${chr} == "chrY" ]; then

        #Get male carrier samples
        echo -e "${samps}" | sed 's/,/\n/g' \
          | fgrep -wf - ${GTDIR}/male.samples.list \
          > ${GTDIR}/carrier_samples.tmp ||true

        #Get non-carrier male controls
        fgrep -wvf ${GTDIR}/carrier_samples.tmp \
          ${GTDIR}/male.samples.list \
          > ${GTDIR}/control_samples.tmp ||true

        #Set default medCN
        medCNdefault=1

      #All other chromosomes: use all samples & assume diploid locus
      else

        #Get all carrier samples
        echo "${samps}" | sed 's/,/\n/g' \
          | sort | uniq > ${GTDIR}/carrier_samples.tmp

        #Get all non-carrier controls
        fgrep -wvf ${GTDIR}/carrier_samples.tmp \
          ${GTDIR}/all.samples.list \
          > ${GTDIR}/control_samples.tmp ||true

        #Set default medCN
        medCNdefault=2
      fi
    fi


    #Gather genotypes corresponding to interval of interest
    if [ $(( ${end} - ${start} )) -lt 1000000 ]; then
      bedtools intersect -wa -r -f 0.95 \
        -a ${GENOTYPES} \
        -b <( echo -e "${chr}\t${start}\t${end}" ) \
        | awk -v OFS="\t" -v VID=${VID} '{ if ($4==VID) print $0 }' \
        > ${GTDIR}/intersected_genotypes.tmp.bed
    else
      zcat ${GENOTYPES} \
        | awk -v OFS="\t" -v VID=${VID} '{ if ($4==VID) print $0 }' \
        | bedtools coverage -f 1 \
          -a - \
          -b <( echo -e "${chr}\t${start}\t${end}" ) \
        | cut -f1-7 \
        > ${GTDIR}/intersected_genotypes.tmp.bed
    fi

    #Get median copy number across predicted non-carrier controls
    if [ $( cat ${GTDIR}/control_samples.tmp | wc -l ) -gt 0 ]; then
      medCN=$( fgrep -wf ${GTDIR}/control_samples.tmp \
                 ${GTDIR}/intersected_genotypes.tmp.bed \
                 | cut -f6 \
                 | sort -nk1,1 \
                 | perl -e '$d=.5;@l=<>;print $l[int($d*$#l)]' ||true)
    fi
    #If no non-carriers exist, assume median CN = default
    if [ -z ${medCN} ]; then
      medCN=${medCNdefault}
    fi


    #For predicted carriers, count number of genotypes lower than,
    # equal to, and greater than the overall median
    if [ $( cat ${GTDIR}/carrier_samples.tmp | wc -l ) -gt 0 ]; then
      fgrep -wf ${GTDIR}/carrier_samples.tmp \
        ${GTDIR}/intersected_genotypes.tmp.bed \
        | awk -v OFS="\t" -v medCN=${medCN} '{ if ($6<medCN) print $6 }' | wc -l ||true
      fgrep -wf ${GTDIR}/carrier_samples.tmp \
        ${GTDIR}/intersected_genotypes.tmp.bed \
        | awk -v OFS="\t" -v medCN=${medCN} '{ if ($6==medCN) print $6 }' | wc -l ||true
      fgrep -wf ${GTDIR}/carrier_samples.tmp \
        ${GTDIR}/intersected_genotypes.tmp.bed \
        | awk -v OFS="\t" -v medCN=${medCN} '{ if ($6>medCN) print $6 }' | wc -l ||true
    else
      echo -e "0\n0\n0"
    fi

    #For predicted non-carriers, count number of genotypes lower than,
    # equal to, and greater than the overall median
    if [ $( cat ${GTDIR}/control_samples.tmp | wc -l ) -gt 0 ]; then
      fgrep -wf ${GTDIR}/control_samples.tmp \
        ${GTDIR}/intersected_genotypes.tmp.bed \
        | awk -v OFS="\t" -v medCN=${medCN} '{ if ($6<medCN) print $6 }' | wc -l ||true
      fgrep -wf ${GTDIR}/control_samples.tmp \
        ${GTDIR}/intersected_genotypes.tmp.bed \
        | awk -v OFS="\t" -v medCN=${medCN} '{ if ($6==medCN) print $6 }' | wc -l ||true
      fgrep -wf ${GTDIR}/control_samples.tmp \
        ${GTDIR}/intersected_genotypes.tmp.bed \
        | awk -v OFS="\t" -v medCN=${medCN} '{ if ($6>medCN) print $6 }' | wc -l ||true
    else
      echo -e "0\n0\n0"
    fi

    #Clean up tmp file
    rm ${GTDIR}/intersected_genotypes.tmp.bed
  done \
    | paste -s \
    | awk -v OFS="\t" \
      '{ ncase=$5+$6+$7; nctrl=$8+$9+$10; if (nctrl==0) nctrl=1; \
         if (ncase>0) print $0, ($5/ncase)-($8/nctrl), ($7/ncase)-($10/nctrl); else print $0, 0, 0 }'
done < <( zcat ${INTERVALS} | fgrep -v "#" ) \
  >> ${GTDIR}/genotype_counts_per_variant.bed


###CLEAN UP INTERVALS
#Write header
echo -e "#chr\tstart\tend\tinterval_size\tinterval_type\tVID\tvariant_svtype\tvariant_cpxtype\
\tcarrier_del\tcarrier_wt\tcarrier_dup\
\tcontrol_del\tcontrol_wt\tcontrol_dup\
\tdiff_case_control_del_frac\tdiff_case_control_dup_frac\
\tCNV_assessment" \
  > ${GTDIR}/genotype_counts_per_variant.cleaned_intervals.bed
#Match intervals per variant with CNV evidence
while read mergedVID; do
  VID=$( echo "${mergedVID}" | sed 's/;/\t/g' | cut -f1 )
  cpxtype=$( echo "${mergedVID}" | sed 's/;/\t/g' | cut -f2 )
  svtype=$( echo "${mergedVID}" | sed 's/;/\t/g' | cut -f3 | head -n1 )
  intervals=$( echo "${mergedVID}" | sed 's/;/\t/g' | cut -f4- | paste -s -d\; )
  if [ ${intervals} == "NA" ]; then
    intervals=$( zcat ${INTERVALS} \
                   | fgrep -w "${mergedVID}" \
                   | awk -v OFS="\t" '{ print "UNK_"$1":"$2"-"$3 }' \
                   | paste -s -d, ||true)
  fi
    echo "${intervals}" | sed -e 's/\_/\t/g' -e 's/\:/\t/g' -e 's/\-/\t/g' -e 's/,/\n/g' \
      | awk -v OFS="\t" -v VID=${VID} -v cpxtype=${cpxtype} -v svtype=${svtype} \
        '{ print $2, $3, $4, $4-$3, $1, VID, svtype, cpxtype }' \
      | bedtools intersect -wb -r -f 0.95 \
          -a - \
          -b ${GTDIR}/genotype_counts_per_variant.bed \
      | cut --complement -f9-12 \
      | awk -v OFS="\t" -v MINSIZE=${MINSIZE} -v MINDIFF=${MINDIFF} \
        '{ if ($(NF-1)>MINDIFF && $NF>MINDIFF) assess="DELDUP"; \
           else if ($(NF-1)>MINDIFF && $NF<=MINDIFF) assess="DEL"; \
           else if ($(NF-1)<=MINDIFF && $NF>MINDIFF) assess="DUP"; \
           else if ($3-$2<MINSIZE) assess="TOO_SMALL"; \
           else assess="WT"; \
           print $0, assess }'
done < <( fgrep -v "#" ${GTDIR}/genotype_counts_per_variant.bed \
            | cut -f4 | sort | uniq ||true) \
  | sort -Vk1,1 -k2,2n -k3,3n -k5,5V -k4,4V \
  >> ${GTDIR}/genotype_counts_per_variant.cleaned_intervals.bed


###CONVERT RELEVANT VARIANTS FROM VCF TO BED
#Start with all but INV single enders
cat <( zcat ${INVCF} | fgrep "#" ) \
    <( zcat ${INVCF} | fgrep -v "#" \
         | fgrep -wf <( fgrep -v "#" \
                          ${GTDIR}/genotype_counts_per_variant.cleaned_intervals.bed \
                          | { fgrep -v "INVERSION_SINGLE_ENDER" || true; } \
                          | cut -f6 | sort | uniq ||true) ) \
  | svtk vcf2bed --no-samples --info ALL \
      - ${GTDIR}/variants_to_reclassify.vcf2bed.bed
#Add inversion single enders and use modified coordinates
cat <( zcat ${INVCF} | fgrep "#" ) \
    <( zcat ${INVCF} | fgrep -v "#" \
         | fgrep -wf <( fgrep -v "#" \
                          ${GTDIR}/genotype_counts_per_variant.cleaned_intervals.bed \
                          | { fgrep "INVERSION_SINGLE_ENDER" || true; } \
                          | cut -f6 | sort | uniq ||true) ) \
  | svtk vcf2bed --no-samples --info ALL \
      - ${GTDIR}/inv_se_vcf2bed.precut.bed
ENDidx=$( head -n1 ${GTDIR}/inv_se_vcf2bed.precut.bed \
            | sed 's/\t/\n/g' \
            | awk '{ if ($1=="END") print NR }' )
awk -v ENDidx=${ENDidx} -v OFS="\t" '{ $3=$ENDidx; print }' \
  ${GTDIR}/inv_se_vcf2bed.precut.bed \
  | fgrep -v "#" || true \
  >> ${GTDIR}/variants_to_reclassify.vcf2bed.bed


###MAKE FINAL ASSESSMENT FOR EACH VARIANT
#Print header
echo -e "#VID\tMODIFICATION\tREASON\tNEW_SVTYPE\tNEW_CPX_TYPE\tNEW_CPX_INTERVALS\tNEW_SVLEN\tNEW_SOURCE\tNEW_START\tNEW_END" \
  > ${GTDIR}/final_variant_reassessment_table.txt
#Get variables required for processing
SOURCEidx=$( head -n1 ${GTDIR}/variants_to_reclassify.vcf2bed.bed \
               | sed 's/\t/\n/g' \
               | awk '{ if ($1=="SOURCE") print NR }' )
#Iterate and evaluate each SV
while read VID; do
  #Get base variant class
  svtype=$( fgrep -w ${VID} ${GTDIR}/genotype_counts_per_variant.cleaned_intervals.bed \
              | cut -f7 | sort | uniq -c | sort -nrk1,1 | head -n1 | awk '{ print $2 }' ||true)

  #Get predicted complex type
  cpxtype=$( fgrep -w ${VID} ${GTDIR}/genotype_counts_per_variant.cleaned_intervals.bed \
               | cut -f8 | sort | uniq -c | sort -nrk1,1 | head -n1 | awk '{ print $2 }' ||true)

  #Logic varies based on variant class
  case ${svtype} in
    #Evaluate existing complex variants
    "CPX")
      #Check if the CPX type has any DEL or DUP intervals predicted
      if [ $( fgrep -w ${VID} ${GTDIR}/genotype_counts_per_variant.cleaned_intervals.bed \
                | cut -f5 | grep -e 'DEL\|DUP' | wc -l ||true) -gt 0 ]; then
        #Check if any CNV interval larger than min size failed to validate
        nfail=$( awk -v VID=${VID} '{ if ($6==VID && ($5=="DEL" || $5=="DUP") && $NF=="WT" ) print $0 }' \
                   ${GTDIR}/genotype_counts_per_variant.cleaned_intervals.bed | wc -l )
        if [ ${nfail} -gt 0 ]; then
          MOD="UNRESOLVED"
          REASON="AT_LEAST_ONE_CNV_INTERVAL_FAILED"

        #Otherwise, check if at least one CNV interval validated
        else
          nval=$( awk -v VID=${VID} '{ if ($6==VID && ($5=="DEL" || $5=="DUP") && ($5==$NF || $NF=="DELDUP") ) print $0 }' \
                    ${GTDIR}/genotype_counts_per_variant.cleaned_intervals.bed | wc -l )
          if [ ${nval} -gt 0 ]; then
            MOD="KEEP"
            REASON="AT_LEAST_ONE_CNV_VALIDATED"

          #Otherwise, not enough information to make a call, so keep it
          else
            MOD="KEEP"
            REASON="NO_CNVS_LARGER_THAN_MIN_SIZE"
          fi
        fi
      else
        MOD="KEEP"
        REASON="NO_PREDICTED_CNV_INTERVALS"
      fi
      ;;

    #Examine CNV evidence for insertion variants
    "INS")
      case ${cpxtype} in
        ###Solve translocational insertions
        #A2B: A interval inserted into B
        #B2A: B interval inserted into A
        #CTX_INS_A2B: A interval inserted into B, interchromosomal, can be inverted
        #CTX_INS_B2A: B interval inserted into A, interchromosomal, can be inverted
        #OTHERWISE: Just test for INS-iDEL
        #Note: CTX_INS_B2A variants are not resolved here, because they are indexed
        # onto a different chromosome, and will be dealt with in that shard
        INS_A2B|INS_B2A|CTX_*INS_A2B|CTX_*INS_B2A)
          #Get inserted sequence coordinates
          sourcecoord=$( fgrep -w "${VID}" ${GTDIR}/variants_to_reclassify.vcf2bed.bed \
                          | cut -f${SOURCEidx} | sed -e 's/_/\t/g' | cut -f2 ||true)
          #Skip variant if no source coordinate is specified
          if ! [ -z ${sourcecoord} ]; then
            #Check if insertion is duplicated
            sourceisdup=$( fgrep -w "${VID}" ${GTDIR}/genotype_counts_per_variant.cleaned_intervals.bed \
                            | bedtools intersect -wa -r -f 0.95 -a - \
                              -b <( echo -e "${sourcecoord}" | sed -e 's/\:/\t/g' -e 's/\-/\t/g' ) \
                            | awk '{ if ($NF=="DUP") print $0 }' | wc -l ||true)
            #Get sink sequence coordinates
            sinkcoord=$( fgrep -w "${VID}" ${GTDIR}/genotype_counts_per_variant.cleaned_intervals.bed \
                            | bedtools intersect -wa -v -r -f 0.95 -a - \
                              -b <( echo -e "${sourcecoord}" | sed -e 's/\:/\t/g' -e 's/\-/\t/g' ) \
                            | awk '{ print $1":"$2"-"$3 }' | head -n1 ||true)
            #Get sink sequence size
            sinksize=$( fgrep -w "${VID}" ${GTDIR}/genotype_counts_per_variant.cleaned_intervals.bed \
                            | bedtools intersect -wa -v -r -f 0.95 -a - \
                              -b <( echo -e "${sourcecoord}" | sed -e 's/\:/\t/g' -e 's/\-/\t/g' ) \
                            | awk '{ print $4 }' | head -n1 ||true)
            #If insertion site size >= $MINSIZEiDEL, check if insertion site is deleted
            if ! [ -z ${sinksize} ] && [ ${sinksize} -ge ${MINSIZEiDEL} ]; then
              sinkisdel=$( fgrep -w "${VID}" ${GTDIR}/genotype_counts_per_variant.cleaned_intervals.bed \
                              | bedtools intersect -wa -v -r -f 0.95 -a - \
                                -b <( echo -e "${sourcecoord}" | sed -e 's/\:/\t/g' -e 's/\-/\t/g' ) \
                              | awk '{ if ($NF=="DEL") print $0 }' | wc -l ||true)
            else
              sinkisdel=0
            fi
            #If source is duplicated and sink is deleted, reclassify as CPX, dDUP_iDEL
            if [ ${sourceisdup} -gt 0 ] && [ ${sinkisdel} -gt 0 ]; then
              svtype="CPX"
              cpxtype="dDUP_iDEL"
              MOD="RECLASSIFY"
              REASON="DISPERSED_DUPLICATION_W_INSERT_SITE_DEL"
              cpxintervals="DUP_${sourcecoord},DEL_${sinkcoord}"
              SOURCE="DUP_${sourcecoord}"
              SVLEN=$( echo "${cpxintervals}" \
                        | sed -e 's/\,/\n/g' -e 's/\:/\t/g' -e 's/\-/\t/g' \
                        | awk '{ sum+=$3-$2 }END{ print sum }' )
            #If source is duplicated and sink is not deleted, reclassify as CPX, dDUP
            elif [ ${sourceisdup} -gt 0 ] && [ ${sinkisdel} -eq 0 ]; then
              svtype="CPX"
              cpxtype="dDUP"
              MOD="RECLASSIFY"
              REASON="DISPERSED_DUPLICATION"
              cpxintervals="DUP_${sourcecoord}"
              SOURCE="DUP_${sourcecoord}"
            #If source is not duplicated but sink is delted, reclassify as CPX, INS_iDEL
            elif [ ${sourceisdup} -eq 0 ] && [ ${sinkisdel} -gt 0 ]; then
              svtype="CPX"
              cpxtype="INS_iDEL"
              MOD="RECLASSIFY"
              REASON="INSERT_SITE_DEL"
              cpxintervals="DEL_${sinkcoord}"
              SOURCE="INS_${sourcecoord}"
              SVLEN=$( echo "${cpxintervals},${sourcecoord}" \
                        | sed -e 's/\,/\n/g' -e 's/\:/\t/g' -e 's/\-/\t/g' \
                        | awk '{ sum+=$3-$2 }END{ print sum }' )
            #Otherwise, leave as canonical insertion
            else
              MOD="KEEP"
              REASON="NON_DUPLICATED_INSERTION"
            fi
          else
            MOD="SKIP"
            REASON="NO_SOURCE_INTERVAL_IN_CURRENT_SHARD"
          fi
        ;;

        ###Solve inverted dispersed duplications vs dupINV / dupINVdel or INVdup / delINVdup
        #DUP5/INS3 or dupINV / dupINVdel
        "DUP5/INS3")
          #Get duplication/insertion interval
          dupinterval=$( fgrep -w "${VID}" ${GTDIR}/variants_to_reclassify.vcf2bed.bed \
                          | cut -f${SOURCEidx} | sed -e 's/_/\t/g' | cut -f2 ||true)
          #Get duplication interval length
          dupsize=$( echo "${dupinterval}" | sed -e 's/\:/\t/g' -e 's/\-/\t/g' \
                      | awk '{ print $3-$2 }' )
          #Test dup interval for evidence of duplication
          dupconfirmed=$( fgrep -w "${VID}" ${GTDIR}/genotype_counts_per_variant.cleaned_intervals.bed \
                            | bedtools intersect -wa -r -f 0.95 -a - \
                              -b <( echo -e "${dupinterval}" | sed -e 's/\:/\t/g' -e 's/\-/\t/g' ) \
                            | awk '{ if ($NF!="WT") print $0 }' | wc -l ||true)
          #As long as dup doesn't explicitly fail depth genotyping, proceed
          if [ ${dupconfirmed} -gt 0 ]; then
            #Get sink interval
            sinkinterval=$( fgrep -w "${VID}" ${GTDIR}/variants_to_reclassify.vcf2bed.bed \
                              | awk '{ print $1":"$2"-"$3 }' ||true)
            #Get sink interval length
            sinksize=$( echo "${sinkinterval}" | sed -e 's/\:/\t/g' -e 's/\-/\t/g' \
                        | awk '{ print $3-$2 }' )
            #If sink length >= $MINSIZEiDEL, test sink interval for evidence of deletion
            if [ ${sinksize} -ge ${MINSIZEiDEL} ]; then
              sinkisdel=$( fgrep -w "${VID}" ${GTDIR}/genotype_counts_per_variant.cleaned_intervals.bed \
                             | bedtools intersect -wa -r -f 0.95 -a - \
                               -b <( echo -e "${sinkinterval}" | sed -e 's/\:/\t/g' -e 's/\-/\t/g' ) \
                             | awk '{ if ($NF=="DEL") print $0 }' | wc -l ||true)
            else
              sinkisdel=0
            fi
            #Get inversion interval. Corresponds to min dup/ins coord and max sink coord
            invinterval=$( paste \
                             <( fgrep -w "${VID}" ${GTDIR}/variants_to_reclassify.vcf2bed.bed | cut -f1 ||true) \
                             <( echo "${dupinterval}" | sed 's/\:/\t/g' | cut -f2 | cut -f1 -d\- ) \
                             <( fgrep -w "${VID}" ${GTDIR}/variants_to_reclassify.vcf2bed.bed | cut -f3 ||true) \
                             | awk '{ print $1":"$2"-"$3 }' )
            #Get inversion size
            invsize=$( echo "${invinterval}" | sed -e 's/\:/\t/g' -e 's/\-/\t/g' \
                         | awk '{ print $3-$2 }' )
            #If inversion length < $MINdDUPTHRESH, classify as dupINV or dupINVdel
            if [ ${invsize} -lt ${MINdDUPTHRESH} ]; then
              #If sink is deleted, classify as dupINVdel
              if [ ${sinkisdel} -gt 0 ]; then
                #Revise inv interval (subtracting del interval)
                invinterval=$( bedtools subtract \
                                 -a <( echo "${invinterval}" | sed -e 's/\:/\t/g' -e 's/\-/\t/g' ) \
                                 -b <( echo "${sinkinterval}" | sed -e 's/\:/\t/g' -e 's/\-/\t/g' ) \
                                 | head -n1 | awk '{ print $1":"$2"-"$3 }' )
                #Update other info
                svtype="CPX"
                cpxtype="dupINVdel"
                MOD="RECLASSIFY"
                REASON="DUP_FLANKED_INVERSION_WITH_DEL"
                cpxintervals="DUP_${dupinterval},INV_${invinterval},DEL_${sinkinterval}"
                SVLEN=${invsize}
                START=$( echo -e "${dupinterval}\n${invinterval}\n${sinkinterval}" \
                           | cut -f2- -d\: | sed 's/\-/\n/g' \
                           | sort -nk1,1 | head -n1 )
                END=$( echo -e "${dupinterval}\n${invinterval}\n${sinkinterval}" \
                           | cut -f2- -d\: | sed 's/\-/\n/g' \
                           | sort -nk1,1 | tail -n1 )
              #If sink is not clearly deleted (or too small), classify as dupINV (no del)
              else
                svtype="CPX"
                cpxtype="dupINV"
                MOD="RECLASSIFY"
                REASON="DUP_FLANKED_INVERSION"
                cpxintervals="DUP_${dupinterval},INV_${invinterval}"
                SVLEN=${invsize}
                START=$( echo -e "${dupinterval}\n${invinterval}" \
                           | cut -f2- -d\: | sed 's/\-/\n/g' \
                           | sort -nk1,1 | head -n1 )
                END=$( echo -e "${dupinterval}\n${invinterval}" \
                           | cut -f2- -d\: | sed 's/\-/\n/g' \
                           | sort -nk1,1 | tail -n1 )
              fi
            #Otherwise, reclassify as dDUP or dDUP_iDEL
            else
              #If sink is deleted, classify as dDUP_iDEL
              if [ ${sinkisdel} -gt 0 ]; then
                svtype="CPX"
                cpxtype="dDUP_iDEL"
                MOD="RECLASSIFY"
                REASON="INVERTED_DISPERSED_DUPLICATION_WITH_DELETION"
                cpxintervals="DUP_${dupinterval},INV_${dupinterval},DEL_${sinkinterval}"
                SOURCE="DUP_${dupinterval}"
                SVLEN=$(( ${dupsize}+${sinksize} ))
              #If sink is not clearly deleted (or too small), classify as dDUP (no iDEL)
              else
                svtype="CPX"
                cpxtype="dDUP"
                MOD="RECLASSIFY"
                REASON="INVERTED_DISPERSED_DUPLICATION"
                cpxintervals="DUP_${dupinterval},INV_${dupinterval}"
                SOURCE="DUP_${dupinterval}"
                SVLEN="${dupsize}"
              fi
            fi
          #If dup explicitly fails depth genotyping, mark as unresolved
          else
            MOD="UNRESOLVED"
            REASON="PREDICTED_DUP_INTERVAL_FAILED_GT"
          fi
        ;;

        #DUP3/INS5 or INVdup / delINVdup
        "DUP3/INS5")
          #Get duplication/insertion interval
          dupinterval=$( fgrep -w "${VID}" ${GTDIR}/variants_to_reclassify.vcf2bed.bed \
                          | cut -f${SOURCEidx} | sed -e 's/_/\t/g' | cut -f2 ||true)
          #Get duplication interval length
          dupsize=$( echo "${dupinterval}" | sed -e 's/\:/\t/g' -e 's/\-/\t/g' \
                      | awk '{ print $3-$2 }' )
          #Test dup interval for evidence of duplication
          dupconfirmed=$( fgrep -w "${VID}" ${GTDIR}/genotype_counts_per_variant.cleaned_intervals.bed \
                            | bedtools intersect -wa -r -f 0.95 -a - \
                              -b <( echo -e "${dupinterval}" | sed -e 's/\:/\t/g' -e 's/\-/\t/g' ) \
                            | awk '{ if ($NF!="WT") print $0 }' | wc -l ||true)
          #As long as dup doesn't explicitly fail depth genotyping, proceed
          if [ ${dupconfirmed} -gt 0 ]; then
            #Get sink interval
            sinkinterval=$( fgrep -w "${VID}" ${GTDIR}/variants_to_reclassify.vcf2bed.bed \
                              | awk '{ print $1":"$2"-"$3 }' ||true)
            #Get sink interval length
            sinksize=$( echo "${sinkinterval}" | sed -e 's/\:/\t/g' -e 's/\-/\t/g' \
                        | awk '{ print $3-$2 }' )
            #If sink length >= $MINSIZEiDEL, test sink interval for evidence of deletion
            if [ ${sinksize} -ge ${MINSIZEiDEL} ]; then
              sinkisdel=$( fgrep -w "${VID}" ${GTDIR}/genotype_counts_per_variant.cleaned_intervals.bed \
                             | bedtools intersect -wa -r -f 0.95 -a - \
                               -b <( echo -e "${sinkinterval}" | sed -e 's/\:/\t/g' -e 's/\-/\t/g' ) \
                             | awk '{ if ($NF=="DEL") print $0 }' | wc -l ||true)
            else
              sinkisdel=0
            fi
            #Get inversion interval. Corresponds to min sink coord and max dup/inv coord
            invinterval=$( paste \
                             <( fgrep -w "${VID}" ${GTDIR}/variants_to_reclassify.vcf2bed.bed | cut -f1 ||true) \
                             <( fgrep -w "${VID}" ${GTDIR}/variants_to_reclassify.vcf2bed.bed | cut -f2 ||true) \
                             <( echo "${dupinterval}" | sed 's/\:/\t/g' | cut -f2 | cut -f2 -d\- ) \
                             | awk '{ print $1":"$2"-"$3 }' )
            #Get inversion size
            invsize=$( echo "${invinterval}" | sed -e 's/\:/\t/g' -e 's/\-/\t/g' \
                         | awk '{ print $3-$2 }' )
            #If inversion length < $MINdDUPTHRESH, classify as INVdup or delINVdup
            if [ ${invsize} -lt ${MINdDUPTHRESH} ]; then
              #If sink is deleted, classify as delINVdup
              if [ ${sinkisdel} -gt 0 ]; then
                #Revise inv interval (subtracting del interval)
                invinterval=$( bedtools subtract \
                                 -a <( echo "${invinterval}" | sed -e 's/\:/\t/g' -e 's/\-/\t/g' ) \
                                 -b <( echo "${sinkinterval}" | sed -e 's/\:/\t/g' -e 's/\-/\t/g' ) \
                                 | head -n1 | awk '{ print $1":"$2"-"$3 }' )
                #Update other info
                svtype="CPX"
                cpxtype="delINVdup"
                MOD="RECLASSIFY"
                REASON="DUP_FLANKED_INVERSION_WITH_DEL"
                cpxintervals="DEL_${sinkinterval},INV_${invinterval},DUP_${dupinterval}"
                SVLEN=${invsize}
                START=$( echo -e "${dupinterval}\n${invinterval}\n${sinkinterval}" \
                           | cut -f2- -d\: | sed 's/\-/\n/g' \
                           | sort -nk1,1 | head -n1 )
                END=$( echo -e "${dupinterval}\n${invinterval}\n${sinkinterval}" \
                           | cut -f2- -d\: | sed 's/\-/\n/g' \
                           | sort -nk1,1 | tail -n1 )
              #If sink is not clearly deleted (or too small), classify as dupINV (no del)
              else
                svtype="CPX"
                cpxtype="INVdup"
                MOD="RECLASSIFY"
                REASON="DUP_FLANKED_INVERSION"
                cpxintervals="INV_${invinterval},DUP_${dupinterval}"
                SVLEN=${invsize}
                START=$( echo -e "${dupinterval}\n${invinterval}" \
                           | cut -f2- -d\: | sed 's/\-/\n/g' \
                           | sort -nk1,1 | head -n1 )
                END=$( echo -e "${dupinterval}\n${invinterval}" \
                           | cut -f2- -d\: | sed 's/\-/\n/g' \
                           | sort -nk1,1 | tail -n1 )
              fi
            #Otherwise, reclassify as dDUP or dDUP_iDEL
            else
              #If sink is deleted, classify as dDUP_iDEL
              if [ ${sinkisdel} -gt 0 ]; then
                svtype="CPX"
                cpxtype="dDUP_iDEL"
                MOD="RECLASSIFY"
                REASON="INVERTED_DISPERSED_DUPLICATION_WITH_DELETION"
                cpxintervals="DUP_${dupinterval},INV_${dupinterval},DEL_${sinkinterval}"
                SOURCE="DUP_${dupinterval}"
                SVLEN=$(( ${dupsize}+${sinksize} ))
              #If sink is not clearly deleted (or too small), classify as dDUP (no iDEL)
              else
                svtype="CPX"
                cpxtype="dDUP"
                MOD="RECLASSIFY"
                REASON="INVERTED_DISPERSED_DUPLICATION"
                cpxintervals="INV_${dupinterval},DUP_${dupinterval}"
                SOURCE="DUP_${dupinterval}"
                SVLEN="${dupsize}"
              fi
            fi
          #If dup explicitly fails depth genotyping, mark as unresolved
          else
            MOD="UNRESOLVED"
            REASON="PREDICTED_DUP_INTERVAL_FAILED_GT"
          fi
        ;;

        #Otherwise, just test for INS-iDEL
        *)
          #Get inserted sequence coordinates
          sourcecoord=$( fgrep -w "${VID}" ${GTDIR}/variants_to_reclassify.vcf2bed.bed \
                          | cut -f${SOURCEidx} | sed -e 's/_/\t/g' | cut -f2 ||true)

          #Get sink sequence coordinates
          sinkcoord=$( fgrep -w "${VID}" ${GTDIR}/variants_to_reclassify.vcf2bed.bed \
                          | awk '{ print $1":"$2"-"$3 }' ||true)
          #Skip variant if no sink coordinate is specified
          if ! [ -z ${sinkcoord} ]; then
            #Get sink sequence size
            sinksize=$( fgrep -w "${VID}" ${GTDIR}/genotype_counts_per_variant.cleaned_intervals.bed \
                            | bedtools intersect -wa -v -r -f 0.95 -a - \
                              -b <( echo -e "${sinkcoord}" | sed -e 's/\:/\t/g' -e 's/\-/\t/g' ) \
                            | awk '{ print $4 }' | head -n1 ||true)
            #If insertion site size >= $MINSIZEiDEL, check if insertion site is deleted
            if ! [ -z ${sinksize} ] && [ ${sinksize} -ge ${MINSIZEiDEL} ]; then
              sinkisdel=$( fgrep -w "${VID}" ${GTDIR}/genotype_counts_per_variant.cleaned_intervals.bed \
                              | bedtools intersect -wa -v -r -f 0.95 -a - \
                                -b <( echo -e "${sinkcoord}" | sed -e 's/\:/\t/g' -e 's/\-/\t/g' ) \
                              | awk '{ if ($NF=="DEL") print $0 }' | wc -l ||true)
              sinkisnotdel=$( fgrep -w "${VID}" ${GTDIR}/genotype_counts_per_variant.cleaned_intervals.bed \
                              | bedtools intersect -wa -v -r -f 0.95 -a - \
                                -b <( echo -e "${sinkcoord}" | sed -e 's/\:/\t/g' -e 's/\-/\t/g' ) \
                              | awk '{ if ($NF=="WT") print $0 }' | wc -l ||true)
            else
              sinkisdel=0
              sinkisnotdel=0
            fi
            #If sink is deleted, reclassify as CPX, INS_iDEL
            if [ ${sourceisdup} -eq 0 ] && [ ${sinkisdel} -gt 0 ]; then
              svtype="CPX"
              cpxtype="INS_iDEL"
              MOD="RECLASSIFY"
              REASON="INSERT_SITE_DEL"
              cpxintervals="DEL_${sinkcoord}"
              SOURCE="INS_${sourcecoord}"
              SVLEN=$( echo "${cpxintervals},${sourcecoord}" \
                        | sed -e 's/\,/\n/g' -e 's/\:/\t/g' -e 's/\-/\t/g' \
                        | awk '{ sum+=$3-$2 }END{ print sum }' )
            #If sink is explicitly *not* deleted and is >= $MINSIZE, reclassify as UNRESOLVED
            elif [ ${sinkisnotdel} -gt 0 ] && [ ${sinksize} -ge ${MINSIZE} ]; then
              MOD="UNRESOLVED"
              REASON="PREDICTED_LARGE_SINK_DEL_INTERVAL_FAILED_GT"
            #Otherwise, leave as canonical insertion
            else
              MOD="KEEP"
              REASON="CANONICAL_INS_NO_SINK_DELETION"
            fi
          else
            MOD="SKIP"
            REASON="NO_SINK_INTERVAL_IN_CURRENT_SHARD"
          fi
        ;;
      esac
      ;;

    #Salvage inverted duplications from inversion single enders
    "BND")
      #Only consider inversion single enders
      if [ ${cpxtype} == "INVERSION_SINGLE_ENDER_--" ] \
         || [ ${cpxtype} == "INVERSION_SINGLE_ENDER_++" ]; then
        #Get span of inversion
        dupinterval=$( fgrep -w "${VID}" ${GTDIR}/variants_to_reclassify.vcf2bed.bed \
                        | awk '{ print $1":"$2"-"$3 }' ||true)
        #Get duplication interval length
        dupsize=$( echo "${dupinterval}" | sed -e 's/\:/\t/g' -e 's/\-/\t/g' \
                    | awk '{ print $3-$2 }' )
        #Test dup interval for explicit confirmation of duplication
        dupconfirmed=$( fgrep -w "${VID}" ${GTDIR}/genotype_counts_per_variant.cleaned_intervals.bed \
                          | bedtools intersect -wa -r -f 0.95 -a - \
                            -b <( echo -e "${dupinterval}" | sed -e 's/\:/\t/g' -e 's/\-/\t/g' ) \
                          | awk '{ if ($NF=="DUP") print $0 }' | wc -l ||true)
        #If dup confirms, resolve as palindromic inverted duplication
        if [ ${dupconfirmed} -gt 0 ]; then
          svtype="CPX"
          MOD="RECLASSIFY"
          REASON="PALINDROMIC_INVERTED_DUPLICATION"
          cpxintervals="DUP_${dupinterval},INV_${dupinterval}"
          SVLEN=${dupsize}
          START=$( echo -e "${dupinterval}" \
                     | cut -f2- -d\: | sed 's/\-/\n/g' \
                     | sort -nk1,1 | head -n1 )
          END=$( echo -e "${dupinterval}" \
                   | cut -f2- -d\: | sed 's/\-/\n/g' \
                   | sort -nk1,1 | tail -n1 )
          #RR single enders = piDUP_FR
          if [ ${cpxtype} == "INVERSION_SINGLE_ENDER_--" ]; then
            cpxtype="piDUP_FR"
          else
            cpxtype="piDUP_RF"
          fi
        #Otherwise, leave as unresolved
        else
          MOD="KEEP"
          REASON="DID_NOT_RESOLVE_AS_piDUP"
        fi
      fi
    ;;

    *)
      MOD="KEEP"
      REASON="IRRELEVANT_SV_TYPE"
      ;;
  esac

  #Leave new intervals and SVLEN as unset if none specified
  if [ -z ${cpxintervals} ]; then
    cpxintervals="."
  fi
  if [ -z ${SVLEN} ]; then
    SVLEN="."
  fi
  if [ -z ${SOURCE} ]; then
    SOURCE="."
  fi
  if [ -z ${START} ]; then
    START="."
  fi
  if [ -z ${END} ]; then
    END="."
  fi

  #Print final result per variant
  echo -e "${VID}\t${MOD}\t${REASON}\t${svtype}\t${cpxtype}\t${cpxintervals}\t${SVLEN}\t${SOURCE}\t${START}\t${END}"

  #Unset all intermediate variables
  unset nval; unset MOD; unset REASON; unset cpxintervals; unset SVLEN; unset SOURCE; unset START; unset END

done < <( fgrep -v "#" ${GTDIR}/genotype_counts_per_variant.cleaned_intervals.bed \
            | cut -f6 | sort | uniq ||true) \
  >> ${GTDIR}/final_variant_reassessment_table.txt
#Copy final reclassification table, if optioned
if ! [ -z ${RTABLE} ]; then
  cp ${GTDIR}/final_variant_reassessment_table.txt ${RTABLE}
fi
#Copy genotyping counts, if optioned
if ! [ -z ${GTCOUNTSTABLE} ]; then
  cp ${GTDIR}/genotype_counts_per_variant.cleaned_intervals.bed ${GTCOUNTSTABLE}
fi


###GENERATE FINAL VCF WITH VARIANT REASSESSMENTS
#Write header to final VCF
zcat ${INVCF} | fgrep "##" \
  > ${GTDIR}/cleaned_output.vcf
#Add lines corresponding to complex types
for wrapper in 1; do
  echo -e "##CPX_TYPE_delINV=\"Complex inversion with 5' flanking deletion.\""
  echo -e "##CPX_TYPE_INVdel=\"Complex inversion with 3' flanking deletion.\""
  echo -e "##CPX_TYPE_dupINV=\"Complex inversion with 5' flanking duplication.\""
  echo -e "##CPX_TYPE_INVdup=\"Complex inversion with 3' flanking duplication.\""
  echo -e "##CPX_TYPE_delINVdel=\"Complex inversion with 5' and 3' flanking deletions.\""
  echo -e "##CPX_TYPE_dupINVdup=\"Complex inversion with 5' and 3' flanking duplications.\""
  echo -e "##CPX_TYPE_delINVdup=\"Complex inversion with 5' flanking deletion and 3' flanking duplication.\""
  echo -e "##CPX_TYPE_dupINVdel=\"Complex inversion with 5' flanking duplication and 3' flanking deletion.\""
  echo -e "##CPX_TYPE_piDUP_FR=\"Palindromic inverted tandem duplication, forward-reverse orientation.\""
  echo -e "##CPX_TYPE_piDUP_RF=\"Palindromic inverted tandem duplication, reverse-forward orientation.\""
  echo -e "##CPX_TYPE_dDUP=\"Dispersed duplication.\""
  echo -e "##CPX_TYPE_dDUP_iDEL=\"Dispersed duplication with deletion at insertion site.\""
  echo -e "##CPX_TYPE_INS_iDEL=\"Insertion with deletion at insertion site.\""
done >> ${GTDIR}/cleaned_output.vcf
#Add samples line to VCF header
zcat ${INVCF} | fgrep "#" | fgrep -v "##" \
  >> ${GTDIR}/cleaned_output.vcf
#Write new vcf without variants in the final reassessment table
zcat ${INVCF} | fgrep -v "#" \
  | fgrep -wvf <( cut -f1 ${GTDIR}/final_variant_reassessment_table.txt | fgrep -v "#" ) \
  >> ${GTDIR}/cleaned_output.vcf ||true
#Get smaller VCF (no header) with all variants in final assessment table
zcat ${INVCF} \
  | fgrep -wf <( cut -f1 ${GTDIR}/final_variant_reassessment_table.txt | fgrep -v "#" ) \
  > ${GTDIR}/variants_to_be_reassessed.vcf||true
#Next, iterate over the final variant reassessment table & take action for each row
while read VID MOD REASON svtype cpxtype cpxintervals SVLEN SOURCE START END; do
  #Perform different actions based on $MOD flag
  case ${MOD} in
    #Emit KEEP or SKIP variants as-is
    "SKIP"|"KEEP")
      fgrep -w ${VID} ${GTDIR}/variants_to_be_reassessed.vcf
      ;;

    #Reassign UNRESOLVED variants to unresolved & scrub previous info
    "UNRESOLVED")

      #Set filter status to UNRESOLVED and add UNRESOLVED
      awk -v vid=${VID} \
         '$3 == vid {OFS="\t"; $7 = "UNRESOLVED"; $8 = $8";UNRESOLVED;UNRESOLVED_TYPE=POSTHOC_RD_GT_REJECTION"; print;}' \
         ${GTDIR}/variants_to_be_reassessed.vcf
      ;;

    #Revise records for RECLASSIFY variants as appropriate
    "RECLASSIFY")
      #Get missing SVLEN, START, and END from existing record if needed
      if [ ${SVLEN} == "." ]; then
        SVLEN=$( fgrep -w ${VID} ${GTDIR}/variants_to_be_reassessed.vcf \
          | cut -f8 \
          | sed 's/\;/\n/g' \
          | fgrep SVLEN= \
          | cut -f2 -d\= )
      fi
      if [ ${START} == "." ]; then
        START=$( fgrep -w ${VID} ${GTDIR}/variants_to_be_reassessed.vcf | cut -f2 )
      fi
      if [ ${END} == "." ]; then
        END=$( fgrep -w ${VID} ${GTDIR}/variants_to_be_reassessed.vcf \
          | cut -f8 \
          | sed 's/\;/\n/g' \
          | fgrep END= \
          | cut -f2 -d\= )
      fi
      #Modify info as needed
      INFO=$( fgrep -w ${VID} ${GTDIR}/variants_to_be_reassessed.vcf \
                | cut -f8 \
                | sed -r -e "s/^END=[^;]*;/END=$END;/" \
                | sed -r -e "s/;END=[^;]*;/;END=$END;/" \
                | sed -r -e "s/;END=[^;]*$/;END=$END/" \
                | sed -r -e "s/^SVTYPE=[^;]*;/SVTYPE=$svtype;/" \
                | sed -r -e "s/;SVTYPE=[^;]*;/;SVTYPE=$svtype;/" \
                | sed -r -e "s/;SVTYPE=[^;]*$/;SVTYPE=$svtype/" \
                | sed -r -e "s/^SVLEN=[^;]*;/SVLEN=$SVLEN;/" \
                | sed -r -e "s/;SVLEN=[^;]*;/;SVLEN=$SVLEN;/" \
                | sed -r -e "s/;SVLEN=[^;]*$/;SVLEN=$SVLEN/" \
                | sed -r -e "s/^CPX_TYPE=[^;]*;/CPX_TYPE=${cpxtype};/" \
                | sed -r -e "s/;CPX_TYPE=[^;]*;/;CPX_TYPE=${cpxtype};/" \
                | sed -r -e "s/;CPX_TYPE=[^;]*$/;CPX_TYPE=${cpxtype}/" \
                | sed -r -e 's/^UNRESOLVED;//' \
                | sed -r -e 's/;UNRESOLVED;/;/' \
                | sed -r -e 's/;UNRESOLVED$//' \
                | sed -r -e 's/^UNRESOLVED_TYPE=[^;]*;//' \
                | sed -r -e 's/;UNRESOLVED_TYPE=[^;]*;/;/' \
                | sed -r -e 's/;UNRESOLVED_TYPE=[^;]*$//' \
                | sed -r -e 's/^EVENT=[^;]*;//' \
                | sed -r -e 's/;EVENT=[^;]*;/;/' \
                | sed -r -e 's/;EVENT=[^;]*$//' )
      #Add/remove/modify CPX_TYPE, if needed
      if [ ${svtype} == "CPX" ]; then
        if [ $( echo "${INFO}" | fgrep "CPX_TYPE" | wc -l ) -eq 0 ]; then
          INFO=$( echo ${INFO} | sed -e "s/\$/;CPX_TYPE=$cpxtype/" )
        fi
      else
        INFO=$( echo "${INFO}" \
                  | sed -r -e 's/;CPX_TYPE=[^;]*;/;/' \
                  | sed -r -e 's/;CPX_TYPE=[^;]*$//' )
      fi
      #Add/remove/correct SOURCE, as needed
      if [ ${SOURCE} == "." ]; then
        INFO=$( echo ${INFO} \
                  | sed -r -e "s/;SOURCE=[^;]*;/;/" \
                  | sed -r -e "s/;SOURCE=[^;]*$//" )
      else
        if [ $( echo ${INFO} | fgrep SOURCE | wc -l ) -gt 0 ]; then
          INFO=$( echo ${INFO} \
                    | sed -r -e "s/;SOURCE=[^;]*;/;SOURCE=${SOURCE};/" \
                    | sed -r -e "s/;SOURCE=[^;]*$/;SOURCE=${SOURCE}/" )
        else
          INFO=$( echo ${INFO} \
                    | sed -r -e "s/$/;SOURCE=${SOURCE}/" )
        fi
      fi
      #Add/remove/correct CPX_INTERVALS, as needed
      if [ ${cpxintervals} == "." ]; then
        INFO=$( echo ${INFO} \
                  | sed -r -e "s/;CPX_INTERVALS=[^;]*;/;/" \
                  | sed -r -e "s/;CPX_INTERVALS=[^;]*$//" )
      else
        if [ $( echo ${INFO} | fgrep CPX_INTERVALS | wc -l ) -gt 0 ]; then
          INFO=$( echo ${INFO} \
                    | sed -r -e "s/;CPX_INTERVALS=[^;]*;/;CPX_INTERVALS=${cpxintervals};/" \
                    | sed -r -e "s/;CPX_INTERVALS=[^;]*$/;CPX_INTERVALS=${cpxintervals}/" )
        else
          INFO=$( echo ${INFO} \
                    | sed -r -e "s/$/;CPX_INTERVALS=${cpxintervals}/" )
        fi
      fi
      #Print record
      paste \
        <( fgrep -w ${VID} ${GTDIR}/variants_to_be_reassessed.vcf | cut -f1-4 \
             | awk -v startpos=${START} -v OFS="\t" '{ print $1, startpos, $3, $4 }' ) \
        <( echo "<${svtype}>" ) \
        <( fgrep -w ${VID} ${GTDIR}/variants_to_be_reassessed.vcf | cut -f6-7 ) \
        <( echo "${INFO}" ) \
        <( fgrep -w ${VID} ${GTDIR}/variants_to_be_reassessed.vcf | cut -f9- )
      ;;
  esac
done < <( fgrep -v "#" ${GTDIR}/final_variant_reassessment_table.txt ) \
  >> ${GTDIR}/cleaned_output.vcf
#Sort & compress final output
#Also strip CPX_TYPE and CPX_INTERVALS from any non-CPX variants
${BIN}/rm_cpx_info.py ${GTDIR}/cleaned_output.vcf /dev/stdout \
  | vcf-sort | bgzip -c \
  > ${OUTVCF}


#Clean up
rm -rf ${GTDIR}

