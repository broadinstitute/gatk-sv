#!/bin/bash
#
# overlapbpchange.sh
#

set -euo pipefail

##Inputs##
vcf=$1
##sr fail##
backgroundlist=$2
##sr support on both sides##
bothendSR=$3

##clean out variants that overlap at one site##
##pull out variants with duplicate bp that are not driven by depth which will be integrated in the clean vcf##
##make sure to flip bed as well so second bp location can be compared with first from other variants##
svtk vcf2bed $vcf stdout -i CHR2 -i STRANDS -i SVLEN -i varGQ -i END -i EVIDENCE -i SVTYPE --split-bnd  \
  | sed "s/+-/+ "$'\t -/g' \
  | sed "s/-+/- "$'\t +/g' \
  | sed "s/++/+ "$'\t +/g' \
  | sed "s/--/- "$'\t -/g' | \
  ##Convert back to 1-based positions##
  awk -v OFS='\t' '{$2=$2+1; print $0}' \
  | awk -v OFS='\t' \
    '{if (!(($NF=="DEL" || $NF=="DUP") && $10>=5000)) print $0 "\n" $7,$12,$2,$4,$5,$6,$1,$9,$8,$10,$11,$2,$13,$14 }' | \
  ###Find duplicated variants that overlap at same bp one side##
  awk 'cnt[$1"_"$2"_"$8]++{if (cnt[$1"_"$2"_"$8]==2) print prev[$1"_"$2"_"$8] "\t" $1"_"$2"_"$8 \
  ; print $0 "\t" $1"_"$2"_"$8} {prev[$1"_"$2"_"$8]=$0}' \
  | awk '!seen[$4"_"$NF]++' \
  | awk 'cnt[$NF]++{if (cnt[$NF]==2) print prev[$NF] \
  ; print $0 } {prev[$NF]=$0}' \
  >dupside1.bed


##Find 50% overlap between samples for overlaps##
join -j 2 <(awk '{print $NF "\t" $6}' dupside1.bed \
            | awk -F'[,\t]' '{for (i=2;i<=NF;i++) print $1 "\t" $i}' \
            | sort \
            | uniq -D \
            | awk '{print $1}'|sort|uniq -c  ) \
          <(awk '{print $NF "\t" $6}' dupside1.bed \
            | awk -F'[,\t]' '{for (i=2;i<=NF;i++) print $1 "\t" $i}' \
            | awk '{print $1}' \
            | sort \
            | uniq -c) \
  | awk '{if ($2 >= 0.5 * $3) print $1}' \
  | (fgrep -wf - dupside1.bed || printf "") \
  > dupside1.freq50.txt

##Add SRfail###
{ fgrep -wf <(awk '{print $NF}' $backgroundlist) dupside1.freq50.txt || true; } \
  | awk '{print $0 "\t" 0}' \
  > dupside1.passSR.txt

{ fgrep -wvf <(awk '{print $NF}' $backgroundlist) dupside1.freq50.txt || true; } \
  | awk '{print $0 "\t" 1}' \
  >> dupside1.passSR.txt

##Attach the % of variants that show SR support at bothends##
join -1 4 -2 1 -e "0" -a 1 -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 1.10 1.11 1.12 1.13 1.14 1.15 1.16 2.2 \
   <(sort -k4,4 dupside1.passSR.txt) \
   <(awk '{print $NF "\t" $1}' $bothendSR | sort -k1,1) \
   | tr ' ' '\t' \
   > dupside1.bothpassfilter.txt
rm dupside1.passSR.txt

##count number of samples and indiciate if size gt 50bp##
join -1 4 -2 1 dupside1.bothpassfilter.txt \
   <(awk '{print $4 "\t" $6}' dupside1.bed \
   | awk -F'[,\t]' '{print $1 "\t" NF-1}' \
   | sort -k1,1) \
   | tr ' ' '\t' \
   | awk '{if ($10>=50) print $0 "\t" 1;else print $0 "\t" 0}' \
   > dupside1.samplecountfilter.txt
rm dupside1.bed dupside1.bothpassfilter.txt

##Convert Evidence column into Integers for scoring and ##
##RD,PE,SR-1,RD,PE-2,PE,SR-3,RD,SR-4,PE-5,RD-6,SR-7##
sed 's/BAF,//g' dupside1.samplecountfilter.txt \
  | awk -v OFS='\t' '
        {
          if ($13=="PE,RD,SR") print $0 "\t" 1
          else if ($13=="PE,RD") print $0 "\t" 2
          else if ($13=="PE,SR") print $0 "\t" 3
          else if ($13=="RD,SR") print $0 "\t" 4
          else if ($13=="PE") print $0 "\t" 5
          else if ($13=="RD") print $0 "\t" 6
          else if ($13=="SR") print $0 "\t" 7
        }' | \
  ##assign BND to bottom
  awk '{if ($14=="BND") print $0 "\t" 0;else print $0 "\t" 1}' \
  > dupside1.allfilter.txt
rm dupside1.samplecountfilter.txt
###DO THIS#####
##


##sort file with overlapping samples LevelofSupport->BothEndsupport->SRfail-> Not BND->Higher varq-> Higher Freq -> Smallest size if gt 5kb##
sort -k20,20n -k17,17nr -nrk16,16 -k21,21nr -k11,11nr -k18,18nr  -k19,19nr -k10,10n dupside1.allfilter.txt \
  | awk '!seen[$15]++' \
  | awk '{print $1}' \
  | (fgrep -wvf - dupside1.freq50.txt || printf "") \
  | awk '{print $4}' \
  > remove.side1.var.txt
rm dupside1.freq50.txt dupside1.allfilter.txt

##remove variants with samebp##
(zgrep -wvf remove.side1.var.txt $vcf || printf "") \
  | bgzip \
  > non_redundant.vcf.gz
