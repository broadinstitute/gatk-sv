#!/bin/bash
#
# clean_vcf_part2.sh
#

set -euo pipefail

##bgzipped combined vcf from clean vcf part1b.sh##
normal_revise_vcf_gz=$1
##whitelist of split ids for parallelization##
whitelist=$2
##list of multiallelic CNVs from 1b##
multi_cnv=$3
##output filename##
outputfile=$4

export LC_ALL=C

##subset vcf to whitelist samples##
bcftools view $normal_revise_vcf_gz -S $whitelist --no-update \
  | gzip \
  > subset.vcf.gz


##create new bed with updated genotypes###
zcat subset.vcf.gz \
  | awk 'BEGIN{FS=OFS="\t"}{if (substr($1,1,1)=="#" || $5=="<DEL>" || $5=="<DUP>") print}' \
  | svtk vcf2bed stdin stdout \
  | sed 1d \
  | sort -k4,4 \
  | gzip \
  > int.afternormalfix.bed.gz

##Find overlapping depth based variants and reassign depth based; note this is necessary because depth call >5kb genotypes are 100% driven by depth ##
##generate a sample list based on depth for depth overlap check below. Necessary because genotype is capped at 1/1 and by direction (i.e no dels in dups)##
##grab all samples per variant with a non normal copy state##
zcat subset.vcf.gz \
  | awk 'BEGIN{FS=OFS="\t"}{if (substr($1,1,1)=="#") print; else if ($5=="<DEL>" || $5=="<DUP>") {$1=$3; print}}' \
  | vcftools --vcf - --stdout --extract-FORMAT-info RD_CN \
  | awk 'BEGIN{FS=OFS="\t"} NR==1{for (i=3;i<=NF;i++) header[i]=$i} NR>1{for(j=3;j<=NF;j++) print $1,header[j],$j }' \
  | sort -k1,1 \
  | awk 'BEGIN{FS=OFS="\t"} \
              {if (($1~"DEL" && $3<2 && $3!=".") || ($1~"DUP" && $3>2 && $3!=".")) a[$1]=a[$1]?a[$1]","$2:$2} \
         END{for (i in a) print i,a[i]}' \
  | sort -k1,1 \
  > afternormal.combined.RD_CN.list.txt

##get a list of samples for actual variants not just those with aberrant copy states##
zcat int.afternormalfix.bed.gz \
  | awk 'BEGIN{FS=OFS="\t"} $6!="" {split($6,samples,","); for (i in samples) print $4"@"samples[i]}' \
  | sort -k1,1 \
  > fullvar.afternormal.list.txt

##create bed with anything that has abnormal copy state##
##do not compress all.bed because of bedtools bug: https://github.com/arq5x/bedtools2/issues/643##
zcat int.afternormalfix.bed.gz \
  | cut -f1-5 \
  | awk 'BEGIN{FS=OFS="\t"}{if($3-$2>=5000)print $0}' \
  | join -1 4 -2 1 -t '	' - afternormal.combined.RD_CN.list.txt \
  | awk 'BEGIN {FS=OFS="\t"} \
         $6!=""{split($6,samples,","); \
                for (i in samples) {s = samples[i]; print $2"_"s,$3,$4,$1,$5,s,$1"@"s}}' \
  | sort -k7,7 \
  | join -t '	' -a 1 -1 7 -2 1 -e "NA" -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,2.1  - fullvar.afternormal.list.txt \
  > all.bed

##intersect variants and always set larger to left##
bedtools intersect -wa -wb -a all.bed -b all.bed \
  | awk 'BEGIN{FS=OFS="\t"} \
              {if ($8=="NA") { if ($16!="NA") print $9,$10,$11,$12,$13,$14,$15,$1,$2,$3,$4,$5,$6,$7; } \
               else if ($4!=$12) {
                if ($3-$2>=$11-$10) print $1,$2,$3,$4,$5,$6,$7,$9,$10,$11,$12,$13,$14,$15; \
                else if ($16!="NA") print $9,$10,$11,$12,$13,$14,$15,$1,$2,$3,$4,$5,$6,$7; \
                else print $1,$2,$3,$4,$5,$6,$7,$9,$10,$11,$12,$13,$14,$15;}}' \
  | sort -k7,7 \
  | uniq \
  | gzip \
  > bed.overlap.txt.gz

zcat subset.vcf.gz \
  | grep '^#' \
  > vcf_header
zcat bed.overlap.txt.gz \
  | awk '{print $4; print $11}' \
  | sort -u \
  > overlap.events
echo '~~~' >> overlap.events
zcat subset.vcf.gz \
  | grep -v '^#' \
  | sort -k3,3 \
  | join -t '	' -1 1 -2 3 overlap.events - \
  | awk 'BEGIN{FS=OFS="\t"}{tmp=$2;$2=$3;$3=tmp;print}' \
  | sort -k2n,2 \
  > overlap.events.vcf

##get info for each variant##
zcat bed.overlap.txt.gz \
  | awk 'BEGIN{FS=OFS="\t"}{print $7; print $14;}' \
  | sort -u \
  > combined.bed
for var in EV RD_CN GT
do
 echo '~~~' >> combined.bed
 cat vcf_header overlap.events.vcf \
  | vcftools --vcf - --stdout --extract-FORMAT-info ${var} \
  | awk 'BEGIN{FS=OFS="\t"} \
         NR==1{for (i=3;i<=NF;i++) header[i]=$i} \
         NR>1 {for (j=3;j<=NF;j++) print $1"@"header[j],$j}' \
  | sort -k1,1 \
  | join -t '	' -1 1 -2 1 combined.bed - \
  > combined.bed.tmp
  mv combined.bed.tmp combined.bed
done

zcat bed.overlap.txt.gz \
  | join -t '	' -1 7 -2 1 - combined.bed \
  | cut -f2- \
  | sort -k13,13 \
  | join -t '	' -1 13 -2 1 - combined.bed \
  | cut -f2- \
  | awk 'BEGIN{FS=OFS="\t"}{print $3-$2,$9-$8,$0}' \
  | sort -k1nr,1 -k2nr,2 \
  | cut -f3- \
  | gzip \
  > all.combined.bed.gz

zcat all.combined.bed.gz \
  | awk -v multiCNVFile=$multi_cnv ' \
     function makeRevision( id, val ) { reviseCN[id] = val; if ( val == 2 ) wasRevisedToNormal[id] = 1; } \
     BEGIN \
       {FS=OFS="\t"; \
        while ( getline < multiCNVFile ) multiCNV[$0] = 1; \
        close(multiCNVFile)} \
       {chr_sample1 = $1; start1 = $2; stop1 = $3; ev1 = $4; svtype1 = $5; sample1 = $6; \
        chr_sample2 = $7; start2 = $8; stop2 = $9; ev2 = $10; svtype2 = $11; sample2 = $12; \
        support1 = $13; RD_CN1 = $14; GT1 = $15; support2 = $16; RD_CN2 = $17; GT2 = $18; \
        id1 = ev1"@"sample1; id2 = ev2"@"sample2; \
        length1 = stop1 - start1; length2 = stop2 - start2; \
        if ( id1 in wasRevisedToNormal ) next; \
        if ( id1 in reviseCN ) RD_CN1 = reviseCN[id1]; \
        if ( id2 in reviseCN ) RD_CN2 = reviseCN[id2]; \
        overlap = (stop1 < stop2 ? stop1 : stop2) - (start1 > start2 ? start1 : start2); \
        smallOverlap50 = overlap/length2 > .5; \
        largeOverlap50 = overlap/length1 > .5; \
        ##Call where smaller depth call is being driven by larger## \
        if ( support1 ~ /RD/ && support1 != "RD" && support2 == "RD" && smallOverlap50 && !(ev1 in multiCNV) ) { \
          if ( RD_CN1 == 0 ) makeRevision(id2, RD_CN2 + 2); \
          else if ( RD_CN1 == 1 ) makeRevision(id2, RD_CN2 + RD_CN1); \
          else if ( RD_CN1 > 1 ) { newCN = RD_CN2 - RD_CN1 + 2; if ( newCN < 0 ) newCN = 0; makeRevision(id2, newCN); } } \
        ##Smaller CNV driving larger CNV genotype## \
        else if ( support1 == "RD" && support2 ~ /RD/ && support2 != "RD" && smallOverlap50 && !(ev2 in multiCNV) && GT2 != "0/0" && largeOverlap50 ) { \
          if ( RD_CN2 == 0 ) makeRevision(id1, RD_CN1 + 2); \
          else if ( RD_CN2 == 1 ) makeRevision(id1, RD_CN1 + RD_CN2); \
          else if ( RD_CN2 > 1 ) { newCN = RD_CN1 - RD_CN2 + 2; if ( newCN < 0 ) newCN = 0; makeRevision(id1, newCN); } } \
        ##Depth only calls where smaller call is being driven by larger## \
        else if ( support1 == "RD" && support2 == "RD" && smallOverlap50 && svtype1 == svtype2 && !(ev1 in multiCNV) ) { \
          if ( RD_CN1 == 0 && RD_CN1 != RD_CN2 ) makeRevision(id2, RD_CN2 + 2); \
          else if ( RD_CN1 == 1 && RD_CN1 > RD_CN2 ) makeRevision(id2, 1); \
          else if ( RD_CN1 > 1 && RD_CN1 < RD_CN2 ) { newCN = RD_CN2 - RD_CN1 + 2; if ( newCN < 0 ) newCN = 0; makeRevision(id2, newCN); } \
          else makeRevision(id2, 2); } \
        ##Any other time a larger call is driving a smaller call## \
        else if ( support1 ~ /RD/ && smallOverlap50 && length2 > 5000 && !(ev1 in multiCNV) ) { \
          if ( RD_CN1 == 0 ) makeRevision(id2, RD_CN2 + 2); \
          else if ( RD_CN1 == 1 ) makeRevision(id2, RD_CN2 + RD_CN1); \
          else if ( RD_CN1 > 1 ) { newCN = RD_CN2 - RD_CN1 + 2; if ( newCN < 0 ) newCN = 0; makeRevision(id2, newCN); } } } \
     END \
       {for ( id in reviseCN ) print id,reviseCN[id]; }' \
  | sed 's/@/	/' \
  | sort \
  | awk 'BEGIN{FS=OFS="\t"}{if ($3<0) $3=0; print $0}' \
  > $outputfile
