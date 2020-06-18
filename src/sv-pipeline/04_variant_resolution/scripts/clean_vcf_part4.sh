#!/bin/bash
#
# clean_VCF_part4.sh
#
# 
#
# Copyright (C) 2018 Harrison Brand<hbrand1@mgh.harvard.edu>
# Distributed under terms of the MIT license.

set -euxo pipefail

##gzipped combined bed file##
##combined output file from clean_vcf_part2.sh##
RD_CN_revise_forgeno=$1
normal_revise_vcf=$2


##seed the vcf lines file which will provide the revisions to vcf file## 
echo "">revise.vcf.lines.txt


##reduce vcf to lines for given shard##
 cat <(zcat $normal_revise_vcf|sed -n '1,1000p' |egrep ^# ) \
  <(zcat $normal_revise_vcf |fgrep -wf <(awk '{print $1}' $RD_CN_revise_forgeno|sort -u)) \
  |bgzip \
  >int.vcf.gz || true


##get column ids##
zcat $normal_revise_vcf \
  |sed -n '1,1000p' \
  |egrep ^# \
  |tail -n 1 \
  |tr '\t' '\n' \
  |cat -n - \
  >col.txt
    

##pull out and revise vcf line that needs to be edited##
while read line
do
 id=$(echo $line|awk '{print $2}' )
 col=$(awk -v id=$id '{if($2==id) print $1}' col.txt)
 variant=$(echo $line|awk '{print $1}')
 cn=$(echo $line|awk '{print $3}')

 zcat int.vcf.gz \
   |{ fgrep -w $variant || true; } \
    >line.txt

 echo $variant $id
 ##Updated genotype and rebuild Format field ##
 GT=$(awk -v col=$col  '{print $col}' line.txt|awk -F":" '{print $1}')
 GQ=$(awk -v col=$col  '{print $col}' line.txt|awk -F":" '{print $2}')
 RD_CN=$(awk -v col=$col  '{print $col}' line.txt|awk -F":" '{print $3}')
 RD_GQ=$(awk -v col=$col  '{print $col}' line.txt|awk -F":" '{print $4}')
 PE_GT=$(awk -v col=$col  '{print $col}' line.txt|awk -F":" '{print $5}')
 PE_GQ=$(awk -v col=$col  '{print $col}' line.txt|awk -F":" '{print $6}')
 SR_GT=$(awk -v col=$col  '{print $col}' line.txt|awk -F":" '{print $7}')
 SR_GQ=$(awk -v col=$col  '{print $col}' line.txt|awk -F":" '{print $8}')
 EV=$(awk -v col=$col  '{print $col}' line.txt|awk -F":" '{print $9}')

 if [ $(cat revise.vcf.lines.txt|fgrep -w $variant|wc -l) -gt 0 ]
 then
  cat revise.vcf.lines.txt \
   |awk -v col=$col -v var=$variant -v GT=$GT -v GQ=$GQ -v RD_CN=$cn -v RD_GQ=$RD_GQ -v PE_GT=$PE_GT -v PE_GQ=$PE_GQ -v SR_GT=$SR_GT -v SR_GQ=$SR_GT -v EV=$EV '{if ($3==var ) $col="0/1:"GQ":"RD_CN":"RD_GQ":"PE_GT":"PE_GQ":"SR_GT":"SR_GQ":"EV ;print}' \
   >int.lines.txt

  cat int.lines.txt > revise.vcf.lines.txt

 else 
  cat line.txt \
  |awk -v col=$col -v var=$variant -v GT=$GT -v GQ=$GQ -v RD_CN=$cn -v RD_GQ=$RD_GQ -v PE_GT=$PE_GT -v PE_GQ=$PE_GQ -v SR_GT=$SR_GT -v SR_GQ=$SR_GT -v EV=$EV '{if ($3==var ) $col="0/1:"GQ":"RD_CN":"RD_GQ":"PE_GT":"PE_GQ":"SR_GT":"SR_GQ":"EV ;print}' \
  >>revise.vcf.lines.txt
  fi

done<$RD_CN_revise_forgeno


bgzip revise.vcf.lines.txt


##get multilallelic genotypes##
##pull out lines for normal vcf for given batch##
total_lines=$(zcat $normal_revise_vcf|egrep -v "^#"|wc -l)
batch=$(ls $RD_CN_revise_forgeno|awk -F'/' '{print $NF}'|awk -F'[._]' '{print $2}'|awk '{if ($1==0) print 1; else print}')
total_batch=$(ls $RD_CN_revise_forgeno|awk -F'/' '{print $NF}'|awk -F'[._]' '{print $3}'|awk '{if ($1==0) print 1; else print}')

segments=$(echo $total_batch $total_lines|awk '{print $2/$1}')

 cat <(zcat $normal_revise_vcf|sed -n '1,1000p' |egrep ^# ) \
  <(zcat $normal_revise_vcf |egrep -v "^#"|awk -v batch=$batch -v segments=$segments '{if (NR<=batch*segments && NR>=((batch-1)*segments) ) print }') \
  |bgzip \
  >split.vcf.gz


for var in PE_GT SR_GT PE_GQ SR_GQ
do
 zcat split.vcf.gz\
  |awk -F'\t' '{if ($1!~"#") $1=$3;print}' OFS="\t" \
  |vcftools --vcf - --stdout --extract-FORMAT-info ${var} \
  |awk -F"\t" 'NR==1{for (i=3;i<=NF;i++) header[i]=$i} NR>1{for(j=3;j<=NF;j++) print $1"@"header[j] "\t" $j }' \
  |sort -k1,1 \
  |gzip \
  >multicheck.${var}.FORMAT.gz
done

##concatenate metrics##
# get each line formatted as SITE@SAMPLE PE_GT PE_GQ SR_GT SR_GQ
join -j 1 <(zcat multicheck.PE_GT.FORMAT.gz) \
  <(zcat multicheck.PE_GQ.FORMAT.gz) \
  |join -j 1  - <(zcat multicheck.SR_GT.FORMAT.gz) \
  |join -j 1  - <(zcat multicheck.SR_GQ.FORMAT.gz) \
  |tr ' ' '\t' \
  |gzip \
  >multi.combined.format.gz

# Set the maximum allowable number of samples with a PE or SR GT > 3 to be 1% or 2, whichever is greater
vf_1=$(zcat split.vcf.gz  |sed -n '1,1000p' |egrep -v "^##"|cut -f10-|awk 'NR==1{print (NF) * 0.01}' |awk '{if ($1 <= 2) {$1 = 2}; print $1}')

# Choose the best of PE and SR genotypes for each site / sample
# Count the number of samples with a GT over 3 for each site
# Add site IDs with sample counts over $vf_1 to the multi.geno.ids.txt.gz file
zcat multi.combined.format.gz \
   |awk '{if ($2>0 && $4==0) print $1"\t" $2; \
   else if ($2==0) print $1 "\t" $4; \
   else if ($3>=$5)print $1"\t" $2; \
   else print $1"\t" $4 }' \
   |tr '@' '\t' \
   |awk '{if ($3>2 && $2!=".") print $1}' \
   |sort \
   |uniq -c \
   |awk -v vf_1=$vf_1 '{if ($1>vf_1)print $2}' \
   |gzip \
   >multi.geno.ids.txt.gz 

