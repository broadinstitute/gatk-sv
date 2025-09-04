#!/bin/bash
#
# IntegrateGQ_depthonly.sh
#

set -euo pipefail

vcf=$1
RD_melted_genotypes=$2 
RD_melted_variants_genotypes=$3

#Final Output for individual genotype file  VID IID RD-Geno, RD-GQ, PE-Geno, PE-GQ, SR-Geno, SR-GQ, Combined-Geno, Combined-GQ#
##Final Output for variant genotype file  VID  PE-GQ,SR-GQ,Combined-GQ##


svtk vcf2bed $vcf pre.int.bed -i ALGORITHMS -i SVLEN -i EVIDENCE -i SVTYPE

##remove header from int.bed##
tail -n +2 pre.int.bed|awk -F'\t' '{if($5=$10); print}' OFS='\t' >int.bed

##depth only CNV##
##Make Matrix with everything combined##
##RD##
##modify to match PE/SR format##
zcat $RD_melted_genotypes \
  |awk '{print $4"@"$5,$4,$5,$6,$7}' OFS='\t' \
  |sort -k1,1 \
  |gzip \
  >rd_indiv_geno.txt.gz
    

##Depth only- Just RD support used for genotype##

awk -F'\t' '{if ($5=="DEL") print $4}' int.bed>depthonly.del.ids.txt
awk -F'\t' '{if ($5=="DUP") print $4}' int.bed>depthonly.dup.ids.txt

# In case of empty inputs
touch genotype.indiv.txt
touch genotype.variant.txt

if [ $(cat depthonly.del.ids.txt|wc -l) -gt 0 ] && [ $(zcat rd_indiv_geno.txt.gz|wc -l) -gt 0 ]
then
 zcat rd_indiv_geno.txt.gz \
  |fgrep -wf depthonly.del.ids.txt \
  |awk '{if ($4>=2) print $0 "\t" 0 "\t" $5"\t" "RD"; \
   else if ($4==1) print $0 "\t" 1 "\t" $5"\t" "RD"; \
   else if ($4==0) print $0 "\t" 2 "\t" $5"\t" "RD" }' \
  |awk '{print $1,$2,$3,$4,$5,".",".",".",".", $6,$7,$8}' \
  |tr ' ' '\t' \
  >genotype.indiv.txt

 cat $RD_melted_variants_genotypes \
  |cut -f4- \
  |fgrep -wf depthonly.del.ids.txt \
  |awk '{print $0,".",".",$2}'\
  |tr ' ' '\t' \
  >genotype.variant.txt
fi

if [ $(cat depthonly.dup.ids.txt|wc -l) -gt 0 ] && [ $(zcat rd_indiv_geno.txt.gz|wc -l) -gt 0 ]
then
 zcat rd_indiv_geno.txt.gz  \
  |fgrep -wf depthonly.dup.ids.txt \
  |awk '{if ($4<=2) print $0 "\t" 0 "\t" $5"\t" "RD"; \
  else if ($4-2>4) print $0 "\t" 2 "\t" $5"\t" "RD" ; \
  else print $0 "\t" $4-2 "\t" $5"\t" "RD"}' \
  |awk '{print $1,$2,$3,$4,$5,".",".",".",".", $6,$7,$8}' \
  |tr ' ' '\t' \
  >>genotype.indiv.txt
 
cat $RD_melted_variants_genotypes |cut -f4- \
  |fgrep -wf depthonly.dup.ids.txt \
  |awk '{print $0,".",".",$2}'\
  |tr ' ' '\t' \
  >>genotype.variant.txt
fi

##any quality score of 0 moved to 1 since assigend genotype is best option score should not be 0###

cat genotype.indiv.txt \
  |awk '{if ($5==0) $5=1;if ($7==0) $7=1; if ($9==0) $9=1; if ($11==0) $11=1; print}' OFS="\t" \
  |gzip \
  >genotype.indiv.depth.txt.gz 

cat genotype.variant.txt \
  |awk '{if ($2==0) $2=1;if ($3==0) $3=1; if ($4==0) $4=1; if ($5==0) $5=1; print}' OFS="\t" \
  |gzip \
  >genotype.variant.depth.txt.gz
