#!/bin/bash
#
# IntegrateGQ_PESR.sh
#

set -euo pipefail

vcf=$1
pegeno_indiv_file=$2
pegeno_variants_file=$3
srgeno_indiv_file=$4
srgeno_variants_file=$5

#Final Output for individual genotype file  VID IID RD-Geno, RD-GQ, PE-Geno, PE-GQ, SR-Geno, SR-GQ, Combined-Geno, Combined-GQ#
##Final Output for variant genotype file  VID  PE-GQ,SR-GQ,Combined-GQ##


svtk vcf2bed $vcf pre.int.bed -i ALGORITHMS -i SVLEN -i EVIDENCE -i SVTYPE

##remove header from int.bed##
tail -n +2 pre.int.bed|awk -F'\t' '{if ($5=$10) print}' OFS='\t' >int.bed

##depth only CNV##
##Make Matrix with everything combined##

##check to make sure PE and SR are same size which they should be##

if [ $(zcat $pegeno_indiv_file|wc -l) !=  $(zcat $srgeno_indiv_file|wc -l)  ] 
then
 echo "ERROR: PE and SR genotype files have different sizes"
 exit
fi

##All PE/SR .'s for samples missing RD##

join -j 1 <(zcat $pegeno_indiv_file|awk '{print $1"@"$2 "\t" $0}'|sort -k1,1) \
   <(zcat $srgeno_indiv_file|awk '{print $1"@"$2 "\t" $3 "\t" $4 "\t" $5}'|sort -k1,1) \
   |tr ' ' '\t' \
   |gzip \
   >PESRall.combined.files.txt.gz

##variant combine##
join -j 1 -o 1.1 1.2 2.2 <(zcat $pegeno_variants_file |sort -k1,1) \
  <(zcat $srgeno_variants_file|sort -k1,1 ) \
  |tr ' ' '\t' \
  |gzip \
  >PESRall.variants.combined.files.txt.gz


##Non-CNV##

zcat PESRall.combined.files.txt.gz \
  | awk '{if ($5==$8 && $6>=$9) print $1,$2,$3,".",".",$5,$6,$8,$9,$5,$6,"PE,SR"; \
  else if ($5==$8 && $6<$9) print $1,$2,$3,".",".",$5,$6,$8,$9,$8,$9,"PE,SR"; \
  else if ($5>0 && $8==0) print $1,$2,$3,".",".",$5,$6,$8,$9,$5,$6,"PE"; \
  else if ($5==0 && $8>0) print $1,$2,$3,".",".",$5,$6,$8,$9,$8,$9,"SR"; \
  else if ($6>=$9)  print $1,$2,$3,".",".",$5,$6,$8,$9,$5,$6,"PE,SR"; \
  else if ($6<$9) print $1,$2,$3,".",".",$5,$6,$8,$9,$8,$9,"PE,SR"}' \
  |tr ' ' '\t' \
  >genotype.indiv.txt

 zcat PESRall.variants.combined.files.txt.gz \
  |awk '{if ($2>=$3)print $1,".",$2,$3,$2 ;else print $1,".",$2,$3,$3 }' \
  |tr ' ' '\t' \
  >genotype.variant.txt


##output- integrated genotypes##


##if you want to cap DEL at 2 for multiallelics but not everyone else##
#if [ $(awk '{if ($5=="DEL") print $4}' int.bed|wc -l) -gt 0 ]
#then

#awk '{if ($5=="DEL") print $4}' int.bed| \
#  fgrep -wf - genotype.indiv.txt| \
#  awk '{if ($10>2) print $1,$2,$3,$4,$5,$6,$7,$8,$9,2,$11,$12;else print}'\
#  >genotype.indiv.del.txt

#awk '{if ($5=="DEL") print $4}' int.bed| \
#  fgrep -wvf - genotype.indiv.txt \
#  >genotype.indiv.nodel.txt
#cat genotype.indiv.del.txt genotype.indiv.nodel.txt|tr ' ' '\t'|gzip>genotype.indiv.txt.gz

#else
#gzip genotype.indiv.txt
#fi

##any quality score of 0 moved to 1 since assigend genotype is best option score should not be 0###

cat genotype.indiv.txt \
  |awk -F'\t' '{if ($10>2) print $1,$2,$3,$4,$5,$6,$7,$8,$9,2,$11,$12;else print}' OFS="\t" \
  |awk '{if ($5==0) $5=1;if ($7==0) $7=1; if ($9==0) $9=1; if ($11==0) $11=1; print}' OFS="\t" \
  |gzip>genotype.indiv.txt.gz

cat genotype.variant.txt \
  |awk '{if ($2==0) $2=1;if ($3==0) $3=1; if ($4==0) $4=1; if ($5==0) $5=1; print}' OFS="\t" \
  |gzip \
  >genotype.variant.txt.gz
