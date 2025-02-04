#!/bin/bash
#
# IntegrateGQ.sh
#

set -euo pipefail

vcf=$1
RD_melted_genotypes=$2 
RD_melted_variants_gentoypes=$3
pegeno_indiv_file=$4
pegeno_variants_file=$5
srgeno_indiv_file=$6
srgeno_variants_file=$7

#Final Output for individual genotype file  VID IID RD-Geno, RD-GQ, PE-Geno, PE-GQ, SR-Geno, SR-GQ, Combined-Geno, Combined-GQ#
##Final Output for variant genotype file  VID  PE-GQ,SR-GQ,Combined-GQ##

svtk vcf2bed $vcf pre.int.bed -i ALGORITHMS -i SVLEN -i EVIDENCE -i SVTYPE

##remove header from int.bed##
tail -n +2 pre.int.bed|awk -F'\t' '{gsub($5,$10); print}'>int.bed

##depth only CNV##
##Make Matrix with everything combined##
##RD##
zcat $RD_melted_genotypes \
  |awk '{print $4"@"$5,$4,$5,$6,$7}' OFS='\t' \
  |sort -k1,1 \
  |gzip \
  >rd_indiv_geno.txt.gz

##PE##
zcat "$pegeno_indiv_file" | \
awk -v OFS="\t" '
  ARGIND==1 {
    if ($5 == "DEL") {
      del[$4]
    } else {
      nodel[$4]
    }
    next
  }
  ARGIND==2 {
    if ($1 in del) {
      ##Deletions, need to PE-SR genotypes to match RD format (2==ref)##
      final_gt = ($4>1 ? 0 : ($4==1 ? 1 : 2))
      print $1"@"$2, $1, $2, $4, final_gt, $5
    } else if ($1 in nodel) {
      ##Duplications and other events, need to PE-SR genotypes to match RD (2==ref)##
      final_gt = ($4>0 ? $4+2 : 2)
      print $1"@"$2, $1, $2, $4, final_gt, $5
    }
  }
' int.bed - \
| awk '!seen[$1]++' \
| sort -k1,1 \
| gzip > pe_indiv_geno.txt.gz

##SR##
zcat "$srgeno_indiv_file" | \
awk -v OFS="\t" '
  ARGIND==1 {
    if ($5 == "DEL") {
      del[$4]
    } else {
      nodel[$4]
    }
    next
  }
  ARGIND==2 {
    if ($1 in del) {
      ##Deletions, need to PE-SR genotypes to match RD format (2==ref)##
      final_gt = ($4>1 ? 0 : ($4==1 ? 1 : 2))
      print $1"@"$2, $1, $2, $4, final_gt, $5
    } else if ($1 in nodel) {
      ##Duplications and other events, need to PE-SR genotysrs to match RD (2==ref)##
      final_gt = ($4>0 ? $4+2 : 2)
      print $1"@"$2, $1, $2, $4, final_gt, $5
    }
  }
' int.bed - \
| awk '!seen[$1]++' \
| sort -k1,1 \
| gzip > sr_indiv_geno.txt.gz


##check to make sure PE and SR are same size which they should be##

if [ $(zcat sr_indiv_geno.txt.gz|wc -l) !=  $(zcat pe_indiv_geno.txt.gz|wc -l)  ] 
then
 echo "ERROR: PE and SR genotype files have different sizes"
 exit
fi

##All PE/SR .'s for samples missing RD##
##Add max PE/SR to use for future GQ integration##
join -j 1 -a 1 -e "." -o 1.1 1.2 1.3 1.4 1.5 2.4 2.5 2.6  \
  <(zcat rd_indiv_geno.txt.gz) <(zcat pe_indiv_geno.txt.gz) \
  |tr ' ' '\t' \
  |join -j 1 -e "." -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 2.4 2.5 2.6 - <(zcat sr_indiv_geno.txt.gz) \
  |tr ' ' '\t' \
  | awk '{if ($6==$9 && $8>=$11) print $0 "\t" $6 "\t" $7 "\t" $8 ; \
  else if ($6==$9 && $8<$11) print $0 "\t" $9 "\t" $10 "\t" $11 ; \
  else if ($6>0 && $9==0) print $0 "\t" $6 "\t" $7 "\t" $8  ; \
  else if ($6==0 && $9>0) print $0 "\t" $9 "\t" $10 "\t" $11 ; \
  else if ($8>=$11)  print $0 "\t"$6 "\t" $7 "\t" $8  ; \
  else if ($11>$8) print $0 "\t" $9 "\t" $10 "\t" $11 }' \
  |gzip \
  >RDall.combined.files.txt.gz

##variant combine##
join -j 1 -a 1 -e "." -o 1.1 1.2 2.2 <(cut -f4- $RD_melted_variants_gentoypes|fgrep -v variant_gq|sort -k1,1 ) \
  <(zcat $pegeno_variants_file|sort -k1,1 ) \
  |join -j 1 -a 1 -e "." -o 1.1 1.2 1.3 2.2 - <(zcat $srgeno_variants_file|sort -k1,1) \
  |tr ' ' '\t' \
  |gzip \
  >RDall.variants.combined.files.txt.gz

##All RD NA's for samples missing RD##
join -j 1 -a 2 -e "." -o 2.1 2.2 2.3 1.4 1.5 2.4 2.6  <(zcat rd_indiv_geno.txt.gz|sort -k1,1) \
   <(zcat pe_indiv_geno.txt.gz|sort -k1,1) \
   |tr ' ' '\t' \
   |join -j 1 -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 2.4 2.6 - <(zcat sr_indiv_geno.txt.gz|sort -k1,1) \
   |tr ' ' '\t' \
   |gzip \
   >PESRall.combined.files.txt.gz

##variant combine##
join -j 1 -a 2 -e "." -o 2.1 1.2 2.2 <(cut -f4- $RD_melted_variants_gentoypes|fgrep -v variant_gq|sort -k1,1) \
  <(zcat $pegeno_variants_file|sort -k1,1 ) \
  |join -j 1 -a 2 -e "." -o 1.1 1.2 1.3 2.2 - <(zcat $srgeno_variants_file|sort -k1,1 ) \
  |tr ' ' '\t' \
  |gzip \
  >PESRall.variants.combined.files.txt.gz


##Depth only- Just RD support used for genotype##

awk -F'\t' '{if ($8=="depth" && $5=="DEL") print $4}' int.bed>depthonly.del.ids.txt
awk -F'\t' '{if ($8=="depth" && $5=="DUP") print $4}' int.bed>depthonly.dup.ids.txt

touch genotype.variant.txt
touch genotype.indiv.txt

if [ $(cat depthonly.del.ids.txt|wc -l) -gt 0 ] 
then
 zcat RDall.combined.files.txt.gz \
  | { fgrep -wf depthonly.del.ids.txt || true; } \
  |awk '{if ($4>=2) print $0 "\t" 0 "\t" $5"\t" "RD"; \
   else if ($4==1) print $0 "\t" 1 "\t" $5"\t" "RD"; \
   else if ($4==0) print $0 "\t" 2 "\t" $5"\t" "RD" }' \
  |cut -f1-6,8-9,11,15- \
  |tr ' ' '\t' \
  >genotype.indiv.txt

 zcat RDall.variants.combined.files.txt.gz \
  | { fgrep -wf depthonly.del.ids.txt || true; }  \
  |awk '{print $0 "\t" $2}'\
  >genotype.variant.txt
fi

if [ $(cat depthonly.dup.ids.txt|wc -l) -gt 0 ]
then
 zcat RDall.combined.files.txt.gz \
  | { fgrep -wf depthonly.dup.ids.txt || true; }  \
  |awk '{if ($4<=2) print $0 "\t" 0 "\t" $5"\t" "RD"; \
  else print $0 "\t" $4-2 "\t" $5"\t" "RD"}' \
  |cut -f1-6,8-9,11,15- \
  |tr ' ' '\t' \
  >>genotype.indiv.txt

 zcat RDall.variants.combined.files.txt.gz \
  | { fgrep -wf depthonly.dup.ids.txt  || true; } \
  |awk '{print $0 "\t" $2}'\
  >>genotype.variant.txt
fi

##Non-CNV##
awk -F'\t' '{if ($5!="DUP" && $5!="DEL") print $4}' int.bed>noncnv.ids.txt
##include any CNV with no depth genotype because of no coverage##
zcat RDall.combined.files.txt.gz \
  |awk '{if ($4==".") print $2}' \
  >>noncnv.ids.txt
if [ $(cat noncnv.ids.txt|wc -l) -gt 0 ]
then
 ##Pull out highest GQ if an event has a genotype and lowest if it doesn't##
 #take max of PE/SR for genotype when they have a genotype## \
 zcat PESRall.combined.files.txt.gz \
  | { fgrep -wf noncnv.ids.txt || true; }  \
  |awk '{if ($6==$8 && $7>=$9) print $0 "\t" $6 "\t" $7 "\t" "PE,SR"; \
  else if ($6==$8 && $7<$9) print $0 "\t" $8 "\t" $9 "\t" "PE,SR"; \
  else if ($6>0 && $8==0) print $0 "\t" $6 "\t" $7 "\t" "PE"; \
  else if ($6==0 && $8>0) print $0 "\t" $8 "\t" $9 "\t" "SR"; \
  else if ($7>=$9)  print $0 "\t" $6 "\t" $7 "\t" "PE,SR" ; \
  else if ($9>$7) print $0 "\t" $8 "\t" $9 "\t" "PE,SR"}' \
  >>genotype.indiv.txt

 zcat PESRall.variants.combined.files.txt.gz \
  | { fgrep -wf noncnv.ids.txt || true; }  \
  |awk '{if ($3>=$4)print $0 "\t" $3 ;else print $0 "\t" $4 }' \
  |tr ' ' '\t' \
  >>genotype.variant.txt
fi

rm PESRall.combined.files.txt.gz

##CNV##
##Recode RD to match PE/SR##
##CNV >5kb and removing any CNV with no depth genotype###
{ fgrep -wvf <(zcat RDall.combined.files.txt.gz  \
    |awk -F'\t' '{if ($4==".") print $2}') int.bed || true; } \
 |awk '{if (($5=="DEL") && $3-$2>=5000 ) print $4}' \
 >gt5kbcnv.del.ids.txt

{ fgrep -wvf <(zcat RDall.combined.files.txt.gz \
   |awk -F'\t' '{if ($4==".") print $2}') int.bed || true; } \
 |awk -F'\t' '{if (($5=="DUP" ) && $3-$2>=5000 ) print $4}'\
 >gt5kbcnv.dup.ids.txt


if [ $(cat gt5kbcnv.del.ids.txt|wc -l) -gt 0 ] 
then

 zcat RDall.combined.files.txt.gz \
  | { fgrep -wf gt5kbcnv.del.ids.txt || true; } \
  |awk  '{if ($5==999) print $0,$4,999; \
   else if ($4==$13) print $0,$4,$5+((999-$5) * ($14/999)*0.5); \
   else print $0,$4,$5-(($5-1) * ($14/999)*0.5) }' OFS='\t' \
  |awk '{if ($4==$7 && $4==$10) print $0 "\t" "RD,PE,SR"; \
  else if ($4==2 && $4==$7 && $4!=$10) print $0 "\t" "RD,PE" ; \
  else if ($4==2 && $4!=$7 && $4==$10) print $0 "\t" "RD,SR" ; \
  else if ($4==2 && $4!=$7 && $4!=$10) print $0 "\t" "RD" ; \
  else if ($4!=2 && $7!=2 && $10!=2) print $0 "\t" "RD,PE,SR" ; \
  else if ($4!=2 && $7!=2 && $10==2) print $0 "\t" "RD,PE" ; \
  else if ($4!=2 && $7==2 && $10!=2) print $0 "\t" "RD,SR" ; \
  else if ($4!=2 && $7==2 && $10==2) print $0 "\t" "RD" }' \
  |cut -f1-6,8-9,11,15- \
  |awk '{if ($4>=2) print $1,$2,$3,$4,$5,$6,$7,$8,$9,0,$11,$12; \
   else if ($4==1) print $1,$2,$3,$4,$5,$6,$7,$8,$9,1,$11,$12; \
   else if ($4==0) print $1,$2,$3,$4,$5,$6,$7,$8,$9,2,$11,$12; }' \
  |tr ' ' '\t' \
  >>genotype.indiv.txt

 zcat PESRall.variants.combined.files.txt.gz \
  | { fgrep -wf gt5kbcnv.del.ids.txt || true; } \
  |awk  '{if ($2==999) print $0,999; \
   else if ($3>=$4) print $0,$2+((999-$2) * ($3/999)*0.5); \
   else if ($4>$3) print $0,$2+((999-$2) * ($4/999)*0.5);  }' \
   |tr ' ' '\t' \
  >>genotype.variant.txt
fi


if [ $(cat gt5kbcnv.dup.ids.txt|wc -l) -gt 0 ] 
then

 zcat RDall.combined.files.txt.gz \
  | { fgrep -wf gt5kbcnv.dup.ids.txt || true; } \
  |awk  '{if ($5==999) print $0,$4,999; \
   else if ($4==$13) print $0,$4,$5+((999-$5) * ($14/999)*0.5); \
   else print $0,$4,$5-(($5-1) * ($14/999)*0.5) }' OFS='\t' \
  |awk '{if ($4==$7 && $4==$10) print $0 "\t" "RD,PE,SR"; \
  else if ($4==2 && $4==$7 && $4!=$10) print $0 "\t" "RD,PE" ; \
  else if ($4==2 && $4!=$7 && $4==$10) print $0 "\t" "RD,SR" ; \
  else if ($4==2 && $4!=$7 && $4!=$10) print $0 "\t" "RD" ; \
  else if ($4!=2 && $7!=2 && $10!=2) print $0 "\t" "RD,PE,SR" ; \
  else if ($4!=2 && $7!=2 && $10==2) print $0 "\t" "RD,PE" ; \
  else if ($4!=2 && $7==2 && $10!=2) print $0 "\t" "RD,SR" ; \
  else if ($4!=2 && $7==2 && $10==2) print $0 "\t" "RD" }' \
  |cut -f1-6,8-9,11,15- \
  |awk '{if ($4<=2) print $1,$2,$3,$4,$5,$6,$7,$8,$9,0,$11,$12; \
   else if ($4>2) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10-2,$11,$12}' \
  |tr ' ' '\t' \
  >>genotype.indiv.txt

 zcat PESRall.variants.combined.files.txt.gz \
  | { fgrep -wf gt5kbcnv.dup.ids.txt || true; } \
  |awk  '{if ($2==999) print $0,999; \
   else if ($3>=$4) print $0,$2+((999-$2) * ($3/999)*0.5); \
   else if ($4>$3) print $0,$2+((999-$2) * ($4/999)*0.5);  }' \
   |tr ' ' '\t' \
  >>genotype.variant.txt
fi

##CNV 1-5kb and removing any CNV with no depth genotype###
{ fgrep -wvf <(zcat RDall.combined.files.txt.gz \
     |awk -F'\t' '{if ($4==".") print $2}') int.bed || true; } \
  |awk -F'\t' '{if (($5=="DUP" || $5=="DEL") && $3-$2<5000 && $3-$2>=1000 ) print $4}' \
  >gt1_5kbcnv.ids.txt


if [ $(cat gt1_5kbcnv.ids.txt|wc -l) -gt 0 ]
then

 zcat RDall.combined.files.txt.gz \
  | { fgrep -wf gt1_5kbcnv.ids.txt || true; } \
  |awk  '{if ($14==999) print $0,$12,999; \
  else if ($4==$13) print $0,$12,$14+((999-$14) * ($5/999)*0.5); \
  else print $0,$12,$14-(($14-1) * ($5/999)*0.5) }' OFS='\t' \
  |awk '{if ($4==$7 && $4==$10) print $0 "\t" "RD,PE,SR"; \
  else if ($13==2 && $4==2 && $7==2 && $10!=2) print $0 "\t" "RD,PE" ; \
  else if ($13==2 && $4==2 && $7!=2 && $10==2) print $0 "\t" "RD,SR" ; \
  else if ($13==2 && $4!=2 && $7==2 && $10!=2) print $0 "\t" "PE" ; \
  else if ($13==2 && $4!=2 && $7!=2 && $10==2) print $0 "\t" "SR" ; \
  else if ($13==2 && $4!=2 && $7==2 && $10==2) print $0 "\t" "PE,SR" ; \
  else if ($13!=2 && $4!=2 && $7!=2 && $10!=2) print $0 "\t" "RD,PE,SR" ; \
  else if ($13!=2 && $4!=2 && $7!=2 && $10==2) print $0 "\t" "RD,PE" ; \
  else if ($13!=2 && $4!=2 && $7==2 && $10!=2) print $0 "\t" "RD,SR" ; \
  else if ($13!=2 && $4==2 && $7==2 && $10!=2) print $0 "\t" "SR" ; \
  else if ($13!=2 && $4==2 && $7!=2 && $10==2) print $0 "\t" "PE"; \
  else if ($13!=2 && $4==2 && $7!=2 && $10!=2) print $0 "\t" "PE,SR" }' \
  |cut -f1-6,8-9,11,15- \
  |tr ' ' '\t' \
  >>genotype.indiv.txt

 zcat PESRall.variants.combined.files.txt.gz \
  | { fgrep -wf gt1_5kbcnv.ids.txt || true; } \
  |awk  '{if ($3==999 || $4==999) print $0,999; \
   else if ($3>=$4) print $0,$3+((999-$3) * ($2/999)*0.5); \
   else if ($4>$3) print $0,$4+((999-$4) * ($2/999)*0.5);  }' \
   |tr ' ' '\t' \
  >>genotype.variant.txt
fi

###CNV <1kb and removing any CNV with no depth genotype####
{ fgrep -wvf <(zcat RDall.combined.files.txt.gz \
     |awk -F'\t' '{if ($4==".") print $2}') int.bed || true; } \
  |awk -F'\t' '{if (($5=="DUP" || $5=="DEL") && $3-$2<1000 ) print $4}' \
  >lt1kbcnv.ids.txt

if [ $(cat lt1kbcnv.ids.txt|wc -l) -gt 0 ] 
then

 zcat RDall.combined.files.txt.gz \
  | { fgrep -wf lt1kbcnv.ids.txt || true; } \
  |awk  '{if ($14==999) print $0,$12,999; \
  else if ($7==$10 && $8>=$11) print $0,$12,$14+((999-$14) * ($11/999)*0.5); \
  else if ($7==$10 && $11>=$8) print $0,$12,$14+((999-$14) * ($8/999)*0.5); \
  else  print $0 "\t" $12 "\t" $14}' OFS='\t' \
  |awk '{if ($4==$7 && $4==$10) print $0 "\t" "RD,PE,SR"; \
  else if ($13==2 && $4==2 && $7==2 && $10!=2) print $0 "\t" "RD,PE" ; \
  else if ($13==2 && $4==2 && $7!=2 && $10==2) print $0 "\t" "RD,SR" ; \
  else if ($13==2 && $4!=2 && $7==2 && $10!=2) print $0 "\t" "PE" ; \
  else if ($13==2 && $4!=2 && $7!=2 && $10==2) print $0 "\t" "SR" ; \
  else if ($13==2 && $4!=2 && $7==2 && $10==2) print $0 "\t" "PE,SR" ; \
  else if ($13!=2 && $4!=2 && $7!=2 && $10!=2) print $0 "\t" "RD,PE,SR" ; \
  else if ($13!=2 && $4==$13 && $7!=2 && $10==2) print $0 "\t" "RD,PE" ; \
  else if ($13!=2 && $4==$13 && $7==2 && $10!=2) print $0 "\t" "RD,SR" ; \
  else if ($13!=2  && $7==2 && $10!=2) print $0 "\t" "SR" ; \
  else if ($13!=2  && $7!=2 && $10==2) print $0 "\t" "PE"; \
  else if ($13!=2  && $7!=2 && $10!=2) print $0 "\t" "PE,SR" }' \
  |cut -f1-6,8-9,11,15- \
  |tr ' ' '\t' \
  >>genotype.indiv.txt

 zcat PESRall.variants.combined.files.txt.gz \
  | { fgrep -wf lt1kbcnv.ids.txt || true; } \
  |awk  '{if ($3==999 || $4==999) print $0,999; \
   else if ($3>=$4) print $0,$3+((999-$3) * ($2/999)*0.5); \
   else if ($4>$3) print $0,$4+((999-$4) * ($2/999)*0.5);  }' \
   |tr ' ' '\t' \
  >>genotype.variant.txt
fi

rm RDall.combined.files.txt.gz

##output- integrated genotypes##

#If Del and other should be treated seperate## 
#awk -F'\t' '{if ($5=="DEL") print $4}' int.bed| \
#  fgrep -wf - genotype.indiv.txt| \
#  awk -F'\t' '{if ($10>2) print $1,$2,$3,$4,$5,$6,$7,$8,$9,2,$11,$12;else print}'\
#  |tr ' ' '\t' \
#  >genotype.indiv.del.txt

#currently DEL and other treated same but past verisions seperate#
#awk -F'\t' '{if ($5=="DEL") print $4}' int.bed| \
#  fgrep -wvf - genotype.indiv.txt \
#  |awk -F'\t' '{if ($10>2) print $1,$2,$3,$4,$5,$6,$7,$8,$9,2,$11,$12;else print}'\
#  |tr ' ' '\t' \
#  >genotype.indiv.nodel.txt

#rm genotype.indiv.txt

#cat genotype.indiv.del.txt genotype.indiv.nodel.txt|gzip>genotype.indiv.txt.gz 

##any quality score of 0 moved to 1 since assigend genotype is best option score should not be 0###

cat genotype.indiv.txt \
  |awk -F'\t' '{if ($10>2) print $1,$2,$3,$4,$5,$6,$7,$8,$9,2,$11,$12;else print}' OFS="\t" \
  |awk '{if ($5==0) $5=1;if ($7==0) $7=1; if ($9==0) $9=1; if ($11==0) $11=1; print}' OFS="\t" \
  |gzip>genotype.indiv.txt.gz

cat genotype.variant.txt \
  |awk '{if ($2==0) $2=1;if ($3==0) $3=1; if ($4==0) $4=1; if ($5==0) $5=1; print}' OFS="\t" \
  |gzip \
  >genotype.variant.txt.gz
