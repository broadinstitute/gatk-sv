#!/bin/bash
#
# SR_genotype.sh
#
# 
#
# Copyright (C) 2018 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.


set -euo pipefail

#whole cohort files that can be split
vcf=$1
SR_counts=$2
SR_sum=$3
metric_file=$4

#sr cutoff taken from metric file for genotyping
sr_count=$(awk '{if($1=="sr_count") print $2}' $metric_file)
median_hom=$(awk '{if($1=="median_hom") print $2}' $metric_file)
sd_het=$(awk '{if($1=="sd_het") print $2}' $metric_file)
rare_min=$(awk '{if($1=="rare_min") print $2}' $metric_file)
rare_max=$(awk '{if($1=="rare_max") print $2}' $metric_file)
common_min=$(awk '{if($1=="common_min") print $2}' $metric_file)
common_max=$(awk '{if($1=="common_max") print $2}' $metric_file)
rare_single=$(awk '{if($1=="rare_single") print $2}' $metric_file)
rare_both=$(awk '{if($1=="rare_both") print $2}' $metric_file)
common_single=$(awk '{if($1=="common_single") print $2}' $metric_file)
common_both=$(awk '{if($1=="common_both") print $2}' $metric_file)


zcat ${SR_counts} \
  | awk -v sr_count=$sr_count '{if ($NF>(sr_count/2)) print $1"@"$3}' \
  | sort \
  | uniq -c \
  | awk '{if ($1==2) print $2}' \
  > two.sided.pass.txt

# grep out variants which pass SR in random forest
svtk vcf2bed $vcf int.bed -i EVIDENCE

awk '{if ($NF~"SR") print $4}' int.bed> pass.srtest.txt

##Genotype SR genotype (0-ref, then estimate copy state based on copy state that is 1 sd from sd_het  )##
if [ $(cat two.sided.pass.txt|wc -l) -gt 0 ]
then
zcat ${SR_sum} \
  | awk '{print $0 "\t" $1"@"$2}' \
  | { fgrep -wf two.sided.pass.txt || [[ $? == 1 ]]; } \
  | cut -f1-3 \
  | awk -v var=$sr_count -v var1=$median_hom -v var2=$sd_het '{if ($3<var) print $1,$2,$3,0;else if ($3<=var1-var2) print $1,$2,$3,1; else print $1,$2,$3,int($3/(var1/2)+0.5)}'  \
  |tr ' ' '\t' \
  > sr.geno.final.txt
fi

zcat ${SR_sum} \
  | awk '{print $0 "\t" $1"@"$2}' \
  | { fgrep -wvf two.sided.pass.txt || [[ $? == 1 ]]; } \
  | cut -f1-3 \
  | awk '{print $1,$2,$3,0}' \
  |tr ' ' '\t' \
  >> sr.geno.final.txt

gzip sr.geno.final.txt

zcat ${SR_sum} \
  | awk '{print $0 "\t" $1"@"$2}' \
  | cut -f1-3 \
  | awk -v var=$sr_count -v var1=$median_hom -v var2=$sd_het '{if ($3<var) print $1,$2,$3,0;else if ($3<=var1-var2) print $1,$2,$3,1; else print $1,$2,$3,int($3/(var1/2)+0.5)}' \
  |tr ' ' '\t' \
  | gzip \
  > sr.geno.final.oneside.txt.gz
  
##Allow just one side##
zcat sr.geno.final.oneside.txt.gz \
  |awk '{if ($3>1) print $1}' \
  |sort|uniq -c|awk '{print $2 "\t" $1}' \
  >background.sr.txt

zcat sr.geno.final.oneside.txt.gz \
  |awk '{if ($NF>0) print $1}' \
  |sort|uniq -c \
  |awk '{print $2 "\t" $1}'|sort -k1,1 \
  |join -j 1 - background.sr.txt \
  |awk '{ print $0 "\t" $2/$3}' \
  >recover.txt

##Require both##
zcat $SR_counts|awk '{if ($NF>0) print $1"@"$3}'|sort|uniq -c|awk '{if ($1==2) print $2}'>two.sided.pass.just1read.txt

join -j 2 <(cut -d"@" -f1 two.sided.pass.txt|sort|uniq -c) \
   <(cut -d"@" -f1 two.sided.pass.just1read.txt|sort|uniq -c) \
   |awk '{print $0"\t" $2/$3}' \
   >recover.bothsides.txt


## Find passing variants based off ROC check##

awk -v var=$rare_max  -v var1=$rare_single '{if ($NF>=var1 && $2<=var) print $1}' recover.txt>single.pass.txt
awk -v var=$rare_max -v var1=$rare_both '{if ($NF>=var1  && $2<=var) print $1}' recover.bothsides.txt>both.pass.txt

awk -v var=$common_min  -v var1=$common_single '{if ($NF>=var1 && $2>=var) print $1}' recover.txt>>single.pass.txt
awk -v var=$common_min -v var1=$common_both '{if ($NF>=var1  && $2>=var) print $1}' recover.bothsides.txt>>both.pass.txt

##SR background variant failures##
if [ $(cat both.pass.txt single.pass.txt \
  |fgrep -wvf - int.bed \
  |awk '{print $4}' \
  |sort -u \
  |fgrep -wf pass.srtest.txt \
  |fgrep -wf <(cat recover.txt recover.bothsides.txt|awk '{print $1}'|sort -u)|wc -l) -gt 0 ]
then
 cat both.pass.txt single.pass.txt \
  |{ fgrep -wvf - int.bed || [[ $? == 1 ]]; } \
  |awk '{print $4}' \
  |sort -u \
  |{ fgrep -wf pass.srtest.txt || [[ $? == 1 ]]; } \
  |{ fgrep -wf <(cat recover.txt recover.bothsides.txt|awk '{print $1}'|sort -u) || [[ $? == 1 ]]; } \
   >background.variant.fail.txt
else 
 echo "">background.variant.fail.txt 
fi

##Pull out variants to genotype##
join -j 1 -a 1 -a 2  -e "." -o 1.1 2.1 1.2 1.3 1.4 1.5 2.2 2.3 2.4 2.5 \
  <(zcat sr.geno.final.txt.gz \
  |fgrep -wf both.pass.txt \
  |awk '{print $1"@"$2"\t"$0}' |sort -k1,1) \
  <(zcat sr.geno.final.oneside.txt.gz \
  |fgrep -wf single.pass.txt|awk '{print $1"@"$2"\t"$0}'|sort -k1,1) \
  |awk '{if ($1!=".") print $1,$3,$4,$5,$6,$9,$10; else print $2,$7,$8,$5,$6,$9,$10}' \
  |gzip>combine.int.txt.gz

join -j 1 -a 1 -a 2 -e "." -o 1.1 2.1 1.2 1.3 1.4 1.5 1.6 1.7 2.2 2.3 2.4 2.5 \
  <(zcat combine.int.txt.gz) <(zcat sr.geno.final.oneside.txt.gz \
  |awk '{print $1"@"$2"\t"$0}' \
  |sort -k1,1)|awk '{ if ($1!=".") print $1,$3,$4,$5,$6,$7,$8,$11,$12; else print $2,$9,$10,$5,$6,$7,$8,$11,$12}' \
  |awk '{if ($5>0) print $2,$3,$4,$5; else if ($7>0) print $2,$3,6,$7;else print $2,$3,$8,$9}' \
  |tr ' ' '\t'|gzip>sr.geno.withfinal.fitler.txt.gz

##Genotype##
##normalization factor##
##Phred score for z score of 0 (p=0.5) normalized (3.0103) to 999 after subtracting nearest alternative genotype like GATK  ##
normalization=$(Rscript -e "print(-10*log10((1-pnorm(($median_hom/2)/$sd_het) )))"  \
                  | tr '\n' '\t' \
                  | awk '{print 999/($NF-3.0103)}')
##null genotype get max quality score##
zcat sr.geno.withfinal.fitler.txt.gz|awk '{if ($NF==0 && $3==0) print $0 "\t" 999}'|gzip>null.geno.txt.gz

##null genotype but has SR reads determined by poisson test##
if [ $(zcat sr.geno.withfinal.fitler.txt.gz|awk '{if ($NF==0 && $3>0) print}'|wc -l) -gt 0 ]
then
zcat sr.geno.withfinal.fitler.txt.gz|awk '{if ($NF==0 && $3>0) print}' \
  | Rscript -e 'd<-read.table("stdin")' \
  -e "z<-cbind(d[1],d[,2],d[,3],d[,4],round(matrix(apply(d[,3,drop=F],1, function (x) -10*log10(1-ppois(0, lambda=x))* $normalization) ,ncol=1)) )" \
  -e 'write.table(z,"null.wreads.geno.txt",col.names=FALSE,quote=FALSE,row.names=FALSE,sep = "\t")'
else 
echo "">null.wreads.geno.txt
fi

##genotype based on z score for observation with genotypes##
if [ $(zcat sr.geno.withfinal.fitler.txt.gz|awk '{if ($NF>0) print}'|wc -l) -gt 0 ]
then
zcat sr.geno.withfinal.fitler.txt.gz|awk '{if ($NF>0) print}' \
  | Rscript -e 'd<-read.table("stdin")' \
  -e "z<-apply(d[,3,drop=F]-(d[,4,drop=F]*$median_hom/2),1, function (x) -10*log10((1-pnorm(abs(x)/$sd_het) )))" \
  -e "z1<-apply(d[,3,drop=F]-((d[,4,drop=F]*$median_hom/2)-($median_hom/2)),1, function (x) -10*log10((1-pnorm(abs(x)/$sd_het) )))" \
  -e "z2<-apply(d[,3,drop=F]-((d[,4,drop=F]*$median_hom/2)+($median_hom/2)),1, function (x) -10*log10((1-pnorm(abs(x)/$sd_het) )))" \
  -e "zfinal<-round((pmin(z1,z2)-z)* $normalization)" \
  -e 'write.table(cbind(d[1],d[,2],d[,3],d[,4],zfinal),"wreads.geno.txt",col.names=FALSE,quote=FALSE,row.names=FALSE,sep = "\t")'
else 
echo "">wreads.geno.txt
fi

cat wreads.geno.txt null.wreads.geno.txt <(zcat null.geno.txt.gz)|sed '/^$/d'|awk '{if ($NF<0) print $1,$2,$3,$4,1;else if ($NF>999) print $1,$2,$3,$4,999 ;else print}' |tr ' ' '\t'|sort -k1,1 -k2,2|gzip>sr.geno.withquality.txt.gz

##Per variant quality scores##
##Assign a normalization factor to scale variant up to 999 if based on expected het median ##
normalization_var=$(Rscript -e "print(-10*log10(ppois(0, lambda=$median_hom/2)))"  \
                  | tr '\n' '\t' \
                  | awk '{print 999/($NF)}')

if [ $(zcat sr.geno.withfinal.fitler.txt.gz|awk '{if ($NF>0) print}'|wc -l) -gt 0 ]
then
zcat sr.geno.withfinal.fitler.txt.gz|awk '{if ($NF>0) print}' \
  | Rscript -e 'd<-read.table("stdin")' \
  -e 'x<-tapply(d[,3],d[,1],median)' \
  -e 'z<-cbind(names(x),(-10*log10(ppois(0, lambda=x))))' \
  -e 'write.table(z,"sr.variant.quality.final.txt",col.names=FALSE,quote=FALSE,row.names=FALSE,sep = "\t")'
else 
echo "">sr.variant.quality.final.txt
fi


##add back variants with no SR support##
if [ $(awk '{print $1}' sr.variant.quality.final.txt |fgrep -wvf - <(zcat $vcf|egrep -v "^#" |awk '{print $3}') |wc -l) -gt 0 ]
then

awk '{print $1}' sr.variant.quality.final.txt \
  |{ fgrep -wvf - <(zcat $vcf|egrep -v "^#" |awk '{print $3}') || [[ $? == 1 ]]; } \
  |awk '{print $1}' \
  |sort -u \
  |awk '{print $1 "\t" 0}' \
  >sr.variant.quality.final.null.txt

else
echo "">sr.variant.quality.final.null.txt
fi

awk -v var=$normalization_var '{if ($2*var>999) print $1 "\t" 999;else print $1 "\t" $2*var}'  sr.variant.quality.final.txt \
  |cat - sr.variant.quality.final.null.txt \
  |awk -F"\t" '{if ($1!="") print}'  \
  |sort -k1,1|gzip \
  >sr.variant.quality.final.txt.gz

