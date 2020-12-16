#!/bin/bash
#
# 
#
# Copyright (C) 2018 Harrison Brand
# Distributed under terms of the MIT license.

set -euo pipefail

PE_counts=$1
metric_file=$2


##load metric file##
pe_count=$(awk '{if($1=="pe_count") print $2}' $metric_file)
median_hom=$(awk '{if($1=="median_hom") print $2}' $metric_file)
sd_het=$(awk '{if($1=="sd_het") print $2}' $metric_file)

zcat ${PE_counts} \
  | { fgrep -v name || [[ $? == 1 ]]; } \
  |awk -v var=$pe_count -v var1=$median_hom -v var2=$sd_het '{if ($3<var) print $1,$2,$3,0;else if ($3<=var1-var2) print $1,$2,$3,1; else print $1,$2,$3,int($3/(var1/2)+0.5)}'  | gzip \
  > pe.geno.final.txt.gz

##Genotype##
##normalization factor##
##Phred score for z score of 0 (p=0.5) normalized (3.0103) to 999 after subtracting nearest alternative genotype like GATK  ##
normalization=$(Rscript -e "print(-10*log10((1-pnorm(($median_hom/2)/$sd_het) )))"  \
                  | tr '\n' '\t' \
                  | awk '{print 999/($NF-3.0103)}')

##null genotype get max quality score##
zcat pe.geno.final.txt.gz|awk '{if ($NF==0 && $3==0) print $0 "\t" 999}'|gzip>null.geno.txt.gz

##null genotype but has SR reads determined by poisson test##
if [ $(zcat pe.geno.final.txt.gz|awk '{if ($NF==0 && $3>0) print}'|wc -l) -gt 0 ] 
then
zcat pe.geno.final.txt.gz|awk '{if ($NF==0 && $3>0) print}' \
  | Rscript -e 'd<-read.table("stdin")' \
  -e "z<-cbind(d[1],d[,2],d[,3],d[,4],round(matrix(apply(d[,3,drop=F],1, function (x) -10*log10(1-ppois(0, lambda=x))* $normalization) ,ncol=1)) )" \
  -e 'write.table(z,"null.wreads.geno.txt",col.names=FALSE,quote=FALSE,row.names=FALSE,sep = "\t")'
else 
echo "">null.wreads.geno.txt
fi

##genotype based on z score for observation with genotypes##
if [ $(zcat pe.geno.final.txt.gz|awk '{if ($NF>0) print}'|wc -l) -gt 0 ]
then
zcat pe.geno.final.txt.gz|awk '{if ($NF>0) print}' \
  | Rscript -e 'd<-read.table("stdin")' \
  -e "z<-apply(d[,3,drop=F]-(d[,4,drop=F]*$median_hom/2),1, function (x) -10*log10((1-pnorm(abs(x)/$sd_het) )))" \
  -e "z1<-apply(d[,3,drop=F]-((d[,4,drop=F]*$median_hom/2)-($median_hom/2)),1, function (x) -10*log10((1-pnorm(abs(x)/$sd_het) )))" \
  -e "z2<-apply(d[,3,drop=F]-((d[,4,drop=F]*$median_hom/2)+($median_hom/2)),1, function (x) -10*log10((1-pnorm(abs(x)/$sd_het) )))" \
  -e "zfinal<-round((pmin(z1,z2)-z)* $normalization)" \
  -e 'write.table(cbind(d[1],d[,2],d[,3],d[,4],zfinal),"wreads.geno.txt",col.names=FALSE,quote=FALSE,row.names=FALSE,sep = "\t")'
else 
echo "">wreads.geno.txt
fi

cat wreads.geno.txt null.wreads.geno.txt <(zcat null.geno.txt.gz)|sed '/^$/d'|awk '{if ($NF<0) print $1,$2,$3,$4,1;else if ($NF>999) print $1,$2,$3,$4,999 ;else print}' |tr ' ' '\t'|sort -k1,1 -k2,2|gzip>pe.geno.withquality.txt.gz

##Per variant quality scores##
##Assign a normalization factor to scale variant up to 999 if based on expected het median ##
normalization_var=$(Rscript -e "print(-10*log10(ppois(0, lambda=$median_hom/2)))"  \
                  | tr '\n' '\t' \
                  | awk '{print 999/($NF)}')

if [ $(zcat pe.geno.final.txt.gz|awk '{if ($NF>0) print}'|wc -l) -gt 0 ]
then
zcat pe.geno.final.txt.gz|awk '{if ($NF>0) print}' \
  | Rscript -e 'd<-read.table("stdin")' \
  -e 'x<-tapply(d[,3],d[,1],median)' \
  -e 'z<-cbind(names(x),(-10*log10(ppois(0, lambda=x))))' \
  -e 'write.table(z,"pe.variant.quality.final.txt",col.names=FALSE,quote=FALSE,row.names=FALSE,sep = "\t")'
else 
 echo "">pe.variant.quality.final.txt
fi


##add back variants with no PE support##
awk '{print $1}' pe.variant.quality.final.txt \
  |{ fgrep -wvf - <(zcat pe.geno.final.txt.gz) || [[ $? == 1 ]]; }  \
  |awk '{print $1}' \
  |sort -u \
  |awk '{print $1 "\t" 0}' \
  >pe.variant.quality.final.null.txt

awk -v var=$normalization_var '{if ($2*var>999) print $1 "\t" 999;else print $1 "\t" $2*var}' pe.variant.quality.final.txt \
  |cat - pe.variant.quality.final.null.txt \
  |sort -k1,1 \
  |awk -F"\t" '{if ($1!="") print}'  |gzip \
  >pe.variant.quality.final.txt.gz


