#!/bin/bash
#

set -euo pipefail

bed=$1
PE_counts=$2
RD_genotypes=$3
RD_melted_genotypes=$4
RF_cutoffs=$5
blacklist=$6
batch=$7

egrep "DEL|DUP" ${bed} > cnv.bed;

# filter out potential overlapping variants##
bedtools intersect -wa -wb -a cnv.bed -b cnv.bed \
  | awk '{if ($4 != $11) print $4"\n"$11}' \
  | sort -u \
  > cnv.exclude.all.bed;

# size filter (>1kb)
awk '{if ($3-$2>=1000) print $4}' cnv.bed > size.pass.txt;

# remove depth only
awk '{if ($NF=="depth") print $4}' ${bed} > depthonly.fail.txt;
  
# only allow up to three copy state (0,1,2 or 1,2) for training 
if [ $(tail -n +2 ${RD_genotypes} \
  | cut -f4- \
  | fgrep -w 1 \
  | fgrep -w 2 \
  | awk '{ for (i = 2; i <= NF; ++i) if ($i>3) print $1 }' \
  | sort -u|wc -l) -gt 0 ]
then
 tail -n +2 ${RD_genotypes} \
  | cut -f4- \
  | fgrep -w 1 \
  | fgrep -w 2 \
  | awk '{ for (i = 2; i <= NF; ++i) if ($i>3) print $1 }' \
  | sort -u \
  > copystatefail.txt;
fi

if [ $(tail -n +2 ${RD_genotypes} \
  | cut -f4- \
  | fgrep -w 2 \
  | fgrep -w 3 \
  | awk '{ for (i = 2; i <= NF; ++i) if ($i>4 || $i<2) print $1 }' \
  | sort -u |wc -l) -gt 0 ]
then
 tail -n +2 ${RD_genotypes} \
  | cut -f4- \
  | fgrep -w 2 \
  | fgrep -w 3 \
  | awk '{ for (i = 2; i <= NF; ++i) if ($i>4 || $i<2) print $1 }' \
  | sort -u \
  >> copystatefail.txt;
fi

if [ $(tail -n +2 ${RD_genotypes} \
  | cut -f4- \
  | fgrep -wvf copystatefail.txt \
  | fgrep -w 1 \
  | fgrep -w 2 \
  | awk '{print $1}'|wc -l) -gt 0 ] 
then
tail -n +2 ${RD_genotypes} \
  | cut -f4- \
  | fgrep -wvf copystatefail.txt \
  | fgrep -w 1 \
  | fgrep -w 2 \
  | awk '{print $1}' \
  > copystate.pass.txt;
fi

if [ $(tail -n +2 ${RD_genotypes} \
  | cut -f4- \
  | fgrep -wvf copystatefail.txt \
  | fgrep -w 2 \
  | fgrep -w 3 \
  | awk '{print $1}'|wc -l) -gt 0 ] 
then  
tail -n +2 ${RD_genotypes} \
  | cut -f4- \
  | fgrep -wvf copystatefail.txt \
  | fgrep -w 2 \
  | fgrep -w 3 \
  | awk '{print $1}' \
  >> copystate.pass.txt;
fi

# check for mulitcopystate
zcat ${RD_melted_genotypes} \
  | awk '!seen[$1"_" $NF]++' \
  | awk '{print $1}' \
  | sort \
  | uniq -c \
  > cnv.copystate.multi.count.txt;
  
if [ $(cat cnv.copystate.multi.count.txt \
   | awk '{if ($1>3) print $2}' \
   | fgrep -wf - cnv.bed \
   | awk '{if ($3-$2<5000 && $3-$2>1000) print $1,$2,$3,$4,$6,$5}' \
   | tr ' ' '\t' \
   | egrep -v "^X|^Y"|wc -l ) -gt 0 ]  
then
 cat cnv.copystate.multi.count.txt \
   | awk '{if ($1>3) print $2}' \
   | fgrep -wf - cnv.bed \
   | awk '{if ($3-$2<5000 && $3-$2>1000) print $1,$2,$3,$4,$6,$5}' \
   | tr ' ' '\t' \
   | egrep -v "^X|^Y" \
   > multi.test.bed;
fi

# remove repetitive breakpoints
cat ${blacklist} \
  | bedtools coverage -a <(awk '{print $1,$2,$2+1,$4 "\n" $1,$3,$3+1,$4}' cnv.bed | tr ' ' '\t')  -b - \
  | awk '{if($NF>0) print $4}' \
  | sort -u \
  > repeat.breakpoint.fail.ids.txt;

# final filter - remove X & Y
# check to make sure you have variants to train PE cutoffs
if [ ! -f copystate.pass.txt ] || [ ! -f size.pass.txt ]
then 
  echo "No variants for training. Increase number of variants." 
  exit
fi

awk '{if ($1!~"X" && $1!~"Y") print $4}' cnv.bed \
  | fgrep -wf copystate.pass.txt \
  | fgrep -wf size.pass.txt \
  | fgrep -wvf <(cat cnv.exclude.all.bed depthonly.fail.txt repeat.breakpoint.fail.ids.txt) \
  > "$batch.pe.train.include.txt";

# select training
pe_pval=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) col[$i]=i; next}
          {if ( $col["metric"]=="PEQ") print $col["cutoff"]}' $RF_cutoffs)
pe_count=$(/opt/sv-pipeline/04_variant_resolution/scripts/convert_poisson_p.py $pe_pval)
zcat ${PE_counts} \
  | awk -v var=$pe_count '{if ($3>=var) print}' \
  | fgrep -v name \
  | gzip -c \
  > pe.geno.all.txt.gz;

# Join RD and PE genotypes
join -j 1 -a 1 -e "2" -o 1.2 1.3 1.4 2.2 \
    <(zcat pe.geno.all.txt.gz \
        | fgrep -wf "$batch.pe.train.include.txt" \
        | awk '{print $1"_"$2 "\t" $0}' \
        | sort -k1,1 ) \
    <(fgrep -wf "$batch.pe.train.include.txt" <(zcat ${RD_melted_genotypes}) \
        | awk '{print $4"_"$5 "\t" $6}' \
        | sort -k1,1) \
  | tr ' ' '\t' \
  > PE.RD.merged.txt; 


# Get cutoffs to filter out incorrectly label hom in R and treat combine het (1 and 3) and hom (0 and 4) copy states 
# throw out any copy state  calls that have reads less than with p=0.05 away from copy state 1 or 3
het_cutoff=$(awk '{print $1"_"$2"\t" $3 "\t" $4}' PE.RD.merged.txt \
  |awk '{if ($NF==3 || $NF==1) print $2}' \
  |Rscript -e 'd<-read.table("stdin")' \
  -e '(z<-1.645*mad(d[,1]))+median(d[,1])' \
  |tr '\n' '\t' \
  |awk '{print $NF}')

# Rerun with excluding 0 or 4 copy states that fail het_cutoff
awk -v var=$het_cutoff '{if (!(($4=="0" || $4=="4") && $3<var)) print}' PE.RD.merged.txt > PE.RD.hetfilter.merged.txt

# get median from copy state 0 and 4
median_hom=$(awk '{if ($NF==0 || $NF==4) print $3}' PE.RD.hetfilter.merged.txt \
              | Rscript -e 'd<-scan("stdin", quiet=TRUE)' \
                        -e 'median(d)'\
              | tr '\n' '\t' \
              | awk '{print $NF}')

##If no median to train on then assign as 2 plus the het cutoff rounded up##
if [ $median_hom == "NA" ]
then
median_hom=$(echo 1 \
   |awk -v het_cutoff=$het_cutoff 'function ceil(x, y){y=int(x); return(x>y?y+1:y)} {print ceil(het_cutoff+2)}')
fi

# get std from 1 && 3 for hom restriction
sd_het=$(awk '{if ($NF==1 || $NF==3) print $3}' PE.RD.hetfilter.merged.txt \
          | Rscript -e 'd<-scan("stdin", quiet=TRUE)' \
                    -e 'mad(d)' \
          | tr '\n' '\t' \
          | awk '{print $NF*1.645}')

##generate PE metric file##
echo -e "pe_count"'\t'$pe_count>"$batch.pe_metric_file.txt"
echo -e "median_hom"'\t'$median_hom>>"$batch.pe_metric_file.txt"
echo -e "sd_het"'\t'$sd_het>>"$batch.pe_metric_file.txt"

zcat ${PE_counts} \
  | fgrep -v name \
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

##null genotype but has PE reads determined by poisson test##
zcat pe.geno.final.txt.gz|awk '{if ($NF==0 && $3>0) print}' > null.wreads.txt
if [ -s null.wreads.txt ]; then
  cat null.wreads.txt \
    | Rscript -e 'd<-read.table("stdin")' \
    -e "z<-cbind(d[1],d[,2],d[,3],d[,4],round(matrix(apply(d[,3,drop=F],1, function (x) -10*log10(1-ppois(0, lambda=x))* $normalization) ,ncol=1)) )" \
    -e 'write.table(z,"null.wreads.geno.txt",col.names=FALSE,quote=FALSE,row.names=FALSE,sep = "\t")'
else
  touch null.wreads.geno.txt
fi

##genotype based on z score for observation with genotypes##
zcat pe.geno.final.txt.gz|awk '{if ($NF>0) print}' \
  | Rscript -e 'd<-read.table("stdin")' \
  -e "z<-apply(d[,3,drop=F]-(d[,4,drop=F]*$median_hom/2),1, function (x) -10*log10((1-pnorm(abs(x)/$sd_het) )))" \
  -e "z1<-apply(d[,3,drop=F]-((d[,4,drop=F]*$median_hom/2)-($median_hom/2)),1, function (x) -10*log10((1-pnorm(abs(x)/$sd_het) )))" \
  -e "z2<-apply(d[,3,drop=F]-((d[,4,drop=F]*$median_hom/2)+($median_hom/2)),1, function (x) -10*log10((1-pnorm(abs(x)/$sd_het) )))" \
  -e "zfinal<-round((pmin(z1,z2)-z)* $normalization)" \
  -e 'write.table(cbind(d[1],d[,2],d[,3],d[,4],zfinal),"wreads.geno.txt",col.names=FALSE,quote=FALSE,row.names=FALSE,sep = "\t")'

cat wreads.geno.txt null.wreads.geno.txt <(zcat null.geno.txt.gz)|awk '{if ($NF<0) print $1,$2,$3,$4,1;else if ($NF>999) print $1,$2,$3,$4,999 ;else print}' |tr ' ' '\t'|sort -k1,1 -k2,2|gzip>"$batch.pe.geno.withquality.txt.gz"

##Per variant quality scores##
##Assign a normalization factor to scale variant up to 999 if based on expected het median ##
normalization_var=$(Rscript -e "print(-10*log10(ppois(0, lambda=$median_hom/2)))"  \
                  | tr '\n' '\t' \
                  | awk '{print 999/($NF)}')

zcat pe.geno.final.txt.gz|awk '{if ($NF>0) print}' \
  | Rscript -e 'd<-read.table("stdin")' \
  -e 'x<-tapply(d[,3],d[,1],median)' \
  -e 'z<-cbind(names(x),(-10*log10(ppois(0, lambda=x))))' \
  -e 'write.table(z,"pe.variant.quality.final.txt",col.names=FALSE,quote=FALSE,row.names=FALSE,sep = "\t")'

##add back variants with no PE support##
awk '{print $1}' pe.variant.quality.final.txt \
  |fgrep -wvf - <(zcat pe.geno.final.txt.gz) \
  |awk '{print $1}' \
  |sort -u \
  |awk '{print $1 "\t" 0}' \
  >pe.variant.quality.final.null.txt

awk -v var=$normalization_var '{if ($2*var>999) print $1 "\t" 999;else print $1 "\t" $2*var}' pe.variant.quality.final.txt \
  |cat - pe.variant.quality.final.null.txt \
  |sort -k1,1|gzip \
  >"$batch.pe.variant.quality.final.txt.gz"

