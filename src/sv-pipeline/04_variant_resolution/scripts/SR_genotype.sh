#!/bin/bash
#
# SR_genotype.sh
#


set -euo pipefail


vcf=$1
SR_counts=$2
SR_sum=$3
RD_melted_genotypes=$4
RF_cutoffs=$5
whitelist=$6
petrainfile=$7
pegenotypes=$8

sr_pval=$( awk -F'\t' '{if ( $5=="SR_sum_log_pval") print $2}' $RF_cutoffs | head -n 1)
sr_count=$(/opt/sv-pipeline/04_variant_resolution/scripts/convert_poisson_p.py $sr_pval)

#Require both sides to have at least half of sr_count for training purposes
zcat ${SR_counts} \
  | awk -v sr_count=$sr_count '{if ($NF>(sr_count/2)) print $1"@"$3}' \
  | sort \
  | uniq -c \
  | awk '{if ($1==2) print $2}' \
  > two.sided.pass.txt

# grep out variants which pass SR in random forest
svtk vcf2bed $vcf int.bed -i EVIDENCE

awk '{if ($NF~"SR") print $4}' int.bed> pass.srtest.txt

echo "step1"

# Join RD and SR genotypes and filter same as PE
cat $petrainfile|fgrep -wf pass.srtest.txt > sr.train.include.txt

join -j 1  -a 1 -e "2" -o 1.2 1.3 1.4 2.2 \
    <(zcat ${SR_sum} \
        | fgrep -wf sr.train.include.txt \
        | awk '{print $1"@"$2 "\t" $0}' \
        | fgrep -wf two.sided.pass.txt \
        | sort -k1,1 ) \
    <(zcat $RD_melted_genotypes|fgrep -wf sr.train.include.txt \
        | awk '{print $4"@"$5 "\t" $6}' \
        | fgrep -wf two.sided.pass.txt \
        | sort -k1,1) \
  | tr ' ' '\t' \
  > SR.RD.merged.txt 

# Get cutoffs to filter out incorrectly label hom in R and treat combine het (1 and 3) and hom (0 and 4) copy states 
# throw out any copy state  calls that have reads less than with p=0.05 away from copy state 1 or 3

het_cutoff=$(awk '{print $1"@"$2"\t" $3 "\t" $4}' SR.RD.merged.txt \
  |awk '{if ($NF==3 || $NF==1) print $2}' \
  |Rscript -e 'd<-read.table("stdin")' \
  -e '(z<-1.645*mad(d[,1]))+median(d[,1])' \
  |tr '\n' '\t' \
  |awk '{print $NF}')

##Rerun without excluding 0 or 4 copy states that fail het_cutoff###
awk -v var=$het_cutoff '{if (!(($4=="0" || $4=="4") && $3<var)) print}' SR.RD.merged.txt > SR.RD.hetfilter.merged.txt

echo "step2"

##median##
##get median from copy state 0 and 4###
median_hom=$(awk '{if ($NF==0 || $NF==4) print $3}'  SR.RD.hetfilter.merged.txt \
               | Rscript -e 'd<-scan("stdin", quiet=TRUE)' \
                         -e 'median(d)' \
               | tr '\n' '\t' \
               | awk '{print $NF}')
##get std from 1 && 3  for hom restriction###          
sd_het=$(awk '{if ($NF==1 || $NF==3) print $3}'  SR.RD.hetfilter.merged.txt \
           | Rscript -e 'd<-scan("stdin", quiet=TRUE)' \
                     -e 'mad(d)' \
           | tr '\n' '\t' \
           | awk '{print $NF*1.645}')

##Genotype SR genotype (0-ref, then estimate copy state based on copy state that is 1 sd from sd_het  )##
zcat ${SR_sum} \
  | awk '{print $0 "\t" $1"@"$2}' \
  | fgrep -wf two.sided.pass.txt \
  | cut -f1-3 \
  | awk -v var=$sr_count -v var1=$median_hom -v var2=$sd_het '{if ($3<var) print $1,$2,$3,0;else if ($3<=var1-var2) print $1,$2,$3,1; else print $1,$2,$3,int($3/(var1/2)+0.5)}'  \
  > sr.geno.final.txt

zcat ${SR_sum} \
  | awk '{print $0 "\t" $1"@"$2}' \
  | fgrep -wvf two.sided.pass.txt \
  | cut -f1-3 \
  | awk '{print $1,$2,$3,0}' \
  >> sr.geno.final.txt


gzip sr.geno.final.txt

zcat ${SR_sum} \
  | awk '{print $0 "\t" $1"@"$2}' \
  | cut -f1-3 \
  | awk -v var=$sr_count -v var1=$median_hom -v var2=$sd_het '{if ($3<var) print $1,$2,$3,0;else if ($3<=var1-var2) print $1,$2,$3,1; else print $1,$2,$3,int($3/(var1/2)+0.5)}' \
  | gzip \
  > sr.geno.final.oneside.txt.gz
  
echo "step3"
##filter by quality of site by looking at % of calls with ##
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


##pull out all observations that have PE_RD genotyping support##

zcat sr.geno.final.oneside.txt.gz|awk '{if ($NF>0) print $1"@"$2}'>sr.final.ids.oneside.txt

echo "step4"

##pull out cnvs gt1kb and not located on x or y##
zcat $RD_melted_genotypes|egrep -v "^X|^Y"|awk '{if ($3-$2>=1000) print $4}'|sort -u>idsgt1kb.txt


zcat $pegenotypes \
  |fgrep -wf <(cut -d '@' -f1 sr.final.ids.oneside.txt|sort -u) \
  |awk '{if ($NF>0) print $1"@"$2}' \
  |cat - <(fgrep -wf <(cut -d '@' -f1 sr.final.ids.oneside.txt|sort -u) \
  <(zcat $RD_melted_genotypes|awk '{if ($6!=2) print $4"@"$5}')) \
  |fgrep -wf idsgt1kb.txt \
  |fgrep -wf pass.srtest.txt \
  |sort -u \
  >pass.pe_rd.txt

##look for optimal cutoffs for SR variants using a 1% freq cutoff##

##passing ids##
##filter the recover file which has variants passing both single and both ends ##
cat recover.txt \
  |sort -k1,1 \
  |join -j 1 - <(zcat sr.geno.final.oneside.txt.gz|awk '{if ($NF>0) print $1 "\t" $1"@"$2 }'|sort -k1,1) \
  |tr ' ' '\t' \
  |fgrep -wf pass.pe_rd.txt \
  >recover.single.txt

cat recover.bothsides.txt \
  |sort -k1,1 \
  |join -j 1 - <(zcat sr.geno.final.oneside.txt.gz|awk '{if ($NF>0) print $1 "\t" $1"@"$2 }'|sort -k1,1) \
  |tr ' ' '\t' \
  |fgrep -wf pass.pe_rd.txt \
  >recover.both.txt

cat recover.txt \
  |sort -k1,1 \
  |join -j 1 - <(zcat sr.geno.final.oneside.txt.gz|awk '{if ($NF>0) print $1 "\t" $1"@"$2 }'|sort -k1,1) \
  |tr ' ' '\t' \
  |fgrep -wvf pass.pe_rd.txt \
  >recover.single.fail.txt

cat recover.bothsides.txt \
  |sort -k1,1 \
  |join -j 1 - <(zcat sr.geno.final.oneside.txt.gz|awk '{if ($NF>0) print $1 "\t" $1"@"$2 }'|sort -k1,1) \
  |tr ' ' '\t' \
  |fgrep -wvf pass.pe_rd.txt \
  >recover.both.fail.txt

echo "step5"

rare_min=0
rare_max=$(awk '{print NR/100}' $whitelist|tail -n 1)
common_min=$(awk '{print NR/100}' $whitelist|tail -n 1)
common_max=$(cat $whitelist|wc -l )

for i in 0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1
do
  for j in 0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1
  do
     /opt/sv-pipeline/04_variant_resolution/scripts/optimalsrcutoff.sh $i $j $rare_min $rare_max
     /opt/sv-pipeline/04_variant_resolution/scripts/optimalsrcutoff.sh $i $j $common_min $common_max
  done
done

##combine values from different freq checks##
cat pe_support.combined.check.*.$rare_min.$rare_max.txt>pe_support.combined.rare.txt
cat pe_support.combined.check.*.$common_min.$common_max.txt>pe_support.combined.common.txt

echo "step6"

##Roc Check##

##Sens Max Value##
baseline_common=$(awk '{if ($1==0 && $2==0) print $3}' pe_support.combined.common.txt)
baseline_rare=$(awk '{if ($1==0 && $2==0) print $3}' pe_support.combined.rare.txt)

rare_single=$(awk -v var=$baseline_rare '{print $0 "\t" ((($3/($3+$4))-1)^2) + ((($3/var)-1)^2)  }' pe_support.combined.rare.txt|sort -nk5,5|head -n 1|awk '{print $1}')
rare_both=$(awk -v var=$baseline_rare '{print $0 "\t" ((($3/($3+$4))-1)^2) + ((($3/var)-1)^2)  }' pe_support.combined.rare.txt|sort -nk5,5|head -n 1|awk '{print $2}')

common_single=$(awk -v var=$baseline_common '{print $0 "\t" ((($3/($3+$4))-1)^2) + ((($3/var)-1)^2)  }' pe_support.combined.common.txt|sort -nk5,5|head -n 1|awk '{print $1}')
common_both=$(awk -v var=$baseline_common '{print $0 "\t" ((($3/($3+$4))-1)^2) + ((($3/var)-1)^2)  }' pe_support.combined.common.txt|sort -nk5,5|head -n 1|awk '{print $2}')

echo "step7"

## Find passing variants based off ROC check##

awk -v var=$rare_max  -v var1=$rare_single '{if ($NF>=var1 && $2<=var) print $1}' recover.txt>single.pass.txt
awk -v var=$rare_max -v var1=$rare_both '{if ($NF>=var1  && $2<=var) print $1}' recover.bothsides.txt>both.pass.txt

awk -v var=$common_min  -v var1=$common_single '{if ($NF>=var1 && $2>=var) print $1}' recover.txt>>single.pass.txt
awk -v var=$common_min -v var1=$common_both '{if ($NF>=var1  && $2>=var) print $1}' recover.bothsides.txt>>both.pass.txt

rm pe_support.combined.check.*.txt

##SR background variant failures##
if [ $(cat both.pass.txt single.pass.txt \
  |fgrep -wvf - int.bed \
  |awk '{print $4}' \
  |sort -u \
  |fgrep -wf pass.srtest.txt|wc -l) -gt 0 ]
then
 cat both.pass.txt single.pass.txt \
  |fgrep -wvf - int.bed \
  |awk '{print $4}' \
  |sort -u \
  |fgrep -wf pass.srtest.txt \
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
  <(zcat combine.int.txt.gz) <(zcat ${SR_sum} \
  |awk '{print $1"@"$2 "\t" $1 "\t" $2 "\t" $3 "\t" 0 }' \
  |sort -k1,1)|awk '{ if ($1!=".") print $1,$3,$4,$5,$6,$7,$8,$11,$12; else print $2,$9,$10,$5,$6,$7,$8,$11,$12}' \
  |awk '{if ($5>0) print $2,$3,$4,$5; else if ($7>0) print $2,$3,6,$7;else print $2,$3,$8,$9}' \
  |tr ' ' '\t'|gzip>sr.geno.withfinal.fitler.txt.gz

##Old slow method on firecloud##
#cat <(zcat sr.geno.final.txt.gz|fgrep -wf both.pass.txt) \
#  <(zcat sr.geno.final.oneside.txt.gz|fgrep -wf single.pass.txt) \
#  <(zcat ${SR_sum}|awk '{print $1 "\t" $2 "\t" $3 "\t" 0 }') \
#  |awk '{if ($NF>0) print 1"\t" $0;else print 0 "\t" $0}' \
#  |sort -nrk1,1 \
#  |cut -f2-|tr ' ' '\t' \
#  |awk '!seen[$1"_"$2]++'|gzip \
#  >sr.geno.withfinal.fitler.txt.gz

##Genotype##
##normalization factor##
##Phred score for z score of 0 (p=0.5) normalized (3.0103) to 999 after subtracting nearest alternative genotype like GATK  ##
normalization=$(Rscript -e "print(-10*log10((1-pnorm(($median_hom/2)/$sd_het) )))"  \
                  | tr '\n' '\t' \
                  | awk '{print 999/($NF-3.0103)}')
##null genotype get max quality score##
zcat sr.geno.withfinal.fitler.txt.gz|awk '{if ($NF==0 && $3==0) print $0 "\t" 999}'|gzip>null.geno.txt.gz

##null genotype but has SR reads determined by poisson test##
zcat sr.geno.withfinal.fitler.txt.gz|awk '{if ($NF==0 && $3>0) print}' \
  | Rscript -e 'd<-read.table("stdin")' \
  -e "z<-cbind(d[1],d[,2],d[,3],d[,4],round(matrix(apply(d[,3,drop=F],1, function (x) -10*log10(1-ppois(0, lambda=x))* $normalization) ,ncol=1)) )" \
  -e 'write.table(z,"null.wreads.geno.txt",col.names=FALSE,quote=FALSE,row.names=FALSE,sep = "\t")'

##genotype based on z score for observation with genotypes##
zcat sr.geno.withfinal.fitler.txt.gz|awk '{if ($NF>0) print}' \
  | Rscript -e 'd<-read.table("stdin")' \
  -e "z<-apply(d[,3,drop=F]-(d[,4,drop=F]*$median_hom/2),1, function (x) -10*log10((1-pnorm(abs(x)/$sd_het) )))" \
  -e "z1<-apply(d[,3,drop=F]-((d[,4,drop=F]*$median_hom/2)-($median_hom/2)),1, function (x) -10*log10((1-pnorm(abs(x)/$sd_het) )))" \
  -e "z2<-apply(d[,3,drop=F]-((d[,4,drop=F]*$median_hom/2)+($median_hom/2)),1, function (x) -10*log10((1-pnorm(abs(x)/$sd_het) )))" \
  -e "zfinal<-round((pmin(z1,z2)-z)* $normalization)" \
  -e 'write.table(cbind(d[1],d[,2],d[,3],d[,4],zfinal),"wreads.geno.txt",col.names=FALSE,quote=FALSE,row.names=FALSE,sep = "\t")'


cat wreads.geno.txt null.wreads.geno.txt <(zcat null.geno.txt.gz)|awk '{if ($NF<0) print $1,$2,$3,$4,1;else if ($NF>999) print $1,$2,$3,$4,999 ;else print}' |tr ' ' '\t'|sort -k1,1 -k2,2|gzip>sr.geno.withquality.txt.gz

##Per variant quality scores##
##Assign a normalization factor to scale variant up to 999 if based on expected het median ##
normalization_var=$(Rscript -e "print(-10*log10(ppois(0, lambda=$median_hom/2)))"  \
                  | tr '\n' '\t' \
                  | awk '{print 999/($NF)}')
echo "step8"
zcat sr.geno.withfinal.fitler.txt.gz|awk '{if ($NF>0) print}' \
  | Rscript -e 'd<-read.table("stdin")' \
  -e 'x<-tapply(d[,3],d[,1],median)' \
  -e 'z<-cbind(names(x),(-10*log10(ppois(0, lambda=x))))' \
  -e 'write.table(z,"sr.variant.quality.final.txt",col.names=FALSE,quote=FALSE,row.names=FALSE,sep = "\t")'

##add back variants with no SR support##
awk '{print $1}' sr.variant.quality.final.txt \
  |fgrep -wvf - <(zcat $vcf|egrep -v "^#" |awk '{print $3}') \
  |awk '{print $1}' \
  |sort -u \
  |awk '{print $1 "\t" 0}' \
  >sr.variant.quality.final.null.txt

awk -v var=$normalization_var '{if ($2*var>999) print $1 "\t" 999;else print $1 "\t" $2*var}'  sr.variant.quality.final.txt \
  |cat - sr.variant.quality.final.null.txt \
  |sort -k1,1|gzip \
  >sr.variant.quality.final.txt.gz
