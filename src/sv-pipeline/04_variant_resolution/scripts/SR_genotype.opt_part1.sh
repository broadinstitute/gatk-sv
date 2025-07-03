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
batch=$9

sr_pval=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) col[$i]=i; next}
          {if ( $col["metric"]=="SR_sum_log_pval") print $col["cutoff"]}' $RF_cutoffs | head -n 1)
echo "sr_pval: $sr_pval"

sr_count=$(/opt/sv-pipeline/04_variant_resolution/scripts/convert_poisson_p.py $sr_pval)
echo "sr_count: $sr_count"

#Require both sides to have at least half of sr_count for training purposes
zcat ${SR_counts} \
  | awk -v sr_count=$sr_count '{if ($NF>(sr_count/2)) print $1"@"$3}' \
  | sort \
  | uniq -c \
  | awk '{if ($1==2) print $2}' \
  > two.sided.pass.txt
echo "len(two.sided.pass.txt): $(wc -l two.sided.pass.txt)"

# grep out variants which pass SR in random forest
svtk vcf2bed $vcf int.bed -i EVIDENCE

awk '{if ($NF~"SR") print $4}' int.bed> pass.srtest.txt
echo "len(pass.srtest.txt): $(wc -l pass.srtest.txt)"
echo "--- step1 ---"

# Join RD and SR genotypes and filter same as PE
cat $petrainfile \
  |awk -F'\t' -v OFS='\t' 'ARGIND==1{ids[$1]; next} ($1 in ids)' pass.srtest.txt - \
  > sr.train.include.txt
echo "len(sr.train.include.txt): $(wc -l sr.train.include.txt)"

zcat ${SR_sum} \
  | awk 'ARGIND==1{ids[$1]; next} ($1 in ids)' sr.train.include.txt - \
  | awk '{print $1"@"$2 "\t" $0}' \
  | awk 'ARGIND==1{ids[$1]; next} ($1 in ids)' two.sided.pass.txt - \
  | sort -k1,1 > sr_temp.txt

zcat $RD_melted_genotypes \
  | awk 'ARGIND==1{ids[$1]; next} ($4 in ids)' sr.train.include.txt - \
  | awk '{print $4"@"$5 "\t" $6}' \
  | awk 'ARGIND==1{ids[$1]; next} ($1 in ids)' two.sided.pass.txt - \
  | sort -k1,1 > rd_temp.txt

join -j 1 -a 1 -e "2" -o 1.2 1.3 1.4 2.2 sr_temp.txt rd_temp.txt \
  | tr ' ' '\t' > SR.RD.merged.txt
echo "len(SR.RD.merged.txt): $(wc -l SR.RD.merged.txt)"

# Get cutoffs to filter out incorrectly label hom in R and treat combine het (1 and 3) and hom (0 and 4) copy states
# throw out any copy state  calls that have reads less than with p=0.05 away from copy state 1 or 3

het_cutoff=$(awk '{print $1"@"$2"\t" $3 "\t" $4}' SR.RD.merged.txt \
  |awk '{if ($NF==3 || $NF==1) print $2}' \
  |Rscript -e 'd<-read.table("stdin")' \
  -e '(z<-1.645*mad(d[,1]))+median(d[,1])' \
  |tr '\n' '\t' \
  |awk '{print $NF}')
echo "het_cutoff: $het_cutoff"

##Rerun without excluding 0 or 4 copy states that fail het_cutoff###
awk -v var=$het_cutoff '{if (!(($4=="0" || $4=="4") && $3<var)) print}' SR.RD.merged.txt > SR.RD.hetfilter.merged.txt
echo "len(SR.RD.hetfilter.merged.txt): $(wc -l SR.RD.hetfilter.merged.txt)"
echo "--- step2 ---"

##median##
##get median from copy state 0 and 4###
median_hom=$(awk '{if ($NF==0 || $NF==4) print $3}'  SR.RD.hetfilter.merged.txt \
               | Rscript -e 'd<-scan("stdin", quiet=TRUE)' \
                         -e 'median(d)' \
               | tr '\n' '\t' \
               | awk '{print $NF}')
echo "median_hom: $median_hom"

##get std from 1 && 3  for hom restriction###
sd_het=$(awk '{if ($NF==1 || $NF==3) print $3}'  SR.RD.hetfilter.merged.txt \
           | Rscript -e 'd<-scan("stdin", quiet=TRUE)' \
                     -e 'mad(d)' \
           | tr '\n' '\t' \
           | awk '{print $NF*1.645}')
echo "sd_het: $sd_het"
##Genotype SR genotype (0-ref, then estimate copy state based on copy state that is 1 sd from sd_het  )##
zcat ${SR_sum} \
  | awk '{print $0 "\t" $1"@"$2}' \
  | awk -F'\t' -v OFS='\t' 'ARGIND==1{ids[$1]; next} ($4 in ids)' two.sided.pass.txt - \
  | cut -f1-3 \
  | awk -v var=$sr_count -v var1=$median_hom -v var2=$sd_het '{if ($3<var) print $1,$2,$3,0;else if ($3<=var1-var2) print $1,$2,$3,1; else print $1,$2,$3,int($3/(var1/2)+0.5)}'  \
  > sr.geno.final.txt
echo "len(sr.geno.final.txt): $(wc -l sr.geno.final.txt)"

zcat ${SR_sum} \
  | awk '{print $0 "\t" $1"@"$2}' \
  | awk -F'\t' -v OFS='\t' 'ARGIND==1{ids[$1]; next} (!($4 in ids))' two.sided.pass.txt - \
  | cut -f1-3 \
  | awk '{print $1,$2,$3,0}' \
  >> sr.geno.final.txt
echo "len(sr.geno.final.txt): $(wc -l sr.geno.final.txt)"


gzip -f sr.geno.final.txt

zcat ${SR_sum} \
  | awk '{print $0 "\t" $1"@"$2}' \
  | cut -f1-3 \
  | awk -v var=$sr_count -v var1=$median_hom -v var2=$sd_het '{if ($3<var) print $1,$2,$3,0;else if ($3<=var1-var2) print $1,$2,$3,1; else print $1,$2,$3,int($3/(var1/2)+0.5)}' \
  | gzip \
  > sr.geno.final.oneside.txt.gz
echo "len(sr.geno.final.oneside.txt.gz): $(wc -l sr.geno.final.oneside.txt.gz)"

echo "--- step3 ---"
##filter by quality of site by looking at % of calls with ##
##Allow just one side##
zcat sr.geno.final.oneside.txt.gz \
  |awk '{if ($3>1) print $1}' \
  |sort|uniq -c|awk '{print $2 "\t" $1}' \
  >background.sr.txt
echo "len(background.sr.txt): $(wc -l background.sr.txt)"

zcat sr.geno.final.oneside.txt.gz \
  |awk '{if ($NF>0) print $1}' \
  |sort|uniq -c \
  |awk '{print $2 "\t" $1}'|sort -k1,1 \
  |join -j 1 - background.sr.txt \
  |awk '{ print $0 "\t" $2/$3}' \
  >recover.txt
echo "len(recover.txt): $(wc -l recover.txt)"

##Require both##
zcat $SR_counts|awk '{if ($NF>0) print $1"@"$3}'|sort|uniq -c|awk '{if ($1==2) print $2}'>two.sided.pass.just1read.txt
echo "len(two.sided.pass.just1read.txt): $(wc -l two.sided.pass.just1read.txt)"

join -j 2 <(cut -d"@" -f1 two.sided.pass.txt|sort|uniq -c) \
   <(cut -d"@" -f1 two.sided.pass.just1read.txt|sort|uniq -c) \
   |awk '{print $0"\t" $2/$3}' \
   >recover.bothsides.txt
echo "len(recover.bothsides.txt): $(wc -l recover.bothsides.txt)"


##pull out all observations that have PE_RD genotyping support##

zcat sr.geno.final.oneside.txt.gz|awk '{if ($NF>0) print $1"@"$2}'>sr.final.ids.oneside.txt

echo "--- step4 ---"

##pull out cnvs gt1kb and not located on x or y##
zcat $RD_melted_genotypes|egrep -v "^X|^Y"|awk '{if ($3-$2>=1000) print $4}'|sort -u>idsgt1kb.txt
echo "len(idsgt1kb.txt): $(wc -l idsgt1kb.txt)"

awk -F'\t' -v OFS='\t' 'ARGIND==1{ids[$1]; next} ($1 in ids)' <(cut -d '@' -f1 sr.final.ids.oneside.txt|sort -u) \
  <(zcat $RD_melted_genotypes|awk -F'\t' -v OFS='\t' '{if ($6!=2) print $4,$5}') \
  > nonref_rd.txt
echo "len(nonref_rd.txt): $(wc -l nonref_rd.txt)"

zcat $pegenotypes \
  |awk -F'\t' -v OFS='\t' 'ARGIND==1{ids[$1]; next} ($1 in ids)' <(cut -d '@' -f1 sr.final.ids.oneside.txt|sort -u) - \
  |awk -F'\t' -v OFS='\t' '{if ($NF>0) print $1,$2}' \
  |cat - nonref_rd.txt \
  |awk -F'\t' -v OFS='\t' 'ARGIND==1{ids[$1]; next} ($1 in ids)' idsgt1kb.txt - \
  |awk -F'\t' -v OFS='\t' 'ARGIND==1{ids[$1]; next} ($1 in ids)' pass.srtest.txt - \
  |sort -u \
  |tr '\t' '@' \
  >pass.pe_rd.txt
echo "len(pass.pe_rd.txt): $(wc -l pass.pe_rd.txt)"

##look for optimal cutoffs for SR variants using a 1% freq cutoff##

##passing ids##
##filter the recover file which has variants passing both single and both ends ##
zcat sr.geno.final.oneside.txt.gz \
  | awk '{if ($NF>0) print $1 "\t" $1"@"$2 }' \
  | sort -k1,1 > sr_oneside_temp.txt

cat recover.txt \
  | sort -k1,1 \
  | join -j 1 - sr_oneside_temp.txt \
  | tr ' ' '\t' \
  | sort -k5,5 \
  | join -1 5 -2 1 - <(sort pass.pe_rd.txt) \
  | awk '{print $2,$3,$4,$5,$1}' \
  | tr ' ' '\t' \
  > recover.single.txt
echo "len(recover.single.txt): $(wc -l recover.single.txt)"

cat recover.bothsides.txt \
  | sort -k1,1 \
  | join -j 1 - sr_oneside_temp.txt \
  | tr ' ' '\t' \
  | sort -k5,5 \
  | join -1 5 -2 1 - <(sort pass.pe_rd.txt) \
  | awk '{print $2,$3,$4,$5,$1}' \
  | tr ' ' '\t' \
  > recover.both.txt
echo "len(recover.both.txt): $(wc -l recover.both.txt)"

cat recover.txt \
  | sort -k1,1 \
  | join -j 1 - sr_oneside_temp.txt \
  | tr ' ' '\t' \
  | sort -k5,5 \
  | join -1 5 -2 1 -v 1 - <(sort pass.pe_rd.txt) \
  | awk '{print $2,$3,$4,$5,$1}' \
  | tr ' ' '\t' \
  > recover.single.fail.txt
echo "len(recover.single.fail.txt): $(wc -l recover.single.fail.txt)"

cat recover.bothsides.txt \
  | sort -k1,1 \
  | join -j 1 - sr_oneside_temp.txt \
  | tr ' ' '\t' \
  | sort -k5,5 \
  | join -1 5 -2 1 -v 1 - <(sort pass.pe_rd.txt) \
  | awk '{print $2,$3,$4,$5,$1}' \
  | tr ' ' '\t' \
  > recover.both.fail.txt
echo "len(recover.both.fail.txt): $(wc -l recover.both.fail.txt)"

echo "--- step5 ---"

##set rare max to be atleast 2 people with a given variant for cohorts with families we want to make sure to include reare private variants to a family##
rare_min=0
rare_max=$(awk '{if (NR/100>=2) print NR/100;else print 2}' $whitelist|tail -n 1)
common_min=$rare_max
common_max=$(cat $whitelist|wc -l )
echo "rare_min: $rare_min"
echo "rare_max: $rare_max"
echo "common_min: $common_min"
echo "common_max: $common_max"

for i in 0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1
do
  for j in 0 .1 .2 .3 .4 .5 .6 .7 .8 .9 1
  do
     /opt/sv-pipeline/04_variant_resolution/scripts/optimalsrcutoff.sh $i $j $rare_min $rare_max
     /opt/sv-pipeline/04_variant_resolution/scripts/optimalsrcutoff.sh $i $j $common_min $common_max
  done
done

echo "rare_min: $rare_min"
echo "rare_max: $rare_max"
echo "common_min: $common_min"
echo "common_max: $common_max"

##combine values from different freq checks##
cat pe_support.combined.check.*.$rare_min.$rare_max.txt>pe_support.combined.rare.txt
cat pe_support.combined.check.*.$common_min.$common_max.txt>pe_support.combined.common.txt

echo "--- step6 ---"

##Roc Check##

##Sens Max Value##
baseline_common=$(awk '{if ($1==0 && $2==0) print $3}' pe_support.combined.common.txt)
baseline_rare=$(awk '{if ($1==0 && $2==0) print $3}' pe_support.combined.rare.txt)
echo "baseline_common: $baseline_common"
echo "baseline_rare: $baseline_rare"

rare_single=$(awk -v var=$baseline_rare '{print $0 "\t" ((($3/($3+$4))-1)^2) + ((($3/var)-1)^2)  }' pe_support.combined.rare.txt|sort -nk5,5|head -n 1|awk '{print $1}')
rare_both=$(awk -v var=$baseline_rare '{print $0 "\t" ((($3/($3+$4))-1)^2) + ((($3/var)-1)^2)  }' pe_support.combined.rare.txt|sort -nk5,5|head -n 1|awk '{print $2}')
echo "rare_single: $rare_single"
echo "rare_both: $rare_both"

common_single=$(awk -v var=$baseline_common '{print $0 "\t" ((($3/($3+$4))-1)^2) + ((($3/var)-1)^2)  }' pe_support.combined.common.txt|sort -nk5,5|head -n 1|awk '{print $1}')
common_both=$(awk -v var=$baseline_common '{print $0 "\t" ((($3/($3+$4))-1)^2) + ((($3/var)-1)^2)  }' pe_support.combined.common.txt|sort -nk5,5|head -n 1|awk '{print $2}')
echo "common_single: $common_single"
echo "common_both: $common_both"

rm pe_support.combined.check.*.txt

echo -e "sr_count"'\t'$sr_count>sr_metric_file.txt
echo -e "median_hom"'\t'$median_hom>>sr_metric_file.txt
echo -e "sd_het"'\t'$sd_het>>sr_metric_file.txt
echo -e "rare_min"'\t'$rare_min>>sr_metric_file.txt
echo -e "rare_max"'\t'$rare_max>>sr_metric_file.txt
echo -e "common_min"'\t'$common_min>>sr_metric_file.txt
echo -e "common_max"'\t'$common_max>>sr_metric_file.txt
echo -e "rare_single"'\t'$rare_single>>sr_metric_file.txt
echo -e "rare_both"'\t'$rare_both>>sr_metric_file.txt
echo -e "common_single"'\t'$common_single>>sr_metric_file.txt
echo -e "common_both"'\t'$common_both>>sr_metric_file.txt
mv sr_metric_file.txt "$batch.sr_metric_file.txt"
