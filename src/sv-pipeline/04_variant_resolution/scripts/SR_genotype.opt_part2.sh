#!/bin/bash
#
# SR_genotype.sh
#

set -euo pipefail

##Genotype##
function assign_genotype_quality () {

  # Arguments
  INFILE=$1
  VIDS_LIST=$2
  MEDIAN_HET=$3
  STDEV=$4
  GQ_FILE=$5
  VARQ_FILE=$6

  ##normalization factor##
  ##Normalize phred scores to 999 if count is the mean  ##
  normalization=$(Rscript -e "print(((pnorm(0) )))"  \
                    | tr '\n' '\t' \
                    | awk '{print 999/($NF)}')
  ##null genotype get max quality score##
  zcat $INFILE \
    | awk '{if ($NF==0 && $3==0) print $0 "\t" 999}' \
    | gzip \
    > null.geno.txt.gz

  ##null genotype but has SR reads determined by poisson test##
  if [ $(zcat $INFILE|awk '{if ($NF==0 && $3>0) print}'|wc -l) -gt 0 ]
  then
    zcat $INFILE|awk '{if ($NF==0 && $3>0) print}' \
      | Rscript -e 'd<-read.table("stdin")' \
      -e "z<-cbind(d[1],d[,2],d[,3],d[,4],round(matrix(apply(d[,3,drop=F],1, function (x) (1-ppois(0, lambda=x))* $normalization) ,ncol=1)) )" \
      -e 'write.table(z,"null.wreads.geno.txt",col.names=FALSE,quote=FALSE,row.names=FALSE,sep = "\t")'
  else
    echo "" > null.wreads.geno.txt
  fi

  ##genotype based on z score for observation with genotypes##
  if [ $(zcat $INFILE|awk '{if ($NF==1) print}'|wc -l) -gt 0 ]
  then
    zcat $INFILE|awk '{if ($NF==1) print}' \
      | Rscript -e 'd<-read.table("stdin")' \
      -e "z<-apply(d[,3,drop=F]-$MEDIAN_HET,1, function (x) ((pnorm(x/$STDEV) )))" \
      -e "zfinal<-round(z* $normalization)" \
      -e 'write.table(cbind(d[1],d[,2],d[,3],d[,4],zfinal),"wreads.geno.het.txt",col.names=FALSE,quote=FALSE,row.names=FALSE,sep = "\t")'
  else
    echo "" > wreads.geno.het.txt
  fi
  if [ $(zcat $INFILE|awk '{if ($NF>1) print}'|wc -l) -gt 0 ]
  then
    zcat $INFILE|awk '{if ($NF>1) print}' \
      | Rscript -e 'd<-read.table("stdin")' \
      -e "z<-apply(d[,3,drop=F]-2*$MEDIAN_HET,1, function (x) ((pnorm(x/$STDEV) )))" \
      -e "zfinal<-round(z* $normalization)" \
      -e 'write.table(cbind(d[1],d[,2],d[,3],d[,4],zfinal),"wreads.geno.hom.txt",col.names=FALSE,quote=FALSE,row.names=FALSE,sep = "\t")'
  else
    echo "" > wreads.geno.hom.txt
  fi

  cat wreads.geno.het.txt wreads.geno.hom.txt null.wreads.geno.txt <(zcat null.geno.txt.gz) \
    | sed '/^$/d' \
    | awk '{if ($NF<0) print $1,$2,$3,$4,1;else if ($NF>999) print $1,$2,$3,$4,999 ;else print}' \
    | tr ' ' '\t' \
    | sort -k1,1 -k2,2 \
    | gzip \
    > $GQ_FILE

  ##Per variant quality scores##
  ##Assign a normalization factor to scale variant up to 999 if based on expected het median ##
  normalization_var=$(Rscript -e "print(-10*log10(ppois(0, lambda=$MEDIAN_HET)))"  \
                    | tr '\n' '\t' \
                    | awk '{print 999/($NF)}')

  if [ $(zcat $INFILE|awk '{if ($NF>0) print}'|wc -l) -gt 0 ]
  then
    zcat $INFILE|awk '{if ($NF>0) print}' \
      | Rscript -e 'd<-read.table("stdin")' \
      -e 'x<-tapply(d[,3],d[,1],median)' \
      -e 'z<-cbind(names(x),(-10*log10(ppois(0, lambda=x))))' \
      -e 'write.table(z,"sr.variant.quality.final.txt",col.names=FALSE,quote=FALSE,row.names=FALSE,sep = "\t")'
  else
    echo "" > sr.variant.quality.final.txt
  fi

  ##add back variants with no SR support##
  awk '{print $1}' sr.variant.quality.final.txt \
    | awk -F'\t' -v OFS='\t' 'ARGIND==1{keys[$1]; next} { if (!($1 in keys)) print }' - $VIDS_LIST \
    | awk '{print $1}' \
    | sort -u \
    | awk '{print $1 "\t" 0}' \
    > sr.variant.quality.final.null.txt

  awk -v var=$normalization_var '{if ($2*var>999) print $1 "\t" 999;else print $1 "\t" $2*var}'  sr.variant.quality.final.txt \
    | cat - sr.variant.quality.final.null.txt \
    | awk -F"\t" '{if ($1!="") print}'  \
    | sort -k1,1 \
    | gzip \
    > $VARQ_FILE
}

#whole cohort files that can be split
vcf=$1
SR_counts=$2
SR_sum=$3
metric_file=$4
svtype=$5
# Recommended 1.6
hom_cutoff_multiplier=$6
# Recommended 105
median_hom_ins=$7

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

if [ $svtype == "INS" ]; then
  # Insertion hom-var counts follow a consistent distribution across cohorts
  # This is a minimal correction to the genotyping model to minimize hom-var over-calling
  median_hom=$median_hom_ins
fi

bothside_median_het=$(echo "print(0.5*$median_hom)" | python)
bothside_hom_cutoff=$(echo "print($hom_cutoff_multiplier*$bothside_median_het)" | python)
bothside_sd_het=$sd_het

oneside_median_het=$(echo "print(0.5*$bothside_median_het)" | python)
oneside_hom_cutoff=$(echo "print($hom_cutoff_multiplier*$oneside_median_het)" | python)
oneside_sd_het=$(echo "print(0.5*$sd_het)" | python)

zcat ${SR_counts} \
  | awk -v sr_count=$sr_count '{if ($NF>(sr_count/2)) print $1"@"$3}' \
  | sort \
  | uniq -c \
  | awk '{if ($1==2) print $2}' \
  > two.sided.pass.txt
# Output columns:
# 1) VID@Sample (with sufficient SR count on both sides)

# get variants to genotype
svtk vcf2bed $vcf int.bed -i EVIDENCE

awk '{if ($NF~"SR") print $4}' int.bed> pass.srtest.txt

##Genotype SR genotype (0-ref, then estimate copy state based on copy state that is 1 sd from sd_het  )##
zcat ${SR_sum} \
  | awk '{print $0 "\t" $1"@"$2}' \
  | awk -F'\t' -v OFS='\t' 'ARGIND==1{keys[$1]; next} { if ($NF in keys) print }' two.sided.pass.txt - \
  | cut -f1-3 \
  | awk -v var=$sr_count -v hom_cutoff=$bothside_hom_cutoff -v median_het=$bothside_median_het '{if ($3<var) print $1,$2,$3,0;else if ($3<=hom_cutoff) print $1,$2,$3,1; else print $1,$2,$3,(int(($3/median_het)+0.5) < 2) ? 2 : int(($3/median_het)+0.5)}'  \
  | tr ' ' '\t' \
  | gzip \
  > sr.geno.final.txt.gz

zcat ${SR_sum} \
  | awk '{print $0 "\t" $1"@"$2}' \
  | cut -f1-3 \
  | awk -v var=$sr_count -v hom_cutoff=$oneside_hom_cutoff -v median_het=$oneside_median_het '{if ($3<var) print $1,$2,$3,0;else if ($3<=hom_cutoff) print $1,$2,$3,1; else print $1,$2,$3,(int(($3/median_het)+0.5) < 2) ? 2 : int(($3/median_het)+0.5)}'  \
  | tr ' ' '\t' \
  | gzip \
  > sr.geno.final.oneside.txt.gz
# Output columns:
# 1) VID
# 2) Sample ID
# 3) SR sum
# 4) Genotype (one-sided)
  
##Allow just one side##
zcat sr.geno.final.oneside.txt.gz \
  | awk '{if ($3>1) print $1}' \
  | sort \
  | uniq -c \
  | awk '{print $2 "\t" $1}' \
  > background.sr.txt
# Output columns:
# 1) VID
# 2) Number of samples with an SR count over 1

zcat sr.geno.final.oneside.txt.gz \
  | awk '{if ($NF>0) print $1}' \
  | sort \
  | uniq -c \
  | awk '{print $2 "\t" $1}'|sort -k1,1 \
  | join -j 1 - background.sr.txt \
  | awk '{ print $0 "\t" $2/$3}' \
  > recover.txt
# Output columns:
# 1) VID
# 2) Number of non-ref genotypes
# 3) Number of samples with more than 1 SR
# 4) Fraction of samples with more than 1 SR that have non-ref genotypes, i.e. (2)/(3)

##Require both##
zcat $SR_counts|awk '{if ($NF>0) print $1"@"$3}'|sort|uniq -c|awk '{if ($1==2) print $2}'>two.sided.pass.just1read.txt
# Output columns:
# 1) VID@Sample (with at least 1 SR on both sides)

join -j 2 <(cut -d"@" -f1 two.sided.pass.txt|sort|uniq -c) \
   <(cut -d"@" -f1 two.sided.pass.just1read.txt|sort|uniq -c) \
   | awk '{print $0"\t" $2/$3}' \
   > recover.bothsides.txt
# Output columns:
# 1) VID (with at least 1 sample having bothside support)
# 2) Number of samples with bothside support
# 3) Number of samples with non-zero counts on both sides
# 4) Fraction of samples with non-zero SR count that have sufficient counts on both sides, i.e. (2)/(3)

## Find passing variants based off ROC check##
awk -v var=$rare_max  -v var1=$rare_single '{if ($NF>=var1 && $2<=var) print $1}' recover.txt>single.pass.txt
awk -v var=$rare_max -v var1=$rare_both '{if ($NF>=var1  && $2<=var) print $1}' recover.bothsides.txt>both.pass.txt

awk -v var=$common_min  -v var1=$common_single '{if ($NF>=var1 && $2>=var) print $1}' recover.txt>>single.pass.txt
# Output columns:
# 1) VID (one-sided variants that don't have too many background samples,
#         i.e. samples with non-zero SR that are genotyped as hom-ref)
awk -v var=$common_min -v var1=$common_both '{if ($NF>=var1  && $2>=var) print $1}' recover.bothsides.txt>>both.pass.txt
# Output columns:
# 1) VID (two-sided variants that don't have too many background samples,
#         i.e. samples with non-zero SR that have two-sided support)

# Note the VIDs in single.pass.txt and both.pass.txt are not mutually exclusive

##SR background variant failures##
cat both.pass.txt single.pass.txt \
  | awk -F'\t' -v OFS='\t' 'ARGIND==1{keys[$1]; next} { if ($4 != "name" && !($4 in keys)) print }' - int.bed \
  | awk '{print $4}' \
  | sort -u \
  | awk -F'\t' -v OFS='\t' 'ARGIND==1{keys[$1]; next} { if ($1 in keys) print }' pass.srtest.txt - \
  | awk -F'\t' -v OFS='\t' 'ARGIND==1{keys[$1]; next} { if ($1 in keys) print }' <(cat recover.txt recover.bothsides.txt|awk '{print $1}'|sort -u) - \
  > background.variant.fail.txt

##Pull out variants to genotype##
join -j 1 -a 1 -a 2  -e "." -o 1.1 2.1 1.2 1.3 1.4 1.5 2.2 2.3 2.4 2.5 \
  <(zcat sr.geno.final.txt.gz \
  |awk '{print $1"@"$2"\t"$0}' |sort -k1,1) \
  <(zcat sr.geno.final.oneside.txt.gz \
  |awk '{print $1"@"$2"\t"$0}'|sort -k1,1) \
  |awk '{if ($1!="." && $5) print $1,$3,$4,$5,$6,$9,$10; else print $2,$7,$8,$5,$6,$9,$10}' \
  |awk '{if ($4>0) print $1,$2,$3,$4,$5; else print $1,$2,$3,$6,$7}' \
  |gzip>combine.int.txt.gz
# Output columns:
# 1) VID@Sample (passing variants only)
# 2) VID
# 3) Sample
# 4) SR count (prioritizing two-sided if its two-sided genotype is non-ref)
# 5) Genotype (prioritizing two-sided if its two-sided genotype is non-ref)

# Reformat
zcat combine.int.txt.gz | awk '{print $2,$3,$4,$5}' |tr ' ' '\t'|gzip>sr.geno.withfinal.filter.txt.gz
# Output columns:
# 1) VID (passing variants only)
# 2) Sample
# 3) SR count
# 4) Genotype (prioritizing two-sided if it is passing)

# Genotype bothsided and onesided variants separately
echo "Genotyping bothside..."
zcat sr.geno.withfinal.filter.txt.gz \
  | awk -F'\t' -v OFS='\t' '{ print $1"@"$2,$0 }' \
  | awk -F'\t' -v OFS='\t' 'ARGIND==1{keys[$1]; next} { if ($1 in keys) print }' two.sided.pass.txt - \
  | awk -F'\t' -v OFS='\t' '{ print $2,$3,$4,$5 }' \
  | gzip \
  > sr.geno.withfinal.filter.bothside.txt.gz
bcftools query -f '%ID\n' $vcf > vids.list
assign_genotype_quality "sr.geno.withfinal.filter.bothside.txt.gz" "vids.list" $bothside_median_het $bothside_sd_het "sr.geno.withquality.bothside.txt.gz" "sr.variant.quality.bothside.txt.gz"

echo "Genotyping oneside..."
zcat sr.geno.withfinal.filter.txt.gz \
  | awk -F'\t' -v OFS='\t' '{ print $1"@"$2,$0 }' \
  | awk -F'\t' -v OFS='\t' 'ARGIND==1{keys[$1]; next} { if (!($1 in keys)) print }' two.sided.pass.txt - \
  | awk -F'\t' -v OFS='\t' '{ print $2,$3,$4,$5 }' \
  | gzip \
  > sr.geno.withfinal.filter.oneside.txt.gz
assign_genotype_quality "sr.geno.withfinal.filter.oneside.txt.gz" "vids.list" $oneside_median_het $oneside_sd_het "sr.geno.withquality.oneside.txt.gz" "sr.variant.quality.oneside.txt.gz"

echo "Concatenating results..."
zcat sr.geno.withquality.bothside.txt.gz sr.geno.withquality.oneside.txt.gz \
  | sort -k1,1 \
  | gzip \
  > sr.geno.withquality.txt.gz

# Computes weighted average of one-sided and two-sided variant qualities, weighted by fraction of two-sided samples
join -j 1 -a 1 -a 2 -e "." -o 1.1 1.2 2.1 2.2 <(zcat sr.variant.quality.bothside.txt.gz) <(zcat sr.variant.quality.oneside.txt.gz) \
  | awk '{ if ($1==".") {print $3,$2,$4} else {print $1,$2,$4} }' \
  | join -j 1 -e "." -a 1 -o 1.1 1.2 1.3 2.4 - recover.bothsides.txt \
  | awk '{ if ($2==".") {print $1,$3} else if ($3==".") {print $1,$2} else {print $1,($2*$4)+($3*(1-$4))} }' \
  | gzip \
  > sr.variant.quality.final.txt.gz
