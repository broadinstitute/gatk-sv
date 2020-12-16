#!/bin/bash
#
# clean_vcf_part1b.sh
#

set -euxo pipefail

##gzipped vcf from clean vcf part1.sh##
int_vcf_gz=$1

##Remove CNVs that are improperly genotyped by depth because they are nested within a real CNV##

##Determine columns of VCF after header##
zcat $int_vcf_gz\
  |sed -n '1,1000p'\
  |egrep ^# \
  |tail -n 1 \
  |tr '\t' '\n' \
  |cat -n - \
  >col.txt

##Only affects CNV so pull those out##
zcat $int_vcf_gz \
  |awk '{if ($5~"DEL" || $5~"DUP" || $1~"#") print}' \
  |svtk vcf2bed stdin stdout \
  |awk -F"\t" '{if ($6=="") $6="blanksample";print $0}' OFS='\t' \
  |gzip>int.bed.gz

##list of potenital overlaps with a normal copy state variant (>5kb variants require depth but nested events could be missed; i.e a duplication with a nest deletion will have a normal copy state for the deletion)##
##flip bed intersect so largest is CNV is always first##
bedtools intersect -wa -wb -a <(zcat int.bed.gz|awk '{if ($3-$2>=5000 ) print}') \
-b <(zcat int.bed.gz|awk '{if ($3-$2>=5000) print}') \
  |awk -F'\t' '{if ($4!=$10 && $3-$2>=$9-$8 && $5!=$11) print ;\
 else if ($4!=$10 && $5!=$11) print $7,$8,$9,$10,$11,$12,$1,$2,$3,$4,$5,$6}' OFS='\t' \
  |awk -F'\t' '{if ($6!="blanksample") print}' \
  |sort -u \
  >normaloverlap.txt
  

##pull out the depth based copy number variant for each normal overlapping variant## 
{ cat <(zcat $int_vcf_gz|awk -F"\t" '{if ($1~"#") print}') \
  <(awk '{print $4 "\n" $10}' normaloverlap.txt|sort -u|fgrep -wf - <(zcat $int_vcf_gz)) || true; }\
  |awk '{if ($1!~"#") $1=$3;print}' OFS="\t" \
  |awk '{if ($1~"#" || $5=="<DEL>" || $5=="<DUP>") print}' \
  |vcftools --vcf - --stdout --extract-FORMAT-info RD_CN \
  |awk -F"\t" 'NR==1{for (i=3;i<=NF;i++) header[i]=$i} NR>1{for(j=3;j<=NF;j++) print $1"@"header[j] "\t" $j }' \
  |sort -k1,1 \
  |gzip \
  >RD_CN.normalcheck.FORMAT.gz


##pull out evidence supporting each normal overlapping variant## 
{ cat <(zcat $int_vcf_gz|awk -F"\t" '{if ($1~"#") print}') \
  <(awk '{print $4 "\n" $10}' normaloverlap.txt|sort -u|fgrep -wf - <(zcat $int_vcf_gz)) || true; }\
  |awk '{if ($1!~"#") $1=$3;print}' OFS="\t"\
  |vcftools --vcf - --stdout --extract-FORMAT-info EV \
  |awk -F"\t" 'NR==1{for (i=3;i<=NF;i++) header[i]=$i} NR>1{for(j=3;j<=NF;j++) print $1"@"header[j] "\t" $j }' \
  |sort -k1,1 \
  |gzip \
  >EV.normalcheck.FORMAT.gz


##check if nested is incorrectly classified as normal##
touch overlap.test.txt
while read bed
do
 echo $bed|tr ' ' '\t'|cut -f1-6 >large.bed
 echo $bed|tr ' ' '\t'|cut -f7-12>small.bed
 ##require at least 50% coverage to consider a variant overlapping##
 overlap=$(bedtools coverage -a small.bed -b large.bed|awk '{if ($NF>=0.50) print "YES";else print "NO"}')

 if [ "$overlap" == "YES" ]
 then
  smallid=$(awk '{print $4}' small.bed)
  
  ##pull out variants that are called a variants for both the smaller and larger CNVs (don't have normal copy state to check for)##
   if [ $(awk '{print $NF}' small.bed \
       |tr ',' '\n' \
       |fgrep -wvf - <(awk -F"[,\t]" -v var=$smallid '{for(i=6;i<=NF;i++) print var"@"$i "\t" $4"@"$i "\t" $5}' large.bed)|wc -l) -gt 0 ]
    then
         awk '{print $NF}' small.bed \
          |tr ',' '\n' \
          |fgrep -wvf - <(awk -F"[,\t]" -v var=$smallid '{for(i=6;i<=NF;i++) print var"@"$i "\t" $4"@"$i "\t" $5}' large.bed) \
          >>overlap.test.txt
   fi
 fi
done<normaloverlap.txt


##determine variants that need to be revised from a normal copy state into a CNV##
cat overlap.test.txt \
  |sort -k1,1 \
  |join -j 1 - <(zcat RD_CN.normalcheck.FORMAT.gz) \
  |join -j 1 - <(zcat EV.normalcheck.FORMAT.gz) \
  |tr ' ' '\t' \
  |sort -k2,2 \
  |join -1 2 -2 1 - <(zcat RD_CN.normalcheck.FORMAT.gz) \
  |awk '{if ($3=="DUP" && $4==2 && $6==3) print $2 "\t" 1; else if ($3=="DEL" && $4==2 && $6==1)  print $2 "\t" 3 }' \
  |tr '@' '\t'\
  |sort -u \
  >geno.normal.revise.txt

##Update genotypes##
{ zfgrep -wf <(awk '{print $1}' geno.normal.revise.txt|sort -u) $int_vcf_gz || true; }\
  |bgzip \
  >subset.vcf.gz || true

##pull out and revise vcf line that needs to be edited##
while read variant
do
 
 echo $variant
 #note no longer change depth from id.txt (column 2)##
 { fgrep -w $variant geno.normal.revise.txt || true; }|awk '{print $2 "\t" $3}'>id.txt
 zcat subset.vcf.gz |{ fgrep -w $variant || true; }>line.txt
 
 cat line.txt  \
   |tr '\t' '\n' \
   |paste col.txt - \
   |tr ':' '\t' \
   |awk 'NR==FNR{inFileA[$1]=$2; next} {if ($2 in inFileA ) $3="0/1"; print }' OFS='\t' id.txt - \
   |awk 'NR==FNR{inFileA[$1]=$2; next} {if ($2 in inFileA ) $4=$6; print }' OFS='\t' id.txt - \
   |cut -f3-|tr '\t' ':' \
   |tr '\n' '\t' \
   |awk '{print $0}' \
   |awk '{ sub(/[ \t]+$/, ""); print }' \
   >>normal.revise.vcf.lines.txt

done< <(awk '{print $1}' geno.normal.revise.txt|sort -u)


##rewrite vcf with updated genotypes##

cat <(zcat $int_vcf_gz|fgrep -wvf <(awk '{print $3}' normal.revise.vcf.lines.txt|sort -u)) \
  <(sed 's/\t$//' normal.revise.vcf.lines.txt) \
  |vcf-sort \
  |bgzip \
  >normal.revise.vcf.gz || true

 bcftools index normal.revise.vcf.gz 

##get copy state per variant##
zcat normal.revise.vcf.gz \
  |awk '{if ($1!~"#") $1=$3;print}' OFS="\t" \
  |vcftools --vcf - --stdout --extract-FORMAT-info RD_CN \
  |gzip \
  >copystate.RD_CN.FORMAT.gz

##get copy state per variant##
zcat copystate.RD_CN.FORMAT.gz \
  |awk 'NR>1{for(i=3;i<=NF;i++) lines[$1 "\t" $i]++ } END{for (x in lines) print x}' \
  |gzip \
  >copystate.per.variant.txt.gz

##Find multi-allelic for del or dup ; CNV >1kb we trust depth ##
##del##
zcat copystate.per.variant.txt.gz \
  |awk '{if ($2!="." && $2>3) print $1}' \
  |sort -u \
  |fgrep -wf <(zcat int.bed.gz|awk -F"\t" '{if ($5=="DEL" && $3-$2>=1000) print $4}' ) \
  >multi.cnvs.txt || true

##dup##
zcat copystate.per.variant.txt.gz \
  |awk '{if ($2!="." && ($2<1 || $2>4)) print $1}' \
  |sort -u \
  |fgrep -wf <(zcat int.bed.gz|awk -F"\t" '{if ($5=="DUP" && $3-$2>=1000) print $4}' ) \
  >>multi.cnvs.txt || true
