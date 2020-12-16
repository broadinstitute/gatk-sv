#!/bin/bash
#
# clean_VCF_part5.sh
#

# Deprecated in favor of CleanVcf5.wdl

set -euo pipefail

##gzipped combined bed file##
##combined output file from clean_vcf_part2.sh##
revise_vcf_lines=$1
normal_revise_vcf=$2
famfile=$3
sexchr_revise=$4
multi_geno_ids_txt=$5
outliers_samples_list=$6

# use BCFTOOLS 1.9
BCFTOOLS=/usr/local/bin/bcftools

cat <(zcat $normal_revise_vcf|fgrep -wvf <(zcat $revise_vcf_lines|awk '{if ($1!="") print $3}'|sort -u)) \
  <(zcat $revise_vcf_lines|awk '{if ($1!="") print}' |tr ' ' '\t') \
  |vcf-sort \
  |bgzip \
  >overlap.revise.vcf.gz || true
  
##create bed of VCF##
svtk vcf2bed overlap.revise.vcf.gz stdout|gzip> overlap.revise.bed.gz 

##multi check##
zcat overlap.revise.vcf.gz \
  |awk '{if ($1!~"#") $1=$3;print}' OFS="\t" \
  |vcftools --vcf - --remove $outliers_samples_list --stdout --extract-FORMAT-info RD_CN \
  |gzip \
  >copystate.RD_CN.FORMAT.gz

zcat overlap.revise.vcf.gz \
  |awk '{if ($1!~"#") $1=$3;print}' OFS="\t" \
  |vcftools --vcf - --remove $outliers_samples_list --stdout --extract-FORMAT-info GT \
  |gzip \
  >genotype.gt.FORMAT.gz  

##New method for determining copy state based on >1% of people having an multi-allelic copy state as define above##
vf_1=$(zcat copystate.RD_CN.FORMAT.gz|awk 'NR==1{print (NF-2) * 0.01}'|awk '{if ($1<=1) print 2; else print }' )

zcat copystate.RD_CN.FORMAT.gz \
   |fgrep -wf <(zcat overlap.revise.bed.gz|awk -F"\t" '{if ($5=="DEL" && $3-$2>=1000) print $4}' ) \
   |awk 'NR>1{for(i=3;i<=NF;i++) if ($i!="." && $i>3) print  $1 }' \
   |sort \
   |uniq -c \
   |awk -v vf_1=$vf_1 '{if ($1>vf_1)print $2}' \
   |gzip \
    >multi.del.ids.txt.gz || true

zcat copystate.RD_CN.FORMAT.gz \
   |fgrep -wf <(zcat overlap.revise.bed.gz|awk -F"\t" '{if ($5=="DUP" && $3-$2>=1000) print $4}' ) \
   |awk 'NR>1{for(i=3;i<=NF;i++) if ($i!="." && $i>4) print  $1 }' \
   |sort \
   |uniq -c \
   |awk -v vf_1=$vf_1 '{if ($1>vf_1)print $2}' \
    >multi.dup.ids.txt || true
##Case with CN 0,1,2,3,4##
zcat copystate.RD_CN.FORMAT.gz \
   |fgrep -wf <(zcat overlap.revise.bed.gz \
   |awk -F"\t" '{if ($5=="DUP" && $3-$2>=1000) print $4}') \
   |awk 'NR>1{for(i=3;i<=NF;i++) if ($i!="." && ($i<1 || $i>4)) print  $1 "\t" $i }'\
   |sort -u \
   |awk  '{print $1}' \
   |sort \
   |uniq -c \
   |awk '{if ($1>4) print $2}'>gt4copystate.txt ||true
zcat copystate.RD_CN.FORMAT.gz \
   |fgrep -wf <(zcat overlap.revise.bed.gz|awk -F"\t" '{if ($5=="DUP" && $3-$2>=1000) print $4}' ) \
   |awk 'NR>1{for(i=3;i<=NF;i++) if ($i!="." && ($i<1 || $i>4)) print  $1 }' \
   |sort \
   |uniq -c \
   |fgrep -wf gt4copystate.txt \
   |awk -v vf_1=$vf_1 '{if ($1>vf_1)print $2}' \
    >>multi.dup.ids.txt || true
sort -u multi.dup.ids.txt|gzip >multi.dup.ids.txt.gz||true

##Regenotype to determine multiallelic; we just change copy state for some nested variants and we need to make sure we get proper genotype for these; also previous stages have different notaion for multiallelic and we need to make this uniform; this is a CN based regenotyping so restricted to >5kb ##
##Genotype big dup##
svtk vcf2bed overlap.revise.vcf.gz stdout \
  |gzip>regeno.bed.gz

##add variants that are <5kb because clustering but have a mutliallelic genotype from before##
zcat genotype.gt.FORMAT.gz \
   |awk '{if ($1~"DUP") print}' \
   |awk '{for (i = 3; i <= NF; ++i) print $1 "\t" $i}' \
   |awk '{if ($2!="1/1" && $2!="0/0" && $2!="0/1" && $2!="./.") print $1}' \
   |fgrep -wvf <(zcat multi.dup.ids.txt.gz) \
   |sort -u>gt5kb.dup.ids.txt || true

zcat genotype.gt.FORMAT.gz \
   |awk '{if ($1~"DEL") print}' \
   |awk '{for (i = 3; i <= NF; ++i) print $1 "\t" $i}' \
   |awk '{if ($2!="1/1" && $2!="0/0" && $2!="0/1" && $2!="./.") print $1}' \
   |fgrep -wvf <(zcat multi.del.ids.txt.gz) \
   |sort -u>gt5kb.del.ids.txt || true

##generate list##
##CNV >5kb, split del and dup ##
if [ -f multi.dup.ids.txt.gz ]
then
 zcat regeno.bed.gz  \
  |awk '{if ($3-$2>=5000 && $5=="DUP")print $4}' \
  |fgrep -wvf <(zcat multi.dup.ids.txt.gz) \
  >>gt5kb.dup.ids.txt || true
else  
 zcat regeno.bed.gz  \
  |awk '{if ($3-$2>=5000 && $5=="DUP")print $4}' \
  >>gt5kb.dup.ids.txt
fi  

if [ -f multi.del.ids.txt.gz ]
then
 zcat regeno.bed.gz \
  |awk '{if ($3-$2>=5000 && $5=="DEL")print $4}' \
  |fgrep -wvf <(zcat multi.del.ids.txt.gz) \
  >>gt5kb.del.ids.txt || true
else 
 zcat regeno.bed.gz \
  |awk '{if ($3-$2>=5000 && $5=="DEL")print $4}' \
  >>gt5kb.del.ids.txt
fi


zcat overlap.revise.vcf.gz \
  |fgrep -wf gt5kb.dup.ids.txt \
  >>dup.int.txt || true

zcat overlap.revise.vcf.gz \
  |fgrep -wf gt5kb.del.ids.txt \
  >>del.int.txt || true

##regenotype VCF##
dellen=$(cat del.int.txt|wc -l)
columnlen=$(less del.int.txt|cut -f10-|tr '\t' '\n' |wc -l)
dellenchange=$(echo $dellen $columnlen|awk '{if ($1 == 0) { print "0" } else { print $2/$1}}')

paste <(less del.int.txt|cut -f1-9) <(less del.int.txt|cut -f10-|tr '\t' '\n' \
 |awk -F':' '{if ($3>=2 && $1!="./.") $1="0/0"; \
  else if ($3==1 && $1!="./.") $1="0/1"; \
  else if ($1!="./.")$1="1/1";print}' OFS=":" \
  |awk -v lenchange=$dellenchange 'NR%lenchange {printf("%s\t", $0); next} \
    {print $0}')>del.revise.txt

duplen=$(cat dup.int.txt|wc -l)
columnlen=$(less dup.int.txt|cut -f10-|tr '\t' '\n' |wc -l)
duplenchange=$(echo $duplen $columnlen|awk '{if ($1 == 0) { print "0" } else { print $2/$1}}')

  
paste <(less dup.int.txt|cut -f1-9) <(less dup.int.txt|cut -f10-|tr '\t' '\n' \
 |awk -F':' '{if ($3<=2 && $1!="./.") $1="0/0"; \
  else if ($3==3 && $1!="./.") $1="0/1"; \
  else if ($1!="./.") $1="1/1";print}' OFS=":" \
  |awk -v lenchange=$duplenchange 'NR%lenchange {printf("%s\t", $0); next} \
    {print $0}') >dup.revise.txt
    

cat <(zcat overlap.revise.vcf.gz|fgrep -wvf <(cat gt5kb.dup.ids.txt gt5kb.del.ids.txt)) \
  <(cat dup.revise.txt del.revise.txt) \
  |vcf-sort \
  |bgzip \
  >newdepth.geno.vcf.gz || true


##Tag multi##
##Add filters to header##
zcat newdepth.geno.vcf.gz \
  |awk -F'\t' 'NR==FNR{inFileA[$1]; next} {if ($3 in inFileA && $1!~"#" && $7!~"PESR_GT_OVERDISPERSION") $7=$7";PESR_GT_OVERDISPERSION"; print }' OFS='\t' <(cat <(zcat $multi_geno_ids_txt) <(printf "\n")) - \
  |awk -F'\t' 'NR==FNR{inFileA[$1]; next} {if ($3 in inFileA && $1!~"#") $7=$7";MULTIALLELIC"; print }' OFS='\t' \
   <(cat <(zcat multi.del.ids.txt.gz multi.dup.ids.txt.gz |sort -u) <(printf "\n")) - \
  |sed 's\PASS;\\g' \
  |awk '{if (NR==2) print $0 "\n" "##FILTER=<ID=PESR_GT_OVERDISPERSION,Description=\"High PESR dispersion count\">" ;else print}' \
  |awk '{if (NR==2) print $0 "\n" "##FILTER=<ID=MULTIALLELIC,Description=\"Multiallelic copy number variant>" ;else print}' \
  |bgzip \
  >multitagged.vcf.gz
tabix multitagged.vcf.gz

touch all.multi.revised.list

touch dup.multi.revise.vcf
if [ $(zcat multi.dup.ids.txt.gz|wc -l) -ge 1  ]
then
  /opt/sv-pipeline/04_variant_resolution/scripts/reset_multiallelic_format_fields.py multitagged.vcf.gz <(zcat multi.dup.ids.txt.gz) > dup.multi.revise.vcf
 ${BCFTOOLS} query -f '%ID\n' dup.multi.revise.vcf >> all.multi.revised.list
fi

touch del.multi.revise.vcf
if [ $(zcat multi.del.ids.txt.gz|wc -l) -ge 1 ]
then
  /opt/sv-pipeline/04_variant_resolution/scripts/reset_multiallelic_format_fields.py multitagged.vcf.gz <(zcat multi.del.ids.txt.gz) > del.multi.revise.vcf
 ${BCFTOOLS} query -f '%ID\n' del.multi.revise.vcf >> all.multi.revised.list
fi

# make sure that the new header includes CN and CNQ format fields if we set any
if [ -s dup.multi.revise.vcf ]
then
  grep '^#' dup.multi.revise.vcf > new_header.vcf
elif [ -s  del.multi.revise.vcf ]
then
  grep '^#' del.multi.revise.vcf > new_header.vcf
else
  zcat multitagged.vcf.gz | grep '^#' > new_header.vcf
fi

# combine the revised variants with the unrevised variants, reheader, resort, and compress
 cat <(zcat multitagged.vcf.gz| \
  fgrep -wvf all.multi.revised.list) \
  <(cat del.multi.revise.vcf dup.multi.revise.vcf \
  | grep -v '^#' \
  |awk '!seen[$3]++') \
  |${BCFTOOLS} reheader -h new_header.vcf \
  |vcf-sort \
  |bgzip \
  >multitagged.geno.vcf.gz || true
         
##remove overlapping multi###
zcat multitagged.vcf.gz \
  |awk -F'\t' '{if ($1~"#" || ($7~"MULTIALLELIC" &&  ($5=="<DEL>" || $5=="<DUP>"))) print}' \
  |svtk vcf2bed stdin stdout  \
  |cut -f1-5 \
  |gzip \
  >multi.bed.gz

##strip out overlapping multiallelics##
bedtools intersect -wa -wb -a  multi.bed.gz -b  multi.bed.gz \
  |awk -F'\t' '{if ($4!=$9 && $3-$2>=$8-$7) print $0; \
  else if ($4!=$9) print $6,$7,$8,$9,$10,$1,$2,$3,$4,$5}' OFS="\t" \
  |sort -u \
  |awk '{print $3-$2,$8-$7,$0}' OFS="\t"  \
  |sort -nrk1,1 -k2,2nr \
  |cut -f3- \
  >multi.bed.overlap.txt

echo "">multi.remove.txt

while read bed
do
  echo "$bed"|cut -d$'\t' -f1-5 >large.bed
  echo "$bed"|cut -d$'\t' -f6-10>small.bed
  overlap=$(bedtools coverage -a small.bed -b large.bed|awk '{if ($NF>0.50) print "YES";else print "NO"}')  
  echo $bed|awk '{print $4}'
  if [ "$overlap" == "YES" ] && [ $(awk '{print $4}' large.bed|fgrep -wf - multi.remove.txt|wc -l) -eq 0 ]
  then
  awk '{print $4}' small.bed >>multi.remove.txt
  fi
done< multi.bed.overlap.txt

##get alt tag for multiallelics##
## produces a file with a row for each distinct multialllic variant ID and copy number combination
${BCFTOOLS} query -i 'FILTER = "MULTIALLELIC"' -f '[%ID\t%CN\n]' multitagged.geno.vcf.gz \
 |sort -u >multi.cn.txt

##strip out variants with no genotypes and overlapping multiallelics##
### Find missing genotype and then add multiallelics that need to be removed###
##change multiallelics svtype into mCNV##
##add CN information to ALT column##
zcat multitagged.geno.vcf.gz \
  |${BCFTOOLS} view -e 'FILTER == "MULTIALLELIC"'  \
  |svtk vcf2bed stdin stdout \
  |awk -F'\t' '{if ($6=="") print $4}' \
  |cat - multi.remove.txt \
  |sed '/^$/d' \
  |fgrep -wvf - <(zcat multitagged.geno.vcf.gz ) \
  |awk -F';' '{if ($1~"MULTIALLELIC" && ( $2~"DEL" || $2~"DUP")) $2="SVTYPE=CNV"; print}' OFS=';' \
  |awk '{OFS="\t"; if ($8~"SVTYPE=CNV;") $5="<CNV>"; print}' \
  |bgzip \
  >cleantagandmulti.vcf.gz || true
  
##add back original CN for sex variants which had to be changed for multiallelic##

if [ $(zcat cleantagandmulti.vcf.gz|awk '{if (($1~"X" || $1~"Y") && $1!~"#") print}'|wc -l) -gt 0 ]
then
##Determine columns male columns##
zcat cleantagandmulti.vcf.gz\
  |egrep ^# \
  |tail -n 1 \
  |tr '\t' '\n' \
  |cat -n - \
  >col.txt
  
awk '{if ($5==1) print $2}' $famfile \
  |fgrep -wf - col.txt \
   >malecols.txt || true

##regenotype male calls on sex chr and add 1 to copy state for multialleic check##
zcat cleantagandmulti.vcf.gz \
   |fgrep -wf <(grep . $sexchr_revise || true) \
   |awk  -v OFS='\t' 'NR == FNR {list[$1]; next} { for (col in list) $col="MALE"$col; print $0 }' malecols.txt - \
   |awk '{print $0 "\t" "ENDOFLINE"}' \
   |tr '\t' '\n' \
   |awk -F':' '{ if ($0!~"SVTYPE" && NF>4 && $1~"MALE" && $1!="GT" && $3-1>=0 && $3!=".") $3=$3-1;print}' OFS=":" \
   |sed 's/^MALE//g' \
   |tr '\n' '\t' \
   |sed 's/ENDOFLINE/\n/g' \
   |sed -e 's/^[ \t]*//' \
   |sed -e 's/[\t]$//g' \
   |bgzip \
   >sexchr.backtoorig.txt.gz || true

cat <(zcat cleantagandmulti.vcf.gz|fgrep -wvf <(zcat sexchr.backtoorig.txt.gz|awk '{print $3}'  )) \
  <(zcat sexchr.backtoorig.txt.gz |awk '{if ($1!="") print}' |tr ' ' '\t') \
  |vcf-sort \
  |bgzip \
  >cleansexCN.vcf.gz || true

else 
cp cleantagandmulti.vcf.gz cleansexCN.vcf.gz

fi

mv cleansexCN.vcf.gz cleanGQ.vcf.gz  

##find blank variants with no samples##
svtk vcf2bed cleanGQ.vcf.gz stdout \
 |awk -F'\t' '{if ($5!~"CN" && $6=="") print $4}' \
 >blankcheck.ids.txt

##Fix header##
##get header to clean##
##add new filters##
zcat cleanGQ.vcf.gz \
  |awk '{if ($1~"##" && NR>1)  print}' \
  |fgrep -v "MULTIALLELIC" \
  |awk '{if (NR==2) print $0 "\n" "##FILTER=<ID=MULTIALLELIC,Description=\"Multiallelic site\">" ;else print}' \
  |awk '{if (NR==2) print $0 "\n" "##ALT=<ID=CNV,Description=\"Copy Number Polymorphism\">" ;else print}' \
  |sort -k1,1 \
  |egrep -v "CIPOS|CIEND|RMSSTD|EVENT|INFO=<ID=UNRESOLVED,|source|varGQ|bcftools|ALT=<ID=UNR" \
  |cat <(zcat cleanGQ.vcf.gz|head -n 1) - <(zcat cleanGQ.vcf.gz|fgrep -wvf blankcheck.ids.txt |awk '{if ($1!~"##")  print}') \
  |bgzip >polished.vcf.gz || true
