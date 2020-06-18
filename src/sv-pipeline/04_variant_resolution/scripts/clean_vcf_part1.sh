#!/bin/bash
#
# clean_VCF.sh
#
# 
#
# Copyright (C) 2018 Harrison Brand<hbrand1@mgh.harvard.edu>
# Distributed under terms of the MIT license.

##requires >= vcftools/0.1.15 ##

set -euo pipefail


##gzipped vcf##
vcf=$1
backgroundlist=$2
famfile=$3

##get sampleids from VCF##
zcat $vcf \
  |sed -n '1,1000p' \
  |egrep "^#" \
  |tail -n -1 \
  |cut -f10- \
  |tr '\t' '\n' \
  >whitelist.txt

##convert EV integer back into string##
zcat $vcf \
      | awk '{print $0 "\t"}' \
      | sed -e 's/:7'"\t"'/:RD,PE,SR'"\t"'/g' \
      | sed -e 's/:6'"\t"'/:PE,SR'"\t"'/g' \
      | sed -e 's/:5'"\t"'/:RD,SR'"\t"'/g' \
      | sed -e 's/:3'"\t"'/:RD,PE'"\t"'/g' \
      | sed -e 's/:2'"\t"'/:PE'"\t"'/g' \
      -e 's/:4'"\t"'/:SR'"\t"'/g' \
      -e 's/:1'"\t"'/:RD'"\t"'/g' \
      |sed 's/'"\t"'$//g' \
      |sed 's/ID=EV,Number=1,Type=Integer/ID=EV,Number=1,Type=String/g' \
      | bgzip > EV.update.vcf.gz

##convert all alt to svtype and alt to N##
svtk vcf2bed EV.update.vcf.gz stdout -i SVTYPE  \
  |awk -F"\t" '{ if ($5!~"ME")$5=$7; print $4"\t" "<"$5 ">"}' \
  |gzip \
  >vcf.convert.svtype.gz

zcat EV.update.vcf.gz \
  |awk 'NR==FNR{inFileA[$1]=$2; next} {if ($3 in inFileA && $1!~"#") $5=inFileA[$3]; print }' OFS='\t'  \
   <(zcat vcf.convert.svtype.gz) - \
   |awk '{if ($1!~"#") $4="N"; print}' OFS='\t' \
   |bgzip \
   >convertsvtype.vcf.gz

##get rid of multiallelic tage in INFO field and add varGQ to QUAL column and Members field##
svtk vcf2bed convertsvtype.vcf.gz stdout -i varGQ \
  |awk -F"\t" '{print $4 "\t" $7}' \
  >vargq.persample

zcat convertsvtype.vcf.gz \
  |sed 's/;MULTIALLELIC//g' \
  |sed 's/UNRESOLVED;//g' \
  |sed 's/;varGQ=[0-9]*//g' \
  |awk 'NR==FNR{inFileA[$1]=$2; next} {if ($3 in inFileA && $1!~"#") $6=inFileA[$3]; print }' OFS='\t' vargq.persample - \
  |bgzip \
   >cleaninfo.vcf.gz

   
##fix sex chr if necessary##
if [ $(zcat cleaninfo.vcf.gz|awk '{if (($1~"X" || $1~"Y") && $1!~"#" ) print}'|wc -l) -gt 0 ]
then


svtk vcf2bed cleaninfo.vcf.gz stdout \
  |awk '{if (($5=="DEL" || $5=="DUP") && $3-$2>=5000 && ($1~"X" || $1~"Y") && $1!~"#") print}' \
  >clean.bed || true

awk '{print $4}' clean.bed>clean.bed.ids.txt


##male##
awk '{if ($5==1) print $2}' $famfile \
   |fgrep -wf <(zcat cleaninfo.vcf.gz|head -n 1000|fgrep "CHROM"|fgrep POS|cut -f10-|tr '\t' '\n') >male.txt

##female##
awk '{if ($5==2) print $2}' $famfile \
   |fgrep -wf <(zcat cleaninfo.vcf.gz|head -n 1000|fgrep "CHROM"|fgrep POS|cut -f10-|tr '\t' '\n') >female.txt

   if [ $(cat clean.bed.ids.txt|wc -l) -gt 0 ]
   then

    zcat cleaninfo.vcf.gz \
     |awk '{if ($1!~"#" && ($1~"X" || $1~"Y")) $1=$3;print}' OFS="\t"\
     |vcftools --vcf - --stdout --extract-FORMAT-info RD_CN \
     |awk -F"\t" 'NR==1{for (i=3;i<=NF;i++) header[i]=$i} NR>1{for(j=3;j<=NF;j++) print $1"\t"header[j] "\t" $j }' \
     |fgrep -wf clean.bed.ids.txt \
     |awk '{if ($3!=".") print}' \
     |gzip \
      >RD_CN.sexcheck.FORMAT.gz
     
    zcat RD_CN.sexcheck.FORMAT.gz|fgrep -wf male.txt| Rscript -e 'd<-read.table("stdin")' \
    -e 'x<-tapply(d[,3],d[,1],median)' \
    -e 'write.table(x,"male.median.value.pervar.txt",col.names=FALSE,quote=FALSE,sep = "\t")'


     zcat RD_CN.sexcheck.FORMAT.gz|fgrep -wf female.txt| Rscript -e 'd<-read.table("stdin")' \
     -e 'x<-tapply(d[,3],d[,1],median)' \
     -e 'write.table(x,"female.median.value.pervar.txt",col.names=FALSE,quote=FALSE,sep = "\t")'
    fi
##Pull out ids where male copy state 1 to normal when female normal and on X##
 echo "">sexchr.revise.txt 

 if [ $(awk '{if (($5=="DEL" || $5=="DUP") && $3-$2>=5000) print }' clean.bed|awk '{if (($1~"X") && $1!~"#" ) print}'|wc -l) -gt 0 ]
 then
 awk '{if ($2==1) print $1}' male.median.value.pervar.txt \
  |fgrep -wf <(awk '{if ($2==2) print $1}' female.median.value.pervar.txt) \
  |fgrep -wf  - <(zcat cleaninfo.vcf.gz|awk '{if ($1~"X" && $1!~"#") print $3}') \
  >sexchr.revise.txt || true
 fi

 if [ $(awk '{if (($5=="DEL" || $5=="DUP") && $3-$2>=5000) print }' clean.bed|awk '{if (($1~"Y") && $1!~"#" ) print}'|wc -l) -gt 0 ]
 then
 awk '{if ($2==1) print $1}' male.median.value.pervar.txt \
  |fgrep -wf <(awk '{if ($2==0) print $1}' female.median.value.pervar.txt) \
  |fgrep -wf - <(zcat cleaninfo.vcf.gz|awk '{if ($1~"Y" && $1!~"#") print $3}') \
  >>sexchr.revise.txt || true
 fi


bcftools index cleaninfo.vcf.gz 

##Pull out male and females sex chr##
bcftools view cleaninfo.vcf.gz -S male.txt -r chrX:1-1000000000,chrY:1-1000000000,X:1-1000000000,Y:1-1000000000 --no-update|bgzip>male.vcf.gz
bcftools view cleaninfo.vcf.gz -S female.txt -r chrX:1-1000000000,chrY:1-1000000000,X:1-1000000000,Y:1-1000000000 --no-update|bgzip>female.vcf.gz

bcftools index male.vcf.gz
bcftools index female.vcf.gz

zcat male.vcf.gz\
   |awk -F'\t' '{if ($5~"DEL" && $1!~"#") print $0 "\t" "ENDOFLINE"}' \
   |fgrep -wf sexchr.revise.txt \
   |tr '\t' '\n' \
   |awk -F':' '{if ($3>=1 && NF>4 && $1!="GT") $1="0/0";else if ($3==0 && NF>4 && $1!="GT" ) $1="0/1"; if (NF>4 && $1!="GT") $3=$3+1;print}' OFS=":" \
   |tr '\n' '\t' \
   |sed 's/ENDOFLINE/\n/g' \
   |sed -e 's/^[ \t]*//' \
   |sed -e 's/[\t]$//g' \
   |bgzip \
   >male_del.revise.txt.gz ||true


zcat male.vcf.gz\
   |awk -F'\t' '{if ($5~"DUP" && $1!~"#") print $0 "\t" "ENDOFLINE"}' \
   |fgrep -wf sexchr.revise.txt \
   |tr '\t' '\n' \
   |awk -F':' '{if ($3<=1 && NF>4 && $1!="GT") $1="0/0";else if ($3==2 && NF>4 && $1!="GT" ) $1="0/1";else if (NF>4 && $1!="GT" ) $1="1/1"; if (NF>4 && $1!="GT" ) $3=$3+1;print}' OFS=":" \
   |tr '\n' '\t' \
   |sed 's/ENDOFLINE/\n/g' \
   |sed -e 's/^[ \t]*//' \
   |sed -e 's/[\t]$//g' \
   |bgzip \
   >male_dup.revise.txt.gz ||true

  if [ $(cat male_dup.revise.txt.gz male_del.revise.txt.gz|wc -l) -gt 0 ]
  then 
   cat <(zcat male.vcf.gz|fgrep -wvf <(zcat male_dup.revise.txt.gz male_del.revise.txt.gz|awk '{print $3}' )) \
     <(zcat male_del.revise.txt.gz male_dup.revise.txt.gz|awk '{if ($1!="") print}'|tr ' ' '\t') \
    |vcf-sort \
    |bgzip \
     >cleanmale.vcf.gz
  else
    cp male.vcf.gz cleanmale.vcf.gz
  fi

 bcftools index cleanmale.vcf.gz 

  ##Modify female only for chrY###
  if [ $(zcat cleaninfo.vcf.gz |awk '{if ($1~"Y" && $1!~"#") print}'|wc -l) -gt 0 ]
  then
   zcat female.vcf.gz\
    |awk -F'\t' '{if ($1!~"#" && $1~"Y") print $0 "\t" "ENDOFLINE"}' \
    |tr '\t' '\n' \
    |awk -F':' '{ if (NF>4 && $1!="GT" ) $1="./."  \
    ;if (NF>4 && $1!="GT" ) $2=$3=$4=$5=$6=$7=$8=$9=".";print}' OFS=":" \
    |tr '\n' '\t' \
    |sed 's/ENDOFLINE/\n/g' \
    |sed -e 's/^[ \t]*//' \
    |sed -e 's/[\t]$//g' \
    |bgzip \
    >female.y.revise.txt.gz

   cat <(zcat female.vcf.gz \
    |fgrep -wvf <(zcat female.y.revise.txt.gz|awk '{print $3}' )) \
    <(zcat female.y.revise.txt.gz) \
    |vcf-sort \
    |bgzip \
    >cleanfemale.vcf.gz 

    bcftools index cleanfemale.vcf.gz 

  else 
   cp female.vcf.gz cleanfemale.vcf.gz
   bcftools index cleanfemale.vcf.gz
  fi


  ##replace genotype to ./. for other sex calls##
  ##sex anueplodies ##

  if [ $(awk '{if ($5!=2 && $5!=1) print $2}' $famfile|wc -l) -gt 0 ]
  then    
    awk '{if ($5!=2 && $5!=1) print $2}' $famfile>other.txt
    bcftools view cleaninfo.vcf.gz -S other.txt -r chrX:1-1000000000,chrY:1-1000000000,X:1-1000000000,Y:1-1000000000 --no-update|bgzip>other.vcf.gz
    bcftools index other.vcf.gz

   zcat other.vcf.gz\
    |awk -F'\t' '{if ($1!~"#") print $0 "\t" "ENDOFLINE"}' \
    |tr '\t' '\n' \
    |awk -F':' '{ if (NF>4 && $1!="GT" ) $1="./.";print}' OFS=":" \
    |tr '\n' '\t' \
    |sed 's/ENDOFLINE/\n/g' \
    |sed -e 's/^[ \t]*//' \
    |sed -e 's/[\t]$//g' \
    |bgzip \
    >other.revise.txt.gz

  cat <(zcat other.vcf.gz \
    |fgrep -wvf <(zcat other.revise.txt.gz|awk '{print $3}' )) \
    <(zcat other.revise.txt.gz) \
    |vcf-sort \
    |bgzip \
    >cleanother.vcf.gz 

  bcftools index cleanother.vcf.gz 
  
   cat <(zcat cleanmale.vcf.gz|egrep "##") \
    <(paste <(zcat cleanmale.vcf.gz|egrep -v "##") <(zcat cleanfemale.vcf.gz|cut -f10-|egrep -v "##") <(zcat cleanother.vcf.gz|cut -f10-|egrep -v "##") ) \
    |bgzip \
    >combinedsex.vcf.gz

else 
   cat <(zcat cleanmale.vcf.gz|egrep "##") \
    <(paste <(zcat cleanmale.vcf.gz|egrep -v "##") <(zcat cleanfemale.vcf.gz|cut -f10-|egrep -v "##"))  \
    |bgzip \
    >combinedsex.vcf.gz
fi



  tabix -p vcf combinedsex.vcf.gz
  tabix -p vcf cleaninfo.vcf.gz

zcat combinedsex.vcf.gz|awk '{if ($1!~"#") print $3}'>modified.ids.txt

##shuffle sex ids backinto place to match original vcf and back to initial vcf##
 vcf-shuffle-cols -t cleaninfo.vcf.gz combinedsex.vcf.gz \
  |awk '{if ($1!~"#") print}' \
  |cat <(zcat cleaninfo.vcf.gz|fgrep -wvf modified.ids.txt ) - \
  |vcf-sort \
  |bgzip \
  >cleanallo.vcf.gz

else 
 cp cleaninfo.vcf.gz cleanallo.vcf.gz
 echo "">sexchr.revise.txt
fi   

# the code below will not print any lines if the background list file is empty, so add a dummy sentinel record at the end
cat $backgroundlist <(echo "XXX_SENTINEL_XXX") > background_list_with_sentinel.list

##change tag for SR background failures and Unresolved##
zcat cleanallo.vcf.gz\
  |awk 'NR==FNR{inFileA[$1]; next} {if ($3 in inFileA && $1!~"#") $7=$7";HIGH_SR_BACKGROUND"; print }' OFS='\t' <(awk '{print $NF}' background_list_with_sentinel.list) - \
  |awk '{if (NR==2) print $0 "\n" "##FILTER=<ID=HIGH_SR_BACKGROUND,Description=\"High number of SR splits in background samples indicating messy region\">" ;else print}' \
  |awk '{if ($8~"UNRESOLVED") $7=$7";UNRESOLVED";print}' OFS='\t' \
  |awk '{if (NR==2) print $0 "\n" "##FILTER=<ID=UNRESOLVED,Description=\"Variant is unresolved\">" ;else print}' \
  |bgzip \
  >int.vcf.gz  
