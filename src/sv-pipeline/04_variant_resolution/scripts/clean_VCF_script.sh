#!/bin/bash
#
# clean_VCF.sh
#
# 
#
# Copyright (C) 2018 Harrison Brand<hbrand1@mgh.harvard.edu>
# Distributed under terms of the MIT license.

##requires >= vcftools/0.1.15 ##

set -e

##gzipped vcf##
vcf=$1
backgroundlist=$2


##get sampleids from VCF##
zcat $vcf \
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
  |awk '{ if ($5!~"ME")$5=$7; print $4"\t" "<"$5 ">"}' \
  |gzip \
  >vcf.convert.svtype.gz

zcat EV.update.vcf.gz \
  |awk 'NR==FNR{inFileA[$1]=$2; next} {if ($3 in inFileA && $1!~"#") $5=inFileA[$3]; print }' OFS='\t'  \
   <(zcat vcf.convert.svtype.gz) - \
   |awk '{if ($1!~"#") $4="N"; print}' OFS='\t' \
   |bgzip \
   >convertsvtype.vcf.gz

##get rid of multiallelic tage in INFO field and add varGQ to QUAL column##
svtk vcf2bed convertsvtype.vcf.gz stdout -i varGQ \
  |awk '{print $4 "\t" $7}' \
  >vargq.persample

zcat convertsvtype.vcf.gz \
  |sed 's/;MULTIALLELIC//g' \
  |sed 's/;varGQ=[0-9]*//g' \
  |awk 'NR==FNR{inFileA[$1]=$2; next} {if ($3 in inFileA && $1!~"#") $6=inFileA[$3]; print }' OFS='\t' vargq.persample - \
  |bgzip \
   >cleaninfo.vcf.gz
   
   
##change tag for SR background failures##
zcat cleaninfo.vcf.gz \
  |awk 'NR==FNR{inFileA[$1]; next} {if ($3 in inFileA && $1!~"#") sub($7,"HIGH_SR_BACKGROUND"); print }' $backgroundlist - \
  |awk '{if (NR==2) print $0 "\n" "##FILTER=<ID=HIGH_SR_BACKGROUND,Description=\"High number of SR splits in background samples indicating messy region\">" ;else print}' \
  |bgzip \
  >int.vcf.gz  
 

##Remove CNVs that are improperly genotyped by depth because they are nested within a real CNV##

##Only affects CNV so pull those out##
svtk vcf2bed int.vcf.gz stdout \
  |awk '{if ($5=="DEL" || $5=="DUP") print}' \
  |gzip>int.bed.gz

##list of potenital overlaps with a normal copy state variant (>5kb variants require depth but nested events could be missed; i.e a duplication with a nest deletion will have a normal copy state for the deletion)##
bedtools intersect -wa -wb -a <(zcat int.bed.gz|awk '{if ($3-$2>=5000 ) print}') \
-b <(zcat int.bed.gz|awk '{if ($3-$2>=5000) print}') \
  |awk -F"\t" '{if ($4!=$10 && $3-$2>=$9-$8 && $5!=$11) print ;\
 else if ($4!=$10 && $5!=$11) print $7,$8,$9,$10,$11,$12,$1,$2,$3,$4,$5,$6}' OFS='\t' \
  |awk -F'\t' '{if ($6!="") print}' \
  |sort -u \
  >normaloverlap.txt


##pull out the depth based copy number variant for each normal overlapping variant## 
zcat int.vcf.gz \
  |awk '{if ($1!~"#") sub($1,$3);print}' \
  |awk '{if ($1~"#" || $5=="<DEL>" || $5=="<DUP>") print}' \
  |vcftools --vcf - --stdout --extract-FORMAT-info RD_CN \
  |awk -F"\t" 'NR==1{for (i=3;i<=NF;i++) header[i]=$i} NR>1{for(j=3;j<=NF;j++) print $1"@"header[j] "\t" $j }' \
  |sort -k1,1 \
  |gzip \
  >RD_CN.normalcheck.FORMAT.gz


##pull out evidence supporting each normal overlapping variant## 
cat <(zcat int.vcf.gz|awk -F"\t" '{if ($1~"#") print}') \
  <(awk '{print $4 "\n" $10}' normaloverlap.txt|sort -u|fgrep -wf - <(zcat int.vcf.gz))\
  |awk '{if ($1!~"#") sub($1,$3);print}' \
  |vcftools --vcf - --stdout --extract-FORMAT-info EV \
  |awk -F"\t" 'NR==1{for (i=3;i<=NF;i++) header[i]=$i} NR>1{for(j=3;j<=NF;j++) print $1"@"header[j] "\t" $j }' \
  |sort -k1,1 \
  |gzip \
  >EV.normalcheck.FORMAT.gz

##check if nested is incorrectly classified as normal##

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
  awk '{print $NF}' small.bed \
    |tr ',' '\n' \
    |fgrep -wvf - <(awk -F"[,\t]" -v var=$smallid '{for(i=6;i<=NF;i++) print var"@"$i "\t" $4"@"$i "\t" $5}' large.bed) \
    >>overlap.test.txt
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
  |awk '{if ($3=="DUP" && $4==2 && $6==3) print $2 "\t" "\t" 1; else if ($3=="DEL" && $4==2 && $6==1)  print $2 "\t" 3 }' \
  |tr '@' '\t'\
  >geno.normal.revise.txt

##Update genotypes##

##Determine columns of VCF after header##
zcat int.vcf.gz \
  |egrep ^# \
  |tail -n 1 \
  |tr '\t' '\n' \
  |cat -n - \
  >col.txt


##seed the vcf lines file which will provide the revisions to vcf file## 
echo "">normal.revise.vcf.lines.txt


##pull out and revise vcf line that needs to be edited##
while read line
do
 id=$(echo $line|awk '{print $2}' )
 col=$(fgrep -w $id col.txt|awk '{print $1}')
 variant=$(echo $line|awk '{print $1}')
 cn=$(echo $line|awk '{print $3}')

 zcat int.vcf.gz |fgrep -w $variant >line.txt

 ##Updated genotype and rebuild Format field ##
 GT=$(awk -v col=$col  '{print $col}' line.txt|awk -F":" '{print $1}')
 GQ=$(awk -v col=$col  '{print $col}' line.txt|awk -F":" '{print $2}')
 RD_CN=$(awk -v col=$col  '{print $col}' line.txt|awk -F":" '{print $3}')
 RD_GQ=$(awk -v col=$col  '{print $col}' line.txt|awk -F":" '{print $4}')
 PE_GT=$(awk -v col=$col  '{print $col}' line.txt|awk -F":" '{print $5}')
 PE_GQ=$(awk -v col=$col  '{print $col}' line.txt|awk -F":" '{print $6}')
 SR_GT=$(awk -v col=$col  '{print $col}' line.txt|awk -F":" '{print $7}')
 SR_GQ=$(awk -v col=$col  '{print $col}' line.txt|awk -F":" '{print $8}')
 EV=$(awk -v col=$col  '{print $col}' line.txt|awk -F":" '{print $9}')

 if [ $(cat normal.revise.vcf.lines.txt|fgrep -w $variant|wc -l) -gt 0 ]
 then
  cat normal.revise.vcf.lines.txt \
   |awk -v col=$col -v var=$variant -v GT=$GT -v GQ=$GQ -v RD_CN=$cn -v RD_GQ=$RD_GQ -v PE_GT=$PE_GT -v PE_GQ=$PE_GQ -v SR_GT=$SR_GT -v SR_GQ=$SR_GT -v EV=$EV '{if ($3==var ) $col="0/1:"GQ":"RD_CN":"RD_GQ":"PE_GT":"PE_GQ":"SR_GT":"SR_GQ":"EV ;print}' \
   >int.lines.txt

  cat int.lines.txt > normal.revise.vcf.lines.txt

 else 
  cat line.txt \
  |awk -v col=$col -v var=$variant -v GT=$GT -v GQ=$GQ -v RD_CN=$cn -v RD_GQ=$RD_GQ -v PE_GT=$PE_GT -v PE_GQ=$PE_GQ -v SR_GT=$SR_GT -v SR_GQ=$SR_GT -v EV=$EV '{if ($3==var ) $col="0/1:"GQ":"RD_CN":"RD_GQ":"PE_GT":"PE_GQ":"SR_GT":"SR_GQ":"EV ;print}' \
  >>normal.revise.vcf.lines.txt
  fi

done<geno.normal.revise.txt


##rewrite vcf with updated genotypes##

cat <(zcat int.vcf.gz|fgrep -wvf <(awk '{print $1}' geno.normal.revise.txt|sort -u)) \
  <(awk '{if ($1!="") print}' normal.revise.vcf.lines.txt|tr ' ' '\t') \
  |vcf-sort \
  |bgzip \
  >normal.revise.vcf.gz

##create new bed with updated genotypes###
svtk vcf2bed normal.revise.vcf.gz stdout \
  |awk '{if ($5=="DEL" || $5=="DUP") print}' \
  |sort -k4,4 \
  |gzip \
  >int.afternormalfix.bed.gz


###Find overlapping depth based variants and reassign depth based; note this is necessary because depth call >5kb genotypes are 100% driven by depth ##

## generate a sample list based on depth for depth overlap check below. Necessary because genotype is capped at 1/1 and by direction (i.e no dels in dups)##
zcat normal.revise.vcf.gz \
  |awk '{if ($1~"#" || ($5=="<DEL>" || $5=="<DUP>")) print}'\
  |awk '{if ($1!~"#") sub($1,$3);print}' \
  |vcftools --vcf - --stdout --extract-FORMAT-info RD_CN \
  |awk -F"\t" 'NR==1{for (i=3;i<=NF;i++) header[i]=$i} NR>1{for(j=3;j<=NF;j++) print $1 "\t" header[j] "\t" $j }' \
  |sort -k1,1 \
  |gzip \
  >RD_CN.afternormalfix.FORMAT.gz
  
##grab all samples per variant with a non normal copy state## 
zcat RD_CN.afternormalfix.FORMAT.gz \
  |awk '{if ($3!="2") print $1 "\t" $2}' \
  |awk '{a[$1]=a[$1]?a[$1]","$2:$2;}END{for (i in a)print i "\t" a[i];}' \
  |sort -k1,1 \
  >afternormal.combined.RD_CN.list.txt

#overlapping##
zcat int.afternormalfix.bed.gz \
  |cut -f1-5 \
  |join -1 4 -2 1 -t $'\t' - afternormal.combined.RD_CN.list.txt \
  |awk -F"\t" '{if ($6!="") print }' \
  |awk -F'[,\t]' '{for(i=6;i<=NF;i++) print $2"_"$i,$3,$4,$1,$5,$i,$1"@"$i}' \
  |tr ' ' '\t' \
  |gzip \
  >all.bed.gz


##intersect variants and always set larger to left##
bedtools intersect -wa -wb -a all.bed.gz -b all.bed.gz \
  |awk '{if ($4!=$11 && $3-$2>=$10-$9) print $0;else if ($4!=$11) print $8,$9,$10,$11,$12,$13,$14,$1,$2,$3,$4,$5,$6,$7}' \
  |tr ' ' '\t' \
  |sort -u \
  |sort -k7,7 \
  |gzip \
  >bed.overlap.txt.gz


##pull out per variant metrics from the INFO field## 

for var in EV RD_CN PE_GT SR_GT PE_GQ SR_GQ
do
 cat <(zcat normal.revise.vcf.gz|awk -F"\t" '{if ($1~"#") print}')  \
  <(zcat bed.overlap.txt.gz|awk '{print $4 "\n" $11}' |sort -u|fgrep -wf - <(zcat normal.revise.vcf.gz)) \
  |awk '{if ($1!~"#") sub($1,$3);print}' \
  |vcftools --vcf - --stdout --extract-FORMAT-info ${var} \
  |awk -F"\t" 'NR==1{for (i=3;i<=NF;i++) header[i]=$i} NR>1{for(j=3;j<=NF;j++) print $1"@"header[j] "\t" $j }' \
  |sort -k1,1 \
  |gzip \
  >${var}.FORMAT.gz
done

##Append info field to bed file##
join -1 7 -2 1 <(zcat bed.overlap.txt.gz) \
  <(zcat EV.FORMAT.gz)|join -j 1 - <(zcat RD_CN.FORMAT.gz) \
  |join -j 1 - <(zcat PE_GT.FORMAT.gz) \
  |join -j 1 - <(zcat PE_GQ.FORMAT.gz) \
  |join -j 1 - <(zcat SR_GT.FORMAT.gz) \
  |join -j 1 - <(zcat SR_GQ.FORMAT.gz) \
  |sort -k14,14 \
  |join -1 14 -2 1 - <(zcat EV.FORMAT.gz) \
  |join -j 1 - <(zcat RD_CN.FORMAT.gz) \
  |join -j 1 - <(zcat PE_GT.FORMAT.gz) \
  |join -j 1 - <(zcat PE_GQ.FORMAT.gz) \
  |join -j 1 - <(zcat SR_GT.FORMAT.gz) \
  |join -j 1 - <(zcat SR_GQ.FORMAT.gz) \
  |tr ' ' '\t' \
  |cut -f3- \
  |awk '{print $3-$2,$10-$9,$0}' \
  |tr ' ' '\t' \
  |sort -nrk1,1 -k2,2nr \
  |cut -f3- \
  |gzip \
  >all.combined.bed.gz


####If Multi-allelic is driving depth difference ignore###

##get copy state per variant##
zcat normal.revise.vcf.gz \
  |awk '{if ($1!~"#") sub($1,$3);print}' \
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
  |awk '{if ($2!="." && $2>2) print $1}' \
  |sort -u \
  |fgrep -wf <(zcat int.bed.gz|awk -F"\t" '{if ($5=="DEL" && $3-$2>=1000) print $4}' ) \
  >multi.cnvs.txt

##dup##
zcat copystate.per.variant.txt.gz \
  |awk '{if ($2!="." && ($2<2 || $2>4)) print $1}' \
  |sort -u \
  |fgrep -wf <(zcat int.bed.gz|awk -F"\t" '{if ($5=="DUP" && $3-$2>=1000) print $4}' ) \
  >>multi.cnvs.txt


##update copy state which will lead to a new genotype when genotyping is rerun towards end of script ##
echo "">RD_CN.revise.txt

while read id
do
 echo $id
 zcat all.combined.bed.gz \
  |awk -v id=$id '{if ($6==id) print $0 "\t" $4"@"$10}' \
  >overlap.bed.ids.txt

 while read bed
 do
  compareID=$(echo $bed |awk '{print $NF}')
  id1=$(echo $bed |awk '{print $4"@"$6}')
  id2=$(echo $bed |awk '{print $10"@"$12}')
  vID1=$(echo $bed |awk '{print $4}')
  vID2=$(echo $bed |awk '{print $10}')
  svtype1=$(echo $bed |awk '{print $5}')
  svtype2=$(echo $bed |awk '{print $11}')
  support1=$(echo $bed |awk '{print $13}')
  support2=$(echo $bed |awk '{print $19}')
  length1=$(echo $bed|awk '{print $3-$2}')
  length2=$(echo $bed|awk '{print $9-$8}')
  RD_CN1=$(echo $bed|awk '{print $14}')
  RD_CN2=$(echo $bed|awk '{print $20}')
  PE_GT1=$(echo $bed|awk '{print $15}')
  PE_GT2=$(echo $bed|awk '{print $21}')
  PE_GQ1=$(echo $bed|awk '{print $16}')
  PE_GQ2=$(echo $bed|awk '{print $22}')
  SR_GT1=$(echo $bed|awk '{print $17}')
  SR_GT2=$(echo $bed|awk '{print $23}')
  SR_GQ1=$(echo $bed|awk '{print $18}')
  SR_GQ2=$(echo $bed|awk '{print $24}')

  echo $bed|tr ' ' '\t'|cut -f1-6 >large.bed
  echo $bed|tr ' ' '\t'|cut -f7-12>small.bed
  overlap=$(bedtools coverage -a small.bed -b large.bed|awk '{if ($NF>0.50) print "YES";else print "NO"}')  
  
  ##remove any large CNV comparisons that have been revised to normal copy state of 2##
  awk '{if ($2==2) print $1}' RD_CN.revise.txt>depthnormal.exclude.txt

  if [ $(fgrep -w $id1 depthnormal.exclude.txt |wc -l) -eq 0 ] 
  then 
  ##classification##
   ##Call where smaller depth call is being driven by larger##
   if [[ $support1 =~ "RD" ]] && [[ $support2 = "RD" ]] && [ "$overlap" == "YES"  ] && [[ $support1 != "RD" ]] && [ $(fgrep -w $vID1 multi.cnvs.txt |wc -l) -eq 0 ] 
   then
    echo $bed \
     |awk -v id2=$id2 -v svtype1=$svtype1 -v RD_CN1=$RD_CN1 -v RD_CN2=$RD_CN2  '{if (RD_CN1==1) print id2   "\t" RD_CN2+RD_CN1 ; \
     else if(RD_CN1>1) print id2  "\t"  RD_CN2-(RD_CN1-2) }' \
     >>RD_CN.revise.txt
   ##Smaller CNV driving larger CNV genotype##
   elif [[ $support1 = "RD"  ]] && [[ $support2 =~ "RD"  ]] && [ "$overlap" == "YES"  ] && [[ $support2 != "RD" ]] && [ $(fgrep -w $vID2 multi.cnvs.txt |wc -l) -eq 0 ] 
   then
    echo $bed \
     |awk -v id1=$id1 -v svtype1=$svtype1 -v RD_CN1=$RD_CN1 -v RD_CN2=$RD_CN2  '{if (RD_CN2==1) print id1   "\t" RD_CN1+RD_CN2 ; \
     else if(RD_CN2>1) print id1  "\t"  RD_CN1-(RD_CN2-2) }' \
     >>RD_CN.revise.txt
   ##Depth only calls where smaller call is being driven by larger##
   elif [[ $support1 = "RD" ]]  && [[ $support2 = "RD" ]] && [ "$overlap" == "YES"  ] && [  "$svtype1" == "$svtype2"  ]  && [ $(fgrep -w $vID1 multi.cnvs.txt |wc -l) -eq 0 ] 
   then
    echo $bed \
     |awk -v id2=$id2 -v svtype1=$svtype1 -v RD_CN1=$RD_CN1 -v RD_CN2=$RD_CN2 '{if (RD_CN1==1 && RD_CN1>RD_CN2  ) print id2 "\t" 1; \
      else if (RD_CN1>1 && RD_CN1<RD_CN2 ) print id2  "\t"  RD_CN2-(RD_CN1-2) ; \
      else print id2   "\t" 2 }' \
      >>RD_CN.revise.txt
   ##Any other time a larger call is driving a smaller call##
   elif [[ $support1 =~ "RD" ]]  && [ "$overlap" == "YES"  ] && [ $length2 -gt 5000  ] && [ $(fgrep -w $vID1 multi.cnvs.txt |wc -l) -eq 0 ] 
   then
    echo $bed \
    |awk -v id2=$id2 -v svtype1=$svtype1 -v RD_CN1=$RD_CN1 -v RD_CN2=$RD_CN2  '{if (RD_CN1==1) print id2   "\t" RD_CN2+RD_CN1 ; \
    else if(RD_CN1>1) print id2  "\t"  RD_CN2-(RD_CN1-2) }' \
    >>RD_CN.revise.txt
   fi
  fi
done<overlap.bed.ids.txt
done<whitelist.txt

##fix negative copy states that will rarely arise and assign them 0 copy state##
sed '/^$/d' RD_CN.revise.txt \
  |tr '@' '\t' \
  |awk '{if ($NF<0) print $1 "\t" $2 "\t" 0; else print}' \
  >RD_CN.revise.forgeno.txt

##Determine columns of VCF after header##
zcat normal.revise.vcf.gz\
  |egrep ^# \
  |tail -n 1 \
  |tr '\t' '\n' \
  |cat -n - \
  >col.txt


##seed the vcf lines file which will provide the revisions to vcf file## 
echo "">revise.vcf.lines.txt


##pull out and revise vcf line that needs to be edited##
while read line
do
 id=$(echo $line|awk '{print $2}' )
 col=$(fgrep -w $id col.txt|awk '{print $1}')
 variant=$(echo $line|awk '{print $1}')
 cn=$(echo $line|awk '{print $3}')

 zcat normal.revise.vcf.gz |fgrep -w $variant >line.txt

 echo $variant $id
 ##Updated genotype and rebuild Format field ##
 GT=$(awk -v col=$col  '{print $col}' line.txt|awk -F":" '{print $1}')
 GQ=$(awk -v col=$col  '{print $col}' line.txt|awk -F":" '{print $2}')
 RD_CN=$(awk -v col=$col  '{print $col}' line.txt|awk -F":" '{print $3}')
 RD_GQ=$(awk -v col=$col  '{print $col}' line.txt|awk -F":" '{print $4}')
 PE_GT=$(awk -v col=$col  '{print $col}' line.txt|awk -F":" '{print $5}')
 PE_GQ=$(awk -v col=$col  '{print $col}' line.txt|awk -F":" '{print $6}')
 SR_GT=$(awk -v col=$col  '{print $col}' line.txt|awk -F":" '{print $7}')
 SR_GQ=$(awk -v col=$col  '{print $col}' line.txt|awk -F":" '{print $8}')
 EV=$(awk -v col=$col  '{print $col}' line.txt|awk -F":" '{print $9}')

 if [ $(cat revise.vcf.lines.txt|fgrep -w $variant|wc -l) -gt 0 ]
 then
  cat revise.vcf.lines.txt \
   |awk -v col=$col -v var=$variant -v GT=$GT -v GQ=$GQ -v RD_CN=$cn -v RD_GQ=$RD_GQ -v PE_GT=$PE_GT -v PE_GQ=$PE_GQ -v SR_GT=$SR_GT -v SR_GQ=$SR_GT -v EV=$EV '{if ($3==var ) $col="0/1:"GQ":"RD_CN":"RD_GQ":"PE_GT":"PE_GQ":"SR_GT":"SR_GQ":"EV ;print}' \
   >int.lines.txt

  cat int.lines.txt > revise.vcf.lines.txt

 else 
  cat line.txt \
  |awk -v col=$col -v var=$variant -v GT=$GT -v GQ=$GQ -v RD_CN=$cn -v RD_GQ=$RD_GQ -v PE_GT=$PE_GT -v PE_GQ=$PE_GQ -v SR_GT=$SR_GT -v SR_GQ=$SR_GT -v EV=$EV '{if ($3==var ) $col="0/1:"GQ":"RD_CN":"RD_GQ":"PE_GT":"PE_GQ":"SR_GT":"SR_GQ":"EV ;print}' \
  >>revise.vcf.lines.txt
  fi

done<RD_CN.revise.forgeno.txt

cat <(zcat normal.revise.vcf.gz|fgrep -wvf <(awk '{print $1}' RD_CN.revise.forgeno.txt|sort -u)) \
  <(awk '{if ($1!="") print}' revise.vcf.lines.txt|tr ' ' '\t') \
  |vcf-sort \
  |bgzip \
  >overlap.revise.vcf.gz

##multi check##
zcat overlap.revise.vcf.gz \
  |awk '{if ($1!~"#") sub($1,$3);print}' \
  |vcftools --vcf - --stdout --extract-FORMAT-info RD_CN \
  |gzip \
  >copystate.RD_CN.FORMAT.gz

zcat copystate.RD_CN.FORMAT.gz \
  |awk 'NR>1{for(i=3;i<=NF;i++) lines[$1 "\t" $i]++ } END{for (x in lines) print x}' \
  |gzip \
  >copystate.per.variant.txt.gz

##Copy state just del and dup ; CNV >1kb we trust depth ##
zcat copystate.per.variant.txt.gz \
  |awk '{if ($2!="." && $2>2) print $1}' \
  |sort -u \
  |fgrep -wf <(zcat int.bed.gz|awk -F"\t" '{if ($5=="DEL" && $3-$2>=1000) print $4}' ) \
  |gzip \
  >multi.del.ids.txt.gz

zcat copystate.per.variant.txt.gz \
  |awk '{if ($2!="." && ($2<2 || $2>4)) print $1}' \
  |sort -u \
  |fgrep -wf <(zcat int.bed.gz|awk -F"\t" '{if ($5=="DUP" && $3-$2>=1000) print $4}' ) \
  |gzip \
  >multi.dup.ids.txt.gz 

##Regenotype to determine multiallelic##
##Genotype big dup##
svtk vcf2bed overlap.revise.vcf.gz stdout \
  |gzip>regeno.bed.gz

##generate list##
##CNV >5kb, split del and dup ##
##  ##
zcat regeno.bed.gz  \
 |awk '{if ($3-$2>=5000 && $5=="DUP")print $4}' \
 |fgrep -wvf <(zcat multi.dup.ids.txt.gz)  \
 >gt5kb.dup.ids.txt
 
zcat regeno.bed.gz \
  |awk '{if ($3-$2>=5000 && $5=="DEL")print $4}' \
  |fgrep -wvf <(zcat multi.del.ids.txt.gz) \
  >gt5kb.del.ids.txt

end=$(zcat overlap.revise.vcf.gz|awk '{if ($1!~"#") print}'|head -n 1 |awk -F'[:\t]' '{print NF}' )

zcat overlap.revise.vcf.gz \
  |fgrep -wf gt5kb.dup.ids.txt \
  >dup.int.txt

zcat overlap.revise.vcf.gz \
  |fgrep -wf gt5kb.del.ids.txt \
  >del.int.txt

##regenotype VCF##
for ((i=18;i<=$end;i+=9))
do
 echo $i
 cat dup.int.txt \
  |awk -F'[:\t]' -v i=$i '{if ($(i+2)==2) sub($i,"0/0"); \
  else if ($(i+2)==3) sub($i,"0/1"); \
  else sub($i,"1/1");print}' \
  >dup.revise.txt
  
 cat del.int.txt \
  |awk -F'[:\t]' -v i=$i '{if ($(i+2)==2) sub($i,"0/0"); \
  else if ($(i+2)==1) sub($i,"0/1"); \
  else sub($i,"1/1");print}' \
  >del.revise.txt
  
  cat dup.revise.txt>dup.int.txt
  cat del.revise.txt>del.int.txt
done

cat <(zcat overlap.revise.vcf.gz|fgrep -wvf <(cat gt5kb.dup.ids.txt gt5kb.del.ids.txt)) \
  <(cat dup.revise.txt del.revise.txt) \
  |vcf-sort \
  |bgzip \
  >newdepth.geno.vcf.gz


##Tag VCF##
##find individual level metrics to determine multi allelic by PE/SR genotypes##

for var in PE_GT SR_GT PE_GQ SR_GQ
do
 zcat newdepth.geno.vcf.gz \
  |awk '{if ($1!~"#") sub($1,$3);print}' \
  |vcftools --vcf - --stdout --extract-FORMAT-info ${var} \
  |awk -F"\t" 'NR==1{for (i=3;i<=NF;i++) header[i]=$i} NR>1{for(j=3;j<=NF;j++) print $1"@"header[j] "\t" $j }' \
  |sort -k1,1 \
  |gzip \
  >multicheck.${var}.FORMAT.gz
done

##concatenate metrics##
join -j 1 <(zcat multicheck.PE_GT.FORMAT.gz) \
  <(zcat multicheck.PE_GQ.FORMAT.gz) \
  |join -j 1  - <(zcat multicheck.SR_GT.FORMAT.gz) \
  |join -j 1  - <(zcat multicheck.SR_GQ.FORMAT.gz) \
  |tr ' ' '\t' \
  |gzip \
  >multi.combined.format.gz


##check by genotype##
zcat multi.combined.format.gz \
  |awk '{if ($2>0 && $4==0) print $1"\t" $2; \
  else if ($2==0) print $1 "\t" $4; \
  else if ($3>=$5)print $1"\t" $2; \
  else print $1"\t" $4 }' \
  |tr '@' '\t' \
  |awk '{if ($3>2 && $2!=".") print $1}' \
  |sort -u \
  |gzip \
  >multi.geno.ids.txt.gz 

##Tag multi##
zcat newdepth.geno.vcf.gz \
  |awk 'NR==FNR{inFileA[$1]; next} {if ($3 in inFileA && $1!~"#") sub($7,"MULTIALLELIC"); print }' OFS='\t' \
   <(zcat multi.del.ids.txt.gz multi.dup.ids.txt.gz multi.geno.ids.txt.gz|sort -u) - \
   |bgzip \
   >multitagged.vcf.gz

###genotype multiallelics##
##pull out multiallelic lines of vcf###
zcat multitagged.vcf.gz \
   |fgrep -wf <(zcat multi.geno.ids.txt.gz) \
   >multi.gt.int.txt

zcat multitagged.vcf.gz \
   |fgrep -wf <(zcat multi.dup.ids.txt.gz) \
   >multi.dup.int.txt

zcat multitagged.vcf.gz \
   |fgrep -wf <(zcat multi.del.ids.txt.gz) \
   >multi.del.int.txt  

end=$(zcat multitagged.vcf.gz|awk '{if ($1!~"#") print}'|head -n 1 |awk -F'[:\t]' '{print NF}' )

for ((i=18;i<=$end;i+=9))
do
 echo $i
 cat multi.dup.int.txt \
  |awk -F'[:\t]' -v i=$i '{sub($i,"./"$(i+2));print}' \
  >dup.multi.revise.txt
  
 cat multi.del.int.txt \
  |awk -F'[:\t]' -v i=$i '{sub($i,"./"$(i+2));print}' \
  >del.multi.revise.txt
  
 cat multi.gt.int.txt \
  |awk -F'[:\t]' -v i=$i '{if ($(i+4)>0 && $(i+6)==0)  sub($i,"./"$(i+4)); \
  else if ($(i+4)==0)  sub($i,"./"$(i+6)); \
  else if ($(i+5)>=$(i+7)) sub($i,"./"$(i+4)); \
  else  sub($i,"./"$(i+6)) ;print }' \
  >gt.multi.revise.txt
  
  cat dup.multi.revise.txt>multi.dup.int.txt
  cat del.multi.revise.txt>multi.del.int.txt
  cat gt.multi.revise.txt>multi.gt.int.txt
done

##remove overlapping multi###
zcat multitagged.vcf.gz \
  |awk '{if ($1~"#" || ($7=="MULTIALLELIC" &&  ($5=="<DEL>" || $5=="<DUP>"))) print}' \
  |svtk vcf2bed stdin stdout  \
  |gzip \
  >multi.bed.gz

##strip out overlapping multiallelics##
bedtools intersect -wa -wb -a  multi.bed.gz -b  multi.bed.gz \
  |awk '{if ($4!=$10 && $3-$2>=$9-$8) print $0; \
  else if ($4!=$10) print $7,$8,$9,$10,$11,$12,$1,$2,$3,$4,$5,$6}' \
  |tr ' ' '\t' \
  |sort -u \
  |awk '{print $3-$2,$9-$8,$0}' \
  |tr ' ' '\t' \
  |sort -nrk1,1 -k2,2nr \
  |cut -f3- \
  |gzip \
  >multi.bed.overlap.txt.gz

echo "">multi.remove.txt

while read bed
do
  echo $bed|tr ' ' '\t'|cut -f1-6 >large.bed
  echo $bed|tr ' ' '\t'|cut -f7-12>small.bed
  overlap=$(bedtools coverage -a small.bed -b large.bed|awk '{if ($NF>0.50) print "YES";else print "NO"}')  
  
  if [ "$overlap" == "YES" ] && [ $(awk '{print $4}' large.bed|fgrep -wf - multi.remove.txt|wc -l) -eq 0 ]
  then
  awk '{print $4}' small.bed >>multi.remove.txt
  fi
done< <(zcat multi.bed.overlap.txt.gz)


##strip out variants with no genotypes and overlapping multiallelics##
### Find missing genotype and then add multiallelics that need to be removed###

svtk vcf2bed multitagged.vcf.gz stdout \
  |awk -F'\t' '{if ($6=="") print $4}' \
  |cat - multi.remove.txt \
  |sed '/^$/d' \
  |fgrep -wvf - <(zcat multitagged.vcf.gz) \
  |gzip \
  >cleantagandmulti.vcf.gz

##Fix header##
##get header to clean##
##add new filters##
zcat cleantagandmulti.vcf.gz \
  |awk '{if (NR==2) print $0 "\n" "##FILTER=<ID=MULTIALLELIC,Description=\"Multiallelic site\">" ;else print}' \
  |awk '{if ($1~"##" && NR>1)  print}' \
  |sort -k1,1 \
  |egrep -v "CIPOS|CIEND|RMSSTD|MEMBERS|UNRESOLVED|source|MULTIALLELIC|varGQ|bcftools|ALT=<ID=UNR" \
  |cat <(zcat cleantagandmulti.vcf.gz|head -n 1) - <(zcat cleantagandmulti.vcf.gz|awk '{if ($1!~"##")  print}') \
  |gzip \
  >polished.vcf.gz

