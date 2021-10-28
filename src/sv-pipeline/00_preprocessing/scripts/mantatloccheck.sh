#!/bin/bash
#
# mantatloc_check.sh
#
# Script to check raw manta calls## 
#

##standardized manta vcf##
std_vcf=$1
discfile=$2
id=$3
meibed=$4
cytobands=$5
set -e
##format manta vcf to work with svtk resolve##
cat <(zcat $std_vcf|egrep ^##) \
   <(echo "##FORMAT=<ID=manta,Number=1,Type=Integer,Description=\"manta genotype\">") \
   <(echo "##INFO=<ID=MEMBERS,Number=.,Type=String,Description=\"IDs of cluster's constituent records.\">") \
   <(echo "##INFO=<ID=EVIDENCE,Number=.,Type=String,Description=\"Classes of random forest support.\">") \
   <(zcat $std_vcf|egrep -v ^##|awk '{if ($1!~"#")$8=$8";EVIDENCE=PE;MEMBERS="$3; print}' OFS='\t') \
    |bgzip -c >manta.vcf.gz 

##INFO=<ID=EVIDENCE,Number=.,Type=String,Description="Classes of random forest support.">

tabix -p vcf manta.vcf.gz

svtk resolve manta.vcf.gz $id.manta.complex.vcf --mei-bed $meibed --cytobands $cytobands --discfile $discfile -u manta.unresolved.vcf


bgzip $id.manta.complex.vcf
