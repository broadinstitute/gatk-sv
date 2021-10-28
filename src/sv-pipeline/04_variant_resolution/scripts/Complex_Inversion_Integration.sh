#!/bin/bash
#
# Complex_Inversion_Integration.sh
#

set -euo pipefail

##Input##
inv_resolved_vcf=$1
full_resolved_vcf=$2
##gzipped##
outputvcf=$3

##make bed of the inversion resolved vcf##
zcat $inv_resolved_vcf \
 |awk '{if ($8!~"UNRESOLVED" || $1~"#") print}' \
 |svtk vcf2bed stdin inv.resolve.bed -i MEMBERS

##make beds of the fully resolved vcf##
zcat $full_resolved_vcf \
 |awk '{if ($8~"UNRESOLVED" || $1~"#") print}' \
 |svtk vcf2bed stdin all.unresolved.inv.bed -i MEMBERS
 
zcat $full_resolved_vcf \
 |awk '{if ($8!~"UNRESOLVED" || $1~"#") print}' \
 |svtk vcf2bed stdin all.resolved.inv.bed -i MEMBERS 
 
##get unresolved variants from full vcf that are resolved in inversion resolved vcf### 
zcat $inv_resolved_vcf \
 | fgrep -v "#" \
 |awk '{if ($8!~"UNRESOLVED") print}' \
 |fgrep -wvf <(awk '{if ($NF!="MEMBERS") print $NF}' all.resolved.inv.bed  \
 |tr ',' '\n') \
 >add.vcf.lines.txt || true

##get unresolved variants id from full vcf to strip since they are resolved in inversion resolved vcf### 
##inversions that cluster were other variants (rare) are kept as unresolved though they will also be part of a resolved variant in add.vcf.lines.txt## 
awk '{if ($NF!="MEMBERS") print $NF}' inv.resolve.bed \
  |tr ',' '\n'\
  |fgrep -wf - all.unresolved.inv.bed \
  |awk '{if ($NF!~",")print $4}' \
  >remove.unresolved.vcf.ids.txt || true
 
zcat $full_resolved_vcf \
  |fgrep -wvf remove.unresolved.vcf.ids.txt \
  |cat - add.vcf.lines.txt \
  |vcf-sort -c \
  |bgzip -c \
  >$outputvcf 
