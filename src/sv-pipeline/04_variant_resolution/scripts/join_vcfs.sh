#!/bin/bash
#
# join_vcfs.sh
#
# Join VCFs which contain the same set of variants in different samples 
#

set -euo pipefail

vcflist=$1
outfile=$2

#Make VCF header
echo -e "join_vcfs.sh :: $(date) :: MAKING VCF HEADER..."
first=$( sed -n '1p' ${vcflist} )
zcat ${first} | sed -n '1,1000p' | fgrep "##" > header.txt
zcat ${first} | sed -n '1,1000p' | fgrep "#" | fgrep -v "##" | cut -f1-9 > header2.txt
while read vcf; do
  zcat ${vcf} | sed -n '1,1000p' | fgrep "#" | fgrep -v "##" | cut -f10- | sed 's/\t/\n/g' >> header2.txt
done < ${vcflist}
cat header2.txt | paste -s - | cat header.txt - > header3.txt

###Transpose all vcfs and cat, then retranspose
#First nine columns from the first file only
echo -e "join_vcfs.sh :: $(date) :: GATHERING VARIANT INFO..."
zcat ${first} | fgrep -v "#" | cut -f1-9 \
| gawk -v FS="\t" '{ for (i=1; i<=NF; i++)  { a[NR,i] = $i } } NF>p { p = NF } END { for(j=1; j<=p; j++) { str=a[1,j]; for(i=2; i<=NR; i++){ str=str" "a[i,j]; } print str } }' \
| sed 's/\ /\t/g' \
> topnine.txt
#Transpose genotypes for each vcf
echo -e "join_vcfs.sh :: $(date) :: TRANSPOSING EACH VCF..."
nvcfs=$( cat ${vcflist} | wc -l )
i=0
while read vcf; do
  i=$(( ${i}+1 ))
  echo -e "join_vcfs.sh :: $(date) :: TRANSPOSING VCF ${i}/${nvcfs} (${vcf})..."
  zcat ${vcf} | fgrep -v "#" | cut -f10- \
  | gawk -v FS="\t" '{ for (i=1; i<=NF; i++)  { a[NR,i] = $i } } NF>p { p = NF } END { for(j=1; j<=p; j++) { str=a[1,j]; for(i=2; i<=NR; i++){ str=str" "a[i,j]; } print str } }' \
  | sed 's/\ /\t/g' \
  | gzip -c -1 \
  > transposed_${i}.txt.gz
done < ${vcflist}
#Cat all files & re-transpose
echo -e "join_vcfs.sh :: $(date) :: MERGING AND RETRANSPOSING ALL VCFS..."
i=0
while read vcf; do
  i=$(( ${i}+1 ))
  zcat transposed_${i}.txt.gz
done < ${vcflist} \
  | cat topnine.txt - \
  | gzip -c -1 \
  > all_transposed.txt.gz
rm transposed_*.txt.gz
#Shard & retranspose by chunks of variants
ncols=$( zcat all_transposed.txt.gz | awk '{ print NF }' | sed -n '1p' )
increment=5000
echo -e "join_vcfs.sh :: $(date) :: SHARDING MERGED TRANSPOSED RECORDS INTO $(( ${ncols} / ${increment} )) SHARDS..."
for start in $( seq 1 ${increment} ${ncols} ); do
  end=$(( ${start} + ${increment} - 1 ))
  echo -e "join_vcfs.sh :: $(date) :: NOW SHARDING COLUMNS ${start}-${end}"
  zcat all_transposed.txt.gz | cut -f${start}-${end} \
  | gawk -v FS="\t" '{ for (i=1; i<=NF; i++)  { a[NR,i] = $i } } NF>p { p = NF } END { for(j=1; j<=p; j++) { str=a[1,j]; for(i=2; i<=NR; i++){ str=str" "a[i,j]; } print str } }' \
  | sed 's/\ /\t/g' \
  | cat header3.txt - \
  | bgzip -c \
  > retransposed_shard_${start}_${end}.vcf.gz
done
rm all_transposed.txt.gz
#Combine all chunks
echo -e "join_vcfs.sh :: $(date) :: COMBINING ALL RETRANSPOSED SHARDS..."
vcf-concat retransposed_shard_*.vcf.gz \
  | vcf-sort \
  | bgzip -c \
  > ${outfile}

