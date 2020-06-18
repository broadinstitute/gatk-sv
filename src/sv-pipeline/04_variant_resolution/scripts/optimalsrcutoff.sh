#!/bin/bash

set -e

i=$1
j=$2
freqmin=$3
freqmax=$4

awk -v var=$i -v freqmin=$freqmin -v freqmax=$freqmax '{if ($4>=var && $2<=freqmax && $2>freqmin) print }' recover.single.txt >one.$i.$j.$freqmin.$freqmax.txt


awk -v var=$j -v freqmin=$freqmin -v freqmax=$freqmax '{if ($4>=var && $2<=freqmax && $2>freqmin) print }' recover.both.txt >both.$i.$j.$freqmin.$freqmax.txt


awk -v var=$i -v freqmin=$freqmin -v freqmax=$freqmax '{if ($4>=var && $2<=freqmax && $2>freqmin) print }' recover.single.fail.txt >one.$i.$j.$freqmin.$freqmax.fail.txt


awk -v var=$j -v freqmin=$freqmin -v freqmax=$freqmax '{if ($4>=var && $2<=freqmax && $2>freqmin) print }' recover.both.fail.txt >both.$i.$j.$freqmin.$freqmax.fail.txt


##passes either single and/or both side cutoff##
combine_pass=$(cat one.$i.$j.$freqmin.$freqmax.txt both.$i.$j.$freqmin.$freqmax.txt|awk '{print $5}'|sort -u|wc -l)


###no depth/pe support but passess cutoff##

combine_fail=$(cat one.$i.$j.$freqmin.$freqmax.fail.txt both.$i.$j.$freqmin.$freqmax.fail.txt|awk '{print $5}'|sort -u|wc -l)

rm one.$i.$j.$freqmin.$freqmax.txt both.$i.$j.$freqmin.$freqmax.txt
rm one.$i.$j.$freqmin.$freqmax.fail.txt both.$i.$j.$freqmin.$freqmax.fail.txt

echo $i $j $combine_pass $combine_fail >pe_support.combined.check.$i.$j.$freqmin.$freqmax.txt
