#!/bin/bash
#
# merge_cutoffs.sh
#

set -e

batch=$1

lower_batch=$(echo $batch | awk '{print tolower($0)}')
workdir=/data/talkowski/Samples/SFARI/DeepSeq/HarrisonCompare/GT_Filtering/${lower_batch}_RF

fout=${batch}.cutoffs.txt

echo "svtype source metric cutoff direction" | sed -e 's/ /\t/g' > $fout

sed -e '1d' ${workdir}/rd.depth.del/rd.gt5000.depth.del.clean.cutoffs \
  | fgrep -v -e "NA" \
  | awk -v OFS="\t" '{print "DEL", "depth", $0}' >> $fout
echo "DEL depth svsize 5000 gt" | sed -e 's/ /\t/g' >> $fout

sed -e '1d' ${workdir}/rd.depth.dup/rd.gt5000.depth.dup.clean.cutoffs \
  | awk -v OFS="\t" '{print "DUP", "depth", $0}' >> $fout
echo "DUP depth svsize 5000 gt" | sed -e 's/ /\t/g' >> $fout

# for svtype in DEL DUP; do
for svtype in CNV; do
  sed -e '1d' ${workdir}/RDgt1kb/rd.gt1kb.nodepth.cutoffs \
    | awk -v OFS="\t" -v svtype=$svtype '{print svtype, "pesr", $0}' >> $fout
  echo "$svtype pesr svsize 1000 gt" | sed -e 's/ /\t/g' >> $fout
done

# for svtype in DEL DUP INV BND; do
for svtype in SV; do
  sed -e '1d' ${workdir}/PE/PE.cutoffs \
    | awk -v OFS="\t" -v svtype=$svtype '{print svtype, "pesr", $0}' >> $fout
  sed -e '1d' ${workdir}/SR/SR.v2.cutoffs \
    | awk -v OFS="\t" -v svtype=$svtype '{print svtype, "pesr", $0}' >> $fout
  sed -e '1d' -e 's/poiss_p/PESR_log_pval/' ${workdir}/PE_SR/PE_SR.cutoffs \
    | awk -v OFS="\t" -v svtype=$svtype '{print svtype, "pesr", $0}' >> $fout
  echo "$svtype pesr svsize 0 gt" | sed -e 's/ /\t/g' >> $fout
done
