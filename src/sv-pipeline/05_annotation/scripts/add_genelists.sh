#!/bin/bash
#
# add_genelists.sh
#

set -e

coding=$1
fout=$2

# Add header
echo -e "name\teffect\tgenelist" > $fout

# First mark all genes
awk -v OFS="\t" '($2 != "INTERGENIC") {print $1, $2, "Any"}' $coding \
  | sort -k1,1 -k2,2 -u \
  | sed -e '1d' >> $fout

while read genelistname; do 
  genelist=genelists/geneSet_${genelistname}.genes.list 
  fgrep -w -f $genelist $coding \
    | awk -v OFS="\t" -v glist=$genelistname '{print $1, $2, glist}' \
    | sort -k1,1 -k2,2 -u \
    >> $fout
done < data/genelists.list
