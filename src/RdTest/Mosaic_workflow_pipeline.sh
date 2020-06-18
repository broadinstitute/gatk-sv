wrkdir=/data/talkowski/Samples/SFARI/DeepSeq/HarrisonCompare/GT_Filtering/phase1_RF/mosaic_CNV/mosaic_pipeline/clean

##Mosaic CNV check##
##We currently are only going to check Mosaic CNV###
##Restrict to events private to individual exclude family structure for this particular module##
##Must come from raw caller output since sep will eliminate most mosaic calls##

##Steps##
##Pull out all unique private variants##
##For SFARI this involves only looking at parents, depthdel, depthdup have seperate cutoff so are pulled out separately from PE supported CNV##

##private to parent##
awk '{if ($1!~"depth" && ($2~"DEL" || $2~"DUP"))  print}' /data/talkowski/Samples/SFARI/DeepSeq/HarrisonCompare/GT_Filtering/phase1_RF/Phase1.all.XY.metrics|awk '{if ( $3<0.0011 && $4==0) print }'  |cat <(head -n 1 /data/talkowski/Samples/SFARI/DeepSeq/HarrisonCompare/GT_Filtering/phase1_RF/Phase1.all.XY.metrics) - >$wrkdir/parentonly.all.nodepth.metrics

awk '{if ($1~"depth" && $(NF-1)>5000 && ($2~"DUP"))  print; else if ($1~"depth" && ($2~"DUP")) print }' /data/talkowski/Samples/SFARI/DeepSeq/HarrisonCompare/GT_Filtering/phase1_RF/Phase1.all.XY.metrics|awk '{if ( $3<0.0011 && $4==0) print }'  |cat <(head -n 1 /data/talkowski/Samples/SFARI/DeepSeq/HarrisonCompare/GT_Filtering/phase1_RF/Phase1.all.XY.metrics) - >$wrkdir/parentonly.all.depthdup.metrics

awk '{if ($1~"depth" && $(NF-1)>5000 && ($2~"DEL"))  print; else if ($1~"depth" && ($2~"DEL")) print }' /data/talkowski/Samples/SFARI/DeepSeq/HarrisonCompare/GT_Filtering/phase1_RF/Phase1.all.XY.metrics|awk '{if ($3<0.0011 && $4==0) print }'  |cat <(head -n 1 /data/talkowski/Samples/SFARI/DeepSeq/HarrisonCompare/GT_Filtering/phase1_RF/Phase1.all.XY.metrics) - >$wrkdir/parentonly.all.depthdel.metrics

##Use cutoffs learned from Random Forest to pull out potential mosaic CNV, note SEP is not included as a metric##
#####must be > 5kb (difficult to confirm BAF smaller than that and removes a bunch of small unlikely calls)####

p_val_col=$(fgrep -wn RD_log_pval  <(head -n 1 $wrkdir/parentonly.all.depthdel.metrics|tr '\t' '\n')|awk -F":" '{print $1}')
secMaxP_col=$(fgrep -wn RD_log_2ndMaxP  <(head -n 1 $wrkdir/parentonly.all.depthdel.metrics|tr '\t' '\n')|awk -F":" '{print $1}')

p_val=$(fgrep -w RD_log_pval /data/talkowski/Samples/SFARI/DeepSeq/HarrisonCompare/GT_Filtering/phase1_RF/RDgt1kb/rd.gt1kb.nodepth.cutoffs|awk '{print $2}'|sed 's/NA/0/g')
secMaxP=$(fgrep -w RD_log_2ndMaxP /data/talkowski/Samples/SFARI/DeepSeq/HarrisonCompare/GT_Filtering/phase1_RF/RDgt1kb/rd.gt1kb.nodepth.cutoffs|awk '{print $2}'|sed 's/NA/0/g')

##gt1kb not depth##
awk -v p_val_col=$p_val_col -v secMaxP_col=$secMaxP_col -v p_val=$p_val -v secMaxP=$secMaxP  '{if ($p_val_col>=p_val && $secMaxP_col>=secMaxP ) print}' $wrkdir/parentonly.all.nodepth.metrics>$wrkdir/parentonly.BAFcheck.nodepth.txt

##depth only del##

p_val=$(fgrep -w RD_log_pval /data/talkowski/Samples/SFARI/DeepSeq/HarrisonCompare/GT_Filtering/phase1_RF/rd.depth.del/rd.gt5000.depth.del.clean.cutoffs|awk '{print $2}'|sed 's/NA/0/g')
secMaxP=$(fgrep -w RD_log_2ndMaxP /data/talkowski/Samples/SFARI/DeepSeq/HarrisonCompare/GT_Filtering/phase1_RF/rd.depth.del/rd.gt5000.depth.del.clean.cutoffs|awk '{print $2}'|sed 's/NA/0/g')

awk  -v p_val_col=$p_val_col -v secMaxP_col=$secMaxP_col  -v p_val=$p_val -v secMaxP=$secMaxP  '{if ($p_val_col>=p_val && $secMaxP_col>=secMaxP ) print}' $wrkdir/parentonly.all.depthdel.metrics>$wrkdir/parentonly.BAFcheck.depthdel.txt


##depth only dup##
p_val=$(fgrep -w RD_log_pval /data/talkowski/Samples/SFARI/DeepSeq/HarrisonCompare/GT_Filtering/phase1_RF/rd.depth.dup/rd.gt5000.depth.dup.clean.cutoffs|awk '{print $2}'|sed 's/NA/0/g')
secMaxP=$(fgrep -w RD_log_2ndMaxP /data/talkowski/Samples/SFARI/DeepSeq/HarrisonCompare/GT_Filtering/phase1_RF/rd.depth.dup/rd.gt5000.depth.dup.clean.cutoffs|awk '{print $2}'|sed 's/NA/0/g')

awk -v p_val_col=$p_val_col -v secMaxP_col=$secMaxP_col -v p_val=$p_val -v secMaxP=$secMaxP  '{if ( $p_val_col>=p_val && $secMaxP_col>=secMaxP ) print}' $wrkdir/parentonly.all.depthdup.metrics>$wrkdir/parentonly.BAFcheck.depthdup.txt


##create mosaic bed for RDtest denovo examination##
###filter out anything less than 5kb that does not have a high enough second max P or is found in a repetitive  region###

 cat $wrkdir/parentonly.BAF*.txt|awk '{print $1}'|fgrep -wf - /data/talkowski/Samples/SFARI/DeepSeq/HarrisonCompare/GT_Filtering/phase1_RF/Phase1.bed|sort -k1,1 -k2,2n|awk '{if ($3-$2>=5000 ) print}'|fgrep -wvf <(cat <(awk '{if ($28<2.05) print $1}' /data/talkowski/Samples/SFARI/DeepSeq/HarrisonCompare/GT_Filtering/phase1_RF/Phase1.all.XY.metrics) <(awk  '{if ($3>0.3) print $1}' /data/talkowski/Samples/SFARI/DeepSeq/HarrisonCompare/GT_Filtering/phase1_RF/filter_region_wSize.bed) ) >$wrkdir/parentonly.pass.bed

python /PHShome/hb875/talkowski/hb875/tools/svtools/svtk/cli/bedcluster.py -m -p parentonly $wrkdir/parentonly.pass.bed|sort -rk6,6|awk '!seen[$4]++'|awk '{print $1,$2,$3,$4,$6,$5}' |tr ' ' '\t' |egrep -v "^X|^Y"|fgrep -v chrom>$wrkdir/parentonly.pass.clustered.bed


##run rdtest with mosaic filter###
mkdir $wrkdir/beds
mkdir $wrkdir/rdtest_results/
cd $wrkdir/beds

shuf $wrkdir/parentonly.pass.clustered.bed|split -l 50 -a 4 -

covfile=/data/talkowski/Samples/SFARI/deep_sv/asc_540/bincov/matrices/ASC540.all.binCov.raw.bed.gz
medianfile=/data/talkowski/Samples/SFARI/deep_sv/asc_540/bincov/matrices/ASC540.all.binCov.median 
famfile=/data/talkowski/Samples/SFARI/lists/SFARI.519.fam
whitelist=/data/talkowski/Samples/SFARI/lists/SFARI.batch1.IDs.txt
outputfolder=$wrkdir/rdtest_results/

for i in x*
do
bed=$wrkdir/beds/$i
outputname=$(echo $i)
bsub -q short -o $outputfolder/$outputname.log "Rscript /PHShome/hb875/talkowski/hb875/code/RdTest_git/RdTest.R -b $bed -c $covfile -m $medianfile -f $famfile  -n $outputname -w $whitelist -o $outputfolder -z TRUE" 
done


##complete##

cat $wrkdir/rdtest_results/x*.metrics|awk '!seen[$0]++'|egrep -v "X|Y"|fgrep -v CNVID>$wrkdir/parentonly.mosaicsep.txt

##Find Sep outlier##
firstIQR=$(wc -l $wrkdir/parentonly.mosaicsep.txt|awk '{print int($1*.25)}')
thirdIQR=$(wc -l $wrkdir/parentonly.mosaicsep.txt|awk '{print int($1*.75)}')

##find cutoff 1.5*3rdIQR-1stIQR##
sepcutoff=$(echo  $(sort -nk 12,12 $wrkdir/parentonly.mosaicsep.txt|head -n $thirdIQR |tail -n 1|awk '{print $NF}') $(sort -nk 12,12   $wrkdir/parentonly.mosaicsep.txt|head -n $firstIQR|tail -n 1|awk '{print $NF}')|awk '{print $2-(($1-$2)*1.5)}')

cat $wrkdir/rdtest_results/x*.metrics|awk '!seen[$0]++'|awk -v var=$sepcutoff '{if ($NF<=var) print}'|cut -f1-6>mosaic.final.bed


