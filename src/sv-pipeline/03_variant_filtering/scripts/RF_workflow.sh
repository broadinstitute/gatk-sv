#!/bin/bash

wrkdir=/PHShome/hb875/SFARI/DeepSeq/HarrisonCompare/GT_Filtering/clean_pipeline
cd $wrkdir

##see earler versions for creating required files##

##Generate initial datasets##
cat /data/talkowski/Samples/SFARI/deep_sv/asc_540/sv-pipeline/rdtest/split_beds/Phase1*>$wrkdir/Phase1.bed

##Remove regions with poor sequencing just for training of data and get size to add to metrics##
##remove X && Y###
module load bedtools2/2.25.0
egrep -hv "^X|^Y" /PHShome/hb875/SFARI/DeepSeq/HarrisonCompare/GT_Filtering/Phase1.bed| sort -k1,1 -k2,2n| coverageBed -a - -b /data/talkowski/rlc47/src/GRch37.segdups_gaps_abParts_heterochrom.lumpy.exclude.bed -sorted |awk '{print  $4 "\t" $3-$2 "\t" $NF}'|cat <(echo -e "name" '\t' "size" '\t'"poor_region_cov") ->filter_region_wSize.bed

##breakpoint check for filtering##

cat /data/talkowski/rlc47/src/GRCh37.*.RMSK.merged.bed /data/talkowski/rlc47/src/GRch37.segdups_gaps_abParts_heterochrom.lumpy.exclude.bed|sort -k1,1 -k2,2n>breakpoint.rmsk.bed

egrep -hv "^X|^Y" /PHShome/hb875/SFARI/DeepSeq/HarrisonCompare/GT_Filtering/Phase1.bed|awk '{print $1,$2,$2+1,$4"\n"$1,$3-1,$3,$4}'| tr ' ' '\t'|sort -k1,1 -k2,2n|bedtools intersect -sorted -a - -b breakpoint.rmsk.bed|awk '{print $4}'|sort -u>one_end_linesine_txt

##get rid of parent only variants when checking  inhrate##

for vcf in /data/talkowski/Samples/SFARI/deep_sv/asc_540/sv-pipeline/01_algorithm_integration/vcfcluster/Phase1.*vcf
do
/data/talkowski/Samples/SFARI/deep_sv/asc_540/sv-pipeline/01_algorithm_integration/scripts/get_called_samples.py $vcf >>SV.samples
done


##create overall metrics file with all metrics and variant types## 
cat /data/talkowski/Samples/SFARI/deep_sv/asc_540/sv-pipeline/03_genotyping/metrics/Phase1.*.metrics| \
    ##remove duplicate lines##	\
    awk '!seen[$0]++'| \
    ##create metric file remove variants where all samples are called to have a CNV## \
    fgrep  -v "All_samples_called_CNV_no_analysis"| \
    ####fix BAF p-value to be on log scale to match others### \
    awk '{if ($34!="NA" && $34!="BAF_KS_pval") {$34=-log($34)/log(10)} print}' | \
    ####flip sign for BAF DEL log likiehood so postive is more significant ### \
    awk '{if ($31!="NA" && $31!="BAF_del_loglik") {$31=-$31} print}' | \
    ##replace inf from -log10 p-value with max -log10(P) of 300## \
    sed 's/inf/300/g' | \
    ##add poor region coverage and size as metrics at the end## \
    sort -k1,1 |join -a 1 -j 1 - filter_region_wSize.bed | \
    ##add size and coverage NA for none CNV sv types## \
    awk '{if (NF==41) print $0,"NA","NA";else print}'| \
    ##remove chr X and Y \
    awk -F"_" '{if ($3!="X" && $3!="Y" && $4!="X" && $4!="Y" ) print}' | \
    ##get rid of straggler header line## \
    tr ' ' '\t'|sort -k1,1|tail -n +2> $wrkdir/Phase1.all.metrics


##remove BAF failures (variants which BAF can not be assessed due to # of snps or ROH)##
cat $wrkdir/Phase1.all.metrics|fgrep DEL|awk '{if ($30 =="NA" || $31 =="NA") print $1 }'>$wrkdir/BAF/Phase1.BAF.noassess

cat $wrkdir/Phase1.all.metrics|fgrep DUP|awk '{if ($33 =="NA" || $34 =="NA") print $1}'>>$wrkdir/BAF/Phase1.BAF.noassess

###Start with BAF##
##Build Training set##
##del restrict to greater than 5 kb for training and in clean regions##
mkdir $wrkdir/BAF

fgrep DEL $wrkdir/Phase1.all.metrics|fgrep -wvf $wrkdir/BAF/Phase1.BAF.noassess|awk '{if ($NF<0.3 && $(NF-1)>=5000 && $26<0.15) print $1,"Fail",$30,$31; else if ($NF<0.3 && $(NF-1)>=5000 && $26>0.4) print $1,"Pass",$30,$31 }'|cat <(awk '{print "name","Status",$30,$31}' $wrkdir/Phase1.all.metrics|head -n 1) - |tr ' ' '\t'>$wrkdir/BAF/BAF.del.metrics

fgrep DEL $wrkdir/Phase1.all.metrics|awk '{if ($(NF-1)>=5000 ) print}'|fgrep -wvf $wrkdir/BAF/Phase1.BAF.noassess|cat <(head -n 1 $wrkdir/Phase1.all.metrics) - >$wrkdir/BAF/BAF.del.all.metrics

fgrep DUP $wrkdir/Phase1.all.metrics|fgrep -wvf $wrkdir/BAF/Phase1.BAF.noassess|awk '{if ($NF<0.3 && $(NF-1)>=5000 && $26<0.15) print $1,"Fail",$33,$34; else if ($NF<0.3 && $(NF-1)>=5000 && $26>0.4) print $1,"Pass",$33,$34 }'|cat <(awk '{print "name","Status",$33,$34}' $wrkdir/Phase1.all.metrics|head -n 1) - |tr ' ' '\t'>$wrkdir/BAF/BAF.dup.metrics

fgrep DUP $wrkdir/Phase1.all.metrics|awk '{if ($(NF-1)>=5000 ) print}'|fgrep -wvf $wrkdir/BAF/Phase1.BAF.noassess|cat <(head -n 1 $wrkdir/Phase1.all.metrics) - >$wrkdir/BAF/BAF.dup.all.metrics

##create file for ROC cut_off structure##
##Dels first filter on BAF_snp_ratio then BAF_del_loglik##
echo "BAF_snp_ratio">$wrkdir/BAF/ROC.del.indep.txt 
echo "BAF_del_loglik">$wrkdir/BAF/ROC.del.dep.txt 

##Dup first filter on BAF_KS_stat then BAF_KS_pval##
echo "BAF_KS_stat">$wrkdir/BAF/ROC.dup.indep.txt 
echo "BAF_KS_pval">$wrkdir/BAF/ROC.dup.dep.txt 

##Run model##
##deletions##
Rscript /PHShome/hb875/talkowski/hb875/code/RF_git/RandomForest.R $wrkdir/BAF/BAF.del.metrics $wrkdir/BAF/BAF.del.all.metrics 1343124 $wrkdir/BAF/BAF.del $wrkdir/BAF/ROC.del.indep.txt $wrkdir/BAF/ROC.del.dep.txt 

##duplications##
Rscript /PHShome/hb875/talkowski/hb875/code/RF_git/RandomForest.R $wrkdir/BAF/BAF.dup.metrics $wrkdir/BAF/BAF.dup.all.metrics 1343124 $wrkdir/BAF/BAF.dup $wrkdir/BAF/ROC.dup.indep.txt $wrkdir/BAF/ROC.dup.dep.txt 

##SR##
mkdir $wrkdir/SR

##require to pass both initial RD and BAF to be considered valid and fail both to be invalid###
##remove events as failing that have potential PE support (p<0.05)##
join -j 1 <(cat $wrkdir/BAF/BAF.del.metrics $wrkdir/BAF/BAF.dup.metrics|cut -f -1,2|sort -k1,1 ) <(cat $wrkdir/BAF/BAF.del.pred $wrkdir/BAF/BAF.dup.pred|sort -k1,1)|egrep -v "name|depth"|awk '{if ($2=="Fail" && $3<0.5) print $1"\t"$3;else if ($2=="Pass" && $3>=0.5) print $1"\t"$3}'|awk '{if ($2>=0.9) print $1 "\t" "Pass";else if ($2<=0.1) print $1 "\t" "Fail"}'|sort -k1,1|join -j 1 - <(awk '{print $1,$11,$14,$17,$6}' $wrkdir/Phase1.all.metrics|sort -k1,1)|awk '{if ($2=="Fail" && $NF<1.30103) print $1,$2,$3,$4,$5 ;else if ($2=="Pass") print $1,$2,$3,$4,$5}'|cat <(awk '{print $1,"Status",$11,$14,$17}' $wrkdir/Phase1.all.metrics|head -n 1) - |tr ' ' '\t'>$wrkdir/SR/SR.metrics

##print all variants to be assessed by SR metrics (exclude depth only and wham INV and BND), NOTE: redundant with PE##

awk '{if ($1!~"depth" && !($1~"wham" && ($2~"INV" || $2~"BND"))) print}' $wrkdir/Phase1.all.metrics>$wrkdir/SR/SR.all.metrics

##Dup first filter on SR_sum_log_pval then SR_sum_bg_median##
echo "SR_sum_log_pval" "\t" "gt" >$wrkdir/SR/ROC.indep.txt 
echo "" >$wrkdir/SR/ROC.dep.txt 

## Build SR random forest model and test all variants##

Rscript /PHShome/hb875/talkowski/hb875/code/RF_git/RandomForest.R $wrkdir/SR/SR.metrics $wrkdir/SR/SR.all.metrics 1343124 $wrkdir/SR/SR $wrkdir/SR/ROC.indep.txt $wrkdir/SR/ROC.dep.txt  

##RD test##
##gt1kb & excluding depth based calls##
mkdir $wrkdir/RDgt1kb

##require PE or SR passing and BAF for training set, restrict data set to >1kb##
cat $wrkdir/SR/SR.pred|cat - $wrkdir/BAF/BAF.del.pred $wrkdir/BAF/BAF.dup.pred|awk '{if ($2>=0.9) print $1 "\t" "Pass";else if ($2<=0.1) print $1 "\t" "Fail"}'|sort|uniq -c|awk '{if ($1==2) print $2 "\t" $3}'|sort -k1,1|join -j 1 - <(cat $wrkdir/Phase1.all.metrics|awk '{if ($(NF-1)>5000  && $NF<0.3) print $1,$26,$27,$28}'|sort -k1,1 )|cat <(awk '{print $1,"Status",$26,$27,$28}' $wrkdir/Phase1.all.metrics|head -n 1) -|tr ' ' '\t'>$wrkdir/RDgt1kb/rd.gt1kb.nodepth.metrics

##print all variants that will be tested excluding depth based calls##
awk '{if ($(NF-1)>1000) print }' $wrkdir/Phase1.all.metrics|fgrep -v depth|egrep -w "DEL|DUP"|cat <(awk '{print}' $wrkdir/Phase1.all.metrics|head -n 1) -|tr ' ' '\t'>$wrkdir/RDgt1kb/rd.all.gt1kb.nodepth.metrics

##Dup first filter on RD_Median_Separation and RD_log_pval then RD_log_2ndMaxP##
echo -e "RD_Median_Separation""\n""RD_log_pval"  >$wrkdir/RDgt1kb/ROC.indep.txt 
echo "RD_log_2ndMaxP">$wrkdir/RDgt1kb/ROC.dep.txt 


Rscript /PHShome/hb875/talkowski/hb875/code/RF_git/RandomForest.R $wrkdir/RDgt1kb/rd.gt1kb.nodepth.metrics $wrkdir/RDgt1kb/rd.all.gt1kb.nodepth.metrics 1343124 $wrkdir/RDgt1kb/rd.gt1kb.nodepth $wrkdir/RDgt1kb/ROC.indep.txt $wrkdir/RDgt1kb/ROC.dep.txt  

##lt1kb##
mkdir $wrkdir/RDlt1kb

##require SR support and restrict to variants greater than 100 bp and less than 1000, restrict to clean regions##

awk '{if ($2<=0.1) print $1 "\t" "Fail";else if ($2>=0.9 ) print $1 "\t" "Pass"}' $wrkdir/SR/SR.pred|fgrep -v name|sort -k1,1|join -j 1 - <(awk '{if ($(NF-1)<=1000 && $(NF-1)>100 && $NF<0.3 ) print $1,$26,$27,$28}' $wrkdir/Phase1.all.metrics|sort -k1,1)|cat <(awk '{print $1,"Status",$26,$27,$28}' $wrkdir/Phase1.all.metrics|head -n 1) -|tr ' ' '\t'>$wrkdir/RDlt1kb/rd.lt1kb.nodepth.metrics

egrep -w "DEL|DUP" $wrkdir/Phase1.all.metrics|fgrep -v depth|awk '{if ($(NF-1)<=1000) print }'| cat <( head -n 1 $wrkdir/Phase1.all.metrics) - >$wrkdir/RDlt1kb/rd.all.lt1kb.nodepth.metrics

##Dup first filter on RD_Median_Separation and RD_log_pval then RD_log_2ndMaxP##
echo -e "RD_Median_Separation""\n""RD_log_pval"  >$wrkdir/RDlt1kb/ROC.indep.txt 
echo "RD_log_2ndMaxP">$wrkdir/RDlt1kb/ROC.dep.txt 

Rscript /PHShome/hb875/talkowski/hb875/code/RF_git/RandomForest.R $wrkdir/RDlt1kb/rd.lt1kb.nodepth.metrics $wrkdir/RDlt1kb/rd.all.lt1kb.nodepth.metrics 1343124 $wrkdir/RDlt1kb/rd.lt1kb.nodepth $wrkdir/RDlt1kb/ROC.indep.txt $wrkdir/RDlt1kb/ROC.dep.txt    

##depth only##

##DEL##

## training set from highest quality BAF variants that are made by a depth caller, variants must be in highly mappable region and greater than 5kb ##
cat $wrkdir/BAF/BAF.del.pred|fgrep depth|awk '{if ($2<0.1) print $1 "\t" "Fail";else if ($2>=0.9) print $1 "\t" "Pass"}'|sort -k1,1|join -j 1 - <(sort -k1,1 $wrkdir/Phase1.all.metrics|awk '{if ($(NF-1)>=5000 && $NF<0.3) print $1,$26,$27,$28 }')|cat <(awk '{print $1,"Status",$26,$27,$28}' $wrkdir/Phase1.all.metrics|head -n 1) -|tr ' ' '\t'>$wrkdir/rd.depth.del/rd.gt5000.depth.del.clean.metrics

## training set from BAF variants that are made by a depth caller, variants must be in poorly mapping region and greater than 5kb ##
cat $wrkdir/BAF/BAF.del.pred|fgrep depth|awk '{if ($2<=0.1) print $1 "\t" "Fail";else if ($2>=0.9) print $1 "\t" "Pass"}'|sort -k1,1|join -j 1 - <(sort -k1,1 $wrkdir/Phase1.all.metrics|awk '{if ($(NF-1)>=5000 && $NF>=0.3) print $1,$26,$27,$28 }')|cat <(awk '{print $1,"Status",$26,$27,$28}' $wrkdir/Phase1.all.metrics|head -n 1) -|tr ' ' '\t'>$wrkdir/rd.depth.del/rd.gt5000.depth.del.poor.metrics

##Cleanly mapping depth variants gt 5kb##

egrep -w "DEL" $wrkdir/Phase1.all.metrics|fgrep depth|awk '{if ($(NF-1)>5000 && $NF<0.3) print }'|cat <(head -n 1 $wrkdir/Phase1.all.metrics) - >$wrkdir/rd.depth.del/rd.all.gt5000.depth.del.clean.metrics

##Poorly mapping depth variants gt 5kb##

egrep -w "DEL" $wrkdir/Phase1.all.metrics|fgrep depth|awk '{if ($(NF-1)>5000 && $NF>=0.3) print }'|cat <(head -n 1 $wrkdir/Phase1.all.metrics) - >$wrkdir/rd.depth.del/rd.all.gt5000.depth.del.poor.metrics


##Del first filter on RD_Median_Separation and RD_log_pval then RD_log_2ndMaxP##
echo -e "RD_Median_Separation""\n""RD_log_pval"  >$wrkdir/rd.depth.del/ROC.indep.clean.txt 
echo "RD_log_2ndMaxP">$wrkdir/rd.depth.del/ROC.dep.clean.txt 


##Clean Random Forest##

Rscript /PHShome/hb875/talkowski/hb875/code/RF_git/RandomForest.R $wrkdir/rd.depth.del/rd.gt5000.depth.del.clean.metrics $wrkdir/rd.depth.del/rd.all.gt5000.depth.del.clean.metrics 1343124 $wrkdir/rd.depth.del/rd.gt5000.depth.del.clean $wrkdir/rd.depth.del/ROC.indep.clean.txt $wrkdir/rd.depth.del/ROC.dep.clean.txt    

##Poor determined by clean metrics##
Rscript /PHShome/hb875/talkowski/hb875/code/RF_git/RandomForest.R $wrkdir/rd.depth.del/rd.gt5000.depth.del.clean.metrics $wrkdir/rd.depth.del/rd.all.gt5000.depth.del.poor.metrics 1343124 $wrkdir/rd.depth.del/rd.gt5000.depth.del.poor $wrkdir/rd.depth.del/ROC.indep.clean.txt $wrkdir/rd.depth.del/ROC.dep.clean.txt $wrkdir/rd.depth.del/rd.gt5000.depth.del.clean.cutoffs   

##DUP##
mkdir $wrkdir/rd.depth.dup

## training set from highest quality BAF variants that are made by a depth caller, variants must be in highly mappable region and greater than 5kb ##
cat $wrkdir/BAF/BAF.dup.pred|fgrep depth|awk '{if ($2<0.1) print $1 "\t" "Fail";else if ($2>=0.9) print $1 "\t" "Pass"}'|sort -k1,1|join -j 1 - <(sort -k1,1 $wrkdir/Phase1.all.metrics|awk '{if ($(NF-1)>=5000 && $NF<0.3) print $1,$26,$27,$28 }')|cat <(awk '{print $1,"Status",$26,$27,$28}' $wrkdir/Phase1.all.metrics|head -n 1) -|tr ' ' '\t'>$wrkdir/rd.depth.dup/rd.gt5000.depth.dup.clean.metrics

## training set from BAF variants that are made by a depth caller, variants must be in poorly mapping region and greater than 5kb ##
cat $wrkdir/BAF/BAF.dup.pred|fgrep depth|awk '{if ($2<=0.1) print $1 "\t" "Fail";else if ($2>=0.9) print $1 "\t" "Pass"}'|sort -k1,1|join -j 1 - <(sort -k1,1 $wrkdir/Phase1.all.metrics|awk '{if ($(NF-1)>=5000 && $NF>=0.3) print $1,$26,$27,$28 }')|cat <(awk '{print $1,"Status",$26,$27,$28}' $wrkdir/Phase1.all.metrics|head -n 1) -|tr ' ' '\t'>$wrkdir/rd.depth.dup/rd.gt5000.depth.dup.poor.metrics

##Cleanly mapping depth variants gt 5kb##

egrep -w "DUP" $wrkdir/Phase1.all.metrics|fgrep depth|awk '{if ($(NF-1)>5000 && $NF<0.3) print }'|cat <(head -n 1 $wrkdir/Phase1.all.metrics) - >$wrkdir/rd.depth.dup/rd.all.gt5000.depth.dup.clean.metrics

##Poorly mapping depth variants gt 5kb##

egrep -w "DUP" $wrkdir/Phase1.all.metrics|fgrep depth|awk '{if ($(NF-1)>5000 && $NF>=0.3) print }'|cat <(head -n 1 $wrkdir/Phase1.all.metrics) - >$wrkdir/rd.depth.dup/rd.all.gt5000.depth.dup.poor.metrics

##Dup first filter on RD_Median_Separation and RD_log_pval then RD_log_2ndMaxP##
echo -e "RD_Median_Separation""\n""RD_log_pval"  >$wrkdir/rd.depth.dup/ROC.indep.clean.txt 
echo "RD_log_2ndMaxP">$wrkdir/rd.depth.dup/ROC.dep.clean.txt 

##Clean Random Forest##

Rscript /PHShome/hb875/talkowski/hb875/code/RF_git/RandomForest.R $wrkdir/rd.depth.dup/rd.gt5000.depth.dup.clean.metrics $wrkdir/rd.depth.dup/rd.all.gt5000.depth.dup.clean.metrics 1343124 $wrkdir/rd.depth.dup/rd.gt5000.depth.dup.clean $wrkdir/rd.depth.dup/ROC.indep.clean.txt $wrkdir/rd.depth.dup/ROC.dep.clean.txt

##Poor Random Forest##

Rscript /PHShome/hb875/talkowski/hb875/code/RF_git/RandomForest.R $wrkdir/rd.depth.dup/rd.gt5000.depth.dup.clean.metrics $wrkdir/rd.depth.dup/rd.all.gt5000.depth.dup.poor.metrics 1343124 $wrkdir/rd.depth.dup/rd.gt5000.depth.dup.poor $wrkdir/rd.depth.dup/ROC.indep.clean.txt $wrkdir/rd.depth.dup/ROC.dep.clean.txt $wrkdir/rd.depth.dup/rd.gt5000.depth.dup.clean.cutoffs  

##PE##
##train PE on BAF results, remove any discordant##
mkdir $wrkdir/PE

##Determine training variants that are concordant with initial depth cutoff and new BAF results## 
##create a training set with PE metrics##
cat $wrkdir/SR/SR.pred $wrkdir/RDgt1kb/rd.gt1kb.nodepth.pred|fgrep -wvf $wrkdir/one_end_linesine_txt|egrep -v "name|depth"|awk '{if ($2>=0.9) print $1 "\t" "Pass";else if ($2<=0.1) print $1 "\t" "Fail"}'|sort|uniq -c|awk '{if ($1==2) print $2 "\t" $3}'|sort -k1,1|join -j 1 - <(sort -k1,1 $wrkdir/Phase1.all.metrics|awk '{if ($NF<0.3) print $1,$6,$7,$8}')|cat <(awk '{print $1,"Status",$6,$7,$8}' $wrkdir/Phase1.all.metrics|head -n 1) -|tr ' ' '\t'>$wrkdir/PE/PE.metrics

##print all variants to be assessed by PE metrics (exclude depth only and wham INV and BND)##

ln -s $wrkdir/SR/SR.all.metrics $wrkdir/PE/PE.all.metrics

ln -s $wrkdir/SR/SR.all.cnv.metrics $wrkdir/PE/PE.all.cnv.metrics


##Dup first filter on PE_log_pval then PE_bg_median##
echo "PE_log_pval">$wrkdir/PE/ROC.indep.txt 
echo "">$wrkdir/PE/ROC.dep.txt 

##Generate PE RF and make predictions##
Rscript /PHShome/hb875/talkowski/hb875/code/RF_git/RandomForest.R $wrkdir/PE/PE.metrics $wrkdir/PE/PE.all.metrics 1343124 $wrkdir/PE/PE $wrkdir/PE/ROC.indep.txt $wrkdir/PE/ROC.dep.txt 

##second iteration necessary for BAF and SR###
##BAF##
##Require pass or fail PE,SR,RD exclude any BAF that fail first round from being passing in training set##

cat $wrkdir/RDgt1kb/rd.gt1kb.nodepth.pred $wrkdir/SR/SR.pred $wrkdir/PE/PE.pred <(awk '{if ($2<=0.1) print }' $wrkdir/BAF/BAF.del.pred)|awk '{if ($2>=0.9) print $1 "\t" "Pass";else if ($2<=0.1) print $1 "\t" "Fail"}'|egrep -v "name|depth"|sort|uniq -c |awk '{if ($1==4 && $3== "Fail") print $2 "\t" "Fail";else if ($1==3 && $3== "Pass") print $2 "\t" "Pass"}'|fgrep -wvf $wrkdir/BAF/Phase1.BAF.noassess|sort -k1,1|join -j 1 - <(fgrep DEL $wrkdir/Phase1.all.metrics|awk '{if ($NF<0.3 && $(NF-1)>5000) print $1,$30,$31}')|cat <(awk '{print "CNVID","Status",$30,$31}' $wrkdir/Phase1.all.metrics|head -n 1) - |tr ' ' '\t'>$wrkdir/BAF/BAF.del.v2.metrics

cat $wrkdir/RDgt1kb/rd.gt1kb.nodepth.pred $wrkdir/SR/SR.pred $wrkdir/PE/PE.pred <(awk '{if ($2<=0.1) print }' $wrkdir/BAF/BAF.dup.pred)|awk '{if ($2>=0.9) print $1 "\t" "Pass";else if ($2<=0.1) print $1 "\t" "Fail"}'|egrep -v "name|depth"|sort|uniq -c |awk '{if ($1==4 && $3== "Fail") print $2 "\t" "Fail";else if ($1==3 && $3== "Pass") print $2 "\t" "Pass"}'|fgrep -wvf $wrkdir/BAF/Phase1.BAF.noassess|sort -k1,1|join -j 1 - <(fgrep dup $wrkdir/Phase1.all.metrics|awk '{if ($NF<0.3 && $(NF-1)>5000) print $1,$33,$34}')|cat <(awk '{print "CNVID","Status",$33,$34}' $wrkdir/Phase1.all.metrics|head -n 1) - |tr ' ' '\t'>$wrkdir/BAF/BAF.dup.v2.metrics


##Run model##
##deletions##
Rscript /PHShome/hb875/talkowski/hb875/code/RF_git/RandomForest.R $wrkdir/BAF/BAF.del.v2.metrics $wrkdir/BAF/BAF.del.all.metrics 1343124 $wrkdir/BAF/BAF.del.v2 $wrkdir/BAF/ROC.del.indep.txt $wrkdir/BAF/ROC.del.dep.txt

##duplications##
Rscript /PHShome/hb875/talkowski/hb875/code/RF_git/RandomForest.R $wrkdir/BAF/BAF.dup.v2.metrics $wrkdir/BAF/BAF.dup.all.metrics 1343124 $wrkdir/BAF/BAF.dup.v2 $wrkdir/BAF/ROC.dup.indep.txt $wrkdir/BAF/ROC.dup.dep.txt 


##SR##

##Require pass or fail PE,SR,RD exclude any SR that fail first round from being passing in training set##

cat $wrkdir/RDgt1kb/rd.gt1kb.nodepth.pred  $wrkdir/PE/PE.pred <(awk '{if ($2<=0.1) print }' $wrkdir/SR/SR.pred)|awk '{if ($2>=0.9) print $1 "\t" "Pass";else if ($2<=0.1) print $1 "\t" "Fail"}'|egrep -v "name|depth"|sort|uniq -c |awk '{if ($1==3 && $3== "Fail") print $2 "\t" "Fail";else if ($1==2 && $3== "Pass") print $2 "\t" "Pass"}'|sort -k1,1|join -j 1 - <(cat $wrkdir/Phase1.all.metrics|awk '{if ($NF<0.3 && $(NF-1)>5000) print $1,$11,$14,$17}')|cat <(awk '{print $1,"Status",$11,$14,$17}' $wrkdir/Phase1.all.metrics|head -n 1) - |tr ' ' '\t'>$wrkdir/SR/SR.v2.metrics


##Run model##

Rscript /PHShome/hb875/talkowski/hb875/code/RF_git/RandomForest.R $wrkdir/SR/SR.v2.metrics $wrkdir/SR/SR.all.metrics 1343124 $wrkdir/SR/SR.v2 $wrkdir/SR/ROC.indep.txt $wrkdir/SR/ROC.dep.txt

##PE/SR##
##create PE_SR metrics##
cat $wrkdir/Phase1.all.metrics|awk '{if ($1!="name")print $1,$7+$14,$8+$17}'|cat <(awk '{print $1,$7"_"$14,$8"_"$17}' $wrkdir/Phase1.all.metrics|head -n 1) - |tr ' ' '\t'>$wrkdir/PE_SR/PE_SR.missingpois.metrics

Rscript $wrkdir/poiss.R /PHShome/hb875/SFARI/DeepSeq/HarrisonCompare/GT_Filtering/clean_pipeline/PE_SR/PE_SR.missingpois.metrics $wrkdir/PE_SR/PE_SR

##train with PE, SR, and RD##
cat $wrkdir/RDgt1kb/rd.gt1kb.nodepth.pred $wrkdir/SR/SR.v2.pred $wrkdir/PE/PE.pred|awk '{if ($2>=0.9) print $1 "\t" "Pass";else if ($2<=0.1) print $1 "\t" "Fail"}'|egrep -v "name|depth"|sort|uniq -c |awk '{if ($1==3 && $3== "Fail") print $2 "\t" "Fail";else if ($1==3 && $3== "Pass") print $2 "\t" "Pass"}'|sort -k1,1|join -j 1 - <(cat $wrkdir/PE_SR/PE_SR.poissP_metrics|fgrep -wf <(awk '{if ($NF<0.3) print $1}' $wrkdir/Phase1.all.metrics))|cat <(awk '{print $1,"Status",$2,$3,$4}' $wrkdir/PE_SR/PE_SR.poissP_metrics|head -n 1) - |tr ' ' '\t'|sed 's/Inf/300/g'>$wrkdir/PE_SR/PE_SR.metrics


awk '{if ($1!~"depth" && !($1~"wham" && ($2~"INV" || $2~"BND"))) print $1}' $wrkdir/Phase1.all.metrics| fgrep -wf - $wrkdir/PE_SR/PE_SR.poissP_metrics|sed 's/Inf/300/g'>$wrkdir/PE_SR/PE_SR.all.metrics

##run model##

##Dup first filter on PE_log_pval then PE_bg_median##
echo "poiss_p">$wrkdir/PE_SR/ROC.indep.txt 
echo "">$wrkdir/PE_SR/ROC.dep.txt 


Rscript /PHShome/hb875/talkowski/hb875/code/RF_git/RandomForest.R $wrkdir/PE_SR/PE_SR.metrics $wrkdir/PE_SR/PE_SR.all.metrics 1343124 $wrkdir/PE_SR/PE_SR $wrkdir/PE_SR/ROC.indep.txt $wrkdir/PE_SR/ROC.dep.txt 


##Put all prediction classes together###

cat $wrkdir/BAF/BAF.del.v2.pred $wrkdir/BAF/BAF.dup.v2.pred|awk '{print $1 "\t" $NF}'|sort -k 1,1 | \
	##join BAF and PE## \
	join -e NA -a 2 -a 1 -j 1 - <(sort -k1,1 $wrkdir/PE/PE.pred|awk '{print $1 "\t" $NF}') -o 1.1,2.1,1.2,2.2| \
	awk '{if ($1=="NA") print $2,$3,$4;else print $1,$3,$4 }'|sort -k1,1|\
	##Join BAF PE with SR## \
	join -e NA -a 2 -a 1 -j 1 -  <(sort -k1,1 $wrkdir/SR/SR.v2.pred|awk '{print $1 "\t" $NF}') -o 1.1,2.1,1.2,1.3,2.2 | \
	awk '{if ($1=="NA") print $2,$3,$4,$5;else print $1,$3,$4,$5 }'|sort -k1,1| \
	##Add PE/SR## \
	join -e NA -a 2 -a 1 -j 1 -  <(cat $wrkdir/PE_SR/PE_SR.pred| \
	awk '{print $1 "\t" $NF}'|sort -k1,1 ) -o 1.1,2.1,1.2,1.3,1.4,2.2| \
	awk '{if ($1=="NA") print $2,$3,$4,$5,$6;else print $1,$3,$4,$5,$6 }'|fgrep -v name | \
	##Add RD## \
	join -e NA -a 2 -a 1 -j 1 -  <(cat $wrkdir/RDlt1kb/rd.lt1kb.nodepth.pred \
	$wrkdir/RDgt1kb/rd.gt1kb.nodepth.pred $wrkdir/rd*/rd.gt5000.depth.*.pred| \
	awk '{print $1 "\t" $NF}'|sort -k1,1 ) -o 1.1,2.1,1.2,1.3,1.4,1.5,2.2| \
	awk '{if ($1=="NA") print $2,$3,$4,$5,$6,$7;else print $1,$3,$4,$5,$6,$7 }'|fgrep -v name | \
	##Add header## \
	awk '{if (NR==1) print "name","BAF","PE","SR","PE/SR","RD" "\n" $0;else print $0}'| \
	tr ' ' '\t'>$wrkdir/metric.table

##combined p-value##

mkdir combined_prob

##PE,SR, probs###
##Pull out max PE/SR score add bonus for the other if > 0.5 and then add an RD bonus if probRD >0.5##
##Then check if PE_SR combined has higher probability, which will be used###

egrep -v "depth|name" $wrkdir/metric.table |awk '{if ($3>=$4 && $4<.5) print $1 "\t" $3 "\t" $5 ; else if ($4>$3 && $3<.5) print $1 "\t" $4 "\t" $5; else if ($3>=$4 && $4>=.5) print $1 "\t" ($3+($4-.5)*(1-$3)) "\t" $5 ;else if ($4>$3 && $3>=.5) print $1 "\t" ($4+($3-.5)*(1-$4)) "\t" $5; }'|awk '{if ($2>=$3) print $1 "\t" $2;else print $1 "\t" $3}' |cat <(echo -e "name" "\t" "PE_SR.prob") ->$wrkdir/combined_prob/PE_SR.all.prob


##RD, BAF##
## if probability fails (p<0.5) use RD value, else if BAF greater than 50% minus 0.5 * 1-probRD to the probRD, which will either increase probability ##

awk '{if ($6!="NA" && $1!="name") print}' $wrkdir/metric.table |awk '{if ($NF<.5 || $2=="NA" || $2<0.5) print $1 "\t" $NF ;else print $1 "\t" $NF+(($2-.5)*(1-$NF))}'|cat <(echo -e "name" "\t" "Depth.prob") ->$wrkdir/combined_prob/RD.all.prob


##PE/SR gt 1kb CNV##
##average PE/SR+RD if RD support or no PE/SR support##
 awk '{print $1}' $wrkdir/RDgt1kb/rd.gt1kb.nodepth.pred|fgrep -wf - $wrkdir/combined_prob/RD.all.prob|sort -k1,1|join -j 1 - <(sort -k1,1 $wrkdir/combined_prob/PE_SR.all.prob) |awk '{if ($2>=0.5 || $3<0.5)print $1 "\t" ($2+$3)/2}'>$wrkdir/combined_prob/CNV.PESR.rdgt1kb_combined.p

##PE/SR support only no Depth gt1kb##
##Pull out anything without read depth support (prob<0.5) and check if PE/SR pass, add bonus ##
 awk '{print $1}' $wrkdir/RDgt1kb/rd.gt1kb.nodepth.pred|fgrep -wf - $wrkdir/combined_prob/RD.all.prob|sort -k1,1|join -j 1 - <(sort -k1,1 $wrkdir/combined_prob/PE_SR.all.prob) |awk '{if ($2<0.5 && 3>=0.5)print $1 "\t" $3}'>$wrkdir/combined_prob/CNV.PESR.nodepth.gt1kb_combined.p
 
#PE/SR lt 1kb CNV###
##Use PE/SR score add bonus for an RD bonus if probRD >0.5##
 awk '{print $1}' $wrkdir/RDlt1kb/rd.lt1kb.nodepth.pred|fgrep -wf - $wrkdir/combined_prob/RD.all.prob|sort -k1,1|join -j 1 - <(sort -k1,1 $wrkdir/combined_prob/PE_SR.all.prob) |awk '{if ($2<0.5 || $3<0.5)print $1 "\t" $3;else print $1 "\t" $3+(($2-.5)*(1-$3))}' >$wrkdir/combined_prob/CNV.PESR.lt1kb_combined.p
 
##depthonly support## 
fgrep depth $wrkdir/combined_prob/RD.all.prob>$wrkdir/combined_prob/CNV.depthvars_combined.p

##BCA##
##Take max PE/SR and then add bonus##
egrep -wv "DEL|DUP|name" $wrkdir/Phase1.all.metrics|awk '{print $1}'|fgrep -wf - PE_SR.all.prob>$wrkdir/combined_prob/BCA.PESR.p

##Find passing variants##

cat $wrkdir/combined_prob/CNV.PESR.rdgt1kb_combined.p $wrkdir/combined_prob/CNV.depthvars_combined.p $wrkdir/combined_prob/CNV.PESR.lt1kb_combined.p |awk '{if ($NF>=.5) print}'|cat <(echo -e "name" '\t' "Probability") - >$wrkdir/combined_prob/combined_CNV.passing.p


cat $wrkdir/combined_prob/BCA.PESR.p $wrkdir/combined_prob/CNV.PESR.nodepth.gt1kb_combined.p|awk '{if ($NF>=.5) print}'|cat <(echo -e "name" '\t' "Probability") - >$wrkdir/combined_prob/combined_BCA.passing.p


##Mosaic CNV check##
##Mosaic CNV check###
mkdir mosaic_CNV
##make a metric file excluding seperation##
##restrict to single family for these or single family and one de novo##

##all CNV gt 1kb except depth calls which must be gt 5kb and called in one family uniparently###
awk '{if ($1!~"depth" && $(NF-1)>1000 && ($2~"DEL" || $2~"DUP"))  print; else if ($1~"depth" && $(NF-1)>=5000) print }' $wrkdir/Phase1.all.metrics|awk '{if ( ($3<0.0011 && $4<0.00314 && $3>0) || ($3==0 && $4<0.0011)) print }'  |cat <(head -n 1 $wrkdir/Phase1.all.metrics) - >$wrkdir/mosaic_CNV/mosaic.all.nodepth.metrics

awk '{if ($1~"depth" && $(NF-1)>5000 && ($2~"DUP"))  print; else if ($1~"depth" && $(NF-1)>=5000) print }' $wrkdir/Phase1.all.metrics|awk '{if ( ($3<0.0011 && $4<0.00314 && $3>0) || ($3==0 && $4<0.0011)) print }'  |cat <(head -n 1 $wrkdir/Phase1.all.metrics) - >$wrkdir/mosaic_CNV/mosaic.all.depthdup.metrics

awk '{if ($1~"depth" && $(NF-1)>5000 && ($2~"DEL"))  print; else if ($1~"depth" && $(NF-1)>=5000) print }' $wrkdir/Phase1.all.metrics|awk '{if ( ($3<0.0011 && $4<0.00314 && $3>0) || ($3==0 && $4<0.0011)) print }'  |cat <(head -n 1 $wrkdir/Phase1.all.metrics) - >$wrkdir/mosaic_CNV/mosaic.all.depthdel.metrics


##Found sep outliers at at 0.437306, conservatively drop to 0.40,this value should be hard coded for future studies###

##Find potential mosaic for BAF testing##
##sep cut-off 0.40##
###RDtest.gt1kb##

sep_col=$(fgrep -wn RD_Median_Separation  <(head -n 1 $wrkdir/mosaic_CNV/mosaic.all.nodepth.metrics|tr '\t' '\n')|awk -F":" '{print $1}')
p_val_col=$(fgrep -wn RD_log_pval  <(head -n 1 $wrkdir/mosaic_CNV/mosaic.all.nodepth.metrics|tr '\t' '\n')|awk -F":" '{print $1}')
secMaxP_col=$(fgrep -wn RD_log_2ndMaxP  <(head -n 1 $wrkdir/mosaic_CNV/mosaic.all.nodepth.metrics|tr '\t' '\n')|awk -F":" '{print $1}')

sep=0.40
p_val=$(fgrep -w RD_log_pval $wrkdir/RDgt1kb/rd.gt1kb.nodepth.cutoffs|awk '{print $2}'|sed 's/NA/0/g')
secMaxP=$(fgrep -w RD_log_2ndMaxP $wrkdir/RDgt1kb/rd.gt1kb.nodepth.cutoffs|awk '{print $2}'|sed 's/NA/0/g')

##gt1kb not depth##
awk -v sep_col=$sep_col -v p_val_col=$p_val_col -v secMaxP_col=$secMaxP_col -v sep=$sep -v p_val=$p_val -v secMaxP=$secMaxP  '{if ($sep_col<sep && $p_val_col>=p_val_col && $secMaxP_col>=secMaxP ) print}' $wrkdir/mosaic_CNV/mosaic.all.nodepth.metrics>$wrkdir/mosaic_CNV/mosaic.BAFcheck.nodepth.txt

##depth only del##
sep=0.40
p_val=$(fgrep -w RD_log_pval $wrkdir/rd.depth.del/rd.gt5000.depth.del.clean.cutoffs|awk '{print $2}'|sed 's/NA/0/g')
secMaxP=$(fgrep -w RD_log_2ndMaxP $wrkdir/rd.depth.del/rd.gt5000.depth.del.clean.cutoffs|awk '{print $2}'|sed 's/NA/0/g')

awk -v sep_col=$sep_col -v p_val_col=$p_val_col -v secMaxP_col=$secMaxP_col -v sep=$sep -v p_val=$p_val -v secMaxP=$secMaxP  '{if ($sep_col<sep && $p_val_col>=p_val_col && $secMaxP_col>=secMaxP ) print}' $wrkdir/mosaic_CNV/mosaic.all.depthdel.metrics>$wrkdir/mosaic_CNV/mosaic.BAFcheck.depthdel.txt

##depth only dup##
sep=0.40
p_val=$(fgrep -w RD_log_pval $wrkdir/rd.depth.dup/rd.gt5000.depth.dup.clean.cutoffs|awk '{print $2}'|sed 's/NA/0/g')
secMaxP=$(fgrep -w RD_log_2ndMaxP $wrkdir/rd.depth.dup/rd.gt5000.depth.dup.clean.cutoffs|awk '{print $2}'|sed 's/NA/0/g')

awk -v sep_col=$sep_col -v p_val_col=$p_val_col -v secMaxP_col=$secMaxP_col -v sep=$sep -v p_val=$p_val -v secMaxP=$secMaxP  '{if ($sep_col<sep && $p_val_col>=p_val_col && $secMaxP_col>=secMaxP ) print}' $wrkdir/mosaic_CNV/mosaic.all.depthdup.metrics>$wrkdir/mosaic_CNV/mosaic.BAFcheck.depthdup.txt


