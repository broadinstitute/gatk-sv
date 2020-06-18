#!/bin/bash

#/data/talkowski/Samples/common-mind/03_variant_filtering/filtered_vcfs##
##define metric file##

set -e

metric_file=$1

##define cutoff columns##
PE_log_pval_col=$(cat $metric_file|head -n 1|tr '\t' '\n'|fgrep -n PE_log_pval|awk -F':' '{print $1}')

SR_log_pval_col=$(cat $metric_file|head -n 1|tr '\t' '\n'|fgrep -n SR_sum_log_pval|awk -F':' '{print $1}')

PESR_log_pval_col=$(cat $metric_file|head -n 1|tr '\t' '\n'|fgrep -n PESR_log_pval|awk -F':' '{print $1}')

RD_Median_Separation_col=$(cat $metric_file|head -n 1|tr '\t' '\n'|fgrep -n RD_Median_Separation|awk -F':' '{print $1}')

RD_log_2ndMaxP_col=$(cat $metric_file|head -n 1|tr '\t' '\n'|fgrep -n RD_log_2ndMaxP|awk -F':' '{print $1}')

RD_log_pval_col=$(cat $metric_file|head -n 1|tr '\t' '\n'|fgrep -n RD_log_pval|awk -F':' '{print $1}')

svlen_col=$(cat $metric_file|head -n 1|tr '\t' '\n'|fgrep -n svlen|awk -F':' '{print $1}')

##define cutoff variables##
PE_log_pval=3.04 
SR_log_pval=2.17 
PESR_log_pval=5.21 
RD_Median_Separation_depth_del=0.40 
RD_log_pval_depth_del=9.41 
RD_log_2ndMaxP_depth_del=0 
RD_Median_Separation_depth_dup=0.40 
RD_log_pval_depth_dup=12.70 
RD_log_2ndMaxP_depth_dup=0.11 
RD_Median_Separation_gt1kb=0.30 
RD_log_pval_gt1kb=3.47 
RD_log_2ndMaxP_gt1kb=0 
  
       
       
##Pull out valid variants for each group##
##CNV gt1kb depth##
awk -v RD_Median_Separation_gt1kb=$RD_Median_Separation_gt1kb\
    -v RD_log_pval_gt1kb=$RD_log_pval_gt1kb\
    -v RD_log_2ndMaxP_gt1kb=$RD_log_2ndMaxP_gt1kb\
    -v RD_Median_Separation_col=$RD_Median_Separation_col\
    -v RD_log_2ndMaxP_col=$RD_log_2ndMaxP_col\
    -v RD_log_pval_col=$RD_log_pval_col\
    -v svlen_col=$svlen_col\
      '{if (RD_Median_Separation_gt1kb<$RD_Median_Separation_col && $RD_Median_Separation_col!="NA" && RD_log_2ndMaxP<$RD_log_2ndMaxP_col && RD_log_pval<$RD_log_pval_col && $svlen_col>=1000 && ($2=="DEL" || $2=="DUP") && $1!~"depth") print $1}' $metric_file>CNV.gt1kb.pass.txt

##PE##
awk -v PE_log_pval=$PE_log_pval\
    -v PE_log_pval_col=$PE_log_pval_col\
    '{if (PE_log_pval<$PE_log_pval_col && $1!~"depth" && $SR_log_pval_col!="NA") print $1 }' $metric_file>PE.pass.txt

##SR##
awk -v SR_log_pval=$SR_log_pval\
    -v SR_log_pval_col=$SR_log_pval_col\
    '{if (SR_log_pval<$SR_log_pval_col && $1!~"depth" && $SR_log_pval_col!="NA"  ) print $1 }' $metric_file>SR.pass.txt


##PE_SR##
awk -v PESR_log_pval=$PESR_log_pval\
    -v PESR_log_pval_col=$PESR_log_pval_col\
    '{if (PESR_log_pval<$PESR_log_pval_col && $1!~"depth" && $PESR_log_pval_col!="NA" ) print $1 }' $metric_file>PESR.pass.txt

##Depth Del##
awk -v RD_Median_Separation_depth_del=$RD_Median_Separation_depth_del\
    -v RD_log_pval_depth_del=$RD_log_pval_depth_del\
    -v RD_log_2ndMaxP_depth_del=$RD_log_2ndMaxP_depth_del\
    -v RD_Median_Separation_col=$RD_Median_Separation_col\
    -v RD_log_2ndMaxP_col=$RD_log_2ndMaxP_col\
    -v RD_log_pval_col=$RD_log_pval_col\
    -v svlen_col=$svlen_col\
      '{if (RD_Median_Separation_depth_del<$RD_Median_Separation_col && $RD_Median_Separation_col!="NA" && RD_log_2ndMaxP<$RD_log_2ndMaxP_col && RD_log_pval<$RD_log_pval_col && $svlen_col>5000 && ($2=="DEL" || $2=="DUP")) print $1}' $metric_file>CNV.depthdel.pass.txt


##Depth Dup##
awk -v RD_Median_Separation_depth_depth_dup=$RD_Median_Separation_depth__depth_dup\
    -v RD_log_pval_depth__depth_dup=$RD_log_pval_depth__depth_dup\
    -v RD_log_2ndMaxP_depth__depth_dup=$RD_log_2ndMaxP_depth__depth_dup\
    -v RD_Median_Separation_col=$RD_Median_Separation_col\
    -v RD_log_2ndMaxP_col=$RD_log_2ndMaxP_col\
    -v RD_log_pval_col=$RD_log_pval_col\
    -v svlen_col=$svlen_col\
      '{if (RD_Median_Separation_depth__depth_dup<$RD_Median_Separation_col && $RD_Median_Separation_col!="NA" && RD_log_2ndMaxP<$RD_log_2ndMaxP_col && RD_log_pval<$RD_log_pval_col && $svlen_col>=5000 && ($2=="_depth_dup" || $2=="DUP")) print $1}' $metric_file>CNV.depthdup.pass.txt

##Make calls##
##PE/SR caller##
##CNV >1kb##
##1) Pass CNV gt1kb RD metrics provided below && pass PE, SR, or PE_SR##
cat CNV.depthdel.pass.txt CNV.depthdup.pass.txt>int.cnv.pass.txt

cat PE.pass.txt SR.pass.txt PESR.pass.txt  |fgrep -wf - CNV.gt1kb.pass.txt>>int.cnv.pass.txt

cat PE.pass.txt SR.pass.txt PESR.pass.txt |fgrep -wf <(cat /data/talkowski/Samples/common-mind/03_variant_filtering/metrics/*.metrics.pos_add|fgrep -v name|awk -v svlen_col=$svlen_col '{if ($svlen_col>=1000) print $1}')|fgrep -wf - CNV.gt1kb.pass.txt>>int.cnv.pass.txt

cat PE.pass.txt SR.pass.txt PESR.pass.txt |fgrep -wf <(cat /data/talkowski/Samples/common-mind/03_variant_filtering/metrics/*.metrics.pos_add|fgrep -v name|awk -v svlen_col=$svlen_col '{if ($svlen_col<1000 && ($2=="DEL" || $2=="DUP")) print $1}')>>int.cnv.pass.txt

##BCA##
cat PE.pass.txt SR.pass.txt PESR.pass.txt |fgrep -wf <(cat /data/talkowski/Samples/common-mind/03_variant_filtering/metrics/*.metrics.pos_add|fgrep -v name|awk -v svlen_col=$svlen_col '{if (!($2=="DEL" || $2=="DUP")) print $1}')>int.BCA.pass.txt

cat PE.pass.txt SR.pass.txt PESR.pass.txt |fgrep -wf <(cat /data/talkowski/Samples/common-mind/03_variant_filtering/metrics/*.metrics.pos_add|fgrep -v name|awk -v svlen_col=$svlen_col '{if ($svlen_col>=1000 && ($2=="DEL" || $2=="DUP")) print $1}')|fgrep -wvf CNV.gt1kb.pass.txt>>int.BCA.pass.txt

sort -u int.BCA.pass.txt|awk '{print $1 "\t" "1"}'>BCA.pass.txt
sort -u int.cnv.pass.txt|awk '{print $1 "\t" "1"}'>CNV.pass.txt

