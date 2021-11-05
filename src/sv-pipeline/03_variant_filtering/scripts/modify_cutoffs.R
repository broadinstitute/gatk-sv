#!R
#!/usr/bin/env Rscript
# This script is used to modify the output scores metrics for common (freq > 50%) SVs.
# Sets scores to 0.51 (ie passing) for sites that passed at in at least one of SR_sum_log_pval, PE_log_pval, or PESR_log_pval

library("optparse")

option_list = list(
         make_option(c("-c","--cutoff"), type="character", default=NULL,	help="name of RF cutoffs file", metavar="character" ),
         make_option(c("-m","--metrics"), type="character", default=NULL,	help="name of metrics file", metavar="character" ),
         make_option(c("-s","--score"), type="character", default=NULL,		help="name of RF scores file", metavar="character" ),
         make_option(c("-o","--output"), type="character", default=NULL,	help="name of output metrics", metavar="character" )
 );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

cutoff_file = opt$cutoff
metrics_file = opt$metrics
scores_file = opt$score
output_file = opt$output

cutoff = read.table(cutoff_file, sep='\t', header=T)
metrics = read.table(metrics_file, sep='\t', header=T)
scores = read.table(scores_file, sep='\t', header=T)

SR_sum_log_pval = cutoff[cutoff$metric=='SR_sum_log_pval',2]
PE_log_pval = cutoff[cutoff$metric=='PE_log_pval',2]
PESR_log_pval = cutoff[cutoff$metric=='PESR_log_pval',2]

SR_sum_log_pval_kept_svid = metrics[!is.na(metrics$SR_sum_log_pval) & metrics$SR_sum_log_pval>SR_sum_log_pval,]
PE_log_pval_kept_svid = metrics[!is.na(metrics$PE_log_pval) & metrics$PE_log_pval>PE_log_pval,]
PESR_log_pval_kept_svid = metrics[!is.na(metrics$PESR_log_pval) & metrics$PESR_log_pval>PESR_log_pval,]
scores[scores[,1]%in%SR_sum_log_pval_kept_svid[,1],]$score = 0.51
scores[scores[,1]%in%PE_log_pval_kept_svid[,1],]$score = 0.51
scores[scores[,1]%in%PESR_log_pval_kept_svid[,1],]$score = 0.51
write.table(scores, output_file, quote = F, sep='\t', col.names=T, row.names=F)

