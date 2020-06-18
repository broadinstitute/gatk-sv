#!/usr/bin/env Rscript
# This script is used to modify the SR bg info for common (freq > 50) SVs with SR support on both breakpoints
# Sets PESR evidence metrics to 0 at common sites where there is SR support on both sides.
library("optparse")

option_list = list(
         make_option(c("-p","--pedfile"), type="character", default=NULL,
              help="name of ped file", metavar="character" ),
         make_option(c("-b","--bedfile"), type="character", default=NULL,
              help="name of input SV in bed format", metavar="character" ),
         make_option(c("-e","--evidence"), type="character", default=NULL,
              help="name of evidence metrics", metavar="character" ),
         make_option(c("-o","--output"), type="character", default=NULL,
              help="name of output metrics", metavar="character" )
 );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

ped_file = opt$pedfile
bed_file = opt$bedfile
evidence_metrics = opt$evidence
out_file = opt$output


ped = read.table(ped_file)
bed = read.table(bed_file)
bed[,ncol(bed)+1] = apply(bed,1,function(x){length(strsplit(as.character(x[6]),',')[[1]])})
evi = read.table(evidence_metrics,sep='\t', header=T)

svid=bed[bed[,ncol(bed)]>nrow(ped)/2-1,4]
evi_modi = evi[evi$name%in%svid & !is.na(evi$SR_posA_called_median) & evi$SR_posA_called_median>0 & !is.na(evi$SR_posB_called_median) & evi$SR_posB_called_median>0 ,]

evi_modi$SR_posA_bg_median = 0
evi_modi$SR_posB_bg_median = 0
evi_modi$SR_sum_bg_median  = 0
evi_modi$SR_posA_bg_frac = 0
evi_modi$SR_posB_bg_frac = 0
evi_modi$SR_sum_bg_frac  = 0
evi_modi$PE_bg_median = 0
evi_modi$PE_bg_frac = 0

evi_oth = evi[!evi$name %in% evi_modi$name,]
evi_new = rbind(evi_modi,evi_oth)

write.table(evi_new,out_file, quote=F, sep='\t', col.names=T, row.names=F)

