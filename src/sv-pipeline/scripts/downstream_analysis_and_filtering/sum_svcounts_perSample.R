#!/usr/bin/env Rscript

# Helper script to merge outputs from count_svtypes task
# in final_outlier_sample_filter.wdl


### Set parameters & read arguments
options(stringsAsFactors=F, scipen=1000)
args <- commandArgs(trailingOnly=TRUE)
INFILE <- args[1]
OUTFILE <- args[2]

### Read input data & reformat
dat <- read.table(INFILE, header=F, check.names=F, sep="\t", comment.char="")
drop.rows <- which(dat[, 1] %in% c("sample", "#sample"))
if(length(drop.rows) > 0){
  dat <- dat[-drop.rows, ]
}
colnames(dat) <- c("sample", "svtype", "count", "chrom")[1:ncol(dat)]
dat$count <- as.numeric(dat$count)
samples <- as.character(unique(dat$sample))
svtypes <- as.character(unique(dat$svtype))

### Get sum of counts per sample per svtype
summed.res <- do.call("rbind", lapply(samples, function(sample){
  do.call("rbind", lapply(svtypes, function(svtype){
    data.frame("sample"=sample,
               "svtype"=svtype,
               "count"=sum(dat[which(dat$sample==sample & dat$svtype==svtype), ]$count,
                           na.rm=T))
  }))
}))

### Write summed results to outfile
write.table(summed.res, OUTFILE, col.names=T, row.names=F, sep="\t", quote=F)
