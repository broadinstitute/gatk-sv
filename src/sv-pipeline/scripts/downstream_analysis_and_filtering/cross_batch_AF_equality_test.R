#!/usr/bin/env Rscript

# Copyright (c) 2022 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# Talkowski SV pipeline downstream analysis helper script

# Test variants for equality of allele frequencies across all batches


###Set global parameters
options(stringsAsFactors=F, scipen=1000)


###################
###HELPER FUNCTIONS
###################
# Test all variants for significant difference in AFs across any batch
# Returns one P-value per variant
test.variants <- function(dat, allpops){
  pvals <- apply(dat[, -c(1:3)], 1, function(vals){
    # Determine optimal population for cross-batch test
    # Optimal pop: AC>0 in at least one batch, then maximize 5th %ile of AN across batches
    AC.gt.0.cols <- names(which(vals[grep("_AC.", names(vals), fixed=T)] > 0))
    AC.gt.0 <- unique(sapply(AC.gt.0.cols, function(x){gsub("_AC", "", unlist(strsplit(x, split=".", fixed=T))[1])}))
    min.AN <- sapply(allpops, function(pop){
      as.numeric(quantile(vals[grep(paste(pop, "AN", sep="_"), names(vals))], prob=0.05))
    })
    if(length(AC.gt.0) > 0){
      min.AN <- min.AN[which(names(min.AN) %in% AC.gt.0)]
    }
    optimal.pop <- names(min.AN)[which(min.AN == max(min.AN))]
    # Compute P-value from proportion test across all batches for optimal.pop
    ACs <- vals[grep(paste(optimal.pop, "AC.", sep="_"), names(vals), fixed=T)]
    ANs <- vals[grep(paste(optimal.pop, "AN.", sep="_"), names(vals), fixed=T)]
    keep.batch.idx <- which(ANs > 0)
    if(length(keep.batch.idx) > 1){
      prop.test(ACs[keep.batch.idx], ANs[keep.batch.idx], correct=F)$p.value
    }else{
      return(NA)
    }
  })
  data.frame("VID"=dat$VID, "prop.test.pvalue"=pvals)
}


########################
###RSCRIPT FUNCTIONALITY
########################
### Read command-line arguments
args <- commandArgs(trailingOnly=T)
infile <- as.character(args[1])
batchlist.in <- as.character(args[2])
OUTFILE <- as.character(args[3])

# #Dev parameters:
# infile <- "~/scratch/gnomAD-SV-v3.1.10pct_NCR.no_outliers.batch_fx.chr22.merged_AF_table.shard_7.txt.gz"
# batchlist.in <- "~/scratch/gnomAD-SV-v3.1.PCRMINUS.batches.list"
# OUTFILE <- "~/scratch/test_batchFx.tsv"

### Process input data and subset to batches present in batchlist.in
dat <- read.table(infile, header=T, sep="\t", comment.char="", check.names=F)
allpops <- unique(sapply(colnames(dat)[4:ncol(dat)],
                         function(x){strsplit(as.character(x),'_')[[1]][1]}))
batches.in.header <- unique(sapply(colnames(dat)[grep("_AN.", colnames(dat), fixed=T)], 
                                   function(x){unlist(strsplit(x, split="_AN.", fixed=T))[2]}))
query.batches <- unique(as.character(read.table(batchlist.in, header=F, sep="\t")[, 1]))
allbatches <- intersect(query.batches, batches.in.header)
# Subset data to batches in query.batches
batch.suffixes <- sapply(colnames(dat)[-c(1:3)], 
                         function(x){unlist(strsplit(x, split="_A[NC]\\." ))[2]})
dat <- dat[, c(1:3, which(batch.suffixes %in% allbatches)+3)]

### Run proportion test for all variants
pvals <- test.variants(dat, allpops)

### Write to outfile
write.table(pvals, OUTFILE, col.names=T, row.names=F, sep="\t", quote=F)
