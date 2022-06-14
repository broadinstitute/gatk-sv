#!/usr/bin/env Rscript

# Copyright (c) 2022 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# Talkowski SV pipeline downstream analysis helper script

# Test for differences in AFs between PCR+ and PCR- batches


###Set global parameters
options(stringsAsFactors=F, scipen=1000)


###################
###HELPER FUNCTIONS
###################
# Sum AC and AN data per population for PCR+ and PCR- samples
collapse.data.by.pcr.status <- function(dat, minus.batches, plus.batches, allpops){
  minus.df <- do.call(cbind, lapply(allpops, function(pop){
    do.call(cbind, lapply(c("AC", "AN"), function(x){
      apply(dat[, paste(pop, "_", x, ".", minus.batches, sep="")], 1, sum)
    }))
  }))
  colnames(minus.df) <- unlist(sapply(allpops, function(pop){
    paste(pop, "_", c("AC", "AN"), ".minus", sep="")}))
  plus.df <- do.call(cbind, lapply(allpops, function(pop){
    do.call(cbind, lapply(c("AC", "AN"), function(x){
      apply(dat[, paste(pop, "_", x, ".", plus.batches, sep="")], 1, sum)
    }))
  }))
  colnames(plus.df) <- unlist(sapply(allpops, function(pop){
    paste(pop, "_", c("AC", "AN"), ".plus", sep="")}))
  cbind(dat[, 1:3], minus.df, plus.df)
}

# Find best population for AF comparisons for a single variant and run chi-squared test
compare.pcr.afs <- function(vals, allpops){
  # Clean data
  VID <- as.character(vals[1])
  vals <- vals[-c(1:3)]
  vnames <- names(vals)
  vals <- as.numeric(vals)
  names(vals) <- vnames
  
  # Find population where AC > 0 that has largest minimum AN between PCR+ and PCR- 
  minus.ACs <- vals[paste(allpops, "AC.minus", sep="_")]
  plus.ACs <- vals[paste(allpops, "AC.plus", sep="_")]
  elig.pops <- allpops[which(apply(cbind(minus.ACs, plus.ACs), 1, sum) > 0)]
  minus.ANs <- vals[paste(elig.pops, "AN.minus", sep="_")]
  plus.ANs <- vals[paste(elig.pops, "AN.plus", sep="_")]
  min.ANs <- apply(cbind(minus.ANs, plus.ANs), 1, min)
  best.pop <- elig.pops[which(min.ANs == max(min.ANs, na.rm=T))]
  
  # Extra data and run chi-squared test
  minus.AN <- as.numeric(vals[paste(best.pop, "AN.minus", sep="_")])
  minus.AC <- as.numeric(vals[paste(best.pop, "AC.minus", sep="_")])
  minus.AF <- minus.AC / minus.AN
  plus.AN <- as.numeric(vals[paste(best.pop, "AN.plus", sep="_")])
  plus.AC <- as.numeric(vals[paste(best.pop, "AC.plus", sep="_")])
  plus.AF <- plus.AC / plus.AN
  chisq.p <- chisq.test(matrix(c(minus.AN-minus.AC, minus.AC,
                                 plus.AN-plus.AC, plus.AC),
                               nrow=2, byrow=F))$p.value
  # Return values
  data.frame("VID"=VID, "pop"=best.pop, "minus.AF"=minus.AF, 
             "plus.AF"=plus.AF, "chisq.p"=chisq.p)
}


### Read command-line arguments
args <- commandArgs(trailingOnly=T)
infile <- as.character(args[1])
minus.batches.in <- as.character(args[2])
plus.batches.in <- as.character(args[3])
OUTFILE <- as.character(args[4])

# #Dev parameters:
# infile <- "~/scratch/test_AF_tables.txt.gz"
# minus.batches.in <- "~/scratch/gnomAD-SV-v3.1.PCRMINUS.batches.list"
# plus.batches.in <- "~/scratch/gnomAD-SV-v3.1.PCRPLUS.batches.list"
# OUTFILE <- "~/scratch/test_batchFx.tsv"

### Process input AF data
dat <- read.table(infile, header=T, sep="\t", comment.char="", check.names=F)
allpops <- unique(sapply(colnames(dat)[4:ncol(dat)],
                         function(x){strsplit(as.character(x),'_')[[1]][1]}))
batches.in.header <- unique(sapply(colnames(dat)[grep("_AN.", colnames(dat), fixed=T)], 
                                   function(x){unlist(strsplit(x, split="_AN.", fixed=T))[2]}))

### Load PCR+ and PCR- batches and match with AF data header
minus.batches <- intersect(unique(read.table(minus.batches.in, header=F)[, 1]),
                           batches.in.header)
plus.batches <- intersect(unique(read.table(plus.batches.in, header=F)[, 1]),
                          batches.in.header)

# Subset data to batches in either PCR+ or PCR- lists
batch.suffixes <- sapply(colnames(dat)[-c(1:3)], 
                         function(x){unlist(strsplit(x, split="_A[NC]\\." ))[2]})
dat <- dat[, c(1:3, which(batch.suffixes %in% c(minus.batches, plus.batches))+3)]

# Collapse AC/AN data by PCR status
dat <- collapse.data.by.pcr.status(dat, minus.batches, plus.batches, allpops)

### Compare AFs between PCR+ and PCR- batches per variant
res <- do.call(rbind, apply(dat, 1, compare.pcr.afs, allpops=allpops))
colnames(res) <- c("VID", "pop", "minus.AF", "plus.AF", "chisq.p")

### Write to outfile
write.table(res, OUTFILE, col.names=T, row.names=F, sep="\t", quote=F)
