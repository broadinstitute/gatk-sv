#!/usr/bin/env Rscript

# Contact: Ryan Collins <rlcollins@g.harvard.edu>

# Helper script to merge benchmarking results from a single sample 
# across multiple versions (filters, technologies, etc)

# Set parameters & read inputs
options(stringsAsFactors=F, scipen=1000)
args <- commandArgs(trailingOnly=TRUE)
inputs <- read.table(args[1], header=F, sep="\t")
vinfo <- read.table(args[2], header=F, sep="\t")
colnames(vinfo) <- c("chr", "start", "end", "svtype", "AF", "EV")

# Helper function to pre-process each input BED file
load.bed <- function(path, dname="overlap"){
  df <- read.table(path, sep="\t", comment.char="", header=T)
  colnames(df)[1] <- "chr"
  df$VID <- NULL
  df$AF <- NULL
  df[, dname] <- 1
  orig.ovr.cidxs <- grep("^ovr", colnames(df))
  df[, "ovr"] <- apply(df[, grep("^ovr[1-2]", colnames(df))], 1, 
                       function(vals){any(vals != "NO_OVR")})
  df[, orig.ovr.cidxs] <- NULL
  return(df[which(!duplicated(df)), ])
}

# Load all input BED files as list
dat <- apply(inputs, 1, function(vals){load.bed(vals[2], vals[1])})
names(dat) <- inputs[, 1]

# Collapse results
key.cols <- c("chr", "start", "end", "svtype", "length", "ovr")
res <- dat[[1]]
if(length(dat) > 1){
  for(i in 2:length(dat)){
    res <- merge(res, dat[[i]], all=T, by=key.cols, sort=F)
  }
}
hit.idxs <- which(colnames(res) %in% names(dat))
res[, hit.idxs] <- apply(res[, hit.idxs], 2, function(vals){
  vals[which(is.na(vals))] <- FALSE
  return(as.numeric(vals))
  })
res <- merge(vinfo, res, by=c("chr", "start", "end", "svtype"), sort=F, all=F)
res <- res[with(res, order(chr, start, end)), ]

# Write to output file
colnames(res)[1] <- "#chr"
write.table(res, args[3], col.names=T, row.names=F, sep="\t", quote=F)
