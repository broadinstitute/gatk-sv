#!/usr/bin/env Rscript

# Contact: Ryan Collins <rlcollins@g.harvard.edu>

# Helper script to merge benchmarking results from a single sample 
# across multiple versions (filters, technologies, etc)

# Set parameters & read inputs
options(stringsAsFactors=F, scipen=1000)
args <- commandArgs(trailingOnly=TRUE)
inputs <- read.table(args[1], header=F, sep="\t")

# Helper function to pre-process each input BED file
load.bed <- function(path, dname="overlap"){
  df <- read.table(path, sep="\t", comment.char="", header=T)
  colnames(df)[1] <- "chr"
  colnames(df)[which(colnames(df) == "AF")] <- paste("AF", dname, sep=".")
  df$VID <- NULL
  orig.ovr.cidxs <- grep("^ovr", colnames(df))
  df[, dname] <- apply(df[, grep("^ovr[1-2]", colnames(df))], 1, 
                       function(vals){any(vals != "NO_OVR")})
  df[, orig.ovr.cidxs] <- NULL
  return(df[which(!duplicated(df)), ])
}

# Load all input BED files as list
dat <- apply(inputs, 1, function(vals){load.bed(vals[2], vals[1])})
names(dat) <- inputs[, 1]

# Collapse results
key.cols <- c("chr", "start", "end", "svtype", "length")
key.col.idxs <- which(colnames(dat[[1]]) %in% key.cols)
key.cols <- colnames(dat[[1]])[key.col.idxs]
res <- dat[[1]]
if(length(dat) > 1){
  for(i in 2:length(dat)){
    res <- merge(res, dat[[i]], all=T, by=key.cols, sort=F)
  }
  # Take average of AFs
  AF.idxs <- grep("^AF\\.", colnames(res))
  mean.AFs <- apply(res[, AF.idxs], 1, mean, na.rm=T)
  res[, AF.idxs] <- NULL
  res$AF <- mean.AFs
}
hit.idxs <- which(colnames(res) %in% names(dat))
res[, hit.idxs] <- apply(res[, hit.idxs], 2, function(vals){
  vals[which(is.na(vals))] <- FALSE
  return(as.numeric(vals))
  })
res <- res[with(res, order(chr, start, end)), ]

# Write to output file
colnames(res)[1] <- "#chr"
write.table(res, args[2], col.names=T, row.names=F, sep="\t", quote=F)
