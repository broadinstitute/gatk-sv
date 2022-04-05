#!/usr/bin/env Rscript

# Contact: Ryan Collins <rlcollins@g.harvard.edu>

# Helper script to summarize benchmarking results from a single sample 
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
