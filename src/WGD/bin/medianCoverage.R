#!/usr/bin/env Rscript

# Script to calculate median bin coverage per sample for all samples in a bincov
# matrix

# Note: loads entire coverage matrix into memory. This may pose a problem for
# large matrices or on small-memory machines. Other workarounds exist from the
# command line, but most are slower. One suggested alternative is to split
# the input coverage matrices by chromosome prior to computing medians.

# Load library
require(optparse)

# Define options
option_list <- list(
  make_option(c("-b","--binwise"), action="store_true", default=FALSE,
              help="compute medians of all samples per bin [default: median of all bins per sample]"),
  make_option(c("-m","--mad"), action="store_true", default=FALSE,
              help="compute median absolute deviation of all bins per sample [default: FALSE]"),
  make_option(c("-H","--header"), action="store_true", default=FALSE,
              help="input coverage matrix has header with sample IDs [default: FALSE]"))

# Get command-line arguments and options
args <- parse_args(OptionParser(usage="%prog [options] covMatrix.bed OUTFILE",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 2)
{cat("Incorrect number of required positional arguments\n\n")
  stop()}

# Read matrix
if(opts$header==T){
  cov <- read.table(args$args[1], header=T, comment.char="", check.names=F)
}else{
  cov <- read.table(args$args[1], header=F, comment.char="#", check.names=F)
}

# Function to compute medians per sample
covPerSample <- function(cov,downsample=1000000,mad=F){
  # Downsample to 1M random rows if nrows > 1M (for computational efficiency)
  if(nrow(cov)>1000000){
    cov <- cov[sample(1:nrow(cov), downsample),]
  }
  # Get medians with and without zero-cov bins
  zerobins <- which(as.integer(apply(as.data.frame(cov[,-c(1:3)]), 1, median, na.rm=T)) == 0)
  withzeros <- as.numeric(apply(as.data.frame(cov[,-c(1:3)]), 2, median, na.rm=T))
  withoutzeros <- as.numeric(apply(as.data.frame(cov[-zerobins,-c(1:3)]), 2, median, na.rm=T))
  #Get SDs with and without zero-cov bins (if optioned)
  if(mad==T){
    withzeros.mad <- as.numeric(apply(as.data.frame(cov[,-c(1:3)]), 2, mad, na.rm=T))
    withoutzeros.mad <- as.numeric(apply(as.data.frame(cov[-zerobins,-c(1:3)]), 2, mad, na.rm=T))
  }
  # Compile results df to return
  if(mad==T){
    res <- data.frame("ID"=paste("Sample",1:(ncol(cov)-3),sep=""),
                      "Med_withZeros"=withzeros,
                      "Med_withoutZeros"=withoutzeros,
                      "MAD_withZeros"=withzeros.mad,
                      "MAD_withoutZeros"=withoutzeros.mad)
  }else{
    res <- data.frame("ID"=paste("Sample",1:(ncol(cov)-3),sep=""),
                      "Med_withZeros"=withzeros,
                      "Med_withoutZeros"=withoutzeros)
  }
  # Replace sample IDs if input matrix has header
  if(opts$header==T){
    res$ID <- names(cov[, -c(1:3), drop = FALSE])
  }
  # Return output df
  return(res)
}

# Function to compute medians per bin
covPerBin <- function(cov,downsample=500,mad=F){
  # Downsample to 500 random samples if nsamples > 500 (for computational efficiency)
  if(ncol(cov)>503){
    cov <- cov[,sample(1:ncol(cov), downsample)]
  }
  # Get medians with and without zero-cov samples
  meds <- t(apply(as.data.frame(cov[,-c(1:3)]), 1, function(vals){
    withzeros <- median(vals, na.rm=T)
    if(any(vals>0)){
      withoutzeros <- median(vals[which(vals>0)], na.rm=T)
    }else{
      withoutzeros <- NA
    }
    return(c(withzeros,withoutzeros))
  }))
  # Get standard deviations (if optioned)
  sds <- t(apply(as.data.frame(cov[,-c(1:3)]), 1, function(vals){
    withzeros <- mad(vals, na.rm=T)
    if(any(vals>0)){
      withoutzeros <- mad(vals[which(vals>0)], na.rm=T)
    }else{
      withoutzeros <- NA
    }
    return(c(withzeros,withoutzeros))
  }))
  # compile results df to return
  if(mad==T){
    res <- data.frame("#chr"=cov[,1],"start"=cov[,2],"end"=cov[,3],
                      "Med_withZeros"=meds[,1],
                      "Med_withoutZeros"=meds[,2],
                      "MAD_withZeros"=sds[,1],
                      "MAD_withoutZeros"=sds[,2])
  }else{
    res <- data.frame("#chr"=cov[,1],"start"=cov[,2],"end"=cov[,3],
                      "Med_withZeros"=meds[,1],
                      "Med_withoutZeros"=meds[,2])
  }
  # Return output df
  return(res)
}

# Compute appropriate medians & write out
if(opts$binwise==TRUE){
  res <- covPerBin(cov,mad=opts$mad)
  names(res)[1] <- "#chr"
  write.table(res,args$args[2], sep="\t", col.names=T, row.names=F, quote=F)
}else{
  res <- covPerSample(cov,mad=opts$mad)
  names(res)[1] <- "#sample_id"
  write.table(res,args$args[2], sep="\t", col.names=T, row.names=F, quote=F)
}
