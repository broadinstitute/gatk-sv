#!/usr/bin/env Rscript

# Script to QC merged per sample bincov files 
# Created by: Chelsea Lowther <clowther@mgh.harvard.edu>

# Inputs: 
# 1) A single bincov file per sample that contains bins from the last 1 Mb per chromosome 
# 2) A sample ID string 

# Usage: bincov_QC.R BINCOV SAMPLE_ID OUTFILE

## Libraries
require(optparse)

## Define Options 
option_list <- list(
  make_option(c("-h", "--help"), action="store_true", default=FALSE,
              help="Script to QC per sample bincov files"))

## Get positional arguments
args <- parse_args(OptionParser(usage="%prog [options] BINCOV SAMPLE_ID OUTFILE",
                                option_list=option_list),
                   positional_arguments = TRUE)

opts <- args$options 

## Calculate number of bins with no values per chromosome 
    bincov <- read.table(args$args[1], header = F, fill=NA)
    bincov_na <- aggregate(V4 ~ V1, data=bincov, function(x) {sum(is.na(x))}, na.action = NULL)
    bincov_na <- as.data.frame(sum(bincov_na$V4))
    
    write.table(bincov_na, args$args[3], row.names=F, sep = "\t", quote = F, col.names =F)                       
