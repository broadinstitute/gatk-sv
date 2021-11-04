#!/usr/bin/env Rscript

# Helper script to count median # of filtered variants per trio
# for minGQ optimization filtering workflow


###Set master parameters
options(stringsAsFactors=F,scipen=1000)



################
###RSCRIPT BLOCK
################
require(optparse)
###List of command-line options
option_list <- list(
  make_option(c("--ID"), type="character", default="condition",
              help="condition ID [default %default]")
)

###Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog INFILE FAMFILE OUTFILE",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

###Checks for appropriate positional arguments
if(length(args$args) != 3){
  stop("Incorrect number of required positional arguments\n")
}

###Writes args & opts to vars
INFILE <- args$args[1]
FAMFILE <- args$args[2]
OUTFILE <- args$args[3]
ID <- opts$ID

###Reads data
dat <- as.character(read.table(INFILE,header=F)[,1])
fams <- unique(as.character(read.table(FAMFILE,header=F)[,1]))

###Computes # of variants per family
counts <- as.integer(sapply(fams,function(fam){
  return(length(which(dat==fam)))
}))

###Reports results
out <- data.frame("condition"=ID,
                  "hetsPerProband_median"=median(counts,na.rm=T),
                  "hetsPerProband_Q1"=quantile(counts,probs=0.25,na.rm=T),
                  "hetsPerProband_Q3"=quantile(counts,probs=0.75,na.rm=T))
colnames(out)[1] <- "#condition"
write.table(out,OUTFILE,col.names=T,row.names=F,sep="\t",quote=F)
