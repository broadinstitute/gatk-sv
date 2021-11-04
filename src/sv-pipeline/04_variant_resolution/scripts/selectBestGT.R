#!/usr/bin/env Rscript

# Simple helper script to iterate over a file containing 
# pipeline-formatted genotype columns and returning the 
# best non-reference genotype


####Set parameters & read args
options(stringsAsFactors=F,scipen=1000)
args <- commandArgs(trailingOnly=T)
if(length(args)!=2){
  stop("Requires two positional arguments: input file and output file.")
}
INFILE <- args[1]
OUTFILE <- args[2]


####Read file
dat <- read.table(args[1],header=F,sep="\t")
dat <- apply(dat,2,as.character)


####Helper function
selectBestGT <- function(GTs){
  melted <- matrix(unlist(lapply(GTs,strsplit,split=":")),
                   ncol=9,byrow=T)
  nonref <- which(!(melted[,1] %in% c("0/0","./.")))
  ref <- which(melted[,1] == "0/0")
  nocall <- which(melted[,1] == "./.")
  #Select best non-ref GT, if any exist
  if(length(nonref)>0){
    best <- tail(order(as.numeric(melted[nonref,2])),n=1)
  #Otherwise, select best reference GT, if any exist
  }else if(length(ref)>0){
    best <- tail(order(as.numeric(melted[ref,2])),n=1)
  #Otherwise, pick first non-ref GT
  }else{
    best <- head(nocall,1)
  }
  return(GTs[best])
}


####Select best genotypes & write to output file
final.GTs <- matrix(apply(dat,2,selectBestGT),nrow=1,byrow=T)
write.table(final.GTs,OUTFILE,col.names=F,row.names=F,sep="\t",quote=F)

