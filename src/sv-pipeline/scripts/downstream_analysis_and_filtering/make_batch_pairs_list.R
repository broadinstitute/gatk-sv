#!/usr/bin/env Rscript

# Talkowski SV pipeline downstream analysis helper script

# Make list of all nonredundant pairs of batches from an input list of batches


###Set global parameters
options(stringsAsFactors=F,scipen=1000)


###Read command-line arguments
args <- commandArgs(trailingOnly=T)
in.list <- as.character(args[1])
OUTFILE <- as.character(args[2])

# #Dev parameters
# in.list <- "~/scratch/af_test/input.list"
# OUTFILE <- "~/scratch/batch_pairs_test.txt"


###Process input data
batches <- as.character(read.table(in.list,header=F,sep="\t")[,1])


###Generate list of all possible pairs
out <- do.call("rbind", lapply(1:length(batches),function(a){
  do.call("rbind", lapply(2:length(batches),function(b){
    if(a<b){
      return(data.frame(batches[a],batches[b]))
    }
  }))
}))


###Write list to outfile
write.table(out,OUTFILE,col.names=F,row.names=F,sep="\t",quote=F)
