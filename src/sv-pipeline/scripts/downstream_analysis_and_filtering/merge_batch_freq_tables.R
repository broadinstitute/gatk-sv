#!/usr/bin/env Rscript

# Talkowski SV pipeline downstream analysis helper script

# Merge AF tables across batches


###Set global parameters
options(stringsAsFactors=F,scipen=1000)
svtypes <- c("DEL","DUP","INS","INV","CPX","BND")
allpops <- c("AFR","ASN","EUR","HSP")


###################
###HELPER FUNCTIONS
###################
#Collapse multiallelic ACs
clean.ACs <- function(ACs){
  ACs.to.clean <- grep(",",ACs,fixed=T)
  cleaned.ACs <- as.vector(as.numeric(sapply(ACs[ACs.to.clean],function(s){
    return(sum(as.numeric(unlist(strsplit(s,split=",")))[-2]))
  })))
  ACs[ACs.to.clean] <- cleaned.ACs
  return(as.numeric(ACs))
}
#Import a single table of freq data
import.freq.table <- function(batch,path){
  #Read data
  dat <- read.table(path,header=T,comment.char="",sep="\t")
  #Clean multiallelic ACs
  dat[,grep("_AC",colnames(dat),fixed=T)] <- apply(dat[,grep("_AC",colnames(dat),fixed=T)],2,clean.ACs)
  #Convert all ANs to numerics
  dat[,grep("_AN",colnames(dat),fixed=T)] <- apply(dat[,grep("_AN",colnames(dat),fixed=T)],2,as.numeric)
  #Add batch name to colnames
  colnames(dat)[-c(1:3)] <- paste(colnames(dat)[-c(1:3)],batch,sep=".")
  colnames(dat)[1:3] <- c("VID","SVLEN","SVTYPE")
  return(dat)
}
#Import a list of freq data tables and merge them
import.freq.data <- function(batches.list){
  merged <- import.freq.table(batch=batches.list[1,1],path=batches.list[1,2])
  for(i in 2:nrow(batches.list)){
    newdat <- import.freq.table(batch=batches.list[i,1],path=batches.list[i,2])
    merged <- merge(x=merged,y=newdat,by=c("VID","SVLEN","SVTYPE"),sort=F,all=T)
    rm(newdat)
  }
  return(merged)
}


###Read command-line arguments
args <- commandArgs(trailingOnly=T)
in.list <- as.character(args[1])
OUTFILE <- as.character(args[2])

# #Dev parameters (local)
# in.list <- "~/scratch/af_test/input.list"
# OUTFILE <- "~/scratch/gnomAD_AF_table.merged.txt"


###Process input data
batches.list <- read.table(in.list,header=F,sep="\t")
colnames(batches.list) <- c("batch","path")
dat <- import.freq.data(batches.list)

###Write output data
write.table(dat,OUTFILE,col.names=T,row.names=F,sep="\t",quote=F)
system(paste("gzip -f ",OUTFILE,sep=""),wait=T,intern=F)
