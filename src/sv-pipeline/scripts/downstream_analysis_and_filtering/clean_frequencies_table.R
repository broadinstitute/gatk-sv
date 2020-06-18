#!/usr/bin/env Rscript

# Copyright (c) 2018 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# Talkowski SV pipeline downstream analysis helper script

# Clean allele frequency table across batches


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
import.freq.table <- function(path,pops){
  #Read data
  dat <- read.table(path,header=T,comment.char="",sep="\t")
  colnames(dat)[1:3] <- c("VID","SVLEN","SVTYPE")
  #Clean multiallelic ACs
  dat[,grep("_AC",colnames(dat),fixed=T)] <- apply(dat[,grep("_AC",colnames(dat),fixed=T)],2,clean.ACs)
  #Convert all ANs to numerics
  dat[,grep("_AN",colnames(dat),fixed=T)] <- apply(dat[,grep("_AN",colnames(dat),fixed=T)],2,as.numeric)
  #Double counts for MCNVs (since first allele for all MCNVs is nulled out)
  MCNV.idx <- grep("_MCNV_",dat[,1],fixed=T)
  dat[MCNV.idx,grep("_AC",colnames(dat),fixed=T)] <- 2*dat[MCNV.idx,grep("_AC",colnames(dat),fixed=T)]
  dat[MCNV.idx,grep("_AN",colnames(dat),fixed=T)] <- 2*dat[MCNV.idx,grep("_AN",colnames(dat),fixed=T)]
  #Adjust calls on sex chromosomes
  for(pop in pops){
    x.idx <- unique(c(grep("_X_",dat[,1],fixed=T),
                      grep("_chrX_",dat[,1],fixed=T)))
    y.idx <- unique(c(grep("_Y_",dat[,1],fixed=T),
                      grep("_chrY_",dat[,1],fixed=T)))
    master.AC.idx <- which(colnames(dat)==paste(pop,"AC",sep="_"))
    master.AN.idx <- which(colnames(dat)==paste(pop,"AN",sep="_"))
    male.AC.idx <- which(colnames(dat)==paste(pop,"MALE","AC",sep="_"))
    male.AN.idx <- which(colnames(dat)==paste(pop,"MALE","AN",sep="_"))
    female.AC.idx <- which(colnames(dat)==paste(pop,"FEMALE","AC",sep="_"))
    female.AN.idx <- which(colnames(dat)==paste(pop,"FEMALE","AN",sep="_"))
    #For variants on chrX, overwrite master counts with female-specific counts
    dat[x.idx,master.AC.idx] <- dat[x.idx,female.AC.idx]
    dat[x.idx,master.AN.idx] <- dat[x.idx,female.AN.idx]
    #For all variants on chrY, overwrite master counts with male-specific counts
    dat[y.idx,master.AC.idx] <- dat[y.idx,male.AC.idx]
    dat[y.idx,master.AN.idx] <- dat[y.idx,male.AN.idx]
  }
  #Drop all male- and female-specific columns from table
  dat <- dat[,-grep("MALE",colnames(dat),fixed=T)]
  #Add batch name to colnames
  return(dat)
}


###Read command-line arguments
args <- commandArgs(trailingOnly=T)
INFILE <- as.character(args[1])
OUTFILE <- as.character(args[2])

# #Dev parameters (local)
# INFILE <- "~/scratch/gnomAD_v2_SV_PCRPLUS_Q1_batch_1.frequencies.preclean.txt.gz"


###Process input data
dat <- import.freq.table(path=INFILE,pops=allpops)


###Write output data
write.table(dat,OUTFILE,col.names=T,row.names=F,sep="\t",quote=F)
system(paste("gzip -f ",OUTFILE,sep=""),wait=T,intern=F)
