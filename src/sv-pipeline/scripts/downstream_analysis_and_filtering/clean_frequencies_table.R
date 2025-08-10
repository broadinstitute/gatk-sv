#!/usr/bin/env Rscript

# Talkowski SV pipeline downstream analysis helper script

# Clean allele frequency table across batches


###Set global parameters
options(stringsAsFactors=F,scipen=1000)


###################
###HELPER FUNCTIONS
###################
#Import a single table of freq data
import.freq.table <- function(path, pops){
  #Read data
  dat <- read.table(path,header=T,comment.char="",sep="\t")
  colnames(dat)[1:3] <- c("VID","SVLEN","SVTYPE")
  #Convert multiallelic CN_NONREF_COUNT and CN_NUMBER to AC and AN, and drop CN_ columns
  mcnv.idxs <- which(dat$SVTYPE %in% c("CNV", "MCNV"))
  ac_cols <- grep("AC_",colnames(dat),fixed=T)
  cn_nonref_cols <- grep("CN_NONREF_COUNT",colnames(dat),fixed=T)
  an_cols <- grep("AN_",colnames(dat),fixed=T)
  cn_number_cols <- grep("CN_NUMBER",colnames(dat),fixed=T)
  
  if(length(mcnv.idxs) > 0 && length(ac_cols) > 0 && length(cn_nonref_cols) > 0) {
    dat[mcnv.idxs, ac_cols] <- dat[mcnv.idxs, cn_nonref_cols]
  }
  if(length(mcnv.idxs) > 0 && length(an_cols) > 0 && length(cn_number_cols) > 0) {
    dat[mcnv.idxs, an_cols] <- dat[mcnv.idxs, cn_number_cols]
  }
  
  cols_to_remove <- union(cn_nonref_cols, cn_number_cols)
  if(length(cols_to_remove) > 0) {
    dat <- dat[, -cols_to_remove]
  }
  #Convert all ANs and ACs to numerics
  dat[, -c(1:3)] <- apply(dat[, -c(1:3)], 2, as.numeric)
  #Adjust calls on sex chromosomes
  for(pop in pops){
    x.idx <- unique(c(grep("_X_",dat[,1],fixed=T),
                      grep("_chrX_",dat[,1],fixed=T)))
    y.idx <- unique(c(grep("_Y_",dat[,1],fixed=T),
                      grep("_chrY_",dat[,1],fixed=T)))
    if(length(c(x.idx, y.idx)) > 0){
      male_ac_colname <- paste(pop,"MALE","AC",sep="_")
      male_an_colname <- paste(pop,"MALE","AN",sep="_")
      female_ac_colname <- paste(pop,"FEMALE","AC",sep="_")
      female_an_colname <- paste(pop,"FEMALE","AN",sep="_")
      if(4 == length(intersect(colnames(dat),
                               c(male_ac_colname, male_an_colname,
                                 female_ac_colname, female_an_colname)))){
        master.AC.idx <- which(colnames(dat)==paste(pop,"AC",sep="_"))
        master.AN.idx <- which(colnames(dat)==paste(pop,"AN",sep="_"))
        male.AC.idx <- which(colnames(dat)==male_ac_colname)
        male.AN.idx <- which(colnames(dat)==male_an_colname)
        female.AC.idx <- which(colnames(dat)==female_ac_colname)
        female.AN.idx <- which(colnames(dat)==female_an_colname)
        #For variants on chrX, overwrite master counts with female-specific counts
        dat[x.idx,master.AC.idx] <- dat[x.idx,female.AC.idx]
        dat[x.idx,master.AN.idx] <- dat[x.idx,female.AN.idx]
        #For all variants on chrY, overwrite master counts with male-specific counts
        dat[y.idx,master.AC.idx] <- dat[y.idx,male.AC.idx]
        dat[y.idx,master.AN.idx] <- dat[y.idx,male.AN.idx]
      }
    }
  }
  #Drop all male- and female-specific columns from table
  dat <- dat[, grep("MALE", colnames(dat), fixed=T, invert=T)]
  #Add batch name to colnames
  return(dat)
}


###Read command-line arguments
args <- commandArgs(trailingOnly=T)
INFILE <- as.character(args[1])
POPFILE <- as.character(args[2])
OUTFILE <- as.character(args[3])


###Process input data
allpops <- sort(unique(as.character(read.table(POPFILE, sep="\t", header=F)[, 2])))
dat <- import.freq.table(path=INFILE, pops=allpops)


###Write output data
write.table(dat,OUTFILE,col.names=T,row.names=F,sep="\t",quote=F)
system(paste("gzip -f ",OUTFILE,sep=""),wait=T,intern=F)
