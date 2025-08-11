#!/usr/bin/env Rscript

# Talkowski SV pipeline downstream analysis helper script

# Make reclassification table for batch effect variants


###Set global parameters
options(stringsAsFactors=F,scipen=1000)


###################
###HELPER FUNCTIONS
###################
#Process master freq table data and compute joint PCR+ & PCR- AFs
import.freqs <- function(freq.table.in){
  dat <- read.table(freq.table.in,header=T,comment.char="",sep="\t")
  plus.ACs <- apply(dat[,intersect(grep("PCRPLUS",colnames(dat),fixed=T),
                        grep("_AC",colnames(dat),fixed=T))],1,sum,na.rm=T)
  plus.ANs <- apply(dat[,intersect(grep("PCRPLUS",colnames(dat),fixed=T),
                                   grep("_AN",colnames(dat),fixed=T))],1,sum,na.rm=T)
  plus.AFs <- plus.ACs/plus.ANs
  plus.AFs[which(is.na(plus.AFs))] <- 0

  minus.ACs <- apply(dat[,grep("_AC",colnames(dat),fixed=T)],1,sum,na.rm=T)
  minus.ANs <- apply(dat[,grep("_AN",colnames(dat),fixed=T)],1,sum,na.rm=T)
  minus.AFs <- minus.ACs/minus.ANs
  minus.AFs[which(is.na(minus.AFs))] <- 0
  out <- data.frame("VID"=dat$VID,"PCRPLUS_AF"=plus.AFs,"PCRMINUS_AF"=minus.AFs)  
  return(out)
}

#Process failure lists for a single comparison type
import.fails <- function(minus.in,prefix){
  # Check file size and return an empty dataframe if no records are marked as batch effects
  file_size <- file.info(minus.in)$size
  if(is.na(file_size) || file_size == 0){
    a <- as.data.frame(matrix(ncol=2))
  }else{
    a <- read.table(minus.in,header=F,sep="\t")
  }
  colnames(a) <- c("VID",paste("fails",prefix,sep="_"))
  a <- a[which(a$VID!="NA" & !is.na(a$VID)), ]
  #m <- read.table(minus.in,header=F,sep="\t")
  m <- a
  colnames(m) <- c("VID",paste("fails",prefix,"minus",sep="_"))
  merged <- merge(a,m,by="VID",sort=F,all=T)
  minus.idx <- which(colnames(merged)==paste("fails",prefix,"minus",sep="_"))
  merged[,minus.idx][which(is.na(merged[,minus.idx]))] <- 0
  merged$frac_plus_fails <- 1-(merged[,which(colnames(merged)==paste("fails",prefix,"minus",sep="_"))]/
                                 merged[,which(colnames(merged)==paste("fails",prefix,sep="_"))])
  merged[,which(colnames(merged)==paste("fails",prefix,"minus",sep="_"))] <- NULL
  colnames(merged)[ncol(merged)] <- paste("frac_plus_fails",prefix,sep="_")
  return(merged)
}

#Categorize failing sites
categorize.failures <- function(dat,pairwise.cutoff,onevsall.cutoff){
  dat$max_plus_frac <- apply(data.frame(dat$frac_plus_fails_pairwise,
                                        dat$frac_plus_fails_onevsall),
                             1,max,na.rm=T)
  pairwise.fail.idx <- which(dat$fails_pairwise>=pairwise.cutoff & dat$fails_onevsall<onevsall.cutoff)
  onevsall.fail.idx <- which(dat$fails_pairwise<pairwise.cutoff & dat$fails_onevsall>=onevsall.cutoff)
  both.fail.idx <- which(dat$fails_pairwise>=pairwise.cutoff & dat$fails_onevsall>=onevsall.cutoff)
  plus.gt.minus <- which(dat$PCRPLUS_AF>=dat$PCRMINUS_AF)
  minus.gt.plus <- which(dat$PCRPLUS_AF<dat$PCRMINUS_AF)
  
  #Build vectors of VIDs for each category
  plus.enriched <- c()
  plus.depleted <- c()
  batch.variable <- c()
  for(i in pairwise.fail.idx){
    row <- dat[i,]
    if(as.numeric(row$frac_plus_fails_pairwise)>=0.11){
      if(as.numeric(row$PCRPLUS_AF)>=as.numeric(row$PCRMINUS_AF)){
        plus.enriched <- c(plus.enriched,i)
      }else{
        plus.depleted <- c(plus.depleted,i)
      }
    }else{
      batch.variable <- c(batch.variable,i)
    }
  }
  for(i in onevsall.fail.idx){
    row <- dat[i,]
    if(as.numeric(row$frac_plus_fails_onevsall)>=0.11){
      if(as.numeric(row$PCRPLUS_AF)>=as.numeric(row$PCRMINUS_AF)){
        plus.enriched <- c(plus.enriched,i)
      }else{
        plus.depleted <- c(plus.depleted,i)
      }
    }else{
      batch.variable <- c(batch.variable,i)
    }
  }
  for(i in both.fail.idx){
    row <- dat[i,]
    if(as.numeric(row$max_plus_frac)>=0.11){
      if(as.numeric(row$PCRPLUS_AF)>=as.numeric(row$PCRMINUS_AF)){
        plus.enriched <- c(plus.enriched,i)
      }else{
        plus.depleted <- c(plus.depleted,i)
      }
    }else{
      batch.variable <- c(batch.variable,i)
    }
  }
  
  #Build classification table
  out.table <- data.frame("VID"=dat$VID[c(plus.enriched,plus.depleted,batch.variable)],
                          "classification"=c(rep("PCRPLUS_ENRICHED",length(plus.enriched)),
                                             rep("PCRPLUS_DEPLETED",length(plus.depleted)),
                                             rep("VARIABLE_ACROSS_BATCHES",length(batch.variable))))
  return(out.table)
}



###Read command-line arguments
args <- commandArgs(trailingOnly=T)
freq.table.in <- as.character(args[1])
onevsall.in <- as.character(args[2])
OUTFILE <- as.character(args[3])
onevsall.cutoff <- as.integer(args[4])
# For PCRMinus_only version, we don't use pairwise comparisons
# Set defaults for missing parameters
pairwise.cutoff <- 999999  # Set very high so no pairwise failures are detected

# #Dev parameters:
# freq.table.in <- "~/scratch/gnomAD_v2_SV_MASTER.merged_AF_table.txt.gz"
# pairwise.minus.in <- "~/scratch/minGQ_test/mod07_mcnv_test.pairwise_comparisons.all.failures.txt"
# onevsall.minus.in <- "~/scratch/minGQ_test/mod07_mcnv_test.one_vs_all_comparisons.all.failures.txt"
# OUTFILE <- "~/scratch/batch_effects.reclassification_table.txt"
# pairwise.cutoff <- 12
# onevsall.cutoff <- 2


###Read data
freq.dat <- import.freqs(freq.table.in)

# For PCRMinus_only version, create empty pairwise.fails since we don't have pairwise comparisons
pairwise.fails <- data.frame(VID=character(0), fails_pairwise=integer(0), frac_plus_fails_pairwise=numeric(0))
onevsall.fails <- import.fails(onevsall.in, prefix="onevsall")
onevsall.fails <- onevsall.fails[which(onevsall.fails$fails_onevsall>=onevsall.cutoff), ]


###Combine data
merged <- merge(pairwise.fails,onevsall.fails,all=T,sort=F,by="VID")
if(nrow(merged) > 0){
  merged[,-1] <- apply(merged[,-1],2,function(vals){
    vals[which(is.na(vals))] <- 0
    return(vals)
  })
}
merged <- merge(merged,freq.dat,by="VID",sort=F,all=F)


##Categorize batch effect failure sites
out.table <- categorize.failures(dat=merged,
                                 pairwise.cutoff=pairwise.cutoff,
                                 onevsall.cutoff=onevsall.cutoff)
write.table(out.table,OUTFILE,col.names=F,row.names=F,sep="\t",quote=F)

  