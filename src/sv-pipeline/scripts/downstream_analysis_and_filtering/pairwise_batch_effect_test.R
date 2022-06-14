#!/usr/bin/env Rscript

# Copyright (c) 2022 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# Talkowski SV pipeline downstream analysis helper script

# Find variants with differences in allele frequencies between pairs of batches


###Set global parameters
options(stringsAsFactors=F, scipen=1000)


###################
###HELPER FUNCTIONS
###################
#For any two batches, find most comparable AFs for each variant and run chi-squared test
compare.batches <- function(dat, batch1, batch2, allpops, min.AN=30, return.all.stats=FALSE){
  #Subset data for each batch (for convenience)
  #1: restrict to sites with >0 AC in at least one batch
  b1.dat <- dat[, c(1:3, grep(paste("\\.", batch1, "$", sep=""), colnames(dat)))]
  b1.maxAC <- apply(b1.dat[, grep("_AC",colnames(b1.dat),fixed=T)],1,max)
  if(batch2 != "ALL_OTHERS"){
    b2.dat <- dat[,c(1:3, grep(paste("\\.", batch2, "$", sep=""), colnames(dat)))]
    b2.maxAC <- apply(b2.dat[,grep("_AC",colnames(b2.dat),fixed=T)],1,max)	
  }else{
    b2.consolidated.dat <- do.call("cbind", lapply(allpops,function(pop){
      ACs <- apply(as.data.frame(dat[,setdiff(grep(paste(pop,"AC",sep="_"),colnames(dat),fixed=T),
                                grep(batch1,colnames(dat),fixed=T))]),
                   1, sum, na.rm=T)
      ANs <- apply(as.data.frame(dat[,setdiff(grep(paste(pop,"AN",sep="_"),colnames(dat),fixed=T),
                                grep(batch1,colnames(dat),fixed=T))]),
                   1, sum, na.rm=T)
      dtmp <- data.frame(ANs,ACs)
      colnames(dtmp) <- c(paste(pop,"_AN.ALL_OTHERS",sep=""),
                          paste(pop,"_AC.ALL_OTHERS",sep=""))
      return(dtmp)
    }))
    b2.dat <- cbind(dat[,1:3],b2.consolidated.dat)
    if(length(grep("_AC",colnames(b2.dat),fixed=T))>1){
      b2.maxAC <- apply(b2.dat[,grep("_AC",colnames(b2.dat),fixed=T)],1,max)
    }else{
      b2.maxAC <- b2.dat[,grep("_AC",colnames(b2.dat),fixed=T)]
    }
  }
  keepers.idx <- which(b1.maxAC > 0 | b2.maxAC > 0)
  b1.dat <- b1.dat[keepers.idx, ]
  b2.dat <- b2.dat[keepers.idx, ]
  
  #Iterate over variants and process each
  res <- do.call("rbind", lapply(as.character(b1.dat$VID), function(VID){
    #Find pop with largest min AN and at least one alternate allele between the two batches
    AN.bypop <- sapply(allpops, function(pop){
      min(b1.dat[which(b1.dat$VID==VID),
                 grep(paste(pop, "AN", sep="_"), colnames(b1.dat))],
          b2.dat[which(b2.dat$VID==VID),
                 grep(paste(pop, "AN", sep="_"), colnames(b2.dat))],
          na.rm=T)
    })
    AC.bypop <- sapply(allpops, function(pop){
      if(batch2 != "ALL_OTHERS"){
        max(b1.dat[which(b1.dat$VID==VID),
                   grep(paste(pop, "AC", sep="_"), colnames(b1.dat))],
            b2.dat[which(b2.dat$VID==VID),
                   grep(paste(pop, "AC", sep="_"), colnames(b2.dat))],
            na.rm=T)
      }else{
        max(b1.dat[which(b1.dat$VID==VID),
                   grep(paste(pop, "AC", sep="_"), colnames(b1.dat))],
            na.rm=T)
      }
    })
    AN.bypop[which(AC.bypop<1)] <- 0
    
    #Only process if at least one pop has min AN > min.AN
    if(any(AN.bypop > min.AN)){
      bestpop <- names(AN.bypop)[which(AN.bypop==max(AN.bypop,na.rm=T))]
      b1.AC <- as.numeric(b1.dat[which(b1.dat$VID==VID),
                                 grep(paste(bestpop,"AC",sep="_"),
                                      colnames(b1.dat),fixed=T)])
      b1.AN <- as.numeric(b1.dat[which(b1.dat$VID==VID),
                                 grep(paste(bestpop,"AN",sep="_"),
                                      colnames(b1.dat),fixed=T)])
      if(b1.AC>b1.AN){
        b1.AC <- b1.AN
      }
      b1.AF <- b1.AC/b1.AN
      b2.AC <- as.numeric(b2.dat[which(b2.dat$VID==VID),
                                 grep(paste(bestpop,"AC",sep="_"),
                                      colnames(b2.dat),fixed=T)])
      b2.AN <- as.numeric(b2.dat[which(b2.dat$VID==VID),
                                 grep(paste(bestpop,"AN",sep="_"),
                                      colnames(b2.dat),fixed=T)])
      if(b2.AC>b2.AN){
        b2.AC <- b2.AN
      }
      b2.AF <- b2.AC/b2.AN
      b1b2.p <- chisq.test(matrix(c(b1.AN-b1.AC, b1.AC,
                                    b2.AN-b2.AC, b2.AC),
                                  nrow=2, byrow=F))$p.value
      #Output row
      out.v <- data.frame("VID"=VID, "batch1"=batch1, "batch2"=batch2,
                          "pop"=bestpop, "b1.AF"=b1.AF, 
                          "b2.AF"=b2.AF, "chisq.p"=b1b2.p)
    }else{
      out.v <- data.frame("VID"=VID, "batch1"=batch1, "batch2"=batch2, 
                          "pop"=NA, "b1.AF"=NA, "b2.AF"=NA, "chisq.p"=NA)
    }
    return(out.v)
  }))
  rownames(res) <- NULL
  res[, -c(1:4)] <- apply(res[, -(1:4)], 2, as.numeric)
  res <- res[which(!is.na(res$pop)), ]
  if(return.all.stats){
    return(res)
  }else{
    return(res[which(res$chisq.p <= 0.05), ])
  }
}


### Read command-line arguments
args <- commandArgs(trailingOnly=T)
infile <- as.character(args[1])
batchlist.in <- as.character(args[2])
OUTFILE <- as.character(args[3])

# #Dev parameters:
# infile <- "~/scratch/test_AF_tables.txt.gz"
# batchlist.in <- "~/scratch/gnomAD-SV-v3.1.10pct_NCR.no_outliers.batch_fx.chr22.nonredundant_batch_pairs.txt"
# OUTFILE <- "~/scratch/test_batchFx.tsv"

### Process input data and subset to batches present in batchlist.in
dat <- read.table(infile, header=T, sep="\t", comment.char="", check.names=F)
allpops <- unique(sapply(colnames(dat)[4:ncol(dat)],
                         function(x){strsplit(as.character(x),'_')[[1]][1]}))
batches.in.header <- unique(sapply(colnames(dat)[grep("_AN.", colnames(dat), fixed=T)], 
                            function(x){unlist(strsplit(x, split="_AN.", fixed=T))[2]}))
allpairs <- read.table(batchlist.in, header=F, sep="\t")
# Subset to pairs containing at least one batch in batches.in.header
keep.pairs <- apply(allpairs, 1, function(pair){
  subpair <- setdiff(pair, c("ALL_OTHERS"))
  length(intersect(subpair, batches.in.header)) == length(subpair)
})
allpairs <- allpairs[keep.pairs, ]
allbatches <- setdiff(unique(unlist(allpairs)), c("ALL_OTHERS"))
# Subset data to batches in allpairs
batch.suffixes <- sapply(colnames(dat)[-c(1:3)], 
                         function(x){unlist(strsplit(x, split="_A[NC]\\." ))[2]})
dat <- dat[, c(1:3, which(batch.suffixes %in% allbatches)+3)]

### Find variant indexes per batch where batch AC > 0
batch.vidxs <- lapply(allbatches, function(batch){
  cidxs <- grep(paste("_AC\\.", batch, "$", sep=""), colnames(dat))
  return(which(apply(dat[, cidxs], 1, sum) > 0))
})
names(batch.vidxs) <- allbatches

### Compute pairwise AF stats for all pairs of batches and all variants
res <- apply(allpairs, 1, function(pair){
  batch1 <- as.character(pair[1])
  batch2 <- as.character(pair[2])
  if(batch2 != "ALL_OTHERS"){
    hits <- intersect(batch.vidxs[[batch1]], batch.vidxs[[batch2]])
  }else{
    hits <- batch.vidxs[[batch1]]
  }
  if(length(hits)>0){
    compare.batches(dat=dat[hits, ], batch1=batch1, batch2=batch2, 
                    allpops=allpops, return.all.stats=T)
  }
})
res <- do.call(rbind, res)
colnames(res) <- c("VID", "batch1", "batch2", "pop", "b1.AF", "b2.AF", "chisq.p")

### Write to outfile
write.table(res, OUTFILE, col.names=T, row.names=F, sep="\t", quote=F)
