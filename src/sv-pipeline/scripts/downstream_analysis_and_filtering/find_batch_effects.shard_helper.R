#!/usr/bin/env Rscript

# Talkowski SV pipeline downstream analysis helper script

# Make list of all nonredundant pairs of batches from an input list of batches


###Set global parameters
options(stringsAsFactors=F,scipen=1000)


###################
###HELPER FUNCTIONS
###################
#For any two batches, find most comparable AFs for each variant and run chi-sqared test
grep_col_by_batch<-function(dat, batch){
  columns=data.frame(colnames(dat))
  columns[,2]=apply(columns, 1, function(x){paste(strsplit(as.character(x[1]),'[.]')[[1]][2:length(strsplit(as.character(x[1]),'[.]')[[1]])], collapse='.')})
  return(which(columns[,2]==batch))
}

compare.batches <- function(dat, batch1, batch2, allpops, min.AN=30){
  #Subset data for each batch (for convenience)
  #1: restrict to sites with >0 AC in at least one batch
  b1.dat <- dat[,c(1:3,grep_col_by_batch(dat, batch1))]
  ac_cols_b1 <- grep("AC_",colnames(b1.dat),fixed=T)
  if(length(ac_cols_b1)>1){
    b1.maxAC <- apply(b1.dat[,ac_cols_b1],1,max)
  }else if(length(ac_cols_b1)==1){
    b1.maxAC <- b1.dat[,ac_cols_b1]
  }else{
    b1.maxAC <- rep(0, nrow(b1.dat))
  }
  if(batch2 != "ALL_OTHERS"){
    b2.dat <- dat[,c(1:3,grep_col_by_batch(dat, batch2))]
    ac_cols_b2 <- grep("AC_",colnames(b2.dat),fixed=T)
    if(length(ac_cols_b2)>1){
      b2.maxAC <- apply(b2.dat[,ac_cols_b2],1,max)	
    }else if(length(ac_cols_b2)==1){
      b2.maxAC <- b2.dat[,ac_cols_b2]
    }else{
      b2.maxAC <- rep(0, nrow(b2.dat))
    }
  }else{
    b2.consolidated.dat <- do.call("cbind", lapply(allpops,function(pop){
      ACs <- apply(as.data.frame(dat[,setdiff(grep(paste("AC",pop,sep="_"),colnames(dat),fixed=T),
                                grep(batch1,colnames(dat),fixed=T))]),
                   1, sum, na.rm=T)
      ANs <- apply(as.data.frame(dat[,setdiff(grep(paste("AN",pop,sep="_"),colnames(dat),fixed=T),
                                grep(batch1,colnames(dat),fixed=T))]),
                   1, sum, na.rm=T)
      dtmp <- data.frame(ANs,ACs)
      colnames(dtmp) <- c(paste("AN",pop,"ALL_OTHERS",sep="_"),
                          paste("AC",pop,"ALL_OTHERS",sep="_"))
      return(dtmp)
    }))
    b2.dat <- cbind(dat[,1:3],b2.consolidated.dat)
    ac_cols_b2_all <- grep("AC",colnames(b2.dat),fixed=T)
    if(length(ac_cols_b2_all)>1){
      b2.maxAC <- apply(b2.dat[,ac_cols_b2_all],1,max)
    }else if(length(ac_cols_b2_all)==1){
      b2.maxAC <- b2.dat[,ac_cols_b2_all]
    }else{
      b2.maxAC <- rep(0, nrow(b2.dat))
    }
  }
  keepers.idx <- which(b1.maxAC > 0 | b2.maxAC > 0)
  b1.dat <- b1.dat[keepers.idx, ]
  b2.dat <- b2.dat[keepers.idx, ]
  
  #Iterate over variants and process each
  res <- do.call("rbind", lapply(as.character(b1.dat$VID),function(VID){
    #Find pop with largest min AN and at least one alternate allele between the two batches
    AN.bypop <- sapply(allpops,function(pop){
      min(b1.dat[which(b1.dat$VID==VID),
                 grep(paste("AN",pop,sep="_"),colnames(b1.dat))],
          b2.dat[which(b2.dat$VID==VID),
                 grep(paste("AN",pop,sep="_"),colnames(b2.dat))],
          na.rm=T)
    })
    AC.bypop <- sapply(allpops,function(pop){
      max(b1.dat[which(b1.dat$VID==VID),
                 grep(paste("AC",pop,sep="_"),colnames(b1.dat))],
          b2.dat[which(b2.dat$VID==VID),
                 grep(paste("AC",pop,sep="_"),colnames(b2.dat))],
          na.rm=T)
    })
    AN.bypop[which(AC.bypop<1)] <- 0
    
    #Only process if at least one pop has min AN > min.AN
    if(any(AN.bypop>min.AN)){
      bestpop <- names(AN.bypop)[which(AN.bypop==max(AN.bypop,na.rm=T))]
      b1.AC <- as.numeric(b1.dat[which(b1.dat$VID==VID),
                                 grep(paste("AC",bestpop,sep="_"),colnames(b1.dat),fixed=T)])
      b1.AN <- as.numeric(b1.dat[which(b1.dat$VID==VID),
                                 grep(paste("AN",bestpop,sep="_"),colnames(b1.dat),fixed=T)])
      if(b1.AC>b1.AN){
        b1.AC <- b1.AN
      }
      b1.AF <- b1.AC/b1.AN
      b2.AC <- as.numeric(b2.dat[which(b2.dat$VID==VID),
                                 grep(paste("AC",bestpop,sep="_"),colnames(b2.dat),fixed=T)])
      b2.AN <- as.numeric(b2.dat[which(b2.dat$VID==VID),
                                 grep(paste("AN",bestpop,sep="_"),colnames(b2.dat),fixed=T)])
      if(b2.AC>b2.AN){
        b2.AC <- b2.AN
      }
      b2.AF <- b2.AC/b2.AN
      b1b2.p <- chisq.test(matrix(c(b1.AN-b1.AC,b1.AC,
                                    b2.AN-b2.AC,b2.AC),
                                  nrow=2,byrow=F))$p.value
      #Output row
      out.v <- data.frame("VID"=VID,"pop"=bestpop,"b1.AF"=b1.AF,"b2.AF"=b2.AF,"chisq.p"=b1b2.p)
    }else{
      out.v <- data.frame("VID"=VID,"pop"=NA,"b1.AF"=NA,"b2.AF"=NA,"chisq.p"=NA)
    }
    return(out.v)
  }))
  rownames(res) <- NULL
  res[,-c(1:2)] <- apply(res[,-(1:2)],2,as.numeric)
  res <- res[which(!is.na(res$pop)), ]
  # res$chisq.bonf <- p.adjust(res$chisq.p,method="bonferroni")
  return(res)
}


###Read command-line arguments
args <- commandArgs(trailingOnly=T)
infile <- as.character(args[1])
batch1 <- as.character(args[2])
batch2 <- as.character(args[3])
OUTFILE <- as.character(args[4])

# #Dev parameters:
# infile <- "~/scratch/minGQ_test/freq_table_shard_00010.txt.gz"
# batch1 <- "synthbatch1"
# batch2 <- "ALL_OTHERS"
# OUTFILE <- "~/scratch/test_batchFx.tsv"

###Process data & write output
dat <- read.table(infile,header=T,sep="\t",comment.char="")
allpops <- unique(sapply(colnames(dat)[4:ncol(dat)],
                         function(x){
                           # Extract population from column names like AN_diverse.Batch1 -> diverse
                           parts <- strsplit(as.character(x),'_')[[1]]
                           if(length(parts) >= 2) {
                             # Get the part after the first underscore and before the dot
                             pop_batch <- parts[2]
                             pop <- strsplit(pop_batch, '\\.')[[1]][1]
                             return(pop)
                           } else {
                             return(parts[1])
                           }
                         }))
res <- compare.batches(dat=dat, batch1=batch1, batch2=batch2, allpops=allpops)
write.table(res,OUTFILE,col.names=T,row.names=F,sep="\t",quote=F)

