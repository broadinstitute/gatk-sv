#!/usr/bin/env Rscript

# Helper script to return list of outlier samples from an svcounts.txt file

# Read command line arguments
options(stringsAsFactors=F)
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3){
  stop("Three positional arguments expected: svcounts.txt, N_IQR_cutoff, and OUTFILE")
}
countsfile <- as.character(args[1])
n.iqr <- as.integer(args[2])
OUTFILE <- as.character(args[3])

# Read counts
counts <- read.table(countsfile,header=T,sep="\t")

# Determine outlier samples
fails <- unique(as.character(unlist(sapply(unique(counts$svtype),function(svtype){
  # Get counts corresponding to svtype
  vals <- counts$count[which(counts$svtype==svtype)]
  
  # Determine cutoff
  quartiles <- as.numeric(quantile(vals,probs=c(0.25,0.75),na.rm=T))
  spacer <- n.iqr*IQR(vals,na.rm=T)
  cutoffs <- c(quartiles[1]-spacer,quartiles[2]+spacer)
  
  # Return list of sample IDs failing cutoffs
  fails <- counts$sample[which(counts$svtype==svtype
                      & (counts$count<cutoffs[1] | counts$count>cutoffs[2]))]
  return(fails)
}))))

#Write list of fail samples to OUTFILE
write.table(fails,OUTFILE,col.names=F,row.names=F,quote=F)
