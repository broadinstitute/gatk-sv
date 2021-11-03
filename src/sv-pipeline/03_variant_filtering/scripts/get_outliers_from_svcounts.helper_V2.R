#!/usr/bin/env Rscript

# Helper script to return list of outlier samples from an svcounts.txt file
# an example of svcounts.txt these columns: svtype count
# CutoffFile should include these columns: algorithm       svtype  lower_cff       higher_cff


# Read command line arguments
options(stringsAsFactors=F)
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=4){
  stop("Four positional arguments expected: svcounts.txt, CutoffFile, OUTFILE and algorithm")
}
CountsFile <- as.character(args[1])
CutoffFile<- as.character(args[2])
OutFile <- as.character(args[3])
Algorithm <- as.character(args[4])

main<-function(CountsFile, CutoffFile,OutFile,  Algorithm){
  # Read counts
  counts <- read.table(CountsFile,header=T,sep="\t")
  cffs <-read.table(CutoffFile, header=T, sep='\t')
  # Determine outlier samples
  fails <- unique(as.character(unlist(sapply(unique(counts$svtype),function(svtype){
    # Get counts corresponding to svtype
    vals <- counts$count[which(counts$svtype==svtype)]
    
    # Determine cutoff
    quartiles <- as.numeric(quantile(vals,probs=c(0.25,0.75),na.rm=T))

    cutoffs <- c(cffs[cffs$algorithm==Algorithm & cffs$svtype == svtype,c('lower_cff','higher_cff')])
    # Return list of sample IDs failing cutoffs
    fails <- counts$sample[which(counts$svtype==svtype
                        & (counts$count<cutoffs[1] | counts$count>cutoffs[2]))]
    return(fails)
  }))))

  #Write list of fail samples to OutFile
  write.table(fails,OutFile,col.names=F,row.names=F,quote=F)
}

main(CountsFile, CutoffFile,OutFile,  Algorithm)