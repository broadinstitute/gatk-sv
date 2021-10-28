#!/usr/bin/env Rscript

# Script to correct for GC content across a list of binCov files
# Note: Designed to be run on 24 chromosome-split files from the same library

#Requires tab-delimmed input file with the following columns:
# 1: full paths to binCov files for each contig
# 2: corresponding GC masks for each contig
# 3: full path to desired output file for each contig

####################################
#####Set parameters & load libraries
####################################
options(scipen=1000,stringsAsFactors=F)

###############################################################################
#####Helper function to read a list of binCov files and return values as a list
###############################################################################
#Note: desiged to iterate over per-chromosome binCov files for a single sample
readBinCovList <- function(paths,labels=c(1:22,"X","Y"),
                           norm=F,exclude=c("X","Y"),
                           returnCoords=F){
  #Iterate over paths and read just the fourth column unless otherwise optioned
  dat <- lapply(paths,function(path){
    dat <- as.data.frame(read.table(path,header=F,sep="\t"))
    colnames(dat) <- c("chr","start","end","cov")
    if(returnCoords==F){
      return(dat$cov)
    }else{
      return(as.data.frame(dat[,1:3]))
    }
  })
  names(dat) <- labels

  #Normalize as median of all bins, if optioned
  if(returnCoords==F & norm==T){
    median.all <- median(unlist(sapply(which(!(names(dat) %in% exclude)),function(i){
      return(dat[[i]])
    })),na.rm=T)
    dat <- lapply(dat,function(vals){
      return(vals/median.all)
    })
    return(dat)
  }

  #Return values
  return(dat)
}

################################################################
#####Helper function to compute mean coverage from a binCov list
################################################################
binCovListMean <- function(binCovList,exclude=c("X","Y")){
  #Convert all non-exclude contigs to single giant vector
  vals <- as.numeric(unlist(binCovList[which(!(names(binCovList) %in% exclude))]))
  #Find middle 95% of data after excluding bins with no coverage
  cutoffs <- quantile(vals[which(vals>0)],probs=c(0.025,0.975))
  #Compute mean on middle 95% of data
  covMean <- mean(vals[which(vals>cutoffs[1] & vals<cutoffs[2])],na.rm=T)
  #Return mean
  return(covMean)
}


#################################################################
#####Helper function to calculate coverage adjustments vs GC bins
#################################################################
GCBinAdjustments <- function(vals,GC,
                             scalingMean=NULL,scalingSD=NULL,
                             do.scaleSD=T){
  #Round GC to integers
  GC <- as.integer(as.character(round(as.numeric(100*GC))))

  #Determine cutoffs for middle 95% of coverage data and GC data
  # after excluding bins with zero coverage
  cov.thresh <- quantile(vals[which(vals>0)],probs=c(0.025,0.975),na.rm=T)
  GC.thresh <- quantile(GC,probs=c(0.025,0.975),na.rm=T)

  #Compute statistics for each GC bin
  GCstats <- as.data.frame(t(sapply(seq(0,100,1),function(GCbin){
    binVals <- vals[which(GC==GCbin)]
    #Only computes binwise mean and SD on the middle 95%
    # of the data to protect against extreme outliers
    binMean <- mean(binVals[which(binVals>=cov.thresh[1] & binVals<=cov.thresh[2])],na.rm=T)
    binSD <- sd(binVals[which(binVals>=cov.thresh[1] & binVals<=cov.thresh[2])],na.rm=T)
    nBins <- length(binVals[which(binVals>=cov.thresh[1] & binVals<=cov.thresh[2])])
    return(c(GCbin,binMean,binSD,nBins))
  })))
  colnames(GCstats) <- c("GCbin","binMean","binSD","nBins")

  #Get mean & SD for scaling, if not specified
  if(is.null(scalingMean)){
    #Only consider middle 95% of GC values to
    # protect against outliers
    scalingMean <- mean(vals[which(GC>=GC.thresh[1] & GC<=GC.thresh[2] & vals>=cov.thresh[1] & vals<=cov.thresh[2])],na.rm=T)
  }
  if(do.scaleSD==T && is.null(scalingSD)){
    #Only consider weighted mean of bin SDs from
    # middle 95% of GC values to protect against extreme outliers
    GC.dat.tmp <- GCstats[which(GCstats[,1]>=GC.thresh[1] & GCstats[,1]<=GC.thresh[2]),c(3,4)]
    scalingSD <- sum(GC.dat.tmp[,1]*GC.dat.tmp[,2],na.rm=T)/sum(GC.dat.tmp[,2],na.rm=T)
  }

  #Iterate over all values and adjust vs. GCstats
  #Adjust per-bin, then reorder based on original order
  newVals <- vals
  for(GCbin in GCstats$GCbin){
    #Get all values corresponding to each bin
    binVals <- newVals[which(GC==GCbin)]
    if(length(binVals)>0){
      #Compute mean after excluding outliers (Median±1.5*IQR)
      cutoff <- c(median(binVals[which(binVals>0)])-1.5*IQR(binVals[which(binVals>0)]),
                  median(binVals[which(binVals>0)])+1.5*IQR(binVals[which(binVals>0)]))
      binMean <- mean(binVals[which(binVals>=cutoff[1] & binVals<=cutoff[2])],na.rm=T)
      #Center binVals
      if(!is.na(binMean)){
        binVals <- binVals-binMean
      }
      #Scale residuals by ratio of scalingSD/binSD, if optioned
      if(do.scaleSD==T){
        binSD <- GCstats[which(GCstats[,1]==GCbin),3]
        if(!is.na(binSD)){
          #Require a positive SD and at least 100 bins before scaling values based on SD
          if(binSD>0 & length(binVals>=100)){
            binVals <- binVals*(scalingSD/binSD)
          }
        }
      }
      #Recenter binVals at scalingMean
      binVals <- binVals+scalingMean
      #Replace old values in vals with binVals
      newVals[which(GC==GCbin)] <- binVals
    }
  }

  #If original value was zero, rewrite as zero
  newVals[which(vals==0)] <- 0

  #Compute adjustments
  adj <- newVals-vals

  #Return new values
  return(adj)
}

############################################################################
#####Helper function to smooth adjustments as 5-window weighted rolling mean
############################################################################
wrmean <- function(vals,weights=c(0.1,0.2,1,0.2,0.1)){
  #Create data frame of values for smoothing
  smooth.df <- data.frame("pre2"=c(NA,NA,vals[1:(length(vals)-2)]),
                          "pre1"=c(NA,vals[1:(length(vals)-1)]),
                          "center"=vals,
                          "post1"=c(vals[2:length(vals)],NA),
                          "post1"=c(vals[3:length(vals)],NA,NA))

  #Smooth each row of smooth.df
  smoothed.vals <- apply(smooth.df,1,function(rowvals){
    rowsum <- sum(rowvals*weights,na.rm=T)
    rowdenom <- sum(weights[which(!is.na(rowvals))])
    smoothed <- rowsum/rowdenom
    return(smoothed)
  })

  #If original value was zero, rewrite as zero
  smoothed.vals[which(vals==0)] <- 0

  #Return smoothed vals
  return(smoothed.vals)
}

##########################
#####Rscript functionality
##########################
require(optparse)
#List of command-line options
option_list <- list(
  make_option(c("-z", "--gzip"),action="store_false",default=TRUE,
              help="gzip output files [default: TRUE]")
)

#Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog [options] input.list",
                                option_list=option_list),
                   positional_arguments=TRUE)
INFILE <- args$args[1]
gzip <- args$options$gzip

#Checks for appropriate positional arguments
if(length(args$args) != 1){
  stop("Must supply an input list\n")
}

#Read input list
INFILE <- read.table(INFILE,header=F,sep="\t")
contigs <- as.character(INFILE[,1])
cov.paths <- as.character(INFILE[,2])
GC.paths <- as.character(INFILE[,3])
out.paths <- as.character(INFILE[,4])

#Remove gz extension from output paths if necessary
out.paths <- gsub(".gz$","",out.paths,fixed=F,ignore.case=F)

#Read coverage files
cov.dat <- readBinCovList(cov.paths,norm=F,labels=contigs,returnCoords=F)
cov.coords <- readBinCovList(cov.paths,labels=contigs,returnCoords=T)
GC.dat <- readBinCovList(GC.paths,labels=contigs,norm=F)

#Calculate library-wide mean coverage for autosomes only
cov.mean <- binCovListMean(cov.dat,exclude=c("X","Y"))

#Run GC correction, smooth, and round to nearest whole integer for all autosomes
cov.dat.adj <- cov.dat
cov.dat.adj[1:length(cov.dat)] <- lapply(1:length(cov.dat),function(i){
  #Compute coverage residuals
  #Allosomes: scale to mean coverage of that contig
  #Otherwise: scale to library-wide autosomal mean
  if(names(cov.dat)[i] %in% c("X","Y")){
    vals.residuals <- GCBinAdjustments(vals=cov.dat[[i]],
                                       GC=GC.dat[[i]],
                                       scalingMean=NULL,
                                       scalingSD=NULL,do.scaleSD=T)
  }else{
    vals.residuals <- GCBinAdjustments(vals=cov.dat[[i]],
                                       GC=GC.dat[[i]],
                                       scalingMean=cov.mean,
                                       scalingSD=NULL,do.scaleSD=T)
  }
  #Smooth coverage residuals
  vals.residuals.smoothed <- wrmean(vals.residuals,weights=c(0.1,0.2,1,0.2,0.1))
  #Adjust coverage
  vals.adj <- cov.dat[[i]]+vals.residuals.smoothed
  #Round values to nearest whole integer
  vals.adj.round <- round(vals.adj)
  #Do not allow values < 0
  vals.adj.round[which(vals.adj.round<0)] <- 0
  #Return adjusted, rounded values
  return(vals.adj.round)
})

#Iterate over contigs, append coordinates to normalized coverage values, and write to file
sapply(1:length(cov.dat),function(i){
  #Prepare output data frame
  out.df <- cbind(cov.coords[[i]],
                  "cov"=cov.dat.adj[[i]])
  #Write to file
  write.table(out.df,out.paths[i],
              col.names=F,row.names=F,
              sep="\t",quote=F)
  #Gzip, if optioned
  if(gzip==T){
    system(paste("gzip -f ",out.paths[i],sep=""))
  }
})
