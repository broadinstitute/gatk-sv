#!/usr/bin/env Rscript

# Script to train WGD dosage bias scoring model

# Requires:
#  An input set of bins to be tested
#  A coverage matrix with all training & testing samples corresponding to the bins to be tested
#  1-5 pairs of training sample sets (pairs = one list each for PCR+/PCR-)
#  A pair of testing sample sets


####################################
#####Set parameters & load libraries
####################################
options(stringsAsFactors=F, 
        scipen=1000)
require(optparse)


# #Local dev parameters
# WRKDIR <- "/Users/rlc/Desktop/Collins/Talkowski/NGS/SV_Projects/gnomAD/"
# DATADIR <- paste(WRKDIR, "WGD_training_localData/", sep="")
# path.to.matrix <- paste(DATADIR, "gnomAD_v2.6F_adjCov.WGD_scoring_masked.100bp.matrix.bed.gz", sep="")
# PLUS.train.1.in <- paste(DATADIR, "PCRPLUS.train_1.samples.list", sep="")
# PLUS.train.2.in <- paste(DATADIR, "PCRPLUS.train_2.samples.list", sep="")
# PLUS.train.3.in <- paste(DATADIR, "PCRPLUS.train_3.samples.list", sep="")
# PLUS.test.in <- paste(DATADIR, "PCRPLUS.test.samples.list", sep="")
# MINUS.train.1.in <- paste(DATADIR, "PCRFREE.train_1.samples.list", sep="")
# MINUS.train.2.in <- paste(DATADIR, "PCRFREE.train_2.samples.list", sep="")
# MINUS.train.3.in <- paste(DATADIR, "PCRFREE.train_3.samples.list", sep="")
# MINUS.test.in <- paste(DATADIR, "PCRFREE.test.samples.list", sep="")
# OUTDIR <- "~/scratch/"
# plot <- T
# args <- list("args"=c(path.to.matrix, OUTDIR), 
#              "options"=list("PCRPLUS_train_1"=PLUS.train.1.in, 
#                             "PCRPLUS_train_2"=PLUS.train.2.in, 
#                             "PCRPLUS_train_3"=PLUS.train.3.in, 
#                             "PCRPLUS_test_1"=PLUS.test.in, 
#                             "PCRMINUS_train_1"=MINUS.train.1.in, 
#                             "PCRMINUS_train_2"=MINUS.train.2.in, 
#                             "PCRMINUS_train_3"=MINUS.train.3.in, 
#                             "PCRMINUS_test_1"=MINUS.test.in, 
#                             "plot"=plot))
# opts <- args$options

# #Local dev parameters (hg38; Harold SFARI)
# WRKDIR <- "~/scratch/hg38_WGD_model_tmp/"
# DATADIR <- "~/scratch/hg38_WGD_model_files/"
# path.to.matrix <- paste(DATADIR, "test.bed.gz", sep="")
# PLUS.train.1.in <- paste(DATADIR, "pcr+training.txt", sep="")
# MINUS.train.1.in <- paste(DATADIR, "pcr-train1.txt", sep="")
# MINUS.train.2.in <- paste(DATADIR, "pcr-train2.txt", sep="")
# MINUS.train.3.in <- paste(DATADIR, "pcr-train3.txt", sep="")
# PLUS.test.in <- paste(DATADIR, "pcr+testing.txt", sep="")
# MINUS.test.in <- paste(DATADIR, "pcr-test.txt", sep="")
# OUTDIR <- WRKDIR
# plot <- T
# args <- list("args"=c(path.to.matrix, OUTDIR), 
#              "options"=list("PCRPLUS_train_1"=PLUS.train.1.in, 
#                             "PCRPLUS_test_1"=PLUS.test.in, 
#                             "PCRMINUS_train_1"=MINUS.train.1.in, 
#                             "PCRMINUS_train_2"=MINUS.train.2.in, 
#                             "PCRMINUS_train_3"=MINUS.train.3.in, 
#                             "PCRMINUS_test_1"=MINUS.test.in, 
#                             "plot"=plot))
# opts <- args$options


#####################
#####Helper functions
#####################
#Read & clean binCov matrix
readCov <- function(path, norm=T, tranche=0.999){
  #Read matrix
  cov <- read.table(path.to.matrix, header=T, comment.char="")
  colnames(cov)[1:3] <- c("chr", "start", "end")

  #Normalize per sample by median, if optioned
  if(norm==T){
    cov[, -c(1:3)] <- apply(cov[, -c(1:3)], 2, function(vals){
      #Do not include bins with zero coverage in normalization
      vals[which(vals>0)] <- vals[which(vals>0)]/median(vals[which(vals>0)], na.rm=T)
      vals <- 2*vals
      return(vals)
    })
  }

  #Exclude outlier bins, if optioned
  if(!is.na(tranche)){
    #Get percentiles to exclude
    tails <- c((1-tranche)/2, mean(c(1, tranche)))

    #Remove bins with median cov = 0
    bin.medians <- apply(cov[, -c(1:3)], 1, median, na.rm=T)
    cov <- cov[which(bin.medians>0), ]

    #Remove bins in tails of mean cov distribution
    bin.means <- apply(cov[, -c(1:3)], 1, mean, na.rm=T)
    cutoffs <- quantile(x=bin.means, probs=tails)
    cov <- cov[which(bin.means>=cutoffs[1] & bin.means<=cutoffs[2]), ]
  }

  #Return cleaned coverage
  return(cov)
}

#Gather means and separation between two sets of coverage values
covStatsSingle <- function(valsA, valsB){
  #Compute means & difference
  mA <- mean(valsA, na.rm=T)
  mB <- mean(valsB, na.rm=T)
  sep <- mA-mB

  #Format & return results
  return(c(mA, mB, sep))
}

#Compare coverage distributions for two lists of samples
covTest <- function(lA, lB, cov, min.sep=0.1){
  #Get indexes
  iA <- which(colnames(cov) %in% lA)
  iB <- which(colnames(cov) %in% lB)

  #Iterate over all bins and get stats
  stats <- as.data.frame(t(sapply(1:nrow(cov), function(i){
    covStatsSingle(as.numeric(cov[i, iA]), 
                   as.numeric(cov[i, iB]))
  })))
  colnames(stats) <- c("meanA", "meanB", "sep")
  stats$pval <- NA

  #Determine which bins meet criteria for t-test
  test.bins <- which(abs(stats$sep)>=min.sep &
                       ((stats$meanA>=2 & stats$meanB<=2) |
                          (stats$meanA<=2 & stats$meanB>=2)))

  #Iterate over all bins in cov & run t-test if conditions are met
  stats$pval[test.bins] <- sapply(test.bins, function(i){
    return(t.test(as.numeric(cov[i, iA]), as.numeric(cov[i, iB]))$p.value)
  })

  #Return stats
  return(stats)
}

#Plot single scatter plot, with options
plotScatterSingle <- function(xvals.all, yvals.all, 
                              xvals.final, yvals.final, 
                              highlight="orangered", 
                              axis.lims=c(1, 3), 
                              spline=F, cor=T, 
                              xaxis.bottom=F, xaxis.top=F, 
                              yaxis.left=F, yaxis.right=F){
  #Prep plot area
  if((xaxis.bottom==T & xaxis.top==T) |
     (yaxis.left==T & yaxis.right==T)){
    par(mar=rep(1.5, 4))
  }else{
    if(xaxis.bottom==T){
      xmar <- c(1.5, 0.5)
    }else{
      if(xaxis.top==T){
        xmar <- c(0.5, 1.5)
      }else{
        xmar <- c(1, 1)
      }
    }
    if(yaxis.left==T){
      ymar <- c(1.5, 0.5)
    }else{
      if(yaxis.right==T){
        ymar <- c(0.5, 1.5)
      }else{
        ymar <- c(1, 1)
      }
    }
    par(mar=c(xmar[1], ymar[1], xmar[2], ymar[2]))
  }
  plot(x=xvals.all, y=yvals.all, type="n", 
       xaxt="n", yaxt="n", xlab="", ylab="", 
       xlim=axis.lims, ylim=axis.lims)
  abline(h=mean(axis.lims), v=mean(axis.lims), col="gray30")

  #Prep gradient based on distance percentiles
  gradientCols <- rev(colorRampPalette(c("#003A94", "#3361A9", "#6689BF", 
                                         "#99B0D4", "#CCD8EA", "white"))(101))
  center <- c(mean(xvals.all, na.rm=T), mean(yvals.all, na.rm=T))
  point.dist <- sapply(1:length(xvals.all), function(i){
    sqrt(((xvals.all[i]-center[1])^2) + ((yvals.all[i]-center[2])^2))
  })
  point.dist.pct <- round(100*rank(point.dist)/length(point.dist))
  point.dist.cols <- unlist(sapply(point.dist.pct, function(k){
    return(gradientCols[k+1])
  }))

  #Add points & crosshairs
  points(xvals.all, yvals.all, col=point.dist.cols, pch=19, cex=0.2)
  if(spline==T){
    points(smooth.spline(xvals.all, yvals.all), type="l")
  }
  points(center[1], center[2], pch=9)
  points(xvals.final, yvals.final, col=highlight, pch=19, cex=0.4)

  #Add axes, if optioned
  if(xaxis.bottom==T){
    axis(1, at=axTicks(1), labels=NA)
    axis(1, at=axTicks(1), tick=F, line=-0.4, cex.axis=0.8)
  }
  if(xaxis.top==T){
    axis(3, at=axTicks(3), labels=NA)
    axis(3, at=axTicks(3), tick=F, line=-0.4, cex.axis=0.8)
  }
  if(yaxis.left==T){
    axis(2, at=axTicks(2), labels=NA)
    axis(2, at=axTicks(2), tick=F, line=-0.4, cex.axis=0.8, las=2)
  }
  if(yaxis.right==T){
    axis(4, at=axTicks(3), labels=NA)
    axis(4, at=axTicks(2), tick=F, line=-0.4, cex.axis=0.8, las=2)
  }

  #Add correlation coefficient to plot, if optioned
  if(cor==T){
    text(x=par("usr")[1], 
         y=par("usr")[3]+(0.875*(par("usr")[4]-par("usr")[3])), 
         pos=4, labels=paste("R=", round(cor.test(xvals.all, yvals.all)$estimate, digits=3), "\n", 
                            "Rho=", round(cor.test(xvals.all, yvals.all, method="spearman")$estimate, digits=3), sep=""), 
         cex=0.8)
  }
}

#Plot evidence for a single bin (series of swarmplots)
plotBinEvidence <- function(chr, start, end, title=NULL, 
                            col.PLUS="blue", col.MINUS="red", 
                            ylims=c(1, 3)){
  #Load required library
  require(beeswarm)
  require(vioplot)

  #Make master list of groups to plot
  plot.groups <- c(lapply(PLUS.train, as.character), lapply(MINUS.train, as.character))
  names(plot.groups) <- c(paste("PLUS.train.", 1:length(PLUS.train), sep=""), 
                          paste("MINUS.train.", 1:length(MINUS.train), sep=""))
  if(!is.null(PLUS.test)){
    plot.groups <- c(plot.groups, list(as.character(PLUS.test)))
    names(plot.groups)[length(plot.groups)] <- "PLUS.test"
  }
  if(!is.null(MINUS.test)){
    plot.groups <- c(plot.groups, list(as.character(MINUS.test)))
    names(plot.groups)[length(plot.groups)] <- "MINUS.test"
  }
  n.groups <- length(plot.groups)

  #Prep plot area
  par(mar=c(2, 3, 2, 1))
  plot(x=c(0, n.groups), y=ylims, type="n", 
       xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i")

  #Add gridlines
  abline(h=axTicks(2), col="gray80")
  abline(h=2)

  #Iterate over groups & plot swarms per group
  sapply(1:n.groups, function(i){
    #Get coverage values
    vals <- as.numeric(cov[which(cov$chr==chr & cov$start==start & cov$end==end), 
                       which(colnames(cov) %in% plot.groups[[i]])])

    #Get color
    if(length(grep("PLUS", names(plot.groups)[i], fixed=T))>0){
      plot.col <- col.PLUS
    }else{
      plot.col <- col.MINUS
    }

    #Plot vioplot & swarm
    vioplot(vals, add=T, at=i-0.5, drawRect=F, col="white", wex=0.5, border=plot.col)
    beeswarm(vals, add=T, at=i-0.5, corral="wrap", corralWidth=0.8, 
             pch=19, cex=0.7, col=plot.col)

    #Add mean per group
    points(x=i-0.5, y=mean(vals, na.rm=T), pch=23, bg="white", lwd=4, cex=2)
    points(x=i-0.5, y=mean(vals, na.rm=T), pch=18, col=plot.col)
  })

  #Add axis labels & title
  axis(1, at=(1:n.groups)-0.5, tick=F, cex.axis=0.8, 
       labels=names(plot.groups), line=-0.9)
  axis(2, at=axTicks(2), labels=NA)
  axis(2, at=axTicks(2), tick=F, line=-0.4, las=2, cex.axis=0.8)
  mtext(2, text="Estimated Copy Number", line=2)
  mtext(3, line=0, text=title, font=2)
}


######################
#####Rscript execution
######################
####List of command-line options
option_list <- list(
  make_option(c("--PCRPLUS_train_1"), type="character", default=NULL, 
              help="path to list of first PCRPLUS training samples [REQUIRED; default %default]", 
              metavar="character"), 
  make_option(c("--PCRPLUS_train_2"), type="character", default=NULL, 
              help="path to list of second PCRPLUS training samples [OPTIONAL; default %default]", 
              metavar="character"), 
  make_option(c("--PCRPLUS_train_3"), type="character", default=NULL, 
              help="path to list of third PCRPLUS training samples [OPTIONAL; default %default]", 
              metavar="character"), 
  make_option(c("--PCRPLUS_train_4"), type="character", default=NULL, 
              help="path to list of fourth PCRPLUS training samples [OPTIONAL; default %default]", 
              metavar="character"), 
  make_option(c("--PCRPLUS_train_5"), type="character", default=NULL, 
              help="path to list of fifth PCRPLUS training samples [OPTIONAL; default %default]", 
              metavar="character"), 
  make_option(c("--PCRPLUS_test"), type="character", default=NULL, 
              help="path to list of PCRPLUS testing samples [OPTIONAL; default %default]", 
              metavar="character"), 
  make_option(c("--PCRMINUS_train_1"), type="character", default=NULL, 
              help="path to list of first PCRMINUS training samples [REQUIRED; default %default]", 
              metavar="character"), 
  make_option(c("--PCRMINUS_train_2"), type="character", default=NULL, 
              help="path to list of second PCRMINUS training samples [OPTIONAL; default %default]", 
              metavar="character"), 
  make_option(c("--PCRMINUS_train_3"), type="character", default=NULL, 
              help="path to list of third PCRMINUS training samples [OPTIONAL; default %default]", 
              metavar="character"), 
  make_option(c("--PCRMINUS_train_4"), type="character", default=NULL, 
              help="path to list of fourth PCRMINUS training samples [OPTIONAL; default %default]", 
              metavar="character"), 
  make_option(c("--PCRMINUS_train_5"), type="character", default=NULL, 
              help="path to list of fifth PCRMINUS training samples [OPTIONAL; default %default]", 
              metavar="character"), 
  make_option(c("--PCRMINUS_test"), type="character", default=NULL, 
              help="path to list of PCRMINUS testing samples [OPTIONAL; default %default]", 
              metavar="character"), 
  make_option(c("--plot"), type="logical", default=TRUE, 
              help="generate scatterplots of coverage correlations [default %default]", 
              metavar="logical")
)
####Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog [options] covMatrix.bed OUTDIR", 
                                option_list=option_list), 
                   positional_arguments=TRUE)
opts <- args$options

####Checks for appropriate required arguments
#Positional args
if(length(args$args) != 2){
  stop("Incorrect number of required positional arguments\n")
}
#At least one training set for PCRPLUS and PCRMINUS
if(any(is.null(c(opts$PCRPLUS_train_1, opts$PCRMINUS_train_1)))){
  stop("At least one training set must be supplied for PCRPLUS and PCRMINUS\n")
}

#####Writes args & opts to vars
path.to.matrix <- args$args[1]
OUTDIR <- args$args[2]
PLUS.train.1.in <- opts$PCRPLUS_train_1
PLUS.train.2.in <- opts$PCRPLUS_train_2
PLUS.train.3.in <- opts$PCRPLUS_train_3
PLUS.train.4.in <- opts$PCRPLUS_train_4
PLUS.train.5.in <- opts$PCRPLUS_train_5
PLUS.test.in <- opts$PCRPLUS_test
MINUS.train.1.in <- opts$PCRMINUS_train_1
MINUS.train.2.in <- opts$PCRMINUS_train_2
MINUS.train.3.in <- opts$PCRMINUS_train_3
MINUS.train.4.in <- opts$PCRMINUS_train_4
MINUS.train.5.in <- opts$PCRMINUS_train_5
MINUS.test.in <- opts$PCRMINUS_test
plot <- opts$plot

#Create OUTDIR, if necessary
if(!dir.exists(OUTDIR)){
  dir.create(OUTDIR)
}

#####Reads data
#Coverage matrix
cov <- readCov(path.to.matrix, norm=T, tranche=0.999)
#PCRPLUS training data
PLUS.train <- lapply(list(PLUS.train.1.in, PLUS.train.2.in, PLUS.train.3.in, PLUS.train.4.in, PLUS.train.5.in), 
                     function(path){
                       if(!is.null(path)){
                         return(read.table(path, header=F)[, 1])
                       }else{
                         return(NULL)
                       }
                     })
PLUS.train <- PLUS.train[which(!unlist(lapply(PLUS.train, is.null)))]
#PCRMINUS training data
MINUS.train <- lapply(list(MINUS.train.1.in, MINUS.train.2.in, MINUS.train.3.in, MINUS.train.4.in, MINUS.train.5.in), 
                      function(path){
                        if(!is.null(path)){
                          return(read.table(path, header=F)[, 1])
                        }else{
                          return(NULL)
                        }
                      })
MINUS.train <- MINUS.train[which(!unlist(lapply(MINUS.train, is.null)))]
#PCRPLUS testing data
if(!is.null(PLUS.test.in)){
  PLUS.test <- read.table(PLUS.test.in, header=F)[, 1]
}
#PCRMINUS testing data
if(!is.null(MINUS.test.in)){
  MINUS.test <- read.table(MINUS.test.in, header=F)[, 1]
}


#####Run pairwise comparisons of training lists
train.res <- lapply(1:length(PLUS.train), function(iPLUS){
  lapply(1:length(MINUS.train), function(iMINUS){
    covTest(lA=PLUS.train[[iPLUS]], 
            lB=MINUS.train[[iMINUS]], 
            cov=cov, 
            min.sep=0.1)
  })
})

#####Extract bins that meet selection criteria
#Get mean coverages per training set
meanCovs.PLUS <- as.data.frame(sapply(1:length(PLUS.train), function(i){
  train.res[[i]][[1]]$meanA
}))
colnames(meanCovs.PLUS) <- paste("PLUS.train.", 1:length(PLUS.train), ".meanCov", sep="")
meanCovs.MINUS <- as.data.frame(sapply(1:length(MINUS.train), function(i){
  train.res[[1]][[i]]$meanB
}))
colnames(meanCovs.MINUS) <- paste("MINUS.train.", 1:length(MINUS.train), ".meanCov", sep="")
meanCovs <- cbind(meanCovs.PLUS, meanCovs.MINUS)
#Determine which bins pass meanCov consistency check
meanCovs.keep.bins <- apply(meanCovs, 1, function(vals){
  #Determine column indexes for PLUS and MINUS training sets
  plus.cols <- grep("PLUS", colnames(meanCovs))
  minus.cols <- grep("MINUS", colnames(meanCovs))
  if((all(vals[plus.cols]>=2) & all(vals[minus.cols<=2])) |
     (all(vals[plus.cols]<=2) & all(vals[minus.cols>=2]))){
    return(TRUE)
  }else{
    return(FALSE)
  }
})
meanCovs.keep.bins <- which(meanCovs.keep.bins==T)

#Get minimum absolute separation between any two training sets
minAbsSep <- sapply(1:nrow(cov), function(i){
  vals <- sapply(1:length(PLUS.train), function(iPLUS){
    sapply(1:length(MINUS.train), function(iMINUS){
      return(train.res[[iPLUS]][[iMINUS]][i, 3])
    })
  })
  vals <- abs(as.numeric(vals))
  return(min(vals))
})
#Determine which bins pass minimum absolute separation
min.sep <- 0.1
minAbsSep.keep.bins <- which(minAbsSep>min.sep)

#Set p-value cutoff
p.cutoff <- 0.01 #Note: this is a cutoff of FDR q<0.01
#Get number of significant tests per bin
sigPcount <- sapply(1:nrow(cov), function(i){
  pvals <- sapply(1:length(PLUS.train), function(iPLUS){
    sapply(1:length(MINUS.train), function(iMINUS){
      return(train.res[[iPLUS]][[iMINUS]][i, 4])
    })
  })
  return(length(which(p.adjust(pvals, method="fdr")<p.cutoff)))
})
#Get maximum number of significant p-values
max.sig.tests <- length(PLUS.train)*length(MINUS.train)
#Determine which bins pass p-value requirements
sigPcount.keep.bins <- which(sigPcount==max.sig.tests)

#Determine final set of bins as strict intersection of all three criteria
final.keep.bins <- intersect(intersect(meanCovs.keep.bins, 
                                       minAbsSep.keep.bins), 
                             sigPcount.keep.bins)


#####Format WGD mask with final bins
#Get coordinates of final bins
WGD.mask <- cov[final.keep.bins, 1:3]
colnames(WGD.mask) <- c("chr", "start", "end")
#Get weights of final bins
PLUS.train.all <- unlist(PLUS.train)
MINUS.train.all <- unlist(MINUS.train)
WGD.mask$weight <- sapply(final.keep.bins, function(i){
  plus.idx <- which(colnames(cov) %in% PLUS.train.all)
  minus.idx <- which(colnames(cov) %in% MINUS.train.all)
  weight <- mean(as.numeric(cov[i, plus.idx]), na.rm=T) - mean(as.numeric(cov[i, minus.idx]), na.rm=T)
  return(weight)
})
#Randomly downsample bins so weights are equal
plus.bins <- which(WGD.mask$weight>0)
minus.bins <- which(WGD.mask$weight<0)
min.weight.sum <- floor(min(sum(abs(WGD.mask$weight[plus.bins])), 
                            sum(abs(WGD.mask$weight[minus.bins]))))
set.seed(123456789)
plus.bins.shuffled <- sample(plus.bins)
plus.bins.shuffled.keep <- plus.bins.shuffled[1:head(which(cumsum(abs(WGD.mask$weight[plus.bins.shuffled]))>min.weight.sum), 1)]
set.seed(123456789)
minus.bins.shuffled <- sample(minus.bins)
minus.bins.shuffled.keep <- minus.bins.shuffled[1:head(which(cumsum(abs(WGD.mask$weight[minus.bins.shuffled]))>min.weight.sum), 1)]
balanced.bins.sorted.keep <- sort(unique(c(plus.bins.shuffled.keep, 
                                           minus.bins.shuffled.keep)))
WGD.mask <- WGD.mask[balanced.bins.sorted.keep, ]

#####Compute difference between testing sets
if(!is.null(PLUS.test) & !is.null(MINUS.test)){
  #Subset coverage file to test sample sets
  test.cov <- cov[final.keep.bins[balanced.bins.sorted.keep], ]

  #Compute mean sep between PCR+ and PCR-
  test.weights <- sapply(1:nrow(test.cov), function(i){
    PLUS.vals <- as.numeric(test.cov[i, which(colnames(test.cov) %in% PLUS.test)])
    MINUS.vals <- as.numeric(test.cov[i, which(colnames(test.cov) %in% MINUS.test)])
    return(mean(PLUS.vals, na.rm=T)-mean(MINUS.vals, na.rm=T))
  })
}


########################
#####Write final outputs
########################
#Write final bins to file
colnames(WGD.mask)[1] <- c("#chr")
write.table(WGD.mask, paste(OUTDIR, "WGD_mask.bed", sep=""), 
            col.names=T, row.names=F, sep="\t", quote=F)

#Write WGD training report to file
logfile <- paste(OUTDIR, "WGD_training_report.txt", sep="")
write("WGD Dosage Bias Model: Training Report\n", 
      file=logfile)
write(paste("Completed on ", date(), "\n", sep=""), 
      file=logfile, append=T)
write(paste("Output directory: ", OUTDIR, "\n", sep=""), 
      file=logfile, append=T)
write(paste("Coverage matrix: ", path.to.matrix, "\n", sep=""), 
      file=logfile, append=T)
write(paste("Training Sets:\n", 
            paste(unlist(lapply(1:length(PLUS.train), function(i){
              return(paste("  -PCR plus training set ", i, 
                           ": n=", prettyNum(length(PLUS.train[[i]]), sep=","), 
                           " samples\n", "     ", 
                           if(i==1 & !is.null(PLUS.train.1.in)){
                             paste(PLUS.train.1.in, "\n", sep="")
                           }, 
                           if(i==2 & !is.null(PLUS.train.2.in)){
                             paste(PLUS.train.2.in, "\n", sep="")
                           }, 
                           if(i==3 & !is.null(PLUS.train.3.in)){
                             paste(PLUS.train.3.in, "\n", sep="")
                           }, 
                           if(i==4 & !is.null(PLUS.train.4.in)){
                             paste(PLUS.train.4.in, "\n", sep="")
                           }, 
                           if(i==5 & !is.null(PLUS.train.5.in)){
                             paste(PLUS.train.5.in, "\n", sep="")
                           }, 
                           sep=""))
            })), collapse=""), 
            paste(unlist(lapply(1:length(MINUS.train), function(i){
              return(paste("  -PCR minus training set ", i, 
                           ": n=", prettyNum(length(MINUS.train[[i]]), sep=","), 
                           " samples\n", "     ", 
                           if(i==1 & !is.null(MINUS.train.1.in)){
                             paste(MINUS.train.1.in, "\n", sep="")
                           }, 
                           if(i==2 & !is.null(MINUS.train.2.in)){
                             paste(MINUS.train.2.in, "\n", sep="")
                           }, 
                           if(i==3 & !is.null(MINUS.train.3.in)){
                             paste(MINUS.train.3.in, "\n", sep="")
                           }, 
                           if(i==4 & !is.null(MINUS.train.4.in)){
                             paste(MINUS.train.4.in, "\n", sep="")
                           }, 
                           if(i==5 & !is.null(MINUS.train.5.in)){
                             paste(MINUS.train.5.in, "\n", sep="")
                           }, 
                           sep=""))
            })), collapse=""), 
            sep=""), 
      file=logfile, append=T)
write(paste("Breakdown of bins:\n", 
            paste("  -Bins in initial matrix: n=", 
                  prettyNum(nrow(cov), sep=","), 
                  " (", "100%", ")\n", sep=""), 
            paste("  -Bins meeting selection criteria: n=", 
                  prettyNum(length(final.keep.bins), sep=","), 
                  " (", round(100*length(final.keep.bins)/nrow(cov), digits=1), "%)\n", sep=""), 
            paste("  -Final bins retained after plus/minus weight balancing: n=", 
                  prettyNum(nrow(WGD.mask), sep=","), 
                  " (", round(100*nrow(WGD.mask)/nrow(cov), digits=1), "%)\n", sep=""), 
            sep=""), 
      file=logfile, append=T)
if(!is.null(PLUS.test) & !is.null(MINUS.test)){
  write(paste("Testing Sets:\n", 
              paste("  -PCR plus testing set: n=", 
                    prettyNum(length(PLUS.test), sep=","), 
                    " samples\n", "     ", 
                    PLUS.test.in, "\n", sep=""), 
              paste("  -PCR minus testing set: n=", 
                    prettyNum(length(MINUS.test), sep=","), 
                    " samples\n", "     ", 
                    MINUS.test.in, "\n", sep=""), 
              sep=""), 
        file=logfile, append=T)
  write(paste("Bin performance on testing sets:\n", 
              paste("  -Pearson correlation of bin weights to observed coverage differences: R=", 
                    round(as.numeric(cor.test(WGD.mask$weight, test.weights)$estimate), digits=4), 
                    " (P=", format(cor.test(WGD.mask$weight, test.weights)$p.value, scientific=T), ")\n", sep=""), 
              paste("  -Spearman correlation of bin weights to observed coverage differences: Rho=", 
                    round(as.numeric(cor.test(WGD.mask$weight, test.weights, method="spearman")$estimate), digits=4), 
                    " (P=", format(cor.test(WGD.mask$weight, test.weights, method="spearman")$p.value, scientific=T), ")\n", sep=""), 
              sep=""), 
        file=logfile, append=T)
}

#Write one file for each training set comparison
lapply(1:length(PLUS.train), function(iPLUS){
  lapply(1:length(MINUS.train), function(iMINUS){
    #Clean table
    dat.out <- cbind(cov[, 1:3], train.res[[iPLUS]][[iMINUS]])
    colnames(dat.out) <- c("#chr", "start", "end", "PLUS_mean", "MINUS_mean", "weight", "tTest_pVal")
    write.table(dat.out, 
                paste(OUTDIR, "WGD_training.PLUS_", iPLUS, "_vs_MINUS_", iMINUS, ".bed", sep=""), 
                col.names=T, row.names=F, sep="\t", quote=F)
  })
})
#Write file with training weights vs testing weights
if(!is.null(PLUS.test) & !is.null(MINUS.test)){
  dat.out <- cbind(WGD.mask, test.weights)
  colnames(dat.out) <- c("#chr", "start", "end", "fitted_training_weight", "observed_testing_weight")
  write.table(dat.out, 
              paste(OUTDIR, "WGD_training_vs_testing_weights.bed", sep=""), 
              col.names=T, row.names=F, sep="\t", quote=F)
}


###################
#####Plotting block
###################
if(plot==T){
  #####Grid of training set correlations
  png(paste(OUTDIR, "WGD_model.training_set_correlations.png", sep=""), 
      width=400*length(MINUS.train), height=400*length(PLUS.train), 
      res=300)
  par(mfrow=c(length(PLUS.train), length(MINUS.train)))
  #Iterate over pairs of training sets
  sapply(1:length(PLUS.train), function(iPLUS){
    sapply(1:length(MINUS.train), function(iMINUS){
      #Set appropriate axes
      xaxis.bottom <- F
      xaxis.top <- F
      yaxis.left <- F
      yaxis.right <- F
      if(iPLUS==1){
        xaxis.top <- T
      }
      # if(iPLUS==length(PLUS.train)){
      #   xaxis.bottom <- T
      # }
      if(iMINUS==1){
        yaxis.left <- T
      }
      # if(iMINUS==length(MINUS.train)){
      #   yaxis.right <- T
      # }
      #Plot scatter
      plotScatterSingle(xvals.all=train.res[[iPLUS]][[iMINUS]]$meanA, 
                        yvals.all=train.res[[iPLUS]][[iMINUS]]$meanB, 
                        xvals.final=train.res[[iPLUS]][[iMINUS]]$meanA[final.keep.bins[balanced.bins.sorted.keep]], 
                        yvals.final=train.res[[iPLUS]][[iMINUS]]$meanB[final.keep.bins[balanced.bins.sorted.keep]], 
                        xaxis.bottom=xaxis.bottom, 
                        xaxis.top=xaxis.top, 
                        yaxis.left=yaxis.left, 
                        yaxis.right=yaxis.right)
    })
  })
  dev.off()
  #####Grid of training weights vs testing weights
  png(paste(OUTDIR, "WGD_model.training_weights_vs_testing_weights.png", sep=""), 
      width=900, height=900, res=300)
  plotScatterSingle(xvals.all=WGD.mask$weight, 
                    yvals.all=test.weights, 
                    xvals.final=WGD.mask$weight, 
                    yvals.final=test.weights, 
                    xaxis.bottom=T, yaxis.left=T, spline=F, 
                    axis.lims=c(-1, 1), highlight="darkorchid")
  dev.off()
  #####Create pdf book of swarmplots per bin
  pdf(paste(OUTDIR, "final_bin_evidence.pdf", sep=""), 
      height=6, width=10)
  sapply(1:nrow(WGD.mask), function(i){
    chr <- WGD.mask[i, 1]
    start <- WGD.mask[i, 2]
    end <- WGD.mask[i, 3]
    plotBinEvidence(chr=chr, start=start, end=end, 
                    title=paste("Final WGD Scoring Bin: ", 
                                chr, ":", prettyNum(start, sep=","), 
                                "-", prettyNum(end, sep=","), "bp", sep=""))
  })
  dev.off()
}



