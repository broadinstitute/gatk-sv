#!/usr/bin/env Rscript

# Script to score dosage bias per sample

#Inputs:
# 1) a binCov matrix subsetted to the testing bins of interest
# 2) a list of bins with signed weights per bin

####################################
#####Set parameters & load libraries
####################################
options(scipen=1000,stringsAsFactors=F)

############################################
#####Helper function to load coverage matrix
############################################
readMatrix <- function(INFILE){
  dat <- read.table(INFILE,check.names=FALSE,comment.char="",header=T)
  colnames(dat)[1] <- "Chr"
  dat[,-1] <- apply(dat[,-1],2,as.double)
  return(dat)
}

#############################################################
#####Helper function to normalize contigs for a single sample
#############################################################
normalizeContigsPerSample <- function(mat,exclude=c("X","Y"),ploidy=2){
  #Convert vals to numeric
  mat[,-c(1:3)] <- apply(mat[,-c(1:3),drop=F],2,as.numeric)

  #Iterate over all samples and scale each sample by median and center/scale to expected ploidy
  mat[,-c(1:3)] <- sapply(4:ncol(mat),function(i){
    #Compute median of values excluding current value & any other specified
    sampVals <- as.vector(mat[which(!(mat[,1] %in% exclude)),i])
    excl.median <- median(sampVals[which(sampVals>0)],na.rm=T)

    #Normalize values by excl.median
    newVals <- mat[,i]/excl.median

    #Scale to expected ploidy
    newVals <- ploidy*newVals

    #Return cleaned values
    return(newVals)
  })

  #Return normalized matrix
  return(mat)
}

###########################################################################
#####Helper function to remove rows where >X% of samples have <=Y% coverage
###########################################################################
filterZeroBins <- function(mat,exclude=c("X","Y"),minSamp=0.05,minCov=0.05){
  #Convert vals to numeric
  mat[,-c(1:3)] <- apply(mat[,-c(1:3),drop=F],2,as.numeric)

  #Find bins where >frac% of samples have coverage=0
  nZeros <- apply(mat[,-c(1:3),drop=F],1,function(vals){
    return(length(which(vals<=minCov)))
  })
  fracZeros <- nZeros/(ncol(mat)-3)
  keep.bins <- which(fracZeros<=minSamp | mat[,1] %in% exclude)

  #Return subsetted matrix
  return(mat[keep.bins,])
}

############################################################
#####Helper function to compute median per contig per sample
############################################################
medianPerContigPerSample <- function(dat){
  #Get list of unique contigs
  contigs <- unique(dat[,1])

  #Iterate over all contigs and compute median bin value per sample
  allMedians <- sapply(contigs,function(contig){
    #Calculate median value per sample
    chrVals <- apply(dat[,-c(1:3),drop=F],2,function(vals){
      #Subset vals to contig
      vals <- vals[which(dat[,1]==contig)]
      #Return median
      return(median(vals,na.rm=T))
    })
    #Return vector of medians
    return(as.vector(chrVals))
  })

  #Compose output data frame
  out.df <- data.frame("ID"=names(dat[,-c(1:3),drop=F]))
  out.df <- cbind(out.df,allMedians)
}

#########################################################################
#####Helper function to normalize contigs for an entire matrix of samples
#########################################################################
normalizeContigsPerMatrix <- function(dat,exclude=NA,scale.exclude=NA,
                                      genome.ploidy=2,contig.ploidy){
  #Iterate over samples & normalize
  suppressWarnings(if(is.na(exclude)){
    dat[,-1] <- t(apply(dat[,-1],1,normalizeContigsPerSample,genome.ploidy))
  }else{
    dat[,-1] <- t(apply(dat[,-1],1,normalizeContigsPerSample,exclude=exclude-1,genome.ploidy))
  })

  #Scale mad to mad of first 12 chromosomes
  mad.others <- mad(unlist(dat[,2:13]),na.rm=T)

  #Iterate over contigs (minus scale.exclude) and scale
  scaledVals <- sapply(setdiff(2:ncol(dat),scale.exclude),function(i){
    #Calculate & apply adjustments
    median.adjust <- median(dat[,i],na.rm=T)-contig.ploidy[i-1]
    newvals <- dat[,i]-median.adjust
    return(newvals)
  })
  suppressWarnings(if(is.na(scale.exclude)){
    dat[,-1] <- scaledVals
  }else{
    dat[,-c(1,scale.exclude)] <- scaledVals
  })

  #Round up values that were normalized below zero
  dat[,-1] <- apply(dat[,-1],2,function(vals){
    vals[which(vals<0 & !is.na(vals))] <- 0
    return(vals)
  })

  #Return transformed data
  return(dat)
}

#####################################################
#####Helper function to compute WGD scores per sample
#####################################################
scoreSamples <- function(cov,bins){
  #Iterate over samples
  scores <- apply(cov[,-c(1:3),drop=F],2,function(vals){
    #Determine plus and minus bins
    plus <- which(bins[,4]>0)
    minus <- which(bins[,4]<0)

    #Compute plus and minus scores
    plus.scores <- abs(bins[plus,4])*(vals[plus]-2)
    minus.scores <- abs(bins[minus,4])*(2-vals[minus])

    #Compute & return composite score
    score <- sum(c(plus.scores,minus.scores),na.rm=T)/length(vals)
    return(score)
  })

  #Sort, name, and return
  names(scores) <- colnames(cov[,-c(1:3),drop=F])
  scores <- sort(scores)
  scores <- data.frame("ID"=names(scores),
                       "score"=as.numeric(scores))
  return(scores)
}

##########################################
#####Helper function to plot dosage scores
##########################################
plotScores <- function(scores,label.outliers=T){
  #Determine range of scores
  max.score <- ceiling(100*max(abs(range(scores$score))))/100
  score.range <- c(-max.score,max.score)

  #Prep plot format
  par(mfrow=c(1,2),bty="n",mar=c(4,4,2,1))

  ###Panel a: histogram of scores
  #Calculate bins, plot dimensions, and colors
  hist.bins <- seq(-max.score,max.score,by=(2*max.score)/50)
  h <- hist(scores$score,plot=F,breaks=hist.bins)
  h.col <- sapply(h$mids,function(val){
    if(val>0){
      return("blue")
    }else{
      return("red")
    }
  })

  #Prep plot area
  plot(x=score.range,y=c(0,1.02*max(h$counts)),type="n",
       xlab="",ylab="",xaxt="n",yaxt="n",yaxs="i")

  #Plot rectangles
  rect(xleft=h$breaks[1:(length(h$breaks)-1)],
       xright=h$breaks[2:length(h$breaks)],
       ybottom=0,ytop=h$counts,col=h.col)
  abline(h=0,v=0,lwd=2)

  #Add axes
  axis(1,at=seq(-max.score,max.score,by=(2*max.score)/50),
       labels=NA,tck=-0.01,col="gray50")
  axis(1,at=seq(-max.score,max.score,by=(2*max.score)/10),
       labels=NA,tck=-0.025)
  axis(1,at=seq(-max.score,max.score,by=(2*max.score)/10),
       labels=round(seq(-max.score,max.score,by=(2*max.score)/10),digits=2),
       tick=F,line=-0.2,las=2,cex.axis=0.8)
  mtext(1,line=2.7,text=expression(paste("Dosage Score (",delta,")",sep="")))
  axis(2,at=seq(0,max(h$counts),by=ceiling(max(h$counts)/5)/5),
       labels=NA,tck=-0.01,col="gray50")
  axis(2,at=seq(0,max(h$counts),by=ceiling(max(h$counts)/5)),
       labels=NA,tck=-0.025)
  axis(2,at=seq(0,max(h$counts),by=ceiling(max(h$counts)/5)),
       labels=prettyNum(seq(0,max(h$counts),by=ceiling(max(h$counts)/5)),big.mark=","),
       tick=F,line=-0.2,las=2,cex.axis=0.8)
  mtext(2,text="Samples",line=2.7)
  mtext(3,font=2,line=0.5,text="Distribution of Dosage Scores")

  #Add legend
  legend("topright",pch=21,pt.bg=c("blue","red"),
         legend=c("Score > 0","Score < 0"),
         border=NA,bty="n",bg="white",cex=0.8,pt.cex=1.7)

  ###Panel b: ordered dotplot with labeled outliers
  #Get point colors
  cols <- sapply(scores$score,function(val){
    if(val>0){
      return("blue")
    }else{
      return("red")
    }
  })

  #Prep plot area
  plot(x=c(0,nrow(scores)),y=c(-max.score,max.score),type="n",
       xlab="",ylab="",xaxt="n",yaxt="n")
  abline(v=seq(0,nrow(scores),by=nrow(scores)/5),
         h=seq(-max.score,max.score,by=(2*max.score)/10),
         col="gray70",lwd=0.7)
  abline(h=0,lwd=2)
  abline(v=nrow(scores)/2,lty=2)

  #Add axes
  axis(1,at=seq(0,nrow(scores),by=nrow(scores)/20),
       tck=-0.01,labels=NA,col="gray50")
  axis(1,at=seq(0,nrow(scores),by=nrow(scores)/5),
       tck=-0.025,labels=NA)
  axis(1,at=seq(0,nrow(scores),by=nrow(scores)/5),
       tick=F,labels=paste(seq(0,100,20),"th",sep=""),
       line=-0.2,cex.axis=0.8,las=2)
  mtext(1,line=2.7,text=expression(paste(delta," Percentile",sep="")))
  axis(2,at=seq(-max.score,max.score,by=(2*max.score)/50),
       labels=NA,tck=-0.01,col="gray50")
  axis(2,at=seq(-max.score,max.score,by=(2*max.score)/10),
       labels=NA,tck=-0.025)
  axis(2,at=seq(-max.score,max.score,by=(2*max.score)/10),
       labels=round(seq(-max.score,max.score,by=(2*max.score)/10),digits=2),
       tick=F,line=-0.2,las=2,cex.axis=0.8)
  mtext(2,line=2.7,text=expression(paste("Dosage Score (",delta,")",sep="")))
  mtext(3,font=2,line=0.5,text="Samples Ranked by Dosage Score")

  #Plot points
  points(c(1:nrow(scores))-0.5,y=scores$score,
         pch=21,bg=cols,lwd=0.05,cex=0.7)

  #Label outlier points
  if(label.outliers==T){
    MAD <- mad(scores$score)
    out.range <- c(summary(scores$score)[2]-1.5*MAD,
                   summary(scores$score)[5]+1.5*MAD)
    outliers <- which(scores$score < min(out.range) |
                        scores$score > max(out.range))
    sapply(outliers,function(i){
      #Determine labeling side & add line
      if(scores$score[i]>0){
        pos <- 2
        segments(x0=i-0.5,x1=i-(0.02*nrow(scores)),
                 y0=scores$score[i],
                 y1=scores$score[i])
      }else{
        pos <- 4
        segments(x0=i-0.5,x1=i+(0.02*nrow(scores)),
                 y0=scores$score[i],
                 y1=scores$score[i])
      }
      #Add label
      text(x=i-1,y=scores$score[i],pos=pos,
           labels=scores$ID[i],cex=0.5,font=2)
    })
  }
}

##########################
#####Rscript functionality
##########################
require(optparse)
#List of command-line options
option_list <- list(
  make_option(c("-O", "--OUTDIR"),type="character",default=NULL,
              help="output directory [default: pwd]",
              metavar="character"),
  make_option(c("--noplot"),action="store_true",default=FALSE,
              help="disable all visualization [default: %default]"),
  make_option(c("--outliers"),action="store_true",default=FALSE,
              help="label outliers on dosage score visualization [default: %default]"),
  make_option(c("-z", "--gzip"),action="store_true",default=FALSE,
              help="gzip output files [default: %default]")
)

#Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog [options] coverage_matrix scoring_bins",
                                option_list=option_list),
                   positional_arguments=TRUE)

#Sanity-check arguments
if(length(args$args) != 2){
  stop("Must supply both an input coverage matrix and list of scoring bins\n")
}

#Clean arguments & options
cov.file <- args$args[1]
WGD.bins.file <- args$args[2]
OUTDIR <- args$options$OUTDIR
plot <- !(args$options$noplot)
outliers <- args$options$outliers
gzip <- args$options$gzip
if(is.null(OUTDIR)){
  OUTDIR <- "./"
}

##DEV TEST RUN (on local machine)
# WRKDIR <- "/Users/rlc/Desktop/Collins/Talkowski/NGS/SV_Projects/gnomAD/gnomad_analysis/WGD_version2_Dec2017/"
# cov.file <- paste(WRKDIR,"WGD_training.6F_adjCov.WGD_scoring_masked.100bp.matrix.bed.gz",sep="")
# WGD.bins.file <- paste(WRKDIR,"WGD_scoring_mask.6F_adjusted.100bp.h37.bed",sep="")
# plot <- T
# batch.ideal <- 100
# OUTDIR <- "~/scratch/WGDscore_testing/"
# gzip <- T

#Create OUTDIR if it doesn't already exist
if(!dir.exists(OUTDIR)){
  dir.create(OUTDIR)
}

#####PART 1: DATA PROCESSING#####
#Read scoring bins
WGD.bins <- read.table(WGD.bins.file,header=T,comment.char="",sep="\t")
if(ncol(WGD.bins) != 4){
  stop("WGD bins file must have four columns; see documentation.")
}
colnames(WGD.bins) <- c("chr","start","end","weight")

#Read, normalize, and clean coverage data
cov <- readMatrix(cov.file)
cov <- normalizeContigsPerSample(cov,exclude=c("X","Y"),ploidy=2)
cov <- filterZeroBins(cov)
colnames(cov)[1:3] <- c("chr","start","end")

#Subset coverage data to only include WGD bins
cov <- merge(WGD.bins[,1:3],cov,by=1:3)

#Subset WGD bins to only include relevant coverage data
WGD.bins <- merge(WGD.bins,cov[,1:3],by=1:3)


#####PART 2: DOSAGE SCORING PER SAMPLE#####
#Compute scores and return as data frame
scores <- scoreSamples(cov,WGD.bins)

#Write scores to file
colnames(scores)[1] <- "sample_id"
write.table(scores,
            paste(OUTDIR,"/WGD_scores.txt",sep=""),
            col.names=T,row.names=F,sep="\t",quote=F)
if(gzip==T){
  system(paste("gzip -f ",OUTDIR,"/WGD_scores.txt",sep=""),intern=F,wait=F)
}
colnames(scores)[1] <- "ID"

#####PART 3: PLOT SCORE DISTRIBUTIONS#####
#Only run if optioned
if(plot==T){
  pdf(paste(OUTDIR,"/WGD_score_distributions.pdf",sep=""),height=4,width=8)
  plotScores(scores,label.outliers=outliers)
  dev.off()
}
