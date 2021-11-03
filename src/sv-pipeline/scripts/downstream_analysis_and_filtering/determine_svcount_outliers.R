#!/usr/bin/env Rscript

# Helper script to return list of outlier samples from an svcounts.txt file,
# their reasons for being outliers, and one plot for each SVTYPE


###Set global options
options(stringsAsFactors=F,scipen=1000)


###################
###HELPER FUNCTIONS
###################
#Get list of outliers for a single svtype
getOutliers <- function(dat,svtype,n.iqr,minPerSVTYPE){
  #Get counts corresponding to svtype
  vals <- dat$count[which(dat$svtype==svtype)]
  
  #Determine cutoffs
  quartiles <- as.numeric(quantile(vals,probs=c(0.25,0.75),na.rm=T))
  spacer <- n.iqr*IQR(vals,na.rm=T)
  if(median(vals,na.rm=T) >= minPerSVTYPE){
    cutoffs <- c(quartiles[1]-spacer,quartiles[2]+spacer)
  }else{
    cutoffs <- c(NA,NA)
  }
  
  #Return list of sample IDs failing cutoffs
  if(any(!is.na(cutoffs))){
    fails <- dat$sample[which(dat$svtype==svtype
                              & (dat$count<cutoffs[1] | dat$count>cutoffs[2]))]
  }else{
    fails <- as.character(c())
  }
  return(fails)
}
#Plot distribution of counts per sample & mark outliers
plotOutliers <- function(dat,svtype,n.iqr,minPerSVTYPE){
  #Get counts corresponding to svtype
  vals <- dat$count[which(dat$svtype==svtype)]
  
  #Determine cutoffs
  quartiles <- as.numeric(quantile(vals,probs=c(0.25,0.75),na.rm=T))
  spacer <- n.iqr*IQR(vals,na.rm=T)
  if(median(vals,na.rm=T) >= minPerSVTYPE){
    cutoffs <- c(quartiles[1]-spacer,quartiles[2]+spacer)
  }else{
    cutoffs <- c(NA,NA)
  }
  
  #Prep plot area
  par(mar=c(1,4,2,4))
  yrange <- c(max(c(0,min(c(cutoffs[1]-IQR(vals,na.rm=T),min(vals,na.rm=T)),na.rm=T))),
              max(c(cutoffs[2]+IQR(vals,na.rm=T),max(vals,na.rm=T)),na.rm=T))
  plot(x=c(0,2),y=yrange,type="n",
       xlab="",ylab="",xaxt="n",yaxt="n")
  rect(xleft=par("usr")[1],xright=par("usr")[2],
       ybottom=quartiles[1],ytop=quartiles[2],
       border=NA,bty="n",col="gray80")
  sapply(1:10,function(i){
    abline(h=quartiles[1]-(i*IQR(vals,na.rm=T)),col="gray70",lty=2)
    abline(h=quartiles[2]+(i*IQR(vals,na.rm=T)),col="gray70",lty=2)
  })
  mtext(3,line=0.5,font=2,text=svtype)
  
  #Add sina plot
  dens <- density(vals,na.rm=T)
  dens$y <- 0.8*dens$y/max(dens$y,na.rm=T)
  getJitter <- function(x,dens){
    jit <- dens$y[head(which(dens$x>=x),1)]
    if(is.na(jit) | is.null(jit)){
      jit <- 0
    }
    return(jit)
  }
  jitterX <- function(y,dens){
    return(jitter(1,amount=getJitter(y,dens)))
  }
  plot.df <- data.frame("x"=sapply(vals,function(y){jitterX(y,dens)}),
                        "y"=vals)
  polygon(x=c(-dens$y,rev(dens$y))+1,y=c(dens$x,rev(dens$x)),
          border=adjustcolor("darkgreen",alpha=0.3),col=NA)
  points(plot.df$x,plot.df$y,pch=21,lwd=1.5,cex=0.1,col="darkgreen")
  points(plot.df$x[which(plot.df$y<cutoffs[1] | plot.df$y>cutoffs[2])],
         plot.df$y[which(plot.df$y<cutoffs[1] | plot.df$y>cutoffs[2])],
         pch=19,cex=0.5,lwd=1.5,col="hotpink")
  
  #Add axes & cutoffs
  axis(2,at=axTicks(2),labels=NA)
  axis(2,at=axTicks(2),tick=F,line=-0.4,cex.axis=0.8,
       labels=prettyNum(axTicks(2),big.mark=","),las=2)
  if(any(!is.na(cutoffs))==T){
    axis(4,at=cutoffs,tick=F,line=-0.8,cex.axis=0.8,
         las=2,labels=c("Lower\nCutoff","Upper\nCutoff"))
    abline(h=cutoffs,lwd=2)
  }
  axis(4,at=mean(quartiles),tick=F,line=-0.8,cex.axis=0.8,las=2,labels="IQR")
  box()
}


################
###RSCRIPT BLOCK
################
require(optparse,quietly=T)
###List of command-line options
option_list <- list(
  make_option(c("-p","--prefix"), type="character", default="outlier_check",
              help="prefix appended to all output files [default %default]"),
  make_option(c("-I","--IQR"), type="numeric", default=3,
              help="number of IQRs above Q3 or below Q1 to use for outlier threshold [default %default]"),
  make_option(c("--minPerSVTYPE"), type="numeric", default=100,
              help="minimum median number of SV per SVTYPE per sample to be considered in outlier analysis [default %default]"),
  make_option(c("--noplot"), type="logical", default=FALSE, action="store_true", 
              help="disable plotting of outliers [default %default]")
)


###Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog INFILE OUTDIR",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options
INFILE <- as.character(args$args[1])
OUTDIR <- paste(as.character(args$args[2]),"/",sep="")
prefix <- as.character(opts$prefix)
n.iqr <- as.numeric(opts$IQR)
minPerSVTYPE <- as.integer(opts$minPerSVTYPE)
plot <- !opts$noplot


###Checks for appropriate positional arguments
if(length(args$args) != 2){
  stop("Incorrect number of required positional arguments\n")
}


###Makes output directory, if needed
if(!dir.exists(OUTDIR)){
  dir.create(OUTDIR)
}


###Reads input data
dat <- read.table(INFILE,header=T)
colnames(dat) <- c("sample","svtype","count")
samples <- unique(as.character(dat$sample))
svtypes <- unique(as.character(dat$svtype))


###Get list of outlier samples & write to outfile
outliers.df <- do.call("rbind", lapply(svtypes,function(svtype){
  out.samples <- getOutliers(dat=dat,svtype=svtype,n.iqr=n.iqr,minPerSVTYPE=minPerSVTYPE)
  out.df <- data.frame("sample"=out.samples,
                       "reason"=rep(paste(svtype,"_count_outlier",sep=""),
                                    times=length(out.samples)))
  return(out.df)
}))
colnames(outliers.df) <- c("#sample","reason")
write.table(outliers.df,paste(OUTDIR,"/",prefix,".SV_count_outlier_samples.txt",sep=""),
            col.names=T,row.names=F,quote=F,sep="\t")


###Plots distributions, if optioned
if(plot==T){
  sapply(svtypes,function(svtype){
    png(paste(OUTDIR,"/",prefix,".",svtype,".counts_per_sample.png",sep=""),
        height=1600,width=1000,res=400)
    plotOutliers(dat=dat,svtype=svtype,n.iqr=n.iqr,minPerSVTYPE=minPerSVTYPE)
    dev.off()
  })
  png(paste(OUTDIR,"/",prefix,".all_SVTYPEs.counts_per_sample.png",sep=""),
      height=1000,width=700*length(svtypes),res=400)
  par(mfrow=c(1,length(svtypes)))
  sapply(svtypes,function(svtype){
    plotOutliers(dat=dat,svtype=svtype,n.iqr=n.iqr,minPerSVTYPE=minPerSVTYPE)
  })
  dev.off()
}

