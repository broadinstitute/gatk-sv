#!/usr/bin/env Rscript

# gnomAD project credits available at http://gnomad.broadinstitute.org/

# Code to perform module 00 (matrix creation) QC across all batches for gnomAD v2


####################################
#####Set parameters & load libraries
####################################
args = commandArgs(trailingOnly=TRUE)
batch=args[1]
options(scipen=1000,stringsAsFactors=F)
contigs <- c(1:22,"X","Y")
evs <- c("PE","SR","RD","BAF")
ev.cols <- c("PE"="darkorchid",
             "SR"="forestgreen",
             "RD"="firebrick",
             "BAF"="darkorange1")


##################
#####Clear PLOTDIR
##################



#####################
#####Helper functions
#####################
#Read & clean data
readDat <- function(batch){
  #Iterate over evidence types & read data
  dat <- lapply(as.list(evs),function(ev){
    #Read data
    dat <- read.table(paste(batch,".",ev,
                            ".QC_matrix.txt",sep=""),
                      header=T,comment.char="")
    colnames(dat)[1] <- "batch"
    return(dat)
  })
  dat <- c(dat,batch)
  names(dat) <- c(evs,"batch")
  return(dat)
}
#Single ev type plotting function
plotEv <- function(dat,ev){
  #Extract data
  x <- dat[[which(names(dat)==ev)]][,-c(1:2)]
  x <- t(apply(x,1,function(vals){
    vals[which(vals==0 | is.na(vals))] <- 1
    return(vals)
  }))
  x <- log10(x)
  
  #Get color
  col <- ev.cols[which(names(ev.cols)==ev)]
  
  #Prep plot area
  par(mar=c(3,4.5,3,0.1),bty="n")
  plot(x=c(1,24),y=c(0,max(x)),type="n",
       xaxt="n",yaxt="n",xlab="",ylab="")
  
  #Add axes & title
  abline(v=1:24,col="gray95",lwd=4)
  axis(1,at=1:24,labels=contigs,line=-0.8,tick=F)
  axis(2,at=0:10,labels=c(0,10^(1:10)),las=2)
  mtext(3,text=paste(ev," evidence in last 1Mb per contig\n",
                            " Batch: ",dat$batch,sep=""))
  
  #Iterate over samples and plot lines
  apply(x,1,function(vals){
    points(x=1:24,y=vals,col=adjustcolor(col,alpha=0.15),type="l")
    points(x=1:24,y=vals,col=col,cex=0.5,pch=19)
  })
  
  #Cleanup box
  rect(xleft=par("usr")[1],xright=par("usr")[2],
       ybottom=par("usr")[3],ytop=par("usr")[4],
       col=NA)
  box(bty="o")
}
#Master data plotting function
plotDat <- function(batch){
  #Read data
  dat <- readDat(batch)
  
  #Plot PE, SR, RD, and BAF
  png(paste(batch,
            ".00_matrix_FC_QC.png",sep=""),res=300,
      height=3000,width=2000)
  par(mfrow=c(4,1))
  sapply(evs,plotEv,dat=dat)
  dev.off()
}

##################
#####Plot all data
##################
#Iterate over PCRPLUS batches
plotDat(batch)


#########################
#####Troubleshooting code
#########################
# dat <- readDat(PCR="PCRMINUS",Q=1,i=8)
# dat$SR[which(apply(dat$SR[-c(1,2)],1,median)==0),2]


