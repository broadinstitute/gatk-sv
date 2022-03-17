#!/usr/bin/env Rscript

# Performs ROC analysis on de novo rate to determine optimal 
# GQ threshold for a single condition in the v2 ROC workflow


###Set master parameters
options(stringsAsFactors=F,scipen=1000)


###################
###HELPER FUNCTIONS
###################
#Read & clean trio variant data
readTrioDat <- function(dn_dat){
  dat <- read.table(dn_dat,header=F)
  colnames(dat) <- c("fam","VID","pro_AC","fa_AC","mo_AC","pro_GQ","fa_GQ","mo_GQ")
  dat <- dat[which(!duplicated(paste(dat$fam,dat$VID,sep="_"))),]
  dat[,(ncol(dat)-2):ncol(dat)] <- apply(dat[,(ncol(dat)-2):ncol(dat)],2,as.integer)
  dat <- dat[which(dat$pro_AC==1),]
  dat$label <- "inh"
  dat$label[which(dat$pro_AC - (dat$fa_AC + dat$mo_AC)>0)] <- "dn"
  dat.out <- data.frame("fam"=dat$fam,"label"=dat$label,
                        "pro_GQ"=dat$pro_GQ,
                        "fa_GQ"=dat$fa_GQ,
                        "mo_GQ"=dat$mo_GQ)
  return(dat.out)
}
#Get counts of all inherited & de novo variants given min GQ
filterGQStats <- function(dat,minGQ=0){
  #Get variant counts
  max.inh <- nrow(dat[which(dat$label=="inh"),])
  ret.inh <- nrow(dat[which(dat$pro_GQ>minGQ & dat$label=="inh"),])
  ret.dn <- nrow(dat[which(dat$label=="dn" &
                             dat$pro_GQ>minGQ & dat$fa_GQ>minGQ & dat$mo_GQ>minGQ),])
  #Get rates
  if(max.inh>0){
    ret.inh.rate <- ret.inh/max.inh
    ret.dn.rate <- ret.dn/(ret.inh+ret.dn)
  }else{
    ret.inh.rate <- NA
    ret.dn.rate <- NA
  }
  #Return line
  c(ret.inh,ret.inh.rate,ret.dn,ret.dn.rate)
}
#Get matrix of sens & spec stats across a range of GQs
getGQStats <- function(dat,trio,rangeGQ){
  trio_dat <- dat[which(dat$fam==trio),]
  iter.res <- as.data.frame(t(sapply(rangeGQ,function(minGQ){
    stats <- filterGQStats(trio_dat,minGQ)
    return(stats)
  })))
  colnames(iter.res) <- c("inh.ret","inh.ret.rate",
                          "dn.ret","dn.rate")
  iter.res <- cbind("minGQ"=rangeGQ,iter.res)
  return(iter.res)
}
#Read data & compute median stats across all trios
processDeNovoDat <- function(dn_dat,famfile_path,minGQ,maxGQ,step){
  #Read & format input files
  dat <- readTrioDat(dn_dat)
  fams <- unique(as.character(read.table(famfile_path,header=F)[,1]))
  #Gather sens & spec stats for all families
  GQrange <- seq(min(c(minGQ,maxGQ)),
                 max(c(minGQ,maxGQ)),
                 step)
  res <- lapply(fams,function(fam){
    return(getGQStats(dat=dat,trio=fam,rangeGQ=GQrange))
  })
  res <- do.call("rbind", res)
  # harmonic.mean <- function(vals){1/mean(1/vals,na.rm=T)}
  # res <- as.data.frame(t(sapply(GQrange,function(i){
  #   apply(dat[which(dat$minGQ==i),],2,harmonic.mean)
  # })))
  res <- as.data.frame(t(sapply(GQrange,function(i){
    apply(res[which(res$minGQ==i),],2,median,na.rm=T)
  })))
  return(res)
}
#ROC analysis to optimize minimum GQ
ROC <- function(minGQ,sens,spec,min.spec=0.05){
  #Compute linear distance
  df <- data.frame("minGQ"=minGQ,"sens"=sens,"spec"=spec)
  df.elig <- df[which(df$spec<=min.spec),]
  if(nrow(df.elig)>0){
    optGQ <- df.elig$minGQ[1]
  }else{
    #If max FDR is never met, take most stringent GQ filter considered
    # where sens > 0
    optGQ <- df[max(which(df$sens>0)),]$minGQ
  }
  return(optGQ)
}
#Generate a ROC curve and mark optimal cutoff
ROCcurve <- function(minGQ,sens,spec,min.spec,opt,
                     xlims=NULL,ylims=NULL,title=NULL,
                     labels=T){
  #Reformat values
  sens <- 1-sens; spec <- 1-spec
  if(is.null(xlims)){
    xlims <- range(sens,na.rm=T)
  }
  if(is.null(ylims)){
    ylims <- range(spec,na.rm=T)
  }
  if(min(ylims,na.rm=T)>0.94){
    ylims[which(ylims==min(ylims))] <- 0.94
  }
  
  #Get colors
  col.pal <- colorRampPalette(c("#440154","#365C8C","#25A584","#FDE725"))(length(minGQ))
  
  #Prep plot area
  par(mar=c(2.5,3.25,1.25,0.5))
  plot(x=sens,y=spec,type="n",xlim=xlims,ylim=ylims,
       xaxt="n",yaxt="n",xlab="",ylab="")
  abline(h=axTicks(2),v=axTicks(1),col="gray90",lwd=0.5)
  abline(h=1-min.spec,lty=2,col="gray40")
  if(labels==T){
    text(x=par("usr")[2],y=1-min.spec-0.04*(par("usr")[4]-par("usr")[3]),
         pos=2,labels=paste("~",round(100*min.spec,1),"% FDR",sep=""),font=3,cex=0.7,col="gray40")
  }
  
  #Add axes
  axis(1,at=axTicks(1),labels=NA)
  axis(1,at=axTicks(1),tick=F,line=-0.6,cex.axis=0.7,
       labels=paste(round(100*(1-axTicks(1)),1),"%",sep=""))
  axis(2,at=axTicks(2),labels=NA)
  axis(2,at=axTicks(2),tick=F,line=-0.4,cex.axis=0.7,las=2,
       labels=paste(round(100*(1-axTicks(2)),1),"%",sep=""))
  if(labels==T){
    mtext(1,line=1.35,text="Inherited Het. SV Retained",cex=0.8)
    mtext(2,line=2,text="Apparent De Novo Rate",cex=0.8)
    mtext(3,line=0.1,font=2,text=title)
  }
  
  #Plot ROC curve
  segments(x0=sens[1:(length(minGQ)-1)],x1=sens[2:length(minGQ)],
           y0=spec[1:(length(minGQ)-1)],y1=spec[2:length(minGQ)],
           col=col.pal,lwd=2)
  points(x=sens,y=spec,pch=19,col=col.pal,cex=0.7)
  
  #Annotate with optimal GQ threshold
  arrow.buffer.x <- 0.01*(par("usr")[2]-par("usr")[1])
  arrow.buffer.y <- 0.01*(par("usr")[4]-par("usr")[3])
  opt.coords <- c(sens[which(minGQ==opt)],spec[which(minGQ==opt)])
  if(opt.coords[2]<(par("usr")[3]+(50*arrow.buffer.y))){
    arrow.buffer.y <- -arrow.buffer.y
  }
  arrow.start <- c(opt.coords[1]+6*arrow.buffer.x,
                   opt.coords[2]-6*arrow.buffer.y)
  arrows(x0=arrow.start[1],x1=opt.coords[1]+2*arrow.buffer.x,
         y0=arrow.start[2],y1=opt.coords[2]-2*arrow.buffer.y,
         length=0.07,lwd=2,angle=30)
  if(labels==T){
    text(x=par("usr")[2],y=par("usr")[3]+(0.1*(par("usr")[4]-par("usr")[3])),pos=2,
         labels=paste("Optimal Threshold:\nGQ > ",opt,"\n",round(100*(1-opt.coords[1]),1),
                      "% of Het. SV Retained\n",round(100*(1-opt.coords[2]),1),
                      "% De Novo",sep=""),cex=0.6,font=2)
  }
}

################
###RSCRIPT BLOCK
################
require(optparse,quietly=T)
###List of command-line options
option_list <- list(
  make_option(c("--prefix"), type="character", default="SV_denovo_analysis",
              help="batch prefix used for naming outfiles [default %default]",
              metavar="character"),
  make_option(c("--fdr"), type="numeric", default=0.05,
              help="maximum tolerated false discovery rate [default %default]",
              metavar="numeric"),
  make_option(c("--minGQ"), type="integer", default=0,
              help="minimum GQ in window to test [default %default]",
              metavar="integer"),
  make_option(c("--maxGQ"), type="integer", default=99,
              help="maximum GQ in window to test [default %default]",
              metavar="integer"),
  make_option(c("--step"), type="integer", default=1,
              help="GQ increment step size [default %default]",
              metavar="integer"),
  make_option(c("--noplot"), type="logical", default=F,
              help="do not generate ROC curve plots [default %default]",
              metavar="logical")
)

###Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog DATA FAMFILE OUTDIR",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

###Checks for appropriate positional arguments
if(length(args$args) != 3){
  stop("Incorrect number of required positional arguments\n")
}

###Writes args & opts to vars
dn_dat <- args$args[1]
famfile_path <- args$args[2]
OUTDIR <- args$args[3]
prefix <- opts$prefix
FDR <- as.numeric(opts$fdr)
minGQ <- opts$minGQ
maxGQ <- opts$maxGQ
step <- opts$step
plot <- !(opts$noplot)

###Create output directory, if necessary
if(!dir.exists(OUTDIR)){
  dir.create(OUTDIR)
}

###Compute ROC stats
res <- processDeNovoDat(dn_dat=dn_dat,famfile_path=famfile_path,minGQ=minGQ,maxGQ=maxGQ,step=step)
optGQ <- ROC(minGQ=res$minGQ,sens=res$inh.ret.rate,spec=res$dn.rate,min.spec=FDR)
res.opt <- cbind(data.frame("condition"=prefix),res[which(res$minGQ==optGQ),])

###Write tables to output files
write.table(res,paste(OUTDIR,"/",prefix,".minGQ_ROC.data.txt",sep=""),
            col.names=T,row.names=F,quote=F,sep="\t")
write.table(res.opt,paste(OUTDIR,"/",prefix,".minGQ_ROC.optimal.txt",sep=""),
            col.names=T,row.names=F,quote=F,sep="\t")

###Generate ROC plot, if optioned
if(plot==T){
  pdf(paste(OUTDIR,"/",prefix,".minGQ_ROC.plot.pdf",sep=""),
      height=3,width=3)
  ROCcurve(minGQ=res$minGQ,sens=res$inh.ret.rate,spec=res$dn.rate,min.spec=FDR,
           opt=optGQ,title=prefix,xlims=c(0,1))
  dev.off()
}



