#!/usr/bin/env Rscript

# Performs ROC curve analysis to determine optimal GQ threshold


###Set master parameters
options(stringsAsFactors=F,scipen=1000)
svtypes <- c("DEL","DUP","INS","INV","CPX","BND")


###########################
###GENERAL HELPER FUNCTIONS
###########################
#Read data & compute median across all trios
readDeNovoDat <- function(dn_dat){
  dat <- read.table(dn_dat,header=T)
  GQrange <- sort(unique(as.numeric(dat$minGQ)))
  # harmonic.mean <- function(vals){1/mean(1/vals,na.rm=T)}
  # res <- as.data.frame(t(sapply(GQrange,function(i){
  #   apply(dat[which(dat$minGQ==i),],2,harmonic.mean)
  # })))
  res <- as.data.frame(t(sapply(GQrange,function(i){
    apply(dat[which(dat$minGQ==i),],2,median,na.rm=T)
  })))
  return(res)
}
#ROC analysis to optimize minimum GQ
ROC <- function(minGQ,sens,spec,min.spec=0.05){
  #Compute linear distance
  df <- data.frame("minGQ"=minGQ,"sens"=sens,"spec"=spec)
  dist <- function(sens.i,spec.i){sqrt( ((1-sens.i)^2) + ((0-spec.i)^2) )}
  df$dist <- apply(df[,-1],1,function(vals){dist(vals[1],vals[2])})
  df.elig <- df[which(df$spec<=min.spec),]
  if(nrow(df.elig)>0){
    optGQ <- df.elig[order(df.elig$dist),]$minGQ[1]
  }else{
    #If max FDR is never met, take most stringent GQ filter considered
    optGQ <- df[nrow(df),]$minGQ
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
  par(mar=c(3,4,2,1))
  plot(x=sens,y=spec,type="n",xlim=xlims,ylim=ylims,
       xaxt="n",yaxt="n",xlab="",ylab="")
  abline(h=axTicks(2),v=axTicks(1),col="gray90",lwd=0.5)
  abline(h=1-min.spec,lty=2,col="gray40")
  if(labels==T){
    text(x=par("usr")[2],y=1-min.spec+0.03*(par("usr")[4]-par("usr")[3]),
         pos=2,labels=paste("~",round(100*min.spec,1),"% FDR",sep=""),font=3,cex=0.7,col="gray40")
  }
  
  #Add axes
  axis(1,at=axTicks(1),labels=NA)
  axis(1,at=axTicks(1),tick=F,line=-0.4,cex.axis=0.8,
       labels=paste(round(100*(1-axTicks(1)),1),"%",sep=""))
  axis(2,at=axTicks(2),labels=NA)
  axis(2,at=axTicks(2),tick=F,line=-0.4,cex.axis=0.8,las=2,
       labels=paste(round(100*(1-axTicks(2)),1),"%",sep=""))
  if(labels==T){
    mtext(1,line=1.75,text="Fraction of Het. SV Retained per Child",cex=0.8)
    mtext(2,line=2.25,text="Apparent De Novo Rate",cex=0.8)
    mtext(3,line=0.1,font=2,text=title)
  }
  
  #Plot ROC curve
  segments(x0=sens[1:(length(minGQ)-1)],x1=sens[2:length(minGQ)],
           y0=spec[1:(length(minGQ)-1)],y1=spec[2:length(minGQ)],
           col=col.pal,lwd=2)
  points(x=sens,y=spec,pch=19,col=col.pal,cex=0.8)
  
  #Annotate with optimal GQ threshold
  arrow.buffer.x <- 0.01*(par("usr")[2]-par("usr")[1])
  arrow.buffer.y <- 0.01*(par("usr")[4]-par("usr")[3])
  opt.coords <- c(sens[which(minGQ==opt)],spec[which(minGQ==opt)])
  if(opt.coords[2]<(par("usr")[3]+(50*arrow.buffer.y))){
    arrow.buffer.y <- -arrow.buffer.y
  }
  arrow.start <- c(opt.coords[1]+5*arrow.buffer.x,
                   opt.coords[2]-5*arrow.buffer.y)
  arrows(x0=arrow.start[1],x1=opt.coords[1]+2*arrow.buffer.x,
         y0=arrow.start[2],y1=opt.coords[2]-2*arrow.buffer.y,
         length=0.07,lwd=2,angle=20)
  if(labels==T){
    text(x=arrow.start[1],y=arrow.start[2]-5*arrow.buffer.y,pos=4,
         labels=paste("Optimal Threshold\nGQ > ",opt-1,"\n",round(100*(1-opt.coords[1]),1),
                      "% of Het. SV Retained\n",round(100*(1-opt.coords[2]),1),
                      "% De Novo",sep=""),cex=0.7,font=2)
  }
}



################
###RSCRIPT BLOCK
################
require(optparse)
###List of command-line options
option_list <- list(
  make_option(c("--prefix"), type="character", default="SV_denovo_analysis",
              help="batch prefix used for naming outfiles [default %default]",
              metavar="character"),
  make_option(c("--fdr"), type="numeric", default=0.05,
              help="maximum tolerated false discovery rate [default %default]",
              metavar="numeric"),
  make_option(c("-S", "--svtypes"), type="character", default=NULL,
              help="tab-delimited file specifying SV types and HEX colors [default %default]",
              metavar="character"),
  make_option(c("--noplot"), type="logical", default=F,
              help="do not generate ROC curve plots [default %default]",
              metavar="logical")
)

###Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog merged_denovo_data OUTDIR",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

###Checks for appropriate positional arguments
if(length(args$args) != 2){
  stop("Incorrect number of required positional arguments\n")
}

###Writes args & opts to vars
dn_dat <- args$args[1]
OUTDIR <- args$args[2]
prefix <- opts$prefix
svtypes.file <- opts$svtypes
plot <- !(opts$noplot)
FDR <- as.numeric(opts$fdr)

#Create output directory, if necessary
if(!dir.exists(OUTDIR)){
  dir.create(OUTDIR)
}

#Sets sv types & colors
if(!is.null(svtypes.file)){
  svtypes <- read.table(svtypes.file,sep="\t",header=F,comment.char="")
  svtypes <- as.data.frame(apply(svtypes,2,as.character))
  colnames(svtypes) <- c("svtype","color")
}else{
  svtypes.v <- svtypes
  svtypes.c <- brewer.pal(length(svtypes.v),"Dark2")
  svtypes <- data.frame("svtype"=svtypes.v,
                        "color"=svtypes.c)
}

###Read & average data
dat <- readDeNovoDat(dn_dat)

###Determine optimal GQ across all variants
opts <- sapply(1:(ncol(dat[,-1])/2),function(i){
  ROC(minGQ=dat[,1],sens=dat[,2*i],spec=dat[,(2*i)+1],min.spec=FDR)
})
names(opts) <- unique(unlist(lapply(strsplit(colnames(dat[,-1]),split=".",fixed=T),
                                    function(l){return(l[[1]])})))

###Generate ROC curves, if optioned
if(plot==T){
  sapply(1:(ncol(dat[,-1])/2),function(i){
    if(!is.na(opts[i])){
      pdf(paste(OUTDIR,"/",prefix,".",names(opts)[i],".minGQ_ROC.wLabels.pdf",sep=""),
          height=4,width=4)
      ROCcurve(minGQ=dat[,1],
               sens=dat[,2*i],
               spec=dat[,(2*i)+1],
               min.spec=FDR,
               opt=opts[i],
               title=names(opts)[i],
               xlims=c(0,1))
      dev.off()
      pdf(paste(OUTDIR,"/",prefix,".",names(opts)[i],".minGQ_ROC.noLabels.pdf",sep=""),
          height=4,width=4)
      ROCcurve(minGQ=dat[,1],
               sens=dat[,2*i],
               spec=dat[,(2*i)+1],
               min.spec=FDR,
               opt=opts[i],
               title=names(opts)[i],
               xlims=c(0,1),
               labels=F)
      dev.off()
    }
  })
}

###Output text file with minGQ cutoffs
df.out <- as.data.frame(t(sapply(1:length(opts),function(i){
  svtype <- names(opts)[i]
  minGQ.opt <- opts[i]
  sens.i <- dat[which(dat$minGQ==minGQ.opt),2*i]
  spec.i <- dat[which(dat$minGQ==minGQ.opt),(2*i)+1]
  return(c(svtype,minGQ.opt,sens.i,spec.i))
})))
colnames(df.out) <- c("SVTYPE","minGQ","het_retention_rate","deNovo_rate")
write.table(df.out,paste(OUTDIR,"/",prefix,".minGQ_ROC_results.txt",sep=""),
            col.names=T,row.names=F,sep="\t",quote=F)

