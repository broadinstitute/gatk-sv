#!/usr/bin/env Rscript

# Copyright (c) 2018 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# Performs ROC analysis on de novo rate to determine optimal 
# GQ threshold for a single condition in the v2 ROC workflow


###Set master parameters
options(stringsAsFactors=F, scipen=1000)


###################
###HELPER FUNCTIONS
###################
#Read & clean trio variant data
readTrioDat <- function(dn_dat, metric){
  dat <- read.table(dn_dat, header=F)
  colnames(dat) <- c("fam", "VID", "pro_AC", "fa_AC", "mo_AC", 
                     paste("pro", metric, sep="_"),
                     paste("fa", metric, sep="_"),
                     paste("mo", metric, sep="_"))
  dat <- dat[which(!duplicated(paste(dat$fam, dat$VID, sep="_"))), ]
  dat[, (ncol(dat)-2):ncol(dat)] <- apply(dat[, (ncol(dat)-2):ncol(dat)], 2, as.integer)
  dat <- dat[which(dat$pro_AC==1), ]
  dat$label <- "inh"
  dat$label[which(dat$pro_AC - (dat$fa_AC + dat$mo_AC)>0)] <- "dn"
  dat.out <- data.frame(dat$fam, dat$label, dat[, paste("pro", metric, sep="_")],
                        dat[, paste("fa", metric, sep="_")],
                        dat[, paste("mo", metric, sep="_")])
  colnames(dat.out) <- c("fam", "label", paste("pro", metric, sep="_"),
                         paste("fa", metric, sep="_"), paste("mo", metric, sep="_"))
  return(dat.out)
}

#Get counts of all inherited & de novo variants given min metric
filterMetricStats <- function(dat, trio, metric, min.metric=0){
  #Get variant counts
  max.inh <- nrow(dat[which(dat$fam==trio & dat$label=="inh"),])
  ret.inh <- nrow(dat[which(dat$fam==trio & dat[, paste("pro", metric, sep="_")] > min.metric & dat$label=="inh"),])
  ret.dn <- nrow(dat[which(dat$fam==trio & dat$label=="dn" &
                             dat[, paste("pro", metric, sep="_")] > min.metric & 
                             dat[, paste("fa", metric, sep="_")] > min.metric & 
                             dat[, paste("mo", metric, sep="_")] > min.metric), ])
  #Get rates
  if(max.inh>0){
    ret.inh.rate <- ret.inh/max.inh
    ret.dn.rate <- ret.dn/(ret.inh+ret.dn)
  }else{
    ret.inh.rate <- NA
    ret.dn.rate <- NA
  }
  #Return line
  c(ret.inh, ret.inh.rate, ret.dn, ret.dn.rate)
}

#Get matrix of sens & spec stats across a range of metric values
getMetricStats <- function(dat, trio, metric, metric.range){
  iter.res <- as.data.frame(t(sapply(metric.range, function(min.metric){
    stats <- filterMetricStats(dat, trio, metric, min.metric)
    return(stats)
  })))
  colnames(iter.res) <- c("inh.ret","inh.ret.rate",
                          "dn.ret","dn.rate")
  iter.res <- cbind("min.value"=metric.range, iter.res)
  return(iter.res)
}

#Read data & compute median stats across all trios
processDeNovoDat <- function(dn_dat, famfile_path, metric, min.metric, max.metric, step){
  #Read & format input files
  dat <- readTrioDat(dn_dat, metric)
  fams <- unique(as.character(read.table(famfile_path, header=F)[, 1]))
  #Gather sens & spec stats for all families
  metric.range <- seq(min(c(min.metric, max.metric)),
                      max(c(min.metric, max.metric)),
                      step)
  res <- lapply(fams, function(fam){
    return(getMetricStats(dat=dat, trio=fam, metric=metric, metric.range=metric.range))
  })
  res <- do.call("rbind", res)
  res <- as.data.frame(t(sapply(metric.range, function(i){
    apply(res[which(res$min.value==i), ], 2, median, na.rm=T)
  })))
  return(res)
}

#ROC analysis to optimize minimum metric value
ROC <- function(min.value, sens, spec, min.spec=0.05){
  #Compute linear distance
  df <- data.frame("min.value"=min.value, "sens"=sens, "spec"=spec)
  df.elig <- df[which(df$spec<=min.spec),]
  if(nrow(df.elig)>0){
    opt.value <- df.elig$min.value[1]
  }else{
    #If max FDR is never met, take most stringent GQ filter considered
    # where sens > 0
    opt.value <- df[max(which(df$sens>0)),]$min.value
  }
  return(opt.value)
}

#Generate a ROC curve and mark optimal cutoff
ROCcurve <- function(metric, min.value, sens, spec, min.spec, opt,
                     xlims=NULL, ylims=NULL, title=NULL,
                     labels=T){
  #Reformat values
  sens <- 1-sens; spec <- 1-spec
  if(is.null(xlims)){
    xlims <- range(sens, na.rm=T)
  }
  if(is.null(ylims)){
    ylims <- range(spec, na.rm=T)
  }
  if(min(ylims, na.rm=T)>0.94){
    ylims[which(ylims==min(ylims))] <- 0.94
  }
  
  #Get colors
  col.pal <- colorRampPalette(c("#440154", "#365C8C", "#25A584", "#FDE725"))(length(min.value))
  
  #Prep plot area
  par(mar=c(2.5, 3.25, 1.25, 0.5))
  plot(x=sens, y=spec, type="n", xlim=xlims, ylim=ylims,
       xaxt="n", yaxt="n", xlab="", ylab="")
  abline(h=axTicks(2), v=axTicks(1), col="gray90", lwd=0.5)
  abline(h=1-min.spec, lty=2, col="gray40")
  if(labels){
    text(x=par("usr")[2], y=1-min.spec-0.04*(par("usr")[4]-par("usr")[3]),
         pos=2, labels=paste("~",round(100*min.spec, 1), "% FDR", sep=""),
         font=3, cex=0.7, col="gray40")
  }
  
  #Add axes
  axis(1, at=axTicks(1), labels=NA)
  axis(1, at=axTicks(1), tick=F, line=-0.6, cex.axis=0.7,
       labels=paste(round(100*(1-axTicks(1)), 1), "%", sep=""))
  axis(2, at=axTicks(2), labels=NA)
  axis(2, at=axTicks(2), tick=F, line=-0.4, cex.axis=0.7, las=2,
       labels=paste(round(100*(1-axTicks(2)),1), "%", sep=""))
  if(labels){
    mtext(1, line=1.35, text="Inherited Het. SV Retained", cex=0.8)
    mtext(2, line=2, text="Apparent De Novo Rate", cex=0.8)
    mtext(3, line=0.1, font=2, text=title)
  }
  
  #Plot ROC curve
  segments(x0=sens[1:(length(min.value)-1)], x1=sens[2:length(min.value)],
           y0=spec[1:(length(min.value)-1)], y1=spec[2:length(min.value)],
           col=col.pal, lwd=2)
  points(x=sens, y=spec, pch=19, col=col.pal, cex=0.7)
  
  #Annotate with optimal GQ threshold
  arrow.buffer.x <- 0.01*(par("usr")[2]-par("usr")[1])
  arrow.buffer.y <- 0.01*(par("usr")[4]-par("usr")[3])
  opt.coords <- c(sens[which(min.value==opt)], spec[which(min.value==opt)])
  if(opt.coords[2] < (par("usr")[3]+(50*arrow.buffer.y))){
    arrow.buffer.y <- -arrow.buffer.y
  }
  arrow.start <- c(opt.coords[1]+6*arrow.buffer.x,
                   opt.coords[2]-6*arrow.buffer.y)
  arrows(x0=arrow.start[1], x1=opt.coords[1]+2*arrow.buffer.x,
         y0=arrow.start[2], y1=opt.coords[2]-2*arrow.buffer.y,
         length=0.07, lwd=2, angle=30)
  if(labels){
    text(x=par("usr")[2], y=par("usr")[3]+(0.1*(par("usr")[4]-par("usr")[3])), 
         pos=2, cex=0.6, font=2,
         labels=paste("Optimal Threshold:\n", metric, " > ", opt, "\n",
                      round(100*(1-opt.coords[1]), 1),
                      "% of Het. SV Retained\n", round(100*(1-opt.coords[2]), 1),
                      "% De Novo", sep=""))
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
  make_option(c("--optimize-metric"), type="character", default="GQ",
              help="FORMAT metric to optimize [default %default]",
              metavar="character"),
  make_option(c("--min-metric-value"), type="numeric", default=0,
              help="minimum metric in window to test [default %default]",
              metavar="integer"),
  make_option(c("--max-metric-value"), type="numeric", default=999,
              help="maximum metric in window to test [default %default]",
              metavar="integer"),
  make_option(c("--step"), type="numeric", default=1,
              help="metric increment step size [default %default]",
              metavar="integer"),
  make_option(c("--noplot"), type="logical", default=FALSE,
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
metric <- opts$`optimize-metric`
min.metric <- opts$`min-metric-value`
max.metric <- opts$`max-metric-value`
step <- opts$step
plot <- !(opts$noplot)

###Create output directory, if necessary
if(!dir.exists(OUTDIR)){
  dir.create(OUTDIR)
}

###Compute ROC stats
res <- processDeNovoDat(dn_dat=dn_dat, famfile_path=famfile_path, metric=metric,
                        min.metric=min.metric, max.metric=max.metric, step=step)
opt.value <- ROC(min.value=res$min.value, sens=res$inh.ret.rate, 
                 spec=res$dn.rate, min.spec=FDR)
res.opt <- cbind(data.frame("condition"=prefix), res[which(res$min.value==opt.value),])

###Write tables to output files
write.table(res, paste(OUTDIR, "/", prefix, ".min", metric, "_ROC.data.txt", sep=""),
            col.names=T, row.names=F, quote=F, sep="\t")
write.table(res.opt, paste(OUTDIR, "/", prefix, ".min", metric, "_ROC.optimal.txt", sep=""),
            col.names=T, row.names=F, quote=F, sep="\t")

###Generate ROC plot, if optioned
if(plot){
  pdf(paste(OUTDIR, "/", prefix, ".min", metric, "_ROC.plot.pdf", sep=""),
      height=3, width=3)
  ROCcurve(metric=metric, min.value=res$min.value, sens=res$inh.ret.rate, 
           spec=res$dn.rate, min.spec=FDR, opt=opt.value, title=prefix, xlims=c(0,1))
  dev.off()
}

