#!/usr/bin/env Rscript

# Helper script to plot VCF summary stats output from clean_vcf2bed_output.R


###Set master parameters
options(stringsAsFactors=F,scipen=1000)
rare.max.freq <- 0.01
uncommon.max.freq <- 0.1
common.max.freq <- 0.5
major.max.freq <- 1
tiny.max.size <- 50
small.max.size <- 100
medium.max.size <- 500
medlarge.max.size <- 5000
large.max.size <- 50000
huge.max.size <- 300000000
sex.chroms <- c(1:22, paste("chr", 1:22, sep=""))


###################
###HELPER FUNCTIONS
###################
#General function to plot stacked bars from a matrix
plotStackedBars <- function(mat,colors,scaled=T,log.y=FALSE,title=NULL){
  #Scale columns, if options
  if(scaled==T){
    mat <- apply(mat,2,function(vals){
      vals/sum(vals,na.rm=T)
    })
  }
  
  #Prepare plot area
  ymax <- max(apply(mat,2,sum,na.rm=T),na.rm=T)
  ymin <- if(log.y) 0.5 else 0
  log.arg <- if(log.y) "y" else ""
  ymax.plot <- if(log.y) 10^ceiling(log10(ymax)) else 1.02*ymax
  par(mar=c(4,4,2,1),bty="n")
  plot(x=c(0,ncol(mat)),y=c(ymin,ymax.plot),type="n",
       xaxt="n",yaxt="n",xaxs="i",yaxs="i",xlab="",ylab="",log=log.arg)
  
  #Add axes & title
  if(log.y){
    ymax.plot <- 10^ceiling(log10(ymax))
    log.ticks <- 10^seq(floor(log10(ymin+0.1)), ceiling(log10(ymax)))
    axis(2,at=log.ticks,labels=NA)
    axis(2,at=log.ticks,tick=F,labels=prettyNum(log.ticks,big.mark=","),las=2,cex.axis=0.8,line=-0.4)
  }else{
    axis(2,at=axTicks(2),labels=NA)
    if(scaled==T){
      ylabs <- paste(round(100*axTicks(2),1),"%",sep="")
    }else{
      ylabs <- prettyNum(axTicks(2),big.mark=",")
    }
    axis(2,at=axTicks(2),tick=F,labels=ylabs,las=2,cex.axis=0.8,line=-0.4)
  }
  mtext(3,line=0.5,text=title,font=2)
  
  #Iterate and plot bars
  sapply(1:ncol(mat),function(i){
    #Get plotting values
    vals <- mat[,i]
    starts <- as.numeric(cumsum(c(0,vals[-length(vals)])))
    ends <- as.numeric(cumsum(vals))
    if(log.y){
      starts <- pmax(starts, 0.5)
      ends <- pmax(ends, 0.5)
    }
    
    #Plot bars
    rect(xleft=i-0.85,xright=i-0.15,ybottom=starts,ytop=ends,
         border=NA,bty="n",col=colors)
    rect(xleft=i-0.85,xright=i-0.15,
         ybottom=min(starts,na.rm=T),
         ytop=max(ends,na.rm=T),
         col=NA)
    
    #Add label
    axis(1,at=i-0.5,las=2,line=-0.8,labels=colnames(mat)[i],cex.axis=0.8,tick=F)
  })
}


###################
#####SV count plots
###################
#Plot single set of bars of total count of SV
plotSVCountBars <- function(dat,svtypes,title=NULL,ylab="Count"){
  #Compute table
  counts <- as.data.frame(t(sapply(svtypes$svtype,function(svtype){
    c(svtype,
      length(which(dat$svtype==svtype)),
      svtypes[which(svtypes$svtype==svtype),2])
  })))
  counts[,2] <- as.numeric(counts[,2])
  
  #Prep plotting area with log-scale y-axis
  par(bty="n",mar=c(6,4.5,2.5,0.5))
  pos.counts <- counts[counts[,2]>0, 2]
  min.y <- if(length(pos.counts)>0) max(0.5, min(pos.counts)*0.5) else 0.5
  max.y <- max(counts[,2], 1, na.rm=T) * 2
  plot(x=c(0,nrow(counts)),y=c(min.y,max.y),type="n",
       xaxt="n",yaxt="n",xlab="",ylab="",xaxs="i",yaxs="i",log="y")
  
  #Add y-axis and title
  log.ticks <- 10^(floor(log10(min.y)):ceiling(log10(max.y)))
  axis(2,at=log.ticks,labels=NA)
  axis(2,at=log.ticks,las=2,tick=F,line=-0.4,cex.axis=0.7,
       labels=prettyNum(log.ticks,big.mark=","))
  mtext(2,text=ylab,line=3)
  mtext(3,line=0.5,text=title,font=2)
  
  #Plot per-svtype information
  sapply(1:nrow(counts),function(i){
    cnt <- counts[i,2]
    if(cnt > 0){
      rect(xleft=i-0.85,xright=i-0.15,
           ybottom=min.y,ytop=cnt,
           lwd=0.7,col=counts[i,3])
      text(x=i-0.5,y=cnt*1.3,col=counts[i,3],
           labels=prettyNum(cnt,big.mark=","),cex=0.7)
    }
    axis(1,at=i-0.5,line=-0.8,tick=F,las=2,cex.axis=0.8,
         labels=counts[i,1],col.axis=counts[i,3])
  })
  
  #Add number of SV to plot
  axis(3,at=mean(par("usr")[1:2]),line=-1.5,tick=F,cex.axis=0.8,
       labels=paste("n=",prettyNum(nrow(dat),big.mark=","),sep=""))
}
#Plot dot for fraction of total SV per chromosome
plotDotsSVperChrom <- function(dat,svtypes,title=NULL,ylab="Fraction"){
  contigs <- sort(unique(dat$chr))
  #Compute table
  mat <- sapply(svtypes$svtype,function(svtype){
    counts <- sapply(contigs, function(contig){
      length(which(dat$chr==contig & dat$svtype==svtype))
    })
    if(sum(counts,na.rm=T)>0){
      return(counts/sum(counts))
    }else{
      return(rep(0, length(contigs)))
    }
  })

  #Ensure mat is always a matrix (single-contig data collapses sapply to a vector)
  if(!is.matrix(mat)){
    mat <- matrix(mat, nrow=length(contigs), ncol=length(svtypes$svtype),
                  dimnames=list(contigs, svtypes$svtype))
  }
  
  #Prep plotting area
  par(bty="n",mar=c(3,4.5,2.5,0.5))
  plot(x=c(0, length(contigs)),y=c(0,1.15*max(mat)),type="n",
       xaxt="n",yaxt="n",xlab="",ylab="",xaxs="i",yaxs="i")
  
  #Add axes and title
  sapply(1:length(contigs),function(i){
    axis(1,at=i-0.5,tick=F,labels=contigs[i],line=-0.8,cex.axis=0.7)
  })
  mtext(1,text="Chromosome",line=1.5)
  axis(2,at=axTicks(2),labels=NA)
  axis(2,at=axTicks(2),las=2,tick=F,line=-0.4,cex.axis=0.7,
       labels=paste(round(100*axTicks(2),digits=1),"%",sep=""))
  mtext(2,text=ylab,line=2.5)
  mtext(3,line=0.5,text=title,font=2)
  
  #Plot per-svtype information
  sapply(1:ncol(mat),function(i){
    points(x=(1:length(contigs))-seq(0.8,0.2,by=-0.6/(nrow(svtypes)-1))[i],
           y=mat[,i],type="l",lwd=0.5,col=adjustcolor(svtypes$color[i],alpha=0.5))
  })
  sapply(1:nrow(mat),function(i){
    #Points
    points(x=seq(i-0.8,i-0.2,by=0.6/(nrow(svtypes)-1)),
           y=mat[i,],pch=19,col=svtypes$color)
  })
  
  #Add number of SV to plot
  axis(3,at=mean(par("usr")[1:2]),line=-1.5,tick=F,cex.axis=0.8,
       labels=paste("n=",prettyNum(nrow(dat),big.mark=","),sep=""))
  
  #Add legend
  legend("topright",legend=svtypes$svtype,
         pch=19,col=svtypes$color,cex=0.7,border=NA,bty="n")
}

#Wrapper to plot all barplots of SV counts
wrapperPlotAllCountBars <- function(){
  #All SV
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/count.all.pdf",sep=""),
      height=4,width=2+(nrow(svtypes)/3))
  plotSVCountBars(dat=dat,svtypes=svtypes,
                  title="Variant Count")
  dev.off()
  #Singletons
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/count.singletons.pdf",sep=""),
      height=4,width=2+(nrow(svtypes)/3))
  plotSVCountBars(dat=dat[which(dat$AC==1),],svtypes=svtypes,
                  title="Variant Count (AC = 1)")
  dev.off()
  #Rare (>1 & <1%)
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/count.rare.pdf",sep=""),
      height=4,width=2+(nrow(svtypes)/3))
  plotSVCountBars(dat=dat[which(dat$AC>1 & dat$AF<rare.max.freq),],svtypes=svtypes,
                  title="Variant Count (AC > 1, AF < 1%)")
  dev.off()
  #Uncommon (≥1% & <10%)
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/count.uncommon.pdf",sep=""),
      height=4,width=2+(nrow(svtypes)/3))
  plotSVCountBars(dat=dat[which(dat$AF>=rare.max.freq & dat$AF<uncommon.max.freq),],svtypes=svtypes,
                  title="Variant Count (AF 1% - 10%)")
  dev.off()
  #Common (≥10% & <50%)
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/count.common.pdf",sep=""),
      height=4,width=2+(nrow(svtypes)/3))
  plotSVCountBars(dat=dat[which(dat$AF>=uncommon.max.freq & dat$AF<common.max.freq),],svtypes=svtypes,
                  title="Variant Count (AF 10% - 50%)")
  dev.off()
  #Major (≥50%%)
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/count.major.pdf",sep=""),
      height=4,width=2+(nrow(svtypes)/3))
  plotSVCountBars(dat=dat[which(dat$AF>=common.max.freq),],svtypes=svtypes,
                  title="Variant Count (AF > 50%)")
  dev.off()
  #Tiny
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/count.tiny.pdf",sep=""),
      height=4,width=2+(nrow(svtypes.merged)/3))
  plotSVCountBars(dat=dat.merged[which(dat.merged$length<tiny.max.size),],svtypes=svtypes.merged,
                  title="Variant Count (< 50bp)")
  dev.off()
  #Small
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/count.small.pdf",sep=""),
      height=4,width=2+(nrow(svtypes.merged)/3))
  plotSVCountBars(dat=dat.merged[which(dat.merged$length>=tiny.max.size & dat.merged$length<small.max.size),],svtypes=svtypes.merged,
                  title="Variant Count (50 - 100bp)")
  dev.off()
  #Medium
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/count.medium.pdf",sep=""),
      height=4,width=2+(nrow(svtypes.merged)/3))
  plotSVCountBars(dat=dat.merged[which(dat.merged$length>=small.max.size & dat.merged$length<medium.max.size),],svtypes=svtypes.merged,
                  title="Variant Count (100bp - 500bp)")
  dev.off()
  #Medium-Large
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/count.medlarge.pdf",sep=""),
      height=4,width=2+(nrow(svtypes.merged)/3))
  plotSVCountBars(dat=dat.merged[which(dat.merged$length>=medium.max.size & dat.merged$length<medlarge.max.size),],svtypes=svtypes.merged,
                  title="Variant Count (500bp - 5kb)")
  dev.off()
  #Large
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/count.large.pdf",sep=""),
      height=4,width=2+(nrow(svtypes.merged)/3))
  plotSVCountBars(dat=dat.merged[which(dat.merged$length>=medlarge.max.size & dat.merged$length<large.max.size),],svtypes=svtypes.merged,
                  title="Variant Count (5 - 50kb)")
  dev.off()
  #Huge
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/count.huge.pdf",sep=""),
      height=4,width=2+(nrow(svtypes.merged)/3))
  plotSVCountBars(dat=dat.merged[which(dat.merged$length>=large.max.size & dat.merged$length<huge.max.size),],svtypes=svtypes.merged,
                  title="Variant Count (> 50kb)")
  dev.off()
  #Dotplot per chromosome
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/count_per_chromosome.all.pdf",sep=""),
      height=4,width=6)
  plotDotsSVperChrom(dat=dat,svtypes=svtypes,
                     title="Variants per Chromosome")
  dev.off()
  #Size x AF Grid
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/counts.size_freq_grid.pdf",sep=""),
      height=6*1.75,width=7*1.75)
  #Iterator data frames
  AF.df <- data.frame("label"=c("AC=1","AF<1%","1-10%",
                                "10-50%",">50%"),
                      "max"=c(1.1/(2*nsamp),rare.max.freq,uncommon.max.freq,
                              common.max.freq,major.max.freq))
  size.df <- data.frame("label"=c("<50bp","50-\n100bp","100bp-\n500bp",
                                  "500bp-\n5kb","5-50kb",">50kb"),
                        "max"=c(tiny.max.size,small.max.size,medium.max.size,
                                medlarge.max.size,large.max.size,huge.max.size))
  #Prepare layout
  par(mfrow=c(nrow(AF.df)+1,nrow(size.df)+1))
  #Iterate over iterator data frames & plot bars
  sapply(0:nrow(AF.df),function(r){
    sapply(0:nrow(size.df),function(c){
      #Get cutoffs
      min.AF <- max(c(0,AF.df[max(c(0,r-1)),2]))
      max.AF <- min(c(AF.df[r,2],1))
      min.size <- max(0,size.df[max(c(0,c-1)),2])
      max.size <- min(c(size.df[c,2],300000000))
      #Get data subset
      plotset <- dat.merged[which(dat.merged$AF>min.AF & dat.merged$AF<=max.AF & 
                                    dat.merged$length>=min.size & dat.merged$length<=max.size),]
      #Set titles
      if(r==0){
        title <- gsub("\n","",c("ALL",size.df[,1])[c+1],fixed=T)
      }else{
        title <- NULL
      }
      if(c==0){
        ylab <- gsub("\n","",c("ALL",AF.df[,1])[r+1],fixed=T)
      }else{
        ylab <- NULL
      }
      #Plot
      if(nrow(plotset)==0){
        par(mar=c(0,0,0,0),bty="n")
        plot(x=c(0,1),y=c(0,1),type="n",xaxt="n",xlab="",yaxt="n",ylab="")
        text(x=0.5,y=0.5,labels="n=0")
      }else{
        plotSVCountBars(dat=plotset,svtypes=svtypes.merged,
                        title=title,ylab=NULL)
        mtext(2,line=3,font=2,text=ylab)
      }
    })
  })
  dev.off()
  #Stacked bars of SV count by frequency
  AF.mat <- as.data.frame(sapply(0:nrow(AF.df),function(i){
    #Get cutoffs
    min.AF <- max(c(0,AF.df[max(c(0,i-1)),2]))
    max.AF <- min(c(AF.df[i,2],1))
    #Get data subset
    plotset <- dat[which(dat$AF>min.AF & dat$AF<=max.AF),]
    #Tabulate counts per cutoff
    sapply(svtypes$svtype,function(svtype){
      length(which(plotset$svtype==svtype))
    })
  }))
  colnames(AF.mat) <- c("ALL",AF.df[,1])
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/counts.by_frequency_scaled.pdf",sep=""),
      height=4,width=4)
  plotStackedBars(mat=AF.mat,colors=svtypes$color,scale=T,
                  title="Count by AF (Scaled)")
  abline(v=1,lty=2,col="gray50")
  axis(1,at=1,labels=NA,col="gray75",tck=-0.22)
  dev.off()  
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/counts.by_frequency_raw.pdf",sep=""),
      height=4,width=4)
  plotStackedBars(mat=AF.mat,colors=svtypes$color,scale=F,
                  title="Count by AF")
  abline(v=1,lty=2,col="gray50")
  axis(1,at=1,labels=NA,col="gray75",tck=-0.22)
  dev.off()
  #Stacked bars of SV count by size
  size.mat <- as.data.frame(sapply(0:nrow(size.df),function(i){
    #Get cutoffs
    min.size <- max(c(0,size.df[max(c(0,i-1)),2]))
    max.size <- min(c(size.df[i,2],300000000))
    #Get data subset
    plotset <- dat[which(dat$length>=min.size & dat$length<=max.size),]
    #Tabulate counts per cutoff
    sapply(svtypes$svtype,function(svtype){
      length(which(plotset$svtype==svtype))
    })
  }))
  colnames(size.mat) <- c("ALL",size.df[,1])
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/counts.by_size_scaled.pdf",sep=""),
      height=4,width=4)
  plotStackedBars(mat=size.mat,colors=svtypes$color,scale=T,
                  title="Count by Size (Scaled)")
  abline(v=1,lty=2,col="gray50")
  axis(1,at=1,labels=NA,col="gray75",tck=-0.22)
  dev.off()  
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/counts.by_size_raw.pdf",sep=""),
      height=4,width=4)
  plotStackedBars(mat=size.mat,colors=svtypes$color,scale=F,
                  title="Count by Size")
  abline(v=1,lty=2,col="gray50")
  axis(1,at=1,labels=NA,col="gray75",tck=-0.22)
  dev.off()
  #Build REGION matrix if column exists
  if("REGION" %in% colnames(dat)){
    regions <- sort(unique(dat$REGION[!is.na(dat$REGION)]))
    region.mat <- as.data.frame(sapply(c("ALL",regions), function(reg){
      plotset <- if(reg=="ALL") dat else dat[which(dat$REGION==reg),]
      sapply(svtypes$svtype, function(svtype) length(which(plotset$svtype==svtype)))
    }))
    colnames(region.mat) <- c("ALL",regions)
  }
  pdf(paste(OUTDIR,"/main_plots/counts_distributions.pdf",sep=""),
      height=7,width=11)
  #Merged
  layout(matrix(c(1,2,3,4,1,5,6,7),byrow=T,nrow=2),
         widths=c(3,2,2,2))
  plotSVCountBars(dat=dat,svtypes=svtypes,
                  title="Variant Count")
  plotStackedBars(mat=AF.mat,colors=svtypes$color,scale=F,log.y=T,
                  title="Count by AF")
  abline(v=1,lty=2,col="gray50")
  axis(1,at=1,labels=NA,col="gray75",tck=-0.22)
  plotStackedBars(mat=size.mat,colors=svtypes$color,scale=F,log.y=T,
                  title="Count by Size")
  abline(v=1,lty=2,col="gray50")
  axis(1,at=1,labels=NA,col="gray75",tck=-0.22)
  if("REGION" %in% colnames(dat)){
    plotStackedBars(mat=region.mat,colors=svtypes$color,scale=F,log.y=T,
                    title="Count by Region")
    abline(v=1,lty=2,col="gray50")
    axis(1,at=1,labels=NA,col="gray75",tck=-0.22)
  }else{ plot.new() }
  plotStackedBars(mat=AF.mat,colors=svtypes$color,scale=T,
                  title="Count by AF (Scaled)")
  abline(v=1,lty=2,col="gray50")
  axis(1,at=1,labels=NA,col="gray75",tck=-0.22)
  plotStackedBars(mat=size.mat,colors=svtypes$color,scale=T,
                  title="Count by Size (Scaled)")
  abline(v=1,lty=2,col="gray50")
  axis(1,at=1,labels=NA,col="gray75",tck=-0.22)
  if("REGION" %in% colnames(dat)){
    plotStackedBars(mat=region.mat,colors=svtypes$color,scale=T,
                    title="Count by Region (Scaled)")
    abline(v=1,lty=2,col="gray50")
    axis(1,at=1,labels=NA,col="gray75",tck=-0.22)
  }else{ plot.new() }
  dev.off()
}


###############
#####Size plots
###############
#Plot single size distribution
plotSizeDistrib <- function(dat, svtypes, n.breaks=150, k=10,
                            min.size=1, max.size=1000000,
                            autosomal=F, biallelic=F,
                            title=NULL, legend=F, lwd.cex=1, text.cex=1){
  #Filter/process sizes & compute range + breaks
  filter.legend <- NULL
  if(autosomal==T){
    dat <- dat[which(dat$chr %in% sex.chroms),]
    filter.legend <- c(filter.legend,"Autosomal variants only")
  }
  if(biallelic==T){
    dat <- dat[which(dat$other_gts==0 & dat$missing_gts<dat$genotyped_samples),]
    filter.legend <- c(filter.legend,"Biallelic variants only")
  }
  sizes <- log10(dat$length)
  if(length(sizes)>0){
    xlims <- range(sizes[which(!is.infinite(sizes))],na.rm=T)
    xlims[1] <- log10(min.size)
    xlims[2] <- min(c(log10(max.size),xlims[2]))
    if(!all(is.finite(xlims)) || xlims[1] >= xlims[2]){
      par(bty="n",mar=c(3.5,3.5,3,0.5))
      plot(x=c(0,1),y=c(0,1),type="n",xaxt="n",yaxt="n",xlab="",ylab="",yaxs="i")
      text(x=0.5,y=0.5,labels="No Data")
      mtext(3,line=1.5,text=title,font=2,cex=text.cex)
      return(invisible(NULL))
    }
    breaks <- seq(xlims[1],xlims[2],by=(xlims[2]-xlims[1])/n.breaks)
    mids <- (breaks[1:(length(breaks)-1)]+breaks[2:length(breaks)])/2
    
    #Gather size densities per class
    dens <- lapply(svtypes$svtype,function(svtype){
      vals <- sizes[which(dat$svtype==svtype)]
      h <- hist(vals[which(!is.infinite(vals) & vals>=xlims[1] & vals<=xlims[2])],plot=F,breaks=breaks)
      h$counts[1] <- h$counts[1]+length(which(!is.infinite(vals) & vals<xlims[1]))
      h$counts[length(h$counts)] <- h$counts[length(h$counts)]+length(which(!is.infinite(vals) & vals>xlims[2]))
      return(h$counts/length(vals))
    })
    names(dens) <- svtypes$svtype
    all.vals <- sizes[which(!is.infinite(sizes) & sizes>=xlims[1] & sizes<=xlims[2])]
    all.h <- hist(all.vals,plot=F,breaks=breaks)
    dens$ALL <- as.numeric(all.h$counts/length(all.vals))
    
    #Prepare plot area
    ylims <- c(0,quantile(unlist(dens),probs=0.99,na.rm=T))
    dens <- lapply(dens,function(vals){
      vals[which(vals>max(ylims))] <- max(ylims)
      return(vals)
    })
    par(bty="n",mar=c(3.5,3.5,3,0.5))
    plot(x=xlims,y=ylims,type="n",
         xaxt="n",yaxt="n",xlab="",ylab="",yaxs="i")
    
    #Add vertical gridlines
    logscale.all <- log10(as.numeric(sapply(0:8,function(i){(1:9)*10^i})))
    logscale.minor <- log10(as.numeric(sapply(0:8,function(i){c(5,10)*10^i})))
    logscale.minor.labs <- as.character(sapply(c("bp","kb","Mb"),function(suf){paste(c(1,5,10,50,100,500),suf,sep="")}))
    logscale.minor.labs <- c(logscale.minor.labs[-1],"1Gb")
    logscale.major <- log10(as.numeric(10^(0:8)))
    abline(v=logscale.all,col="gray97")
    abline(v=logscale.minor,col="gray92")
    abline(v=logscale.major,col="gray85")
    
    #Add axes, title, and Alu/SVA/L1 ticks
    axis(1,at=logscale.all,tck=-0.015,col="gray50",labels=NA)
    axis(1,at=logscale.minor,tck=-0.0225,col="gray20",labels=NA)
    axis(1,at=logscale.major,tck=-0.03,labels=NA)
    axis(1,at=logscale.minor,tick=F,cex.axis=0.8,line=-0.4,las=2,
         labels=logscale.minor.labs)
    mtext(1,text="Size",line=2.25,cex=text.cex)
    axis(2,at=axTicks(2),tck=-0.025,labels=NA)
    axis(2,at=axTicks(2),tick=F,line=-0.4,cex.axis=0.8,las=2,
         labels=paste(round(100*axTicks(2),1),"%",sep=""))
    mtext(2,text="Fraction of Variants",line=2,cex=text.cex)
    sapply(1:2,function(i){
      axis(3,at=log10(c(300,6000))[i],labels=NA,tck=-0.01)
      axis(3,at=log10(c(300,6000))[i],tick=F,line=-0.9,cex.axis=0.8,
           labels=c("Alu","L1")[i],font=3)
    })
    axis(3,at=log10(c(1000,2000)),labels=NA,tck=-0.01)
    axis(3,at=mean(log10(c(1000,2000))),tick=F,line=-0.9,cex.axis=0.8,labels="SVA",font=3)
    mtext(3,line=1.5,text=title,font=2,cex=text.cex)
    
    #Add points per SV type
    sapply(1:length(dens),function(i){
      svtype <- names(dens)[i]
      if(svtype=="ALL"){
        color <- "gray15"
        lwd <- 3
      }else{
        color <- svtypes[which(svtypes$svtype==svtype),2]
        lwd <- 1.5
      }
      #Points per individual bin
      points(x=mids,y=dens[[i]],pch=19,cex=0.25,col=color)
    })  
    
    #Add rolling mean lines per SV type
    sapply(1:length(dens),function(i){
      svtype <- names(dens)[i]
      if(svtype=="ALL"){
        color <- "gray15"
        lwd <- 3
      }else{
        color <- svtypes[which(svtypes$svtype==svtype),2]
        lwd <- 1.5
      }
      #Rolling mean for line
      points(x=mids,y=rollapply(dens[[i]],width=k,mean,partial=T),type="l",lwd=lwd.cex*lwd,col=color)
    })  
    
    #Add sv type legend
    if(legend==T){
      idx.for.legend <- which(unlist(lapply(dens,function(vals){any(!is.na(vals) & !is.infinite(vals) & vals>0)})))
      counts.for.legend <- sapply(names(idx.for.legend), function(svtype){
        if(svtype == "ALL"){
          nrow(dat)
        }else{
          length(which(dat$svtype==svtype))
        }
      })
      legend("right",bg=NA,bty="n",pch=NA,cex=text.cex*0.7,lwd=3,
             legend=paste(rbind(svtypes, c("ALL","gray15"))$svtype[idx.for.legend],
                          " (N=", prettyNum(counts.for.legend, big.mark=","),
                          ")", sep=""),
             col=rbind(svtypes,c("ALL","gray15"))$color[idx.for.legend])
    }
  }else{
    par(bty="n",mar=c(3.5,3.5,3,0.5))
    plot(x=c(0,1),y=c(0,1),type="n",
         xaxt="n",yaxt="n",xlab="",ylab="",yaxs="i")
    text(x=0.5,y=0.5,labels="No Data")
    mtext(3,line=1.5,text=title,font=2,cex=text.cex)
  }
  
  #Add number of SV to plot
  axis(3,at=par("usr")[2],line=-0.9,hadj=1,tick=F,
       labels=paste("n=",prettyNum(length(sizes),big.mark=","),sep=""))
  
  #Add filter labels
  if(!is.null(filter.legend)){
    legend("topright",bg=NA,bty="n",pch=NA,legend=filter.legend,cex=text.cex)
  }
}

#Plot comparative size distributions for a series of AC & AF restrictions
plotSizeDistribSeries <- function(dat, svtypes, max.AFs, legend.labs,
                                  n.breaks=100, min.size=1, max.size=1000000,
                                  autosomal=F, biallelic=T, title=NULL, lwd.cex=1){
  #Process sizes & compute range + breaks
  filter.legend <- NULL
  if(autosomal==T){
    dat <- dat[which(dat$chr %in% sex.chroms),]
    filter.legend <- c(filter.legend,"Autosomal variants only")
  }
  if(biallelic==T){
    dat <- dat[which(dat$other_gts==0 & dat$missing_gts<dat$genotyped_samples),]
    filter.legend <- c(filter.legend,"Biallelic variants only")
  }
  sizes <- log10(dat$length)
  if(length(sizes) > 0){
    sizes <- lapply(1:length(max.AFs),function(i){
      if(i==1){
        return(sizes[which(dat$AF<=max.AFs[i])])
      }else{
        return(sizes[which(dat$AF>max.AFs[i-1] & dat$AF<=max.AFs[i])])
      }
    })
    finite_sizes <- unlist(sizes)[!is.infinite(unlist(sizes))]
    xlims <- range(finite_sizes, na.rm=T)
    xlims[1] <- max(c(log10(min.size),xlims[1]))
    xlims[2] <- min(c(log10(max.size),xlims[2]))
    if(!all(is.finite(xlims)) || xlims[1] >= xlims[2]){
      par(bty="n",mar=c(3.5,3.5,3,0.5))
      plot(x=c(0,1),y=c(0,1),type="n",xaxt="n",yaxt="n",xlab="",ylab="",yaxs="i")
      text(x=0.5,y=0.5,labels="No Data")
      mtext(3,line=1.5,text=title,font=2,cex=lwd.cex)
      return(invisible(NULL))
    }
    breaks <- seq(xlims[1],xlims[2],by=(xlims[2]-xlims[1])/n.breaks)
    mids <- (breaks[1:(length(breaks)-1)]+breaks[2:length(breaks)])/2
    
    #Gather size densities per AF tranche
    dens <- lapply(sizes,function(vals){
      h <- hist(vals[which(!is.infinite(vals) & vals>=xlims[1] & vals<=xlims[2])],plot=F,breaks=breaks)
      h$counts[1] <- h$counts[1]+length(which(!is.infinite(vals) & vals<xlims[1]))
      h$counts[length(h$counts)] <- h$counts[length(h$counts)]+length(which(!is.infinite(vals) & vals>xlims[2]))
      return(h$counts/length(vals))
    })
    
    #Prepare plot area
    ylims <- c(0,quantile(unlist(dens),probs=0.99,na.rm=T))
    dens <- lapply(dens,function(vals){
      vals[which(vals>max(ylims))] <- max(ylims)
      return(vals)
    })
    par(bty="n",mar=c(3.5,3.5,3,0.5))
    plot(x=xlims,y=ylims,type="n",
         xaxt="n",yaxt="n",xlab="",ylab="",yaxs="i")
    
    #Add vertical gridlines
    logscale.all <- log10(as.numeric(sapply(0:8,function(i){(1:9)*10^i})))
    logscale.minor <- log10(as.numeric(sapply(0:8,function(i){c(5,10)*10^i})))
    logscale.minor.labs <- as.character(sapply(c("bp","kb","Mb"),function(suf){paste(c(1,5,10,50,100,500),suf,sep="")}))
    logscale.minor.labs <- c(logscale.minor.labs[-1],"1Gb")
    logscale.major <- log10(as.numeric(10^(0:8)))
    abline(v=logscale.all,col="gray97")
    abline(v=logscale.minor,col="gray92")
    abline(v=logscale.major,col="gray85")
    
    #Add axes, title, and Alu/SVA/L1 ticks
    axis(1,at=logscale.all,tck=-0.015,col="gray50",labels=NA)
    axis(1,at=logscale.minor,tck=-0.0225,col="gray20",labels=NA)
    axis(1,at=logscale.major,tck=-0.03,labels=NA)
    axis(1,at=logscale.minor,tick=F,cex.axis=0.8,line=-0.4,las=2,
         labels=logscale.minor.labs)
    mtext(1,text="Size",line=2.25)
    axis(2,at=axTicks(2),tck=-0.025,labels=NA)
    axis(2,at=axTicks(2),tick=F,line=-0.4,cex.axis=0.8,las=2,
         labels=paste(round(100*axTicks(2),1),"%",sep=""))
    mtext(2,text="Fraction of Variants",line=2)
    sapply(1:2,function(i){
      axis(3,at=log10(c(300,6000))[i],labels=NA,tck=-0.01)
      axis(3,at=log10(c(300,6000))[i],tick=F,line=-0.9,cex.axis=0.8,
           labels=c("Alu","L1")[i],font=3)
    })
    axis(3,at=log10(c(1000,2000)),labels=NA,tck=-0.01)
    axis(3,at=mean(log10(c(1000,2000))),tick=F,line=-0.9,cex.axis=0.8,labels="SVA",font=3)
    mtext(3,line=1.5,text=title,font=2)
    
    #Add points & rolling mean per AF tranche
    col.pal <- rev(colorRampPalette(c("#440154","#365C8C","#25A584","#FDE725"))(length(sizes)))
    sapply(1:length(dens), function(i){
      #Points per individual bin
      points(x=mids, y=dens[[i]], pch=19, cex=0.25, col=col.pal[i])
      #Rolling mean for line
      points(x=mids, y=rollapply(dens[[i]], width=5, mean, partial=T),
             type="l", lwd=lwd.cex, col=col.pal[i])
    })  
  }else{
    par(bty="n",mar=c(3.5,3.5,3,0.5))
    plot(x=c(0,1),y=c(0,1),type="n",
         xaxt="n",yaxt="n",xlab="",ylab="",yaxs="i")
    text(x=0.5,y=0.5,labels="No Data")
    mtext(3,line=1.5,text=title,font=2,cex=lwd.cex)
  }
  
  #Add filter labels
  if(!is.null(filter.legend)){
    legend("topright",bg=NA,bty="n",pch=NA,legend=filter.legend)
  }
  
  #Add freq legend
  legend("right",bg="white",bty="n",lwd=3,col=col.pal,legend=legend.labs,cex=0.8)
}

#Wrapper to plot all size distributions
wrapperPlotAllSizeDistribs <- function(){
  
  #All SV
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/size_distribution.all.pdf",sep=""),
      height=4,width=6)
  plotSizeDistrib(dat=dat,svtypes=svtypes,
                  title="Size Distribution",
                  legend=T)
  dev.off()
  
  #Singletons
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/size_distribution.singletons.pdf",sep=""),
      height=4,width=6)
  plotSizeDistrib(dat=dat[which(dat$AC==1),],svtypes=svtypes,
                  autosomal=F,biallelic=T,
                  title="Size Distribution (Singletons; AC = 1)",
                  legend=T)
  dev.off()
  
  #Rare (>1 & <1%)
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/size_distribution.rare.pdf",sep=""),
      height=4,width=6)
  plotSizeDistrib(dat=dat[which(dat$AC>1 & dat$AF<rare.max.freq),],svtypes=svtypes,
                  autosomal=F, biallelic=T,
                  title="Size Distribution (AC > 1, AF < 1%)",
                  legend=T)
  dev.off()
  
  #Uncommon (≥1% & <10%)
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/size_distribution.uncommon.pdf",sep=""),
      height=4,width=6)
  plotSizeDistrib(dat=dat[which(dat$AF>=rare.max.freq & dat$AF<uncommon.max.freq),],svtypes=svtypes,
                  autosomal=F, biallelic=T,
                  title="Size Distribution (AF 1% - 10%)",
                  legend=T)
  dev.off()
  
  #Common (≥10% & <50%)
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/size_distribution.common.pdf",sep=""),
      height=4,width=6)
  plotSizeDistrib(dat=dat[which(dat$AF>=uncommon.max.freq & dat$AF<common.max.freq),],svtypes=svtypes,
                  autosomal=F, biallelic=T,
                  title="Size Distribution (AF 10% - 50%)",
                  legend=T)
  dev.off()
  
  #Major (≥50%%)
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/size_distribution.major.pdf",sep=""),
      height=4,width=6)
  plotSizeDistrib(dat=dat[which(dat$AF>=common.max.freq),],svtypes=svtypes,
                  autosomal=F, biallelic=T,
                  title="Size Distribution (AF > 50%)",
                  legend=T)
  dev.off()
  
  #Frequency series
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/size_distribution.across_freqs.pdf",sep=""),
      height=4,width=6)
  plotSizeDistribSeries(dat=dat,svtypes=svtypes,
                        max.AFs=c(1.1/(2*nsamp),rare.max.freq,uncommon.max.freq,
                                  common.max.freq,major.max.freq),
                        legend.labs=c("Singleton","<1%","1-10%","10-50%",">50%"),
                        title="Size Distributions by AF")
  dev.off()
  
  #Merged
  pdf(paste(OUTDIR,"/main_plots/size_distributions.pdf",sep=""),
      height=6,width=10)
  layout(matrix(c(1,1,1,2,2,3,4,5,6,7),byrow=T,nrow=2),
         heights=c(4,2))
  plotSizeDistrib(dat=dat,svtypes=svtypes,
                  title="Size Distribution",
                  legend=T, lwd.cex=1.5)
  plotSizeDistribSeries(dat=dat,svtypes=svtypes,
                        max.AFs=c(1.1/(2*nsamp),rare.max.freq,uncommon.max.freq,
                                  common.max.freq,major.max.freq),
                        legend.labs=c("Singleton","<1%","1-10%","10-50%",">50%"),
                        title="Size Distributions by AF",
                        lwd.cex=2)
  plotSizeDistrib(dat=dat[which(dat$AC==1),],svtypes=svtypes,
                  autosomal=F, biallelic=T, title="AC = 1", text.cex=0.75)
  plotSizeDistrib(dat=dat[which(dat$AC>1 & dat$AF<rare.max.freq),],svtypes=svtypes,
                  autosomal=F, biallelic=T, title="n > 1 & AF < 1%", text.cex=0.75)
  plotSizeDistrib(dat=dat[which(dat$AF>=rare.max.freq & dat$AF<uncommon.max.freq),],svtypes=svtypes,
                  autosomal=F,biallelic=T, title="1% - 10%", text.cex=0.75)
  plotSizeDistrib(dat=dat[which(dat$AF>=uncommon.max.freq & dat$AF<common.max.freq),],svtypes=svtypes,
                  autosomal=F, biallelic=T, title="10% - 50%", text.cex=0.75)
  plotSizeDistrib(dat=dat[which(dat$AF>=common.max.freq),],svtypes=svtypes,
                  autosomal=F, biallelic=T, title="> 50%", text.cex=0.75)
  dev.off()
}


###########################
#####Allele frequency plots
###########################
#Plot single AF spectrum
plotFreqDistrib <- function(dat, svtypes,
                            autosomal=F, biallelic=T,
                            title=NULL, lwd.cex=1, legend=F){
  #Process freqs & compute range + breaks
  filter.legend <- NULL
  if(autosomal==T){
    dat <- dat[which(dat$chr %in% sex.chroms),]
    filter.legend <- c(filter.legend,"Autosomal variants only")
  }
  if(biallelic==T){
    dat <- dat[which(dat$other_gts==0 & dat$missing_gts<dat$genotyped_samples),]
    filter.legend <- c(filter.legend,"Biallelic variants only")
  }
  freqs <- log10(dat$AF)
  if(length(freqs)>0){
    xlims <- range(freqs[which(!is.infinite(freqs))],na.rm=T)
    breaks <- seq(xlims[1],xlims[2],by=(xlims[2]-xlims[1])/(25*abs(floor(xlims)[1])))
    mids <- (breaks[1:(length(breaks)-1)]+breaks[2:length(breaks)])/2
    
    #Gather freq densities per class
    dens <- lapply(svtypes$svtype,function(svtype){
      vals <- freqs[which(dat$svtype==svtype)]
      h <- hist(vals[which(!is.infinite(vals) & vals>=xlims[1] & vals<=xlims[2])],plot=F,breaks=breaks)
      h$counts[1] <- h$counts[1]+length(which(!is.infinite(vals) & vals<xlims[1]))
      h$counts[length(h$counts)] <- h$counts[length(h$counts)]+length(which(!is.infinite(vals) & vals>xlims[2]))
      return(h$counts/length(vals))
    })
    names(dens) <- svtypes$svtype
    all.vals <- freqs[which(!is.infinite(freqs) & freqs>=xlims[1] & freqs<=xlims[2])]
    all.h <- hist(all.vals,plot=F,breaks=breaks)
    dens$ALL <- as.numeric(all.h$counts/length(all.vals))
    
    #Prepare plot area
    ylims <- c(0, quantile(unlist(dens), probs=0.99, na.rm=T))
    par(bty="n",mar=c(4.5,3.5,3,0.5))
    plot(x=xlims,y=ylims,type="n",
         xaxt="n",yaxt="n",xlab="",ylab="",yaxs="i")
    
    #Add vertical gridlines
    logscale.all <- log10(as.numeric(sapply(min(floor(xlims)):0,function(i){(1:9)*10^i})))
    logscale.minor <- log10(as.numeric(sapply(min(floor(xlims)):0,function(i){c(5,10)*10^i})))
    logscale.minor.labs <- as.character(paste(100*round(10^logscale.minor,10),"%",sep=""))
    logscale.major <- log10(as.numeric(10^(min(floor(xlims)):0)))
    abline(v=logscale.all,col="gray97")
    abline(v=logscale.minor,col="gray92")
    abline(v=logscale.major,col="gray85")
    
    #Add axes & title
    axis(1,at=logscale.all,tck=-0.015,col="gray50",labels=NA)
    axis(1,at=logscale.minor,tck=-0.0225,col="gray20",labels=NA)
    axis(1,at=logscale.major,tck=-0.03,labels=NA)
    axis(1,at=logscale.minor,tick=F,cex.axis=0.8,line=-0.4,las=2,
         labels=logscale.minor.labs)
    mtext(1,text="AF",line=3,cex=lwd.cex)
    axis(2,at=axTicks(2),tck=-0.025,labels=NA)
    axis(2,at=axTicks(2),tick=F,line=-0.4,cex.axis=0.8,las=2,
         labels=paste(round(100*axTicks(2),1),"%",sep=""))
    mtext(2,text="Fraction of Variants",line=2.2,cex=lwd.cex)
    mtext(3,line=1.5,text=title,font=2,cex=lwd.cex)
    
    #Add dots & smoothed lines per SV type
    sapply(1:length(dens),function(i){
      svtype <- names(dens)[i]
      if(svtype=="ALL"){
        color <- "gray15"
        lwd <- 4
      }else{
        color <- svtypes[which(svtypes$svtype==svtype),2]
        lwd <- 2
      }
      if(all(!is.nan(dens[[i]]))){
        if(any(dens[[i]]>0)){
          points(x=mids[which(dens[[i]]>0)],y=dens[[i]][which(dens[[i]]>0)],col=color,pch=19,cex=0.3*lwd.cex)
          points(x=mids[which(dens[[i]]>0)],
                 y=rollapply(dens[[i]][which(dens[[i]]>0)],width=10,mean,partial=T),
                 type="l",lwd=lwd.cex*lwd,col=color)
        }
      }
    })
    
    #Add legend of svtypes, if optioned
    if(legend==T){
      idx.for.legend <- which(unlist(lapply(dens,function(vals){any(!is.na(vals) & !is.infinite(vals) & vals>0)})))
      legend("right",bg=NA,bty="n",pch=NA,cex=lwd.cex*0.7,lwd=3,
             legend=rbind(svtypes,c("ALL","gray15"))$svtype[idx.for.legend],
             col=rbind(svtypes,c("ALL","gray15"))$color[idx.for.legend])
    }
  }else{
    par(bty="n",mar=c(4.5,3.5,3,0.5))
    plot(x=c(0,1),y=c(0,1),type="n",
         xaxt="n",yaxt="n",xlab="",ylab="",yaxs="i")
    text(x=0.5,y=0.5,labels="No Data")
    mtext(3,line=1.5,text=title,font=2,cex=lwd.cex)
  }
  
  #Add number of SV to plot
  axis(3,at=mean(par("usr")[1:2]),line=-0.9,tick=F,
       labels=paste("n=",prettyNum(length(freqs),big.mark=","),sep=""))
  
  #Add filter labels
  if(!is.null(filter.legend)){
    legend("topright",bg=NA,bty="n",pch=NA,legend=filter.legend,cex=lwd.cex*0.8)
  }
}

#Plot AF spectrum series by sizes
plotFreqDistribSeries <- function(dat, svtypes, max.sizes, legend.labs,
                                  autosomal=F, biallelic=T, title=NULL){
  #Process freqs & compute range + breaks
  filter.legend <- NULL
  if(autosomal==T){
    dat <- dat[which(dat$chr %in% sex.chroms),]
    filter.legend <- c(filter.legend,"Autosomal variants only")
  }
  if(biallelic==T){
    dat <- dat[which(dat$other_gts==0 & dat$missing_gts<dat$genotyped_samples),]
    filter.legend <- c(filter.legend,"Biallelic variants only")
  }
  freqs <- log10(dat$AF)
  if(length(freqs) > 0){
    freqs <- lapply(1:length(max.sizes),function(i){
      if(i==1){
        return(freqs[which(dat$length<=max.sizes[i])])
      }else{
        return(freqs[which(dat$length>max.sizes[i-1] & dat$length<=max.sizes[i])])
      }
    })
    xlims <- range(freqs[which(!is.infinite(unlist(freqs)))],na.rm=T)
    breaks <- seq(xlims[1],xlims[2],by=(xlims[2]-xlims[1])/(20*abs(floor(xlims)[1])))
    mids <- (breaks[1:(length(breaks)-1)]+breaks[2:length(breaks)])/2
    
    #Gather freq densities per class
    dens <- lapply(freqs,function(vals){
      h <- hist(vals[which(!is.infinite(vals) & vals>=xlims[1] & vals<=xlims[2])],plot=F,breaks=breaks)
      h$counts[1] <- h$counts[1]+length(which(!is.infinite(vals) & vals<xlims[1]))
      h$counts[length(h$counts)] <- h$counts[length(h$counts)]+length(which(!is.infinite(vals) & vals>xlims[2]))
      return(h$counts/length(vals))
    })
    
    #Prepare plot area
    ylims <- c(0, quantile(unlist(dens), probs=0.99, na.rm=T))
    par(bty="n",mar=c(4.5,3.5,3,0.5))
    plot(x=xlims,y=ylims,type="n",
         xaxt="n",yaxt="n",xlab="",ylab="",yaxs="i")
    
    #Add vertical gridlines
    logscale.all <- log10(as.numeric(sapply(min(floor(xlims)):0,function(i){(1:9)*10^i})))
    logscale.minor <- log10(as.numeric(sapply(min(floor(xlims)):0,function(i){c(5,10)*10^i})))
    logscale.minor.labs <- as.character(paste(100*round(10^logscale.minor,10),"%",sep=""))
    logscale.major <- log10(as.numeric(10^(min(floor(xlims)):0)))
    abline(v=logscale.all,col="gray97")
    abline(v=logscale.minor,col="gray92")
    abline(v=logscale.major,col="gray85")
    
    #Add axes & title
    axis(1,at=logscale.all,tck=-0.015,col="gray50",labels=NA)
    axis(1,at=logscale.minor,tck=-0.0225,col="gray20",labels=NA)
    axis(1,at=logscale.major,tck=-0.03,labels=NA)
    axis(1,at=logscale.minor,tick=F,cex.axis=0.8,line=-0.4,las=2,
         labels=logscale.minor.labs)
    mtext(1,text="AF",line=3)
    axis(2,at=axTicks(2),tck=-0.025,labels=NA)
    axis(2,at=axTicks(2),tick=F,line=-0.4,cex.axis=0.8,las=2,
         labels=paste(round(100*axTicks(2),1),"%",sep=""))
    mtext(2,text="Fraction of Variants",line=2.2)
    mtext(3,line=1.5,text=title,font=2)
    
    #Add points & rolling mean lines per size tranche
    col.pal <- colorRampPalette(c("#440154","#365C8C","#25A584","#FDE725"))(length(freqs))
    sapply(1:length(dens),function(i){
      if(all(!is.nan(dens[[i]]))){
        if(any(dens[[i]]>0)){
          points(x=mids[which(dens[[i]]>0)],y=dens[[i]][which(dens[[i]]>0)],col=col.pal[i],pch=19,cex=0.3)
          points(x=mids[which(dens[[i]]>0)],
                 y=rollapply(dens[[i]][which(dens[[i]]>0)],3,mean,partial=T),
                 type="l",lwd=2,col=col.pal[i])
        }
      }
    })
  }else{
    par(bty="n",mar=c(4.5,3.5,3,0.5))
    plot(x=c(0,1),y=c(0,1),type="n",
         xaxt="n",yaxt="n",xlab="",ylab="",yaxs="i")
    text(x=0.5,y=0.5,labels="No Data")
    mtext(3,line=1.5,text=title,font=2,cex=lwd.cex)
  }
  
  #Add filter labels
  if(!is.null(filter.legend)){
    legend("topright",bg=NA,bty="n",pch=NA,legend=filter.legend,cex=0.8)
  }
  
  #Add size legend
  legend("right",bg="white",bty="n",lwd=3,col=col.pal,cex=0.7,
         legend=gsub("\n","",legend.labs,fixed=T))
}

#Wrapper to plot all AF distributions
wrapperPlotAllFreqDistribs <- function(){
  #All SV
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/freq_distribution.all.pdf",sep=""),
      height=4,width=4)
  plotFreqDistrib(dat=dat,svtypes=svtypes,
                  title="AF Distribution",
                  legend=T)
  dev.off()
  
  #Tiny (<50bp)
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/freq_distribution.tiny.pdf",sep=""),
      height=4,width=4)
  plotFreqDistrib(dat=dat[which(dat$length<tiny.max.size),],svtypes=svtypes,
                  title="AF Distribution (< 50bp)",
                  legend=T)
  dev.off()
  
  #Small (50-100bp)
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/freq_distribution.small.pdf",sep=""),
      height=4,width=4)
  plotFreqDistrib(dat=dat[which(dat$length>=tiny.max.size & dat$length<small.max.size),],svtypes=svtypes,
                  title="AF Distribution (50 - 100bp)",
                  legend=T)
  dev.off()
  
  #Medium (100-500bp)
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/freq_distribution.medium.pdf",sep=""),
      height=4,width=4)
  plotFreqDistrib(dat=dat[which(dat$length>=small.max.size & dat$length<medium.max.size),],svtypes=svtypes,
                  title="AF Distribution (100bp - 500bp)",
                  legend=T)
  dev.off()
  
  #Med-Large (500bp-5kb)
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/freq_distribution.medlarge.pdf",sep=""),
      height=4,width=4)
  plotFreqDistrib(dat=dat[which(dat$length>=medium.max.size & dat$length<medlarge.max.size),],svtypes=svtypes,
                  title="AF Distribution (500bp - 5kb)",
                  legend=T)
  dev.off()
  
  #Large (5-50kb)
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/freq_distribution.large.pdf",sep=""),
      height=4,width=4)
  plotFreqDistrib(dat=dat[which(dat$length>=medlarge.max.size & dat$length<large.max.size),],svtypes=svtypes,
                  title="AF Distribution (5 - 50kb)",
                  legend=T)
  dev.off()
  
  #Huge (>50kb)
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/freq_distribution.huge.pdf",sep=""),
      height=4,width=4)
  plotFreqDistrib(dat=dat[which(dat$length>=large.max.size),],svtypes=svtypes,
                  title="AF Distribution (> 50kb)",
                  legend=T)
  dev.off()
  
  #Size series
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/freq_distribution.across_sizes.pdf",sep=""),
      height=4,width=4)
  plotFreqDistribSeries(dat=dat,svtypes=svtypes,
                        max.sizes=c(tiny.max.size,small.max.size,medium.max.size,
                                    medlarge.max.size,large.max.size,huge.max.size),
                        legend.labs=c("<50bp","50-\n100bp","100bp-\n500bp",
                                      "500bp-\n5kb","5-50kb",">50kb"),
                        title="AF Distributions by Size")
  dev.off()
  
  #Merged
  pdf(paste(OUTDIR,"/main_plots/af_distributions.pdf",sep=""),
      height=6,width=10)
  layout(matrix(c(1,1,1,2,2,2,
                  3,4,5,6,7,8),
                byrow=T,nrow=2),
         heights=c(4,2))
  plotFreqDistrib(dat=dat,svtypes=svtypes,
                  title="AF Distribution",
                  legend=T)
  plotFreqDistribSeries(dat=dat,svtypes=svtypes,
                        max.sizes=c(tiny.max.size,small.max.size,medium.max.size,
                                    medlarge.max.size,large.max.size,huge.max.size),
                        legend.labs=c("<50bp","50-\n100bp","100bp-\n500bp",
                                      "500bp-\n5kb","5-50kb",">50kb"),
                        title="AF Distributions by Size")
  plotFreqDistrib(dat=dat[which(dat$length<tiny.max.size),],svtypes=svtypes,
                  title="< 50bp",lwd.cex=0.7)
  plotFreqDistrib(dat=dat[which(dat$length>=tiny.max.size & dat$length<small.max.size),],svtypes=svtypes,
                  title="50 - 100bp",lwd.cex=0.7)
  plotFreqDistrib(dat=dat[which(dat$length>=small.max.size & dat$length<medium.max.size),],svtypes=svtypes,
                  title="100bp - 500bp",lwd.cex=0.7)
  plotFreqDistrib(dat=dat[which(dat$length>=medium.max.size & dat$length<medlarge.max.size),],svtypes=svtypes,
                  title="500bp - 5kb",lwd.cex=0.7)
  plotFreqDistrib(dat=dat[which(dat$length>=medlarge.max.size & dat$length<large.max.size),],svtypes=svtypes,
                  title="5 - 50kb",lwd.cex=0.7)
  plotFreqDistrib(dat=dat[which(dat$length>=large.max.size),],svtypes=svtypes,
                  title="> 50kb",lwd.cex=0.7)
  dev.off()
}


######################################
#####Transition/Transversion Ti/Tv plot
######################################
# Returns Ti:Tv ratio for a data.frame with REF and ALT columns (SNV rows only)
calcTiTv <- function(d){
  transitions <- c("AG","GA","CT","TC")
  is.snv <- !is.na(d$REF) & !is.na(d$ALT) &
             nchar(as.character(d$REF))==1 & nchar(as.character(d$ALT))==1
  d <- d[is.snv,]
  pair <- paste(toupper(as.character(d$REF)), toupper(as.character(d$ALT)), sep="")
  n.ti <- sum(pair %in% transitions, na.rm=T)
  n.tv <- sum(!pair %in% transitions, na.rm=T)
  if(n.tv == 0) return(NA)
  return(n.ti / n.tv)
}

# Heatmap of Ti:Tv ratio with REGION on x-axis and AF bucket on y-axis
plotTiTvHeatmap <- function(snv.dat, af.labels, af.mins, af.maxs, regions, title="Ti:Tv by Region & AF"){
  mat <- sapply(regions, function(reg){
    sapply(seq_along(af.labels), function(i){
      sub <- snv.dat[which(!is.na(snv.dat$REGION) & snv.dat$REGION==reg &
                             snv.dat$AF>af.mins[i] & snv.dat$AF<=af.maxs[i]),]
      calcTiTv(sub)
    })
  })
  rownames(mat) <- af.labels
  colnames(mat) <- regions
  if(all(is.na(mat))){ plot.new(); mtext(3,text=title,font=2); return(invisible(NULL)) }
  zlim <- c(1.5, 2.5)
  col.pal <- colorRampPalette(c("#440154","#365C8C","#25A584","#FDE725"))(101)
  mat.capped <- pmin(pmax(mat, zlim[1]), zlim[2])
  par(mar=c(5,5,3,8), bty="n")
  image(x=1:ncol(mat), y=1:nrow(mat), z=t(mat.capped), col=col.pal,
        xaxt="n", yaxt="n", xlab="", ylab="", zlim=zlim)
  axis(1, at=1:ncol(mat), labels=colnames(mat), las=2, cex.axis=0.75, tick=F, line=-0.5)
  axis(2, at=1:nrow(mat), labels=rownames(mat), las=2, cex.axis=0.75, tick=F, line=-0.5)
  mtext(1, text="Region", line=3.5, cex=0.85)
  mtext(2, text="AF Bucket", line=3.5, cex=0.85)
  mtext(3, text=title, font=2, line=1)
  # Color bar
  usr <- par("usr")
  color.bar.x <- usr[2] + (usr[2]-usr[1])*0.10
  color.bar.w <- (usr[2]-usr[1])*0.06
  color.bar.y <- seq(usr[3], usr[4], length.out=102)
  rect(xleft=color.bar.x, xright=color.bar.x+color.bar.w,
       ybottom=color.bar.y[-length(color.bar.y)], ytop=color.bar.y[-1],
       col=col.pal, border=NA, xpd=T)
  axis(4, at=seq(usr[3],usr[4],length.out=5), tick=F, las=2, cex.axis=0.65, line=3.5,
       labels=round(seq(zlim[1],zlim[2],length.out=5),2), xpd=T)
}

wrapperPlotTiTv <- function(){
  if(!all(c("REF","ALT") %in% colnames(dat))) return(invisible(NULL))
  snv.dat <- dat[which(dat$svtype=="SNV"),]
  if(nrow(snv.dat)==0) return(invisible(NULL))

  af.labels <- c("AC=1","AF<1%","1-10%","10-50%",">50%")
  af.mins <- c(0, 0, rare.max.freq, uncommon.max.freq, common.max.freq)
  af.maxs <- c(1.1/(2*nsamp), rare.max.freq, uncommon.max.freq, common.max.freq, 1)
  size.labels <- c("<50bp","50-100bp","100bp-500bp","500bp-5kb","5-50kb",">50kb")
  size.mins <- c(0, tiny.max.size, small.max.size, medium.max.size, medlarge.max.size, large.max.size)
  size.maxs <- c(tiny.max.size, small.max.size, medium.max.size, medlarge.max.size, large.max.size, huge.max.size)

  has.region <- "REGION" %in% colnames(snv.dat) && any(!is.na(snv.dat$REGION))
  regions <- if(has.region) sort(unique(snv.dat$REGION[!is.na(snv.dat$REGION)])) else character(0)

  # Number of columns: AF panel + Size panel + optional heatmap
  n.panels <- 2 + as.integer(has.region && length(regions)>0)
  png(paste(OUTDIR,"/main_plots/ti_tv_distributions.png",sep=""),
      res=300, height=1800, width=n.panels*1800)
  layout(matrix(1:n.panels, nrow=1))

  # Panel 1: Ti:Tv by AF bucket
  titv.af <- sapply(seq_along(af.labels), function(i){
    sub <- snv.dat[which(snv.dat$AF>af.mins[i] & snv.dat$AF<=af.maxs[i]),]
    calcTiTv(sub)
  })
  par(mar=c(7,4.5,3,0.5), bty="n")
  plot(x=c(0,length(af.labels)+1), y=c(0, max(titv.af,na.rm=T)*1.15),
       type="n", xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i")
  col.af <- colorRampPalette(c("#440154","#365C8C","#25A584","#FDE725"))(length(af.labels))
  rect(xleft=1:length(af.labels)-0.35, xright=1:length(af.labels)+0.35,
       ybottom=0, ytop=titv.af, col=col.af, border=NA)
  axis(1, at=1:length(af.labels), labels=af.labels, las=2, cex.axis=0.8, tick=F, line=0.5)
  axis(2, at=axTicks(2), labels=NA); axis(2, at=axTicks(2), tick=F, las=2, cex.axis=0.8, line=-0.4)
  mtext(1, text="AF Bucket", line=5.5, cex=0.9)
  mtext(2, text="Ti:Tv Ratio", line=2.5, cex=0.9)
  mtext(3, text="Ti:Tv by AF", font=2, line=0.5)

  # Panel 2: Ti:Tv by REGION (or show "No REGION data")
  if(has.region && length(regions)>0){
    titv.reg <- sapply(regions, function(reg){
      sub <- snv.dat[which(!is.na(snv.dat$REGION) & snv.dat$REGION==reg),]
      calcTiTv(sub)
    })
    col.reg <- colorRampPalette(c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02"))(length(regions))
    par(mar=c(7,4.5,3,0.5), bty="n")
    plot(x=c(0,length(regions)+1), y=c(0, max(titv.reg,na.rm=T)*1.15),
         type="n", xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i")
    rect(xleft=1:length(regions)-0.35, xright=1:length(regions)+0.35,
         ybottom=0, ytop=titv.reg, col=col.reg, border=NA)
    axis(1, at=1:length(regions), labels=regions, las=2, cex.axis=0.8, tick=F, line=0.5)
    axis(2, at=axTicks(2), labels=NA); axis(2, at=axTicks(2), tick=F, las=2, cex.axis=0.8, line=-0.4)
    mtext(1, text="Region", line=5.5, cex=0.9)
    mtext(2, text="Ti:Tv Ratio", line=2.5, cex=0.9)
    mtext(3, text="Ti:Tv by Region", font=2, line=0.5)

    # Panel 3: Ti:Tv heatmap (Region x AF)
    plotTiTvHeatmap(snv.dat=snv.dat, af.labels=af.labels, af.mins=af.mins, af.maxs=af.maxs,
                    regions=regions, title="Ti:Tv by Region & AF")
  }else{
    par(mar=c(7,4.5,3,0.5), bty="n")
    plot(x=c(0,1),y=c(0,1),type="n",xaxt="n",yaxt="n",xlab="",ylab="")
    text(x=0.5,y=0.5,labels="No REGION data")
    mtext(3,text="Ti:Tv by Region",font=2,line=0.5)
  }
  dev.off()
}


######################################
#####Quality score distribution plot
######################################
plotDistribOverlaid <- function(sub.list, sub.labels, sub.colors, main.val=NULL,
                                xlab="QUAL", main.label="All Variants",
                                title=NULL, xlim=NULL, alpha=0.5){
  # Combine all values to get x range
  all.vals <- unlist(lapply(sub.list, function(x) x[!is.na(x) & is.finite(x)]))
  if(length(all.vals)==0){ par(bty="n",mar=c(4.5,4,3,1)); plot.new(); mtext(3,text=title,font=2,line=1); return(invisible(NULL)) }
  if(is.null(xlim)) xlim <- range(all.vals, na.rm=T)
  breaks <- seq(xlim[1], xlim[2], length.out=51)

  # Compute densities
  dens.list <- lapply(sub.list, function(vals){
    vals <- vals[!is.na(vals) & is.finite(vals) & vals>=xlim[1] & vals<=xlim[2]]
    if(length(vals)<2) return(NULL)
    h <- hist(vals, breaks=breaks, plot=F)
    h$density
  })
  main.dens <- NULL
  if(!is.null(main.val)){
    mv <- main.val[!is.na(main.val) & is.finite(main.val) & main.val>=xlim[1] & main.val<=xlim[2]]
    if(length(mv)>=2) main.dens <- hist(mv, breaks=breaks, plot=F)$density
  }
  ylim <- c(0, max(unlist(c(lapply(dens.list,max), list(max(main.dens,na.rm=T)))), na.rm=T)*1.1)
  mids <- (breaks[-length(breaks)] + breaks[-1])/2

  par(bty="n", mar=c(4.5,4,3,1))
  plot(x=xlim, y=ylim, type="n", xaxt="n", yaxt="n", xlab="", ylab="", yaxs="i")
  axis(1, at=pretty(xlim), labels=NA); axis(1, at=pretty(xlim), tick=F, line=-0.4, cex.axis=0.8)
  axis(2, at=axTicks(2), labels=NA); axis(2, at=axTicks(2), tick=F, las=2, cex.axis=0.8, line=-0.4)
  mtext(1, text=xlab, line=3, cex=0.9); mtext(2, text="Density", line=2.5, cex=0.9)
  mtext(3, text=title, font=2, line=1)

  # Draw main (all variants) density first as dark line
  if(!is.null(main.dens)){
    polygon(x=c(mids[1],mids,mids[length(mids)]), y=c(0,main.dens,0),
            col=adjustcolor("gray25",alpha.f=0.3), border="gray25", lwd=1.5)
  }
  # Draw stratified overlays
  for(i in seq_along(sub.list)){
    if(is.null(dens.list[[i]])) next
    polygon(x=c(mids[1],mids,mids[length(mids)]), y=c(0,dens.list[[i]],0),
            col=adjustcolor(sub.colors[i],alpha.f=alpha), border=sub.colors[i], lwd=1.2)
  }
  legend("topright", bty="n", cex=0.7, lwd=2,
         col=c("gray25",sub.colors),
         legend=c(main.label, sub.labels))
}

wrapperPlotQualDistrib <- function(){
  if(!"QUAL" %in% colnames(dat)) return(invisible(NULL))
  qual.vals <- dat$QUAL[!is.na(dat$QUAL)]
  if(length(qual.vals)==0) return(invisible(NULL))

  has.region <- "REGION" %in% colnames(dat) && any(!is.na(dat$REGION))
  regions <- if(has.region) sort(unique(dat$REGION[!is.na(dat$REGION)])) else character(0)
  af.labels <- c("AC=1","AF<1%","1-10%","10-50%",">50%")
  af.mins <- c(0, 0, rare.max.freq, uncommon.max.freq, common.max.freq)
  af.maxs <- c(1.1/(2*nsamp), rare.max.freq, uncommon.max.freq, common.max.freq, 1)
  size.labels <- c("<50bp","50-100bp","100bp-500bp","500bp-5kb","5-50kb",">50kb")
  size.mins <- c(0, tiny.max.size, small.max.size, medium.max.size, medlarge.max.size, large.max.size)
  size.maxs <- c(tiny.max.size, small.max.size, medium.max.size, medlarge.max.size, large.max.size, huge.max.size)

  col.svtype <- svtypes$color[match(svtypes$svtype, svtypes$svtype)]
  col.af <- colorRampPalette(c("#440154","#365C8C","#25A584","#FDE725"))(length(af.labels))
  col.size <- colorRampPalette(c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02"))(length(size.labels))
  col.reg <- colorRampPalette(c("#1F78B4","#33A02C","#E31A1C","#FF7F00","#6A3D9A","#B15928"))(length(regions))

  n.panels <- 4
  png(paste(OUTDIR,"/main_plots/qual_distributions.png",sep=""),
      res=300, height=1800, width=n.panels*1800)
  layout(matrix(1:n.panels, nrow=1))

  # Panel 1: All variants
  par(bty="n", mar=c(4.5,4,3,1))
  q.all <- pmin(dat$QUAL[!is.na(dat$QUAL) & is.finite(dat$QUAL)], 99)
  xlim <- c(min(q.all, na.rm=T), 99)
  breaks <- seq(xlim[1],xlim[2],length.out=51)
  h <- hist(q.all, breaks=breaks, plot=F)
  plot(x=xlim, y=c(0,max(h$density)*1.1), type="n", xaxt="n", yaxt="n", xlab="", ylab="", yaxs="i")
  polygon(c(h$mids[1],h$mids,h$mids[length(h$mids)]), c(0,h$density,0), col="#4393C3", border="#2166AC")
  axis(1,at=pretty(xlim),labels=NA); axis(1,at=pretty(xlim),tick=F,line=-0.4,cex.axis=0.8)
  axis(2,at=axTicks(2),labels=NA); axis(2,at=axTicks(2),tick=F,las=2,cex.axis=0.8,line=-0.4)
  mtext(1,text="QUAL",line=3,cex=0.9); mtext(2,text="Density",line=2.5,cex=0.9)
  mtext(3,text="QUAL Distribution",font=2,line=1)
  axis(3,at=mean(xlim),tick=F,line=-0.9,labels=paste("n=",prettyNum(length(q.all),big.mark=","),sep=""))

  # Panel 2: By region
  reg.vals <- if(has.region && length(regions)>0) lapply(regions, function(r) pmin(dat$QUAL[!is.na(dat$REGION) & dat$REGION==r & !is.na(dat$QUAL)], 99)) else list()
  plotDistribOverlaid(sub.list=reg.vals, sub.labels=regions, sub.colors=col.reg,
                      main.val=q.all, xlab="QUAL", main.label="All",
                      title="QUAL by Region", xlim=xlim)

  # Panel 3: By size bucket
  size.vals <- lapply(seq_along(size.labels), function(i) pmin(dat$QUAL[!is.na(dat$QUAL) & dat$length>=size.mins[i] & dat$length<=size.maxs[i]], 99))
  plotDistribOverlaid(sub.list=size.vals, sub.labels=size.labels, sub.colors=col.size,
                      main.val=q.all, xlab="QUAL", main.label="All",
                      title="QUAL by Size", xlim=xlim)

  # Panel 4: By AF bucket
  af.vals <- lapply(seq_along(af.labels), function(i) pmin(dat$QUAL[!is.na(dat$QUAL) & dat$AF>af.mins[i] & dat$AF<=af.maxs[i]], 99))
  plotDistribOverlaid(sub.list=af.vals, sub.labels=af.labels, sub.colors=col.af,
                      main.val=q.all, xlab="QUAL", main.label="All",
                      title="QUAL by AF", xlim=xlim)
  dev.off()
}


######################################
#####NCR distribution plot
######################################
wrapperPlotNcrDistrib <- function(){
  if(!"NCR" %in% colnames(dat)) return(invisible(NULL))
  ncr.vals <- dat$NCR[!is.na(dat$NCR)]
  if(length(ncr.vals)==0) return(invisible(NULL))

  has.region <- "REGION" %in% colnames(dat) && any(!is.na(dat$REGION))
  regions <- if(has.region) sort(unique(dat$REGION[!is.na(dat$REGION)])) else character(0)
  af.labels <- c("AC=1","AF<1%","1-10%","10-50%",">50%")
  af.mins <- c(0, 0, rare.max.freq, uncommon.max.freq, common.max.freq)
  af.maxs <- c(1.1/(2*nsamp), rare.max.freq, uncommon.max.freq, common.max.freq, 1)
  size.labels <- c("<50bp","50-100bp","100bp-500bp","500bp-5kb","5-50kb",">50kb")
  size.mins <- c(0, tiny.max.size, small.max.size, medium.max.size, medlarge.max.size, large.max.size)
  size.maxs <- c(tiny.max.size, small.max.size, medium.max.size, medlarge.max.size, large.max.size, huge.max.size)

  col.af <- colorRampPalette(c("#440154","#365C8C","#25A584","#FDE725"))(length(af.labels))
  col.size <- colorRampPalette(c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02"))(length(size.labels))
  col.reg <- colorRampPalette(c("#1F78B4","#33A02C","#E31A1C","#FF7F00","#6A3D9A","#B15928"))(length(regions))
  xlim <- c(0,1)

  n.panels <- 4
  png(paste(OUTDIR,"/main_plots/ncr_distributions.png",sep=""),
      res=300, height=1800, width=n.panels*1800)
  layout(matrix(1:n.panels, nrow=1))

  # Panel 1: All variants
  par(bty="n", mar=c(4.5,4,3,1))
  n.all <- dat$NCR[!is.na(dat$NCR)]
  breaks <- seq(0,1,by=0.02)
  h <- hist(n.all, breaks=breaks, plot=F)
  plot(x=xlim, y=c(0,max(h$density)*1.1), type="n", xaxt="n", yaxt="n", xlab="", ylab="", yaxs="i")
  polygon(c(h$mids[1],h$mids,h$mids[length(h$mids)]), c(0,h$density,0), col="#4393C3", border="#2166AC")
  axis(1,at=seq(0,1,0.2),labels=NA); axis(1,at=seq(0,1,0.2),tick=F,line=-0.4,cex.axis=0.8,
       labels=paste(seq(0,100,20),"%",sep=""))
  axis(2,at=axTicks(2),labels=NA); axis(2,at=axTicks(2),tick=F,las=2,cex.axis=0.8,line=-0.4)
  mtext(1,text="NCR",line=3,cex=0.9); mtext(2,text="Density",line=2.5,cex=0.9)
  mtext(3,text="NCR Distribution",font=2,line=1)
  axis(3,at=0.5,tick=F,line=-0.9,labels=paste("n=",prettyNum(length(n.all),big.mark=","),sep=""))

  # Panel 2: By region
  reg.vals <- if(has.region && length(regions)>0) lapply(regions, function(r) dat$NCR[!is.na(dat$REGION) & dat$REGION==r & !is.na(dat$NCR)]) else list()
  plotDistribOverlaid(sub.list=reg.vals, sub.labels=regions, sub.colors=col.reg,
                      main.val=dat$NCR, xlab="NCR", main.label="All",
                      title="NCR by Region", xlim=xlim)

  # Panel 3: By size bucket
  size.vals <- lapply(seq_along(size.labels), function(i) dat$NCR[!is.na(dat$NCR) & dat$length>=size.mins[i] & dat$length<=size.maxs[i]])
  plotDistribOverlaid(sub.list=size.vals, sub.labels=size.labels, sub.colors=col.size,
                      main.val=dat$NCR, xlab="NCR", main.label="All",
                      title="NCR by Size", xlim=xlim)

  # Panel 4: By AF bucket
  af.vals <- lapply(seq_along(af.labels), function(i) dat$NCR[!is.na(dat$NCR) & dat$AF>af.mins[i] & dat$AF<=af.maxs[i]])
  plotDistribOverlaid(sub.list=af.vals, sub.labels=af.labels, sub.colors=col.af,
                      main.val=dat$NCR, xlab="NCR", main.label="All",
                      title="NCR by AF", xlim=xlim)
  dev.off()
}


######################################
#####gnomAD match rate distribution plot
######################################
# Returns gnomAD match rate (proportion of variants with non-empty gnomAD_V4_match_ID)
calcGnomadMatchRate <- function(d){
  gm <- d$gnomAD_V4_match_ID
  if(length(gm)==0) return(NA)
  n.match <- sum(!is.na(gm) & gm != "" & gm != ".", na.rm=T)
  n.total <- length(gm)
  if(n.total==0) return(NA)
  return(n.match/n.total)
}

wrapperPlotGnomadMatchDistrib <- function(){
  if(!"gnomAD_V4_match_ID" %in% colnames(dat)) return(invisible(NULL))

  has.region <- "REGION" %in% colnames(dat) && any(!is.na(dat$REGION))
  regions <- if(has.region) sort(unique(dat$REGION[!is.na(dat$REGION)])) else character(0)
  af.labels <- c("AC=1","AF<1%","1-10%","10-50%",">50%")
  af.mins <- c(0, 0, rare.max.freq, uncommon.max.freq, common.max.freq)
  af.maxs <- c(1.1/(2*nsamp), rare.max.freq, uncommon.max.freq, common.max.freq, 1)
  size.labels <- c("<50bp","50-100bp","100bp-500bp","500bp-5kb","5-50kb",">50kb")
  size.mins <- c(0, tiny.max.size, small.max.size, medium.max.size, medlarge.max.size, large.max.size)
  size.maxs <- c(tiny.max.size, small.max.size, medium.max.size, medlarge.max.size, large.max.size, huge.max.size)

  col.af <- colorRampPalette(c("#440154","#365C8C","#25A584","#FDE725"))(length(af.labels))
  col.size <- colorRampPalette(c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02"))(length(size.labels))
  col.reg <- colorRampPalette(c("#1F78B4","#33A02C","#E31A1C","#FF7F00","#6A3D9A","#B15928"))(length(regions))

  plotMatchRateBars <- function(vals, labs, cols, xlabel, title){
    vals[is.nan(vals)] <- NA
    ylim <- c(0, max(vals, na.rm=T)*1.15)
    par(bty="n", mar=c(7,4.5,3,0.5))
    plot(x=c(0,length(labs)+1), y=ylim, type="n", xaxt="n", yaxt="n",
         xlab="", ylab="", xaxs="i", yaxs="i")
    rect(xleft=1:length(labs)-0.35, xright=1:length(labs)+0.35,
         ybottom=0, ytop=vals, col=cols, border=NA)
    axis(1,at=1:length(labs),labels=labs,las=2,cex.axis=0.8,tick=F,line=0.5)
    axis(2,at=axTicks(2),labels=NA); axis(2,at=axTicks(2),tick=F,las=2,cex.axis=0.8,line=-0.4,
         labels=paste(round(100*axTicks(2),1),"%",sep=""))
    mtext(1,text=xlabel,line=5.5,cex=0.9); mtext(2,text="gnomAD Match Rate",line=3,cex=0.9)
    mtext(3,text=title,font=2,line=0.5)
  }

  # Determine number of panels
  n.panels <- 1 + as.integer(has.region && length(regions)>0) + 2  # overall + region(optional) + size + AF
  if(!has.region || length(regions)==0) n.panels <- 3
  png(paste(OUTDIR,"/main_plots/gnomad_match_distributions.png",sep=""),
      res=300, height=1800, width=n.panels*1800)
  layout(matrix(1:n.panels, nrow=1))

  # Panel 1: By svtype (overall bars)
  svtype.rates <- sapply(svtypes$svtype, function(st) calcGnomadMatchRate(dat[dat$svtype==st,]))
  plotMatchRateBars(vals=svtype.rates, labs=svtypes$svtype, cols=svtypes$color,
                    xlabel="Variant Type", title="gnomAD Match Rate by Variant Type")

  # Optional Panel 2: By region
  if(has.region && length(regions)>0){
    reg.rates <- sapply(regions, function(r) calcGnomadMatchRate(dat[!is.na(dat$REGION) & dat$REGION==r,]))
    plotMatchRateBars(vals=reg.rates, labs=regions, cols=col.reg,
                      xlabel="Region", title="gnomAD Match Rate by Region")
  }

  # Panel: By size bucket
  size.rates <- sapply(seq_along(size.labels), function(i) calcGnomadMatchRate(dat[!is.na(dat$length) & dat$length>=size.mins[i] & dat$length<=size.maxs[i],]))
  plotMatchRateBars(vals=size.rates, labs=size.labels, cols=col.size,
                    xlabel="Size Bucket", title="gnomAD Match Rate by Size")

  # Panel: By AF bucket
  af.rates <- sapply(seq_along(af.labels), function(i) calcGnomadMatchRate(dat[!is.na(dat$AF) & dat$AF>af.mins[i] & dat$AF<=af.maxs[i],]))
  plotMatchRateBars(vals=af.rates, labs=af.labels, cols=col.af,
                    xlabel="AF Bucket", title="gnomAD Match Rate by AF")

  dev.off()
}


#########################
#####Hardy-Weinberg plots
#########################
#Plot single HW ternary comparison
plotHWSingle <- function(dat,svtypes,title=NULL,full.legend=T,lab.cex=1){
  #Restrict data to biallelic, autosomal sites
  HW.dat <- dat[which(dat$chr %in% sex.chroms & !is.na(dat$AN)),]
  
  #Only run if there's data
  if(nrow(HW.dat)>0){
    #Prep HW matrix
    nsamps <- HW.dat$reference_gts+HW.dat$het_gts+HW.dat$hom_gts
    HW.mat <- data.frame("AA"=HW.dat$reference_gts/nsamps,
                         "AB"=HW.dat$het_gts/nsamps,
                         "BB"=HW.dat$hom_gts/nsamps)
    
    #Gather HW p-values & colors
    HW.p <- HWChisqStats(X=HW.mat*nsamp,x.linked=F,pvalues=T)
    HW.cols <- rep("#4DAC26",times=length(HW.p))
    HW.cols[which(HW.p<0.05)] <- "#81F850"
    HW.cols[which(HW.p<0.05/length(HW.p))] <- "#AC26A1"
    
    #Subsample points if plot is too dense
    max.plot.pts <- 50000
    if(nrow(HW.mat) > max.plot.pts){
      sub.idx <- sort(sample(1:nrow(HW.mat), max.plot.pts))
      HW.mat.plot <- HW.mat[sub.idx,]
      HW.cols.plot <- HW.cols[sub.idx]
    }else{
      HW.mat.plot <- HW.mat
      HW.cols.plot <- HW.cols
    }
    
    #Generate HW plot frame
    par(mar=c(1,3.5,3,0.5),bty="n")
    plot(x=1.15*c(-1/sqrt(3),1/sqrt(3)),y=c(-0.15,1.15),type="n",
         xaxt="n",yaxt="n",xlab="",ylab="",xaxs="i",yaxs="i")
    segments(x0=c(-1/sqrt(3),0,1/sqrt(3)),
             x1=c(0,1/sqrt(3),-1/sqrt(3)),
             y0=c(0,1,0),y1=c(1,0,0))
    HWTernaryPlot(X=HW.mat.plot,n=nsamp,newframe=F,
                  vbounds=F,mafbounds=F,
                  region=1,vertexlab=NA,
                  alpha=0.05,
                  curvecols=c("#4DAC26","#81F850",NA,NA),pch=NA)
    
    #Add axes
    text(x=c(-1/sqrt(3),1/sqrt(3)),y=0,labels=c("Ref.","Hom."),pos=1,font=4)
    text(x=0,y=1,labels="Het.",pos=3,font=4)
    axis(2,at=seq(0,1,0.2),tck=-0.01,labels=NA,line=0.5)
    axis(2,at=seq(0,1,0.2),tck=0.01,labels=NA,line=0.5)
    axis(2,at=seq(0,1,0.2),line=0.1,tick=F,las=2,cex.axis=0.7,
         labels=paste(seq(0,100,20),"%",sep=""))
    mtext(2,text="Fraction of Genotypes",line=2.5,cex=lab.cex)
    
    #Finish HW plot
    HWTernaryPlot(X=HW.mat.plot,n=nsamp,newframe=F,
                  vbounds=F,mafbounds=F,
                  region=1,vertexlab=NA,
                  alpha=0.05/nrow(HW.mat),
                  curvecols=c("#4DAC26","#AC26A1",NA,NA),
                  pch=21,cex=0.5,signifcolour=F,markercol=HW.cols.plot,
                  markerbgcol=adjustcolor(HW.cols.plot,alpha=0.25))
    
    #Add legend
    n.pass <- length(which(HW.p>=0.05))
    n.nom <- length(which(HW.p<0.05 & HW.p>=0.05/nrow(HW.mat)))
    n.bonf <- length(which(HW.p<0.05/nrow(HW.mat)))
    if(full.legend==T){
      legend("topright",pch=19,col=c("#4DAC26","#81F850","#AC26A1"),pt.cex=2,
             legend=c(paste("In H-W equilibrium\n(n=",
                            prettyNum(n.pass,big.mark=","),"; ",
                            round(100*(n.pass/nrow(HW.mat)),2),"%)\n",sep=""),
                      paste("Not in H-W equilibrium\n(Nominal; n=",
                            prettyNum(n.nom,big.mark=","),"; ",
                            round(100*(n.nom/nrow(HW.mat)),2),"%)\n",sep=""),
                      paste("Not in H-W equilibrium\n(Bonferroni; n=",
                            prettyNum(n.bonf,big.mark=","),"; ",
                            round(100*(n.bonf/nrow(HW.mat)),2),"%)\n",sep="")),
             bty="n",bg=NA,cex=0.7)
    }else{
      legend("topright",pch=19,col=c("#4DAC26","#81F850","#AC26A1"),pt.cex=2,
             legend=c(paste(round(100*(n.pass/nrow(HW.mat)),0),"%",sep=""),
                      paste(round(100*(n.nom/nrow(HW.mat)),0),"%",sep=""),
                      paste(round(100*(n.bonf/nrow(HW.mat)),0),"%",sep="")),
             bty="n",bg=NA,cex=0.7)
    }
    
    #Add number of SV on plot
    axis(3,at=mean(par("usr")[1:2]),line=-0.9,tick=F,cex.axis=0.8,
         labels=paste("n=",prettyNum(nrow(HW.mat),big.mark=","),sep=""))
  }else{
    par(mar=c(1,3.5,3,0.5),bty="n")
    plot(x=c(0,1),y=c(0,1),type="n",
         xaxt="n",yaxt="n",xlab="",ylab="",yaxs="i")
    text(x=0.5,y=0.5,labels="No Data")
  }
  
  #Add title
  mtext(3,line=1,text=title,font=2,cex=lab.cex)
  
  #Add filter labels
  if(full.legend==T){
    legend("right",bg=NA,bty="n",pch=NA,cex=0.7,
           legend=c("Biallelic variants only","Autosomal variants only"))
  }
}
#Correlation of carrier frequency & AF
plotAlleleCarrierCorrelation <- function(dat,autosomal=T,biallelic=T,
                                         title="Carrier Frequency vs. AF"){
  #Process freqs
  filter.legend <- NULL
  if(autosomal==T){
    dat <- dat[which(dat$chr %in% sex.chroms),]
    filter.legend <- c(filter.legend,"Autosomal SV only")
  }
  if(biallelic==T){
    dat <- dat[which(dat$other_gts==0 & dat$missing_gts<dat$genotyped_samples),]
    filter.legend <- c(filter.legend,"Biallelic SV only")
  }
  CF <- sort(dat$carrierFreq)
  AF <- dat$AF[order(dat$carrierFreq)]
  keep <- which(!is.na(CF) & !is.na(AF))
  CF <- CF[keep]
  AF <- AF[keep]
  
  #Prep plot area
  par(mar=c(4,4,2,2))
  plot(x=0:1,y=0:1,type="n",
       xlab="",ylab="",xaxt="n",yaxt="n",xaxs="i",yaxs="i")
  abline(h=seq(0,1,0.2),v=seq(0,1,0.2),col="gray85")
  abline(0,1,col="gray85")
  abline(0,0.5,col="gray85")
  
  #Add axes & title
  axis(1,at=seq(0,1,0.2),labels=NA)
  axis(1,at=seq(0,1,0.2),tick=F,line=-0.4,cex.axis=0.8,
       labels=paste(seq(0,100,20),"%",sep=""))
  mtext(1,text="Carrier Frequency",line=2)
  axis(2,at=seq(0,1,0.2),labels=NA)
  axis(2,at=seq(0,1,0.2),tick=F,line=-0.4,cex.axis=0.8,las=2,
       labels=paste(seq(0,100,20),"%",sep=""))
  mtext(2,text="AF",line=2.5)
  mtext(3,text=title,line=0.5,font=2)
  axis(4,at=c(0.5,1),las=2,line=-0.8,tick=F,col.axis="gray50",cex.axis=0.7,
       labels=c("All\nHet.","All\nHom."))
  
  #Plot points & rolling mean
  points(x=CF,y=AF,cex=0.2)
  points(x=rollmean(CF,k=100),
         y=rollmean(AF,k=100),
         type="l",col="red",lwd=1.25)
  
  #Add legend
  legend("bottomright",bg=NA,legend=c("Site","Rolling Mean"),cex=0.8,bty="n",
         lwd=c(NA,2),col=c("black","red"),pch=c(1,NA),pt.cex=c(0.4,NA))
  
  #Add filter labels
  if(!is.null(filter.legend)){
    legend("topleft",bg=NA,bty="n",pch=NA,legend=filter.legend,cex=0.8)
  }
}

#Wrapper to plot all HW distributions
wrapperPlotAllHWDistribs <- function(){
  #All SV
  png(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/gt_distribution.all.png",sep=""),
      res=300,height=1800,width=1800)
  plotHWSingle(dat=dat,svtypes=svtypes,
               title="Genotype Distribution")
  dev.off()
  #Tiny
  png(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/gt_distribution.tiny.png",sep=""),
      res=300,height=1800,width=1800)
  plotHWSingle(dat=dat[which(dat$length<tiny.max.size),],svtypes=svtypes,
               title="Genotype Distribution (< 50bp)")
  dev.off()
  #Small
  png(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/gt_distribution.small.png",sep=""),
      res=300,height=1800,width=1800)
  plotHWSingle(dat=dat[which(dat$length>=tiny.max.size & dat$length<small.max.size),],svtypes=svtypes,
               title="Genotype Distribution (50 - 100bp)")
  dev.off()
  #Medium
  png(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/gt_distribution.medium.png",sep=""),
      res=300,height=1800,width=1800)
  plotHWSingle(dat=dat[which(dat$length>=small.max.size & dat$length<medium.max.size),],svtypes=svtypes,
               title="Genotype Distribution (100bp - 500bp)")
  dev.off()
  #Med-Large
  png(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/gt_distribution.medlarge.png",sep=""),
      res=300,height=1800,width=1800)
  plotHWSingle(dat=dat[which(dat$length>=medium.max.size & dat$length<medlarge.max.size),],svtypes=svtypes,
               title="Genotype Distribution (500bp - 5kb)")
  dev.off()
  #Large
  png(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/gt_distribution.large.png",sep=""),
      res=300,height=1800,width=1800)
  plotHWSingle(dat=dat[which(dat$length>=medlarge.max.size & dat$length<large.max.size),],svtypes=svtypes,
               title="Genotype Distribution (5 - 50kb)")
  dev.off()
  #Huge
  png(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/gt_distribution.huge.png",sep=""),
      res=300,height=1800,width=1800)
  plotHWSingle(dat=dat[which(dat$length>=large.max.size & dat$length<huge.max.size),],svtypes=svtypes,
               title="Genotype Distribution (> 50kb)")
  dev.off()
  #AF vs carrier frequency
  png(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/carrier_freq_vs_allele_freq.png",sep=""),res=300,
      height=1200,width=1200)
  plotAlleleCarrierCorrelation(dat=dat)
  dev.off()
  #Merged
  png(paste(OUTDIR,"main_plots/genotype_distributions.png",sep=""),res=300,
      height=1200,width=3.5*1200)
  layout(matrix(c(1,2,3,4,5,
                  1,2,6,7,8),
                byrow=T,nrow=2),
         widths=c(2,2,1,1,1))
  plotAlleleCarrierCorrelation(dat=dat)
  plotHWSingle(dat=dat,svtypes=svtypes,
               title="Genotype Distribution")
  plotHWSingle(dat=dat[which(dat$length<tiny.max.size),],svtypes=svtypes,
               title="< 50bp",full.legend=F,lab.cex=0.7)
  plotHWSingle(dat=dat[which(dat$length>=tiny.max.size & dat$length<small.max.size),],svtypes=svtypes,
               title="50 - 100bp",full.legend=F,lab.cex=0.7)
  plotHWSingle(dat=dat[which(dat$length>=small.max.size & dat$length<medium.max.size),],svtypes=svtypes,
               title="100bp - 500bp",full.legend=F,lab.cex=0.7)
  plotHWSingle(dat=dat[which(dat$length>=medium.max.size & dat$length<medlarge.max.size),],svtypes=svtypes,
               title="500bp - 5kb",full.legend=F,lab.cex=0.7)
  plotHWSingle(dat=dat[which(dat$length>=medlarge.max.size & dat$length<large.max.size),],svtypes=svtypes,
               title="5 - 50kb",full.legend=F,lab.cex=0.7)
  plotHWSingle(dat=dat[which(dat$length>=large.max.size & dat$length<huge.max.size),],svtypes=svtypes,
               title="> 50kb",full.legend=F,lab.cex=0.7)
  dev.off()
}


########################
###RSCRIPT FUNCTIONALITY
########################
###Load libraries as needed
require(optparse, quietly=T)
require(RColorBrewer, quietly=T)
require(zoo, quietly=T)
require(HardyWeinberg, quietly=T)

###List of command-line options
option_list <- list(
  make_option(c("-N", "--nsamp"), type="integer", default=NULL,
              help="number of samples to be used for allele frequency calculations [default %default]",
              metavar="integer"),
  make_option(c("-S", "--svtypes"), type="character", default=NULL,
              help="tab-delimited file specifying SV types and HEX colors [default %default]",
              metavar="character")
)

###Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog vcf2bed_cleaned.simple.bed OUTDIR",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

###Checks for appropriate positional arguments
if(length(args$args) != 2){
  stop("Incorrect number of required positional arguments\n")
}

###Writes args & opts to vars
INFILE <- args$args[1]
OUTDIR <- args$args[2]
nsamp <- opts$nsamp
svtypes.file <- opts$svtypes

###Prepares I/O files
#Read & clean data
dat <- read.table(INFILE,comment.char="",sep="\t",header=T,check.names=F)
colnames(dat)[1] <- "chr"
#Create output directory structure, if necessary
if(!dir.exists(OUTDIR)){
  dir.create(OUTDIR)
}
if(!dir.exists(paste(OUTDIR,"/main_plots/",sep=""))){
  dir.create(paste(OUTDIR,"/main_plots/",sep=""))
}
if(!dir.exists(paste(OUTDIR,"/supporting_plots/",sep=""))){
  dir.create(paste(OUTDIR,"/supporting_plots/",sep=""))
}
if(!dir.exists(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/",sep=""))){
  dir.create(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/",sep=""))
}
#Sets sv types & colors
if(!is.null(svtypes.file)){
  svtypes <- read.table(svtypes.file,sep="\t",header=F,comment.char="",check.names=F)
  svtypes <- as.data.frame(apply(svtypes,2,as.character))
  colnames(svtypes) <- c("svtype","color")
}else{
  svtypes.v <- unique(dat$svtype)
  svtypes.c <- brewer.pal(length(svtypes.v),"Dark2")
  svtypes <- data.frame("svtype"=svtypes.v,
                        "color"=svtypes.c)
}

# Merged svtypes and dat for size-stratified plots (DEL_SHORT+DEL_SV → DEL, INS_SHORT+INS_SV → INS, DUP+DUP_SV → DUP)
svtypes.merged <- data.frame(
  svtype = c("DEL","INS","DUP",svtypes$svtype[!svtypes$svtype %in% c("DEL_SHORT","DEL_SV","INS_SHORT","INS_SV","DUP","DUP_SV")]),
  color = c(svtypes$color[svtypes$svtype=="DEL_SHORT"],
            svtypes$color[svtypes$svtype=="INS_SHORT"],
            svtypes$color[svtypes$svtype=="DUP"],
            svtypes$color[!svtypes$svtype %in% c("DEL_SHORT","DEL_SV","INS_SHORT","INS_SV","DUP","DUP_SV")]))
dat.merged <- dat
dat.merged$svtype[dat.merged$svtype %in% c("DEL_SHORT","DEL_SV")] <- "DEL"
dat.merged$svtype[dat.merged$svtype %in% c("INS_SHORT","INS_SV")] <- "INS"
dat.merged$svtype[dat.merged$svtype %in% c("DUP","DUP_SV")] <- "DUP"

######################################
#####VEP consequence distribution plot
######################################
wrapperPlotVepDistrib <- function(){
  if(!"VEP_consequences" %in% colnames(dat)) return(invisible(NULL))

  # Expand: one row per (variant × consequence)
  exp.cols <- intersect(c("svtype","length","AF","REGION"), colnames(dat))
  dat.exp <- do.call(rbind, lapply(seq_len(nrow(dat)), function(i){
    v <- dat$VEP_consequences[i]
    if(is.na(v) || v == "") return(NULL)
    csqs <- strsplit(v, ";", fixed=TRUE)[[1]]
    cbind(dat[rep(i, length(csqs)), exp.cols, drop=FALSE],
          consequence=csqs, stringsAsFactors=FALSE)
  }))
  if(is.null(dat.exp) || nrow(dat.exp) == 0) return(invisible(NULL))

  all.csqs <- names(sort(table(dat.exp$consequence), decreasing=TRUE))
  n.csq <- length(all.csqs)
  col.csq <- setNames(colorRampPalette(c("#440154","#31688E","#35B779","#FDE725"))(n.csq), all.csqs)

  has.region <- "REGION" %in% colnames(dat.exp) && any(!is.na(dat.exp$REGION))

  # Merge svtypes in expanded data
  dat.exp$svtype_m <- dat.exp$svtype
  dat.exp$svtype_m[dat.exp$svtype_m %in% c("DEL_SHORT","DEL_SV")] <- "DEL"
  dat.exp$svtype_m[dat.exp$svtype_m %in% c("INS_SHORT","INS_SV")] <- "INS"
  dat.exp$svtype_m[dat.exp$svtype_m %in% c("DUP","DUP_SV")] <- "DUP"

  # Build consequence × stratum matrix
  makeConseqMat <- function(grp.vals, grp.labs){
    mat <- sapply(grp.labs, function(g){
      idx <- which(grp.vals == g)
      sapply(all.csqs, function(csq) sum(dat.exp$consequence[idx] == csq, na.rm=T))
    })
    matrix(mat, nrow=n.csq, ncol=length(grp.labs), dimnames=list(all.csqs, grp.labs))
  }

  af.cuts  <- c("AC=1","AF<1%","1-10%","10-50%",">50%")
  af.breaks <- c(-Inf, 1.1/(2*nsamp), rare.max.freq, uncommon.max.freq, common.max.freq, 1)
  dat.exp$af_grp   <- as.character(cut(dat.exp$AF, breaks=af.breaks, labels=af.cuts, include.lowest=TRUE))
  dat.exp$size_grp <- as.character(cut(dat.exp$length,
    breaks=c(-Inf, tiny.max.size, small.max.size, medium.max.size, medlarge.max.size, large.max.size, Inf),
    labels=c("<50bp","50-100bp","100bp-500bp","500bp-5kb","5-50kb",">50kb"), include.lowest=TRUE))

  panel.w <- max(3600, n.csq * 120)
  png(paste(OUTDIR, "/main_plots/vep_distributions.png", sep=""),
      res=300, height=2*1800, width=2*panel.w)
  layout(matrix(1:4, nrow=2, byrow=TRUE))

  # Panel 1: All variants — x=consequence, y=count
  counts.all <- sapply(all.csqs, function(csq) sum(dat.exp$consequence == csq))
  par(bty="n", mar=c(12, 4.5, 3, 0.5))
  bp <- barplot(counts.all, names.arg=rep("", n.csq), col=col.csq[all.csqs],
                border=NA, ylim=c(0, max(counts.all)*1.15), yaxt="n")
  axis(2, at=axTicks(2), labels=NA)
  axis(2, at=axTicks(2), tick=F, las=2, cex.axis=0.8, line=-0.4, labels=prettyNum(axTicks(2), big.mark=","))
  mtext(2, text="Count", line=3, cex=0.9)
  mtext(3, text="VEP Consequences (All Variants)", font=2, line=0.5)
  axis(3, at=mean(bp), tick=F, line=-0.9, labels=paste("n=", prettyNum(nrow(dat.exp), big.mark=",")))
  text(x=bp, y=par("usr")[3] - diff(par("usr")[3:4])*0.025,
       labels=all.csqs, srt=45, adj=1, xpd=TRUE, cex=0.7)

  # Panel 2: by variant type — x=svtype, stacked by consequence
  st.labs <- svtypes.merged$svtype
  st.labs <- st.labs[st.labs %in% dat.exp$svtype_m]
  mat.st <- makeConseqMat(dat.exp$svtype_m, st.labs)
  plotStackedBars(mat=mat.st, colors=col.csq[all.csqs], scaled=F,
                  title="VEP Consequences by Variant Type")

  # Panel 3: by AF bucket — x=AF, stacked by consequence
  mat.af <- makeConseqMat(dat.exp$af_grp[!is.na(dat.exp$af_grp)],
                          af.cuts[af.cuts %in% dat.exp$af_grp])
  if(ncol(mat.af) > 0)
    plotStackedBars(mat=mat.af, colors=col.csq[all.csqs], scaled=F,
                    title="VEP Consequences by AF")
  else plot.new()

  # Panel 4: by region (if present) or size
  if(has.region && length(unique(dat.exp$REGION[!is.na(dat.exp$REGION)])) > 1){
    reg.labs <- sort(unique(dat.exp$REGION[!is.na(dat.exp$REGION)]))
    mat.reg <- makeConseqMat(dat.exp$REGION[!is.na(dat.exp$REGION)], reg.labs)
    plotStackedBars(mat=mat.reg, colors=col.csq[all.csqs], scaled=F,
                    title="VEP Consequences by Region")
  } else {
    sz.labs <- c("<50bp","50-100bp","100bp-500bp","500bp-5kb","5-50kb",">50kb")
    mat.sz <- makeConseqMat(dat.exp$size_grp[!is.na(dat.exp$size_grp)],
                            sz.labs[sz.labs %in% dat.exp$size_grp])
    plotStackedBars(mat=mat.sz, colors=col.csq[all.csqs], scaled=F,
                    title="VEP Consequences by Size")
  }
  dev.off()
}


######################################
#####SVAnnotate prediction distribution plot
######################################
wrapperPlotSvAnnotateDistrib <- function(){
  pred.cols <- grep("^PREDICTED_", colnames(dat), value=TRUE)
  if(length(pred.cols) == 0) return(invisible(NULL))

  has.pred <- rowSums(as.data.frame(lapply(pred.cols, function(p) as.integer(dat[[p]]))), na.rm=TRUE) > 0
  if(sum(has.pred) == 0) return(invisible(NULL))

  pred.labels <- sub("^PREDICTED_", "", pred.cols)
  n.pred <- length(pred.cols)
  col.pred <- setNames(colorRampPalette(c("#A50026","#F46D43","#FEE090",
                                          "#74ADD1","#4575B4","#762A83",
                                          "#1B7837","#5AAE61","#F7F7F7"))(n.pred), pred.cols)

  has.region <- "REGION" %in% colnames(dat) && any(!is.na(dat$REGION))

  svtype_m <- dat$svtype
  svtype_m[svtype_m %in% c("DEL_SHORT","DEL_SV")] <- "DEL"
  svtype_m[svtype_m %in% c("INS_SHORT","INS_SV")] <- "INS"
  svtype_m[svtype_m %in% c("DUP","DUP_SV")] <- "DUP"

  af.cuts   <- c("AC=1","AF<1%","1-10%","10-50%",">50%")
  af.breaks <- c(-Inf, 1.1/(2*nsamp), rare.max.freq, uncommon.max.freq, common.max.freq, 1)
  af_grp    <- as.character(cut(dat$AF, breaks=af.breaks, labels=af.cuts, include.lowest=TRUE))

  # Build prediction × stratum matrix
  makePredMat <- function(grp.vals, grp.labs){
    mat <- sapply(grp.labs, function(g){
      idx <- which(grp.vals == g)
      sapply(pred.cols, function(p) sum(as.logical(dat[[p]][idx]), na.rm=T))
    })
    matrix(mat, nrow=n.pred, ncol=length(grp.labs), dimnames=list(pred.labels, grp.labs))
  }

  panel.w <- max(3600, n.pred * 140)
  png(paste(OUTDIR, "/main_plots/svannotate_distributions.png", sep=""),
      res=300, height=2*1800, width=2*panel.w)
  layout(matrix(1:4, nrow=2, byrow=TRUE))

  # Panel 1: All variants — x=PREDICTED_ category, y=count
  counts.all <- sapply(pred.cols, function(p) sum(as.logical(dat[[p]]), na.rm=T))
  par(bty="n", mar=c(12, 4.5, 3, 0.5))
  bp <- barplot(counts.all, names.arg=rep("", n.pred), col=col.pred[pred.cols],
                border=NA, ylim=c(0, max(counts.all)*1.15), yaxt="n")
  axis(2, at=axTicks(2), labels=NA)
  axis(2, at=axTicks(2), tick=F, las=2, cex.axis=0.8, line=-0.4, labels=prettyNum(axTicks(2), big.mark=","))
  mtext(2, text="Count", line=3, cex=0.9)
  mtext(3, text="SVAnnotate Predictions (All Variants)", font=2, line=0.5)
  axis(3, at=mean(bp), tick=F, line=-0.9, labels=paste("n=", prettyNum(sum(has.pred), big.mark=",")))
  text(x=bp, y=par("usr")[3] - diff(par("usr")[3:4])*0.025,
       labels=pred.labels, srt=45, adj=1, xpd=TRUE, cex=0.7)

  # Panel 2: by variant type — x=svtype, stacked by PREDICTED_
  st.labs <- svtypes.merged$svtype
  st.labs <- st.labs[st.labs %in% svtype_m]
  mat.st <- makePredMat(svtype_m, st.labs)
  plotStackedBars(mat=mat.st, colors=col.pred[pred.cols], scaled=F,
                  title="SVAnnotate Predictions by Variant Type")

  # Panel 3: by AF bucket
  af.present <- af.cuts[af.cuts %in% af_grp]
  mat.af <- makePredMat(af_grp[!is.na(af_grp)], af.present)
  if(ncol(mat.af) > 0)
    plotStackedBars(mat=mat.af, colors=col.pred[pred.cols], scaled=F,
                    title="SVAnnotate Predictions by AF")
  else plot.new()

  # Panel 4: by region (if present) or size
  if(has.region && length(unique(dat$REGION[!is.na(dat$REGION)])) > 1){
    reg.labs <- sort(unique(dat$REGION[!is.na(dat$REGION)]))
    mat.reg <- makePredMat(dat$REGION[!is.na(dat$REGION)], reg.labs)
    plotStackedBars(mat=mat.reg, colors=col.pred[pred.cols], scaled=F,
                    title="SVAnnotate Predictions by Region")
  } else {
    sz.labs   <- c("<50bp","50-100bp","100bp-500bp","500bp-5kb","5-50kb",">50kb")
    sz.breaks <- c(-Inf, tiny.max.size, small.max.size, medium.max.size, medlarge.max.size, large.max.size, Inf)
    sz_grp    <- as.character(cut(dat$length, breaks=sz.breaks, labels=sz.labs, include.lowest=TRUE))
    mat.sz    <- makePredMat(sz_grp[!is.na(sz_grp)], sz.labs[sz.labs %in% sz_grp])
    plotStackedBars(mat=mat.sz, colors=col.pred[pred.cols], scaled=F,
                    title="SVAnnotate Predictions by Size")
  }
  dev.off()
}


###Plotting block
#SV counts
wrapperPlotAllCountBars()

#SV sizes
wrapperPlotAllSizeDistribs()

#SV frequencies
wrapperPlotAllFreqDistribs()

#Genotype frequencies
wrapperPlotAllHWDistribs()

#Transition/Transversion distribution
wrapperPlotTiTv()

#QUAL distribution
wrapperPlotQualDistrib()

#NCR distribution
wrapperPlotNcrDistrib()

#gnomAD match rate distribution
wrapperPlotGnomadMatchDistrib()

#VEP consequence distribution
wrapperPlotVepDistrib()

#SVAnnotate prediction distribution
wrapperPlotSvAnnotateDistrib()
