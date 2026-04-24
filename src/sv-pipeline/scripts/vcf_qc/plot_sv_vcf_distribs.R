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
# Order REGION labels in canonical order (US, RM, SD, SR, then others alphabetically)
orderRegions <- function(regions){
  canonical <- c("US","RM","SD","SR")
  c(canonical[canonical %in% regions], sort(regions[!regions %in% canonical]))
}

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
         col=NA, border=NA)
    
    #Add label
    axis(1,at=i-0.5,las=2,line=-0.8,labels=colnames(mat)[i],cex.axis=0.8,tick=F)
  })
}


###################
#####SV count plots
###################
# Format a count as "X.XXXK" or "X.XXXM" with 4 significant figures
formatCount <- function(n){
  if(is.na(n) || !is.finite(n)) return(as.character(n))
  if(abs(n) >= 1e6) return(paste0(signif(n/1e6, 4), "M"))
  if(abs(n) >= 1e3) return(paste0(signif(n/1e3, 4), "K"))
  return(as.character(round(n)))
}

#Plot single set of bars of total count of SV
plotSVCountBars <- function(dat,svtypes,title=NULL,ylab="Count",tr.env.dat=NULL){
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
  
  #Compute TR_ENVELOPED overlay counts per svtype (excluded for TR types)
  tr.types <- c("TR_SNV","TR_INS","TR_DEL")
  te.counts <- if(!is.null(tr.env.dat)){
    setNames(sapply(svtypes$svtype, function(svtype){
      if(svtype %in% tr.types) return(NA_real_)
      as.numeric(length(which(tr.env.dat$svtype==svtype)))
    }), svtypes$svtype)
  } else NULL
  
  #Plot per-svtype information
  sapply(1:nrow(counts),function(i){
    cnt <- counts[i,2]
    bar.col <- counts[i,3]
    if(cnt > 0){
      te.cnt <- if(!is.null(te.counts) && !is.na(te.counts[i]) && te.counts[i] > 0)
        max(te.counts[i], min.y) else NULL
      # Draw stacked bars: bottom = TR_ENVELOPED (striped), top = remainder (solid)
      if(!is.null(te.cnt)){
        # Bottom stripe segment (color + white diagonal hatching)
        rect(xleft=i-0.85,xright=i-0.15, ybottom=min.y,ytop=te.cnt,
             lwd=0.7, col=bar.col, border=bar.col)
        rect(xleft=i-0.85,xright=i-0.15, ybottom=min.y,ytop=te.cnt,
             density=15, angle=45, col="white", border=NA)
        # Top solid segment
        rect(xleft=i-0.85,xright=i-0.15, ybottom=te.cnt,ytop=cnt,
             lwd=0.7, col=bar.col, border=bar.col)
      } else {
        rect(xleft=i-0.85,xright=i-0.15,
             ybottom=min.y,ytop=cnt,
             lwd=0.7,col=bar.col)
      }
      # Count label: main count, then bracketed TR_ENVELOPED count below
      text(x=i-0.5,y=cnt*1.3,col=bar.col,
           labels=formatCount(cnt),cex=0.7)
      if(!is.null(te.cnt)){
        text(x=i-0.5,y=cnt*1.1,col=bar.col,cex=0.6,
             labels=paste0("(",formatCount(te.counts[i]),")"))
      }
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
    regions <- orderRegions(unique(dat$REGION[!is.na(dat$REGION)]))
    region.mat <- as.data.frame(sapply(c("ALL",regions), function(reg){
      plotset <- if(reg=="ALL") dat else dat[which(dat$REGION==reg),]
      sapply(svtypes$svtype, function(svtype) length(which(plotset$svtype==svtype)))
    }))
    colnames(region.mat) <- c("ALL",regions)
  }
  pdf(paste(OUTDIR,"/main_plots/counts_distributions.pdf",sep=""),
      height=7,width=11)
  # Compute TR_ENVELOPED overlay data for leftmost count bar (use count-specific data)
  te.dat.count <- if("TR_ENVELOPED" %in% colnames(dat.count)){
    te.flag <- dat.count$TR_ENVELOPED
    dat.count[!is.na(te.flag) & te.flag %in% c("TRUE","1",TRUE,1,"true","T","t"), ]
  } else NULL
  #Merged
  layout(matrix(c(1,2,3,4,1,5,6,7),byrow=T,nrow=2),
         widths=c(3,2,2,2))
  plotSVCountBars(dat=dat.count,svtypes=svtypes.count,
                  title="Variant Count",tr.env.dat=te.dat.count)
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
  dat <- dat.merged; svtypes <- svtypes.merged
  
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
prepFreqHistogram <- function(freqs, bins.per.order){
  finite.freqs <- freqs[which(is.finite(freqs))]
  if(length(finite.freqs) == 0) return(NULL)

  xlims <- range(finite.freqs, na.rm=T)
  if(xlims[1] == xlims[2]){
    pad <- max(0.1, abs(xlims[1]) * 0.05)
    xlims <- c(xlims[1] - pad, min(0, xlims[2] + pad))
  }

  n.breaks <- bins.per.order * max(1, abs(floor(xlims[1])))
  breaks <- seq(xlims[1], xlims[2], length.out=n.breaks + 1)
  mids <- (breaks[1:(length(breaks)-1)] + breaks[2:length(breaks)]) / 2

  list(xlims=xlims, breaks=breaks, mids=mids)
}

#Plot single AF spectrum
plotFreqDistrib <- function(dat, svtypes,
                            autosomal=F, biallelic=T,
                            title=NULL, lwd.cex=1, legend=F, show.dropped=F){
  #Process freqs & compute range + breaks
  filter.legend <- NULL
  if(autosomal==T){
    dat <- dat[which(dat$chr %in% sex.chroms),]
    filter.legend <- c(filter.legend,"Autosomal variants only")
  }
  n.pre.filter <- nrow(dat)
  if(biallelic==T){
    dat <- dat[which(dat$other_gts==0 & dat$missing_gts<dat$genotyped_samples),]
  }
  freqs <- log10(dat$AF)
  freq.setup <- prepFreqHistogram(freqs, bins.per.order=25)
  if(length(freqs)>0 && !is.null(freq.setup)){
    xlims <- freq.setup$xlims
    breaks <- freq.setup$breaks
    mids <- freq.setup$mids
    
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
  if(show.dropped && biallelic){
    n.dropped <- n.pre.filter - nrow(dat)
    if(n.dropped > 0){
      axis(3,at=mean(par("usr")[1:2]),line=-1.9,tick=F,cex.axis=0.75,
           labels=paste0("(",prettyNum(n.dropped,big.mark=",")," dropped - multi-allelic sites)"))
    }
  }

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
    freq.setup <- prepFreqHistogram(unlist(freqs), bins.per.order=20)
    if(is.null(freq.setup)){
      freqs <- list()
    } else {
      xlims <- freq.setup$xlims
      breaks <- freq.setup$breaks
      mids <- freq.setup$mids
    }
  }
  if(length(freqs) > 0){
    
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
    mtext(3,line=1.5,text=title,font=2)
  }
  
  #Add filter labels
  if(!is.null(filter.legend)){
    legend("topright",bg=NA,bty="n",pch=NA,legend=filter.legend,cex=0.8)
  }
  
  #Add size legend
  if(length(freqs) > 0){
    legend("right",bg="white",bty="n",lwd=3,col=col.pal,cex=0.7,
           legend=gsub("\n","",legend.labs,fixed=T))
  }
}

#Wrapper to plot all AF distributions
wrapperPlotAllFreqDistribs <- function(){
  dat <- dat.merged; svtypes <- svtypes.merged
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
                  legend=T, show.dropped=T)
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
  mat.valid <- mat[is.finite(mat) & !is.na(mat)]
  zlim <- c(min(1.5, if(length(mat.valid)>0) min(mat.valid) else 1.5),
            max(2.5, if(length(mat.valid)>0) max(mat.valid) else 2.5))
  col.pal <- colorRampPalette(c("#440154","#365C8C","#25A584","#FDE725"))(101)
  mat.capped <- pmin(pmax(mat, zlim[1]), zlim[2])
  par(mar=c(5,5,3,28), bty="n")
  image(x=1:ncol(mat), y=1:nrow(mat), z=t(mat.capped), col=col.pal,
        xaxt="n", yaxt="n", xlab="", ylab="", zlim=zlim)
  axis(1, at=1:ncol(mat), labels=colnames(mat), las=2, cex.axis=0.75, tick=F, line=-0.5)
  axis(2, at=1:nrow(mat), labels=rownames(mat), las=2, cex.axis=0.75, tick=F, line=-0.5)
  mtext(1, text="Region", line=3.5, cex=0.85)
  mtext(2, text="AF Bucket", line=3.5, cex=0.85)
  mtext(3, text=title, font=2, line=1)
  # Color bar in right margin using par(xpd=NA) globally to ensure labels are not clipped
  old.xpd <- par(xpd=NA)
  usr <- par("usr")
  xw  <- usr[2] - usr[1]
  cb.xl <- usr[2] + xw * 0.20
  cb.xr <- usr[2] + xw * 0.33
  cb.y  <- seq(usr[3], usr[4], length.out=102)
  rect(xleft=cb.xl, xright=cb.xr, ybottom=cb.y[-length(cb.y)], ytop=cb.y[-1], col=col.pal, border=NA)
  rect(xleft=cb.xl, xright=cb.xr, ybottom=usr[3], ytop=usr[4], col=NA, border="gray50")
  label.vals <- round(seq(zlim[1],zlim[2],length.out=5),2)
  label.y    <- seq(usr[3],usr[4],length.out=5)
  segments(x0=cb.xr, x1=cb.xr + xw*0.05, y0=label.y, y1=label.y, lwd=0.8)
  text(x=cb.xr + xw*0.07, y=label.y, labels=sprintf("%.2f", label.vals), cex=0.75, adj=0)
  par(old.xpd)
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
  regions <- if(has.region) orderRegions(unique(snv.dat$REGION[!is.na(snv.dat$REGION)])) else character(0)

  # Number of columns: AF panel + Size panel + optional heatmap
  n.panels <- 2 + as.integer(has.region && length(regions)>0)
  has.heatmap <- has.region && length(regions)>0
  png(paste(OUTDIR,"/main_plots/ti_tv_distributions.png",sep=""),
      res=300, height=1800, width=n.panels*1800 + if(has.heatmap) 1800 else 0)
  layout(matrix(1:n.panels, nrow=1),
         widths=c(rep(1, n.panels - as.integer(has.heatmap)), if(has.heatmap) 1.6 else 1))

  # Panel 1: Ti:Tv by AF bucket
  titv.af <- sapply(seq_along(af.labels), function(i){
    sub <- snv.dat[which(snv.dat$AF>af.mins[i] & snv.dat$AF<=af.maxs[i]),]
    calcTiTv(sub)
  })
  par(mar=c(7,4.5,3,0.5), bty="n")
  plot(x=c(0,length(af.labels)+1), y=c(0, max(titv.af,na.rm=T)*1.15),
       type="n", xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i")
  abline(h=axTicks(2), col="gray85", lwd=0.5)
  col.af <- colorRampPalette(c("#440154","#365C8C","#25A584","#FDE725"))(length(af.labels))
  rect(xleft=1:length(af.labels)-0.35, xright=1:length(af.labels)+0.35,
       ybottom=0, ytop=titv.af, col=col.af, border=NA)
  axis(1, at=1:length(af.labels), labels=af.labels, las=2, cex.axis=0.8, tick=F, line=0.5)
  axis(2, at=axTicks(2), labels=NA); axis(2, at=axTicks(2), tick=F, las=2, cex.axis=0.8, line=-0.4)
  mtext(1, text="AF Bucket", line=5.5, cex=0.9)
  mtext(2, text="Ti:Tv Ratio", line=2.5, cex=0.9)
  mtext(3, text="Ti:Tv by AF", font=2, line=0.5)

  # Panel 2: Ti:Tv by REGION (or show "No region data.")
  if(has.region && length(regions)>0){
    titv.reg <- sapply(regions, function(reg){
      sub <- snv.dat[which(!is.na(snv.dat$REGION) & snv.dat$REGION==reg),]
      calcTiTv(sub)
    })
    col.reg <- colorRampPalette(c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02"))(length(regions))
    par(mar=c(7,4.5,3,0.5), bty="n")
    plot(x=c(0,length(regions)+1), y=c(0, max(titv.reg,na.rm=T)*1.15),
         type="n", xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i")
    abline(h=axTicks(2), col="gray85", lwd=0.5)
    regions.ordered <- orderRegions(regions)
    titv.reg <- titv.reg[regions.ordered]
    col.reg.ord <- col.reg[match(regions.ordered, regions)]
    rect(xleft=1:length(regions.ordered)-0.35, xright=1:length(regions.ordered)+0.35,
         ybottom=0, ytop=titv.reg, col=col.reg.ord, border=NA)
    axis(1, at=1:length(regions.ordered), labels=regions.ordered, las=2, cex.axis=0.8, tick=F, line=0.5)
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
    text(x=0.5,y=0.5,labels="No region data.")
    mtext(3,text="Ti:Tv by Region",font=2,line=0.5)
  }
  dev.off()
}


######################################
#####Quality score distribution plot
######################################
plotDistribOverlaid <- function(sub.list, sub.labels, sub.colors, main.val=NULL,
                                xlab="QUAL", main.label="All Variants",
                                title=NULL, xlim=NULL, alpha=0.5, log.y=FALSE, log.x=FALSE){
  # Combine all values to get x range
  all.vals <- unlist(lapply(sub.list, function(x) x[!is.na(x) & is.finite(x)]))
  if(length(all.vals)==0){
    par(bty="n",mar=c(4.5,4,3,1))
    plot(x=c(0,1),y=c(0,1),type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
    text(x=0.5,y=0.5,labels="No region data.")
    mtext(3,text=title,font=2,line=1)
    return(invisible(NULL))
  }
  if(is.null(xlim)) xlim <- range(all.vals, na.rm=T)
  if(log.x && xlim[1] <= 0){
    pos.vals <- all.vals[all.vals > 0]
    xlim[1] <- if(length(pos.vals)>0) min(pos.vals) else 0.001
  }
  if(log.x){
    breaks <- 10^seq(log10(xlim[1]), log10(xlim[2]), length.out=51)
  } else {
    breaks <- seq(xlim[1], xlim[2], length.out=51)
  }

  # Compute densities
  dens.list <- lapply(sub.list, function(vals){
    vals <- vals[!is.na(vals) & is.finite(vals) & vals>=xlim[1] & vals<=xlim[2]]
    if(log.x) vals <- vals[vals > 0]
    if(length(vals)<2) return(NULL)
    h <- hist(vals, breaks=breaks, plot=F)
    h$density
  })
  main.dens <- NULL
  if(!is.null(main.val)){
    mv <- main.val[!is.na(main.val) & is.finite(main.val) & main.val>=xlim[1] & main.val<=xlim[2]]
    if(log.x) mv <- mv[mv > 0]
    if(length(mv)>=2) main.dens <- hist(mv, breaks=breaks, plot=F)$density
  }
  all.dens <- unlist(c(lapply(dens.list, function(d) if(!is.null(d)) d else NULL), list(main.dens)))
  all.dens.pos <- all.dens[!is.na(all.dens) & all.dens > 0]
  ymin.lin <- 0
  ymax.top <- max(all.dens, na.rm=T)*1.1
  if(log.y){
    ymin.lin <- if(length(all.dens.pos)>0) max(min(all.dens.pos)*0.1, 1e-9) else 1e-9
    ymax.top <- max(all.dens, na.rm=T)*3
  }
  ylim <- c(ymin.lin, ymax.top)
  mids <- (breaks[-length(breaks)] + breaks[-1])/2

  par(bty="n", mar=c(4.5,4,3,1))
  log.arg <- paste(c(if(log.x) "x" else "", if(log.y) "y" else ""), collapse="")
  plot(x=xlim, y=ylim, type="n", xaxt="n", yaxt="n", xlab="", ylab="", yaxs="i", log=log.arg)
  if(log.x){
    x.grid <- 10^seq(floor(log10(xlim[1])), ceiling(log10(xlim[2])))
    abline(v=x.grid, col="gray85", lwd=0.5)
    x.ticks <- x.grid
    axis(1, at=x.ticks, labels=NA)
    axis(1, at=x.ticks, tick=F, line=-0.4, cex.axis=0.8,
         labels=paste0(x.ticks*100, "%"))
  } else {
    abline(v=pretty(xlim), col="gray85", lwd=0.5)
    axis(1, at=pretty(xlim), labels=NA); axis(1, at=pretty(xlim), tick=F, line=-0.4, cex.axis=0.8)
  }
  if(log.y){
    log.ticks <- 10^seq(floor(log10(ymin.lin)), ceiling(log10(ymax.top)))
    axis(2, at=log.ticks, labels=NA)
    axis(2, at=log.ticks, tick=F, las=2, cex.axis=0.8, line=-0.4)
  }else{
    axis(2, at=axTicks(2), labels=NA); axis(2, at=axTicks(2), tick=F, las=2, cex.axis=0.8, line=-0.4)
  }
  mtext(1, text=xlab, line=3, cex=0.9); mtext(2, text="Density", line=2.5, cex=0.9)
  mtext(3, text=title, font=2, line=1)

  # Draw main (all variants) density first as dark line
  if(!is.null(main.dens)){
    poly.y <- if(log.y) pmax(c(0,main.dens,0), ymin.lin) else c(0,main.dens,0)
    polygon(x=c(mids[1],mids,mids[length(mids)]), y=poly.y,
            col=adjustcolor("gray25",alpha.f=0.3), border="gray25", lwd=1.5)
  }
  # Draw stratified overlays
  for(i in seq_along(sub.list)){
    if(is.null(dens.list[[i]])) next
    poly.y <- if(log.y) pmax(c(0,dens.list[[i]],0), ymin.lin) else c(0,dens.list[[i]],0)
    polygon(x=c(mids[1],mids,mids[length(mids)]), y=poly.y,
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
  regions <- if(has.region) orderRegions(unique(dat$REGION[!is.na(dat$REGION)])) else character(0)
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

  n.panels <- 5
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
  abline(v=pretty(xlim), col="gray85", lwd=0.5)
  polygon(c(h$mids[1],h$mids,h$mids[length(h$mids)]), c(0,h$density,0), col="#4393C3", border="#2166AC")
  axis(1,at=pretty(xlim),labels=NA); axis(1,at=pretty(xlim),tick=F,line=-0.4,cex.axis=0.8)
  axis(2,at=axTicks(2),labels=NA); axis(2,at=axTicks(2),tick=F,las=2,cex.axis=0.8,line=-0.4)
  mtext(1,text="QUAL",line=3,cex=0.9); mtext(2,text="Density",line=2.5,cex=0.9)
  mtext(3,text="QUAL Distribution",font=2,line=1)
  axis(3,at=mean(xlim),tick=F,line=-0.9,labels=paste("n=",prettyNum(length(q.all),big.mark=","),sep=""))
  n.qual.dropped <- nrow(dat) - length(dat$QUAL[!is.na(dat$QUAL) & is.finite(dat$QUAL)])
  if(n.qual.dropped > 0){
    axis(3,at=mean(xlim),tick=F,line=-1.9,cex.axis=0.75,
         labels=paste0("(",prettyNum(n.qual.dropped,big.mark=",")," dropped - missing QUAL score)"))
  }

  # Panel 2: By svtype
  type.vals <- lapply(svtypes$svtype, function(st) pmin(dat$QUAL[!is.na(dat$QUAL) & dat$svtype==st], 99))
  plotDistribOverlaid(sub.list=type.vals, sub.labels=svtypes$svtype, sub.colors=col.svtype,
                      main.val=q.all, xlab="QUAL", main.label="All",
                      title="QUAL by Type", xlim=xlim)

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

  # Panel 5: By region
  reg.vals <- if(has.region && length(regions)>0) lapply(regions, function(r) pmin(dat$QUAL[!is.na(dat$REGION) & dat$REGION==r & !is.na(dat$QUAL)], 99)) else list()
  plotDistribOverlaid(sub.list=reg.vals, sub.labels=regions, sub.colors=col.reg,
                      main.val=q.all, xlab="QUAL", main.label="All",
                      title="QUAL by Region", xlim=xlim)
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
  regions <- if(has.region) orderRegions(unique(dat$REGION[!is.na(dat$REGION)])) else character(0)
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

  # Panel 1: All variants - log x-axis
  par(bty="n", mar=c(4.5,4,3,1))
  n.all <- dat$NCR[!is.na(dat$NCR)]
  ncr.floor <- 0.001
  n.clipped <- pmax(n.all[n.all > 0], ncr.floor)
  breaks.log <- 10^seq(log10(ncr.floor), 0, length.out=51)
  h <- hist(n.clipped, breaks=breaks.log, plot=F)
  plot(x=c(ncr.floor, 1), y=c(0, max(h$density)*1.15), type="n", xaxt="n", yaxt="n",
       xlab="", ylab="", yaxs="i", log="x")
  abline(v=10^seq(-3, 0), col="gray85", lwd=0.5)
  polygon(c(h$mids[1],h$mids,h$mids[length(h$mids)]), c(0,h$density,0),
          col="#4393C3", border="#2166AC")
  x.ticks.ncr <- 10^seq(-3, 0)
  axis(1,at=x.ticks.ncr,labels=NA)
  axis(1,at=x.ticks.ncr,tick=F,line=-0.4,cex.axis=0.8,labels=paste0(x.ticks.ncr*100,"%"))
  axis(2,at=axTicks(2),labels=NA); axis(2,at=axTicks(2),tick=F,las=2,cex.axis=0.8,line=-0.4)
  mtext(1,text="NCR",line=3,cex=0.9); mtext(2,text="Density",line=2.5,cex=0.9)
  mtext(3,text="NCR Distribution",font=2,line=1)
  mtext(3,text=paste("n=",prettyNum(length(n.all),big.mark=","),sep=""),line=-0.1,cex=0.8)
  n.ncr.dropped <- nrow(dat) - length(n.all)
  if(n.ncr.dropped > 0){
    mtext(3,text=paste0("(",prettyNum(n.ncr.dropped,big.mark=",")," dropped - missing NCR value)"),line=-1.1,cex=0.7)
  }

  # Panel 2: By region
  reg.vals <- if(has.region && length(regions)>0) lapply(regions, function(r) dat$NCR[!is.na(dat$REGION) & dat$REGION==r & !is.na(dat$NCR)]) else list()
  plotDistribOverlaid(sub.list=reg.vals, sub.labels=regions, sub.colors=col.reg,
                      xlab="NCR", title="NCR by Region", xlim=c(0.001,1), log.x=TRUE)

  # Panel 3: By size bucket
  size.vals <- lapply(seq_along(size.labels), function(i) dat$NCR[!is.na(dat$NCR) & dat$length>=size.mins[i] & dat$length<=size.maxs[i]])
  plotDistribOverlaid(sub.list=size.vals, sub.labels=size.labels, sub.colors=col.size,
                      xlab="NCR", title="NCR by Size", xlim=c(0.001,1), log.x=TRUE)

  # Panel 4: By AF bucket
  af.vals <- lapply(seq_along(af.labels), function(i) dat$NCR[!is.na(dat$NCR) & dat$AF>af.mins[i] & dat$AF<=af.maxs[i]])
  plotDistribOverlaid(sub.list=af.vals, sub.labels=af.labels, sub.colors=col.af,
                      xlab="NCR", title="NCR by AF", xlim=c(0.001,1), log.x=TRUE)
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
  regions <- if(has.region) orderRegions(unique(dat$REGION[!is.na(dat$REGION)])) else character(0)
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
    abline(h=seq(0.2, 1.0, by=0.2), col="gray85")
    abline(h=seq(0.1, 0.9, by=0.2), col="gray92", lwd=0.5)
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
             bty="n",bg=NA,cex=0.35,xpd=TRUE)
    }else{
      legend("topright",pch=19,col=c("#4DAC26","#81F850","#AC26A1"),pt.cex=2,
             legend=c(paste(round(100*(n.pass/nrow(HW.mat)),2),"%",sep=""),
                      paste(round(100*(n.nom/nrow(HW.mat)),2),"%",sep=""),
                      paste(round(100*(n.bonf/nrow(HW.mat)),2),"%",sep="")),
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
    legend("right",bg=NA,bty="n",pch=NA,cex=0.35,
           legend=c("Autosomal variants only"))
  }
}
#Correlation of carrier frequency & AF
plotAlleleCarrierCorrelation <- function(dat,autosomal=T,biallelic=T,
                                         title="Carrier Frequency vs. AF"){
  #Process freqs
  filter.legend <- NULL
  if(autosomal==T){
    dat <- dat[which(dat$chr %in% sex.chroms),]
    filter.legend <- c(filter.legend,"Autosomal variants only")
  }
  if(biallelic==T){
    dat <- dat[which(dat$other_gts==0 & dat$missing_gts<dat$genotyped_samples),]
  }
  CF <- sort(dat$carrierFreq)
  AF <- dat$AF[order(dat$carrierFreq)]
  keep <- which(!is.na(CF) & !is.na(AF))
  CF <- CF[keep]
  AF <- AF[keep]

  # Subsample points for plotting to avoid over-dense scatter
  max.plot.pts <- 50000
  if(length(CF) > max.plot.pts){
    sub.idx <- sort(sample(seq_along(CF), max.plot.pts))
    CF.plot <- CF[sub.idx]; AF.plot <- AF[sub.idx]
  } else {
    CF.plot <- CF; AF.plot <- AF
  }

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
  
  #Plot points & rolling mean (subsampled for display; full series used for rolling mean)
  points(x=CF.plot,y=AF.plot,cex=0.2)
  points(x=rollmean(CF,k=100),
         y=rollmean(AF,k=100),
         type="l",col="red",lwd=1.25)
  
  #Add legend
  legend("bottomright",bg=NA,legend=c("Site","Rolling Mean"),cex=0.8,bty="n",
         lwd=c(NA,2),col=c("black","red"),pch=c(1,NA),pt.cex=c(0.4,NA))
  
  #Add filter labels
  if(!is.null(filter.legend)){
    legend("topleft",bg=NA,bty="n",pch=NA,legend=filter.legend,cex=0.55)
  }
}

#Wrapper to plot all HW distributions
wrapperPlotAllHWDistribs <- function(){
  # Supporting individual plots by size
  supp.dir <- paste(OUTDIR,"/supporting_plots/vcf_summary_plots/",sep="")
  for(item in list(
    list(f="gt_distribution.all.png",     d=dat,                                                                           t="Genotype Distribution"),
    list(f="gt_distribution.tiny.png",    d=dat[which(dat$length<tiny.max.size),],                                          t="Genotype Distribution (< 50bp)"),
    list(f="gt_distribution.small.png",   d=dat[which(dat$length>=tiny.max.size & dat$length<small.max.size),],              t="Genotype Distribution (50 - 100bp)"),
    list(f="gt_distribution.medium.png",  d=dat[which(dat$length>=small.max.size & dat$length<medium.max.size),],            t="Genotype Distribution (100bp - 500bp)"),
    list(f="gt_distribution.medlarge.png",d=dat[which(dat$length>=medium.max.size & dat$length<medlarge.max.size),],         t="Genotype Distribution (500bp - 5kb)"),
    list(f="gt_distribution.large.png",   d=dat[which(dat$length>=medlarge.max.size & dat$length<large.max.size),],          t="Genotype Distribution (5 - 50kb)"),
    list(f="gt_distribution.huge.png",    d=dat[which(dat$length>=large.max.size & dat$length<huge.max.size),],              t="Genotype Distribution (> 50kb)")
  )){
    png(paste(supp.dir,item$f,sep=""), res=300, height=1800, width=1800)
    plotHWSingle(dat=item$d, svtypes=svtypes, title=item$t)
    dev.off()
  }
  png(paste(supp.dir,"carrier_freq_vs_allele_freq.png",sep=""), res=300, height=1200, width=1200)
  plotAlleleCarrierCorrelation(dat=dat)
  dev.off()

  cell.px <- 1200

  # gt_overall_distributions.png: overall HWE + carrier freq vs AF
  png(paste(OUTDIR,"main_plots/gt_overall_distributions.png",sep=""),
      res=300, height=cell.px, width=2*cell.px)
  layout(matrix(1:2, nrow=1))
  plotHWSingle(dat=dat, svtypes=svtypes, title="Genotype Distribution")
  plotAlleleCarrierCorrelation(dat=dat)
  dev.off()

  # gt_type_distributions.png: specific 2x5 grid
  # Row 1: SNV, INS_SHORT, DEL_SHORT, DUP_SHORT, TR_INS
  # Row 2: TR_SNV, INS_SV, DEL_SV, DUP_SV, TR_DEL
  type.grid <- c("SNV","INS_SHORT","DEL_SHORT","DUP_SHORT","TR_INS",
                  "TR_SNV","INS_SV","DEL_SV","DUP_SV","TR_DEL")
  type.grid <- type.grid[type.grid %in% svtypes$svtype]
  n.type.cols <- ceiling(length(type.grid) / 2)
  png(paste(OUTDIR,"main_plots/gt_type_distributions.png",sep=""),
      res=300, height=2*cell.px, width=n.type.cols*cell.px)
  mat.type <- matrix(seq_len(length(type.grid)), nrow=2, byrow=TRUE)
  if(ncol(mat.type) * 2 > length(type.grid)){
    mat.type[length(type.grid)+1] <- length(type.grid)+1
  }
  layout(mat.type)
  sapply(seq_along(type.grid), function(i){
    st <- type.grid[i]
    plotHWSingle(dat=dat[which(dat$svtype==st),], svtypes=svtypes,
                 title=st, full.legend=F, lab.cex=0.75)
  })
  if(length(type.grid) %% 2 != 0) plot.new()
  dev.off()

  # gt_distributions_size.png: per-size-bucket HWE grid (2 rows)
  sz.labels.hwe <- c("<50bp","50-100bp","100bp-500bp","500bp-5kb","5-50kb",">50kb")
  sz.mins.hwe   <- c(0, tiny.max.size, small.max.size, medium.max.size, medlarge.max.size, large.max.size)
  sz.maxs.hwe   <- c(tiny.max.size, small.max.size, medium.max.size, medlarge.max.size, large.max.size, huge.max.size)
  n.sz      <- length(sz.labels.hwe)
  n.sz.cols <- ceiling(n.sz / 2)
  mat.sz    <- matrix(NA, nrow=2, ncol=n.sz.cols)
  for(i in seq_len(n.sz)) mat.sz[((i-1)%%2)+1, ceiling(i/2)] <- i
  mat.sz[is.na(mat.sz)] <- n.sz+1
  png(paste(OUTDIR,"main_plots/gt_size_distributions.png",sep=""),
      res=300, height=2*cell.px, width=n.sz.cols*cell.px)
  layout(mat.sz)
  sapply(seq_len(n.sz), function(i){
    sub <- dat[which(dat$length >= sz.mins.hwe[i] & dat$length < sz.maxs.hwe[i]),]
    plotHWSingle(dat=sub, svtypes=svtypes, title=sz.labels.hwe[i], full.legend=F, lab.cex=0.75)
  })
  if(n.sz %% 2 != 0) plot.new()
  dev.off()

  # gt_distributions_region.png: always produced; per-region HWE grid if REGION column present
  has.region <- "REGION" %in% colnames(dat) && any(!is.na(dat$REGION))
  if(has.region){
    regions <- orderRegions(unique(dat$REGION[!is.na(dat$REGION)]))
    n.reg      <- length(regions)
    n.reg.cols <- ceiling(n.reg / 2)
    mat.reg    <- matrix(NA, nrow=2, ncol=n.reg.cols)
    for(i in seq_len(n.reg)) mat.reg[((i-1)%%2)+1, ceiling(i/2)] <- i
    mat.reg[is.na(mat.reg)] <- n.reg+1
    png(paste(OUTDIR,"main_plots/gt_region_distributions.png",sep=""),
        res=300, height=2*cell.px, width=n.reg.cols*cell.px)
    layout(mat.reg)
    sapply(seq_len(n.reg), function(i){
      sub <- dat[which(!is.na(dat$REGION) & dat$REGION==regions[i]),]
      plotHWSingle(dat=sub, svtypes=svtypes, title=regions[i], full.legend=F, lab.cex=0.75)
    })
    if(n.reg %% 2 != 0) plot.new()
    dev.off()
  } else {
    png(paste(OUTDIR,"main_plots/gt_region_distributions.png",sep=""),
        res=300, height=cell.px, width=cell.px)
    plot(x=c(0,1),y=c(0,1),type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
    text(x=0.5,y=0.5,labels="No region data.",cex=1.5)
    dev.off()
  }
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
# Enforce canonical class ordering
.svtype.order <- c("SNV","INS_SHORT","DEL_SHORT","DUP_SHORT","INS_SV","DEL_SV","DUP_SV","TR_SNV","TR_INS","TR_DEL")
.order.idx <- c(match(.svtype.order[.svtype.order %in% svtypes$svtype], svtypes$svtype),
                which(!svtypes$svtype %in% .svtype.order))
svtypes <- svtypes[.order.idx, ]

# Build count-specific svtypes (TR/VNTR for unique locus counting)
# These use deduplicated original TRV IDs classified by MOTIFS length
svtypes.count <- rbind(
  svtypes[!svtypes$svtype %in% c("TR_SNV","TR_INS","TR_DEL"),],
  data.frame(svtype=c("TR","VNTR"), color=c("#FA931E","#B56A00"), stringsAsFactors=FALSE))
.count.order <- c("SNV","INS_SHORT","DEL_SHORT","DUP_SHORT","INS_SV","DEL_SV","DUP_SV","TR","VNTR")
.count.idx <- c(match(.count.order[.count.order %in% svtypes.count$svtype], svtypes.count$svtype),
                which(!svtypes.count$svtype %in% .count.order))
svtypes.count <- svtypes.count[.count.idx, ]

# Build count-specific data (deduplicate TR biallelics back to TR/VNTR loci)
tr.biallelic.types <- c("TR_SNV","TR_INS","TR_DEL")
dat.count <- dat[!dat$svtype %in% tr.biallelic.types,]
if(any(dat$svtype %in% tr.biallelic.types)){
  tr.dat.tmp <- dat[dat$svtype %in% tr.biallelic.types,]
  tr.dat.tmp$base_vid <- sub("_[0-9]+$", "", tr.dat.tmp$VID)
  first.idx <- !duplicated(tr.dat.tmp$base_vid)
  tr.dedup <- tr.dat.tmp[first.idx,]
  if("max_motif_length" %in% colnames(tr.dedup)){
    tr.dedup$svtype <- ifelse(!is.na(tr.dedup$max_motif_length) & tr.dedup$max_motif_length > 6, "VNTR", "TR")
  } else {
    tr.dedup$svtype <- "TR"
  }
  tr.dedup$base_vid <- NULL
  tr.dat.tmp$base_vid <- NULL
  dat.count <- rbind(dat.count, tr.dedup)
}

# Merged svtypes and dat (DEL_SHORT+DEL_SV->DEL, INS_SHORT+INS_SV->INS, DUP_SHORT+DUP_SV->DUP)
.mg.del  <- c("DEL_SHORT","DEL_SV")
.mg.ins  <- c("INS_SHORT","INS_SV")
.mg.dup  <- c("DUP_SHORT","DUP_SV")
.del.col <- svtypes$color[match(intersect(c("DEL_SHORT","DEL"), svtypes$svtype)[1], svtypes$svtype)]
.ins.col <- svtypes$color[match(intersect(c("INS_SHORT","INS"), svtypes$svtype)[1], svtypes$svtype)]
.dup.col <- svtypes$color[match(intersect(c("DUP_SHORT","DUP"), svtypes$svtype)[1], svtypes$svtype)]
.other.types <- svtypes$svtype[!svtypes$svtype %in% c(.mg.del,.mg.ins,.mg.dup)]
svtypes.merged <- data.frame(
  svtype = c("DEL","INS","DUP", .other.types),
  color = c(.del.col, .ins.col, .dup.col,
            svtypes$color[match(.other.types, svtypes$svtype)]),
  stringsAsFactors=FALSE)
.merged.order <- c("SNV","INS","DEL","DUP","TR_SNV","TR_INS","TR_DEL")
.merged.idx <- c(match(.merged.order[.merged.order %in% svtypes.merged$svtype], svtypes.merged$svtype),
                 which(!svtypes.merged$svtype %in% .merged.order))
svtypes.merged <- svtypes.merged[.merged.idx, ]
dat.merged <- dat
dat.merged$svtype[dat.merged$svtype %in% .mg.del] <- "DEL"
dat.merged$svtype[dat.merged$svtype %in% .mg.ins] <- "INS"
dat.merged$svtype[dat.merged$svtype %in% .mg.dup] <- "DUP"

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

  # Filter to non-zero consequences, cap at 9 + 'All other'
  all.csqs.raw <- names(sort(table(dat.exp$consequence), decreasing=TRUE))
  all.csqs.raw <- all.csqs.raw[sapply(all.csqs.raw, function(csq) sum(dat.exp$consequence==csq))>0]
  max.show <- 9
  other.col <- "#AAAAAA"
  if(length(all.csqs.raw) > max.show){
    top.csqs <- all.csqs.raw[seq_len(max.show)]
    dat.exp$consequence[!dat.exp$consequence %in% top.csqs] <- "All other"
    all.csqs <- c(top.csqs, "All other")
    col.csq.base <- setNames(colorRampPalette(c("#440154","#31688E","#35B779","#FDE725"))(max.show), top.csqs)
    col.csq <- c(col.csq.base, "All other"=other.col)
  } else {
    all.csqs <- all.csqs.raw
    col.csq <- setNames(colorRampPalette(c("#440154","#31688E","#35B779","#FDE725"))(length(all.csqs)), all.csqs)
  }
  n.csq <- length(all.csqs)

  has.region <- "REGION" %in% colnames(dat.exp) && any(!is.na(dat.exp$REGION))

  # Merge svtypes in expanded data
  dat.exp$svtype_m <- dat.exp$svtype
  dat.exp$svtype_m[dat.exp$svtype_m %in% c("DEL_SHORT","DEL_SV")] <- "DEL"
  dat.exp$svtype_m[dat.exp$svtype_m %in% c("INS_SHORT","INS_SV")] <- "INS"
  dat.exp$svtype_m[dat.exp$svtype_m %in% c("DUP_SHORT","DUP_SV")] <- "DUP"
  dat.exp$svtype_m[dat.exp$svtype_m %in% c("TR_SNV","TR_INS","TR_DEL")] <- "TR"

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

  panel.w <- max(3600, n.csq * 150)
  png(paste(OUTDIR, "/main_plots/vep_distributions.png", sep=""),
      res=300, height=2*1800, width=2*panel.w)
  layout(matrix(1:4, nrow=2, byrow=TRUE))

  # Panel 1: All variants — x=consequence, y=count
  counts.all <- sapply(all.csqs, function(csq) sum(dat.exp$consequence == csq))
  par(bty="n", mar=c(12, 4.5, 4.5, 0.5))
  bp <- barplot(counts.all, names.arg=rep("", n.csq), col=col.csq[all.csqs],
                border=NA, ylim=c(0, max(counts.all)*1.15), yaxt="n")
  axis(2, at=axTicks(2), labels=NA)
  axis(2, at=axTicks(2), tick=F, las=2, cex.axis=0.8, line=-0.4, labels=prettyNum(axTicks(2), big.mark=","))
  mtext(2, text="Count", line=3, cex=0.9)
  mtext(3, text="VEP Consequences (All Variants)", font=2, line=2.0)
  mtext(3, text=paste("n=", prettyNum(nrow(dat.exp), big.mark=","), sep=""), line=0.6, cex=0.75)
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
    reg.labs <- orderRegions(unique(dat.exp$REGION[!is.na(dat.exp$REGION)]))
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

  # Filter to non-zero PREDICTED_ cols, cap at 9 + 'All other'
  counts.raw <- sapply(pred.cols, function(p) sum(as.logical(dat[[p]]), na.rm=TRUE))
  nonzero.pred <- pred.cols[counts.raw > 0]
  if(length(nonzero.pred) == 0) return(invisible(NULL))
  max.show.pred <- 9
  other.pred.col <- "#AAAAAA"
  if(length(nonzero.pred) > max.show.pred){
    top.preds <- nonzero.pred[seq_len(max.show.pred)]
    other.preds <- nonzero.pred[-(seq_len(max.show.pred))]
    # Add an "All other" boolean column: TRUE if any of other.preds is TRUE
    dat[["PRED_OTHER"]] <- rowSums(as.data.frame(lapply(other.preds, function(p) as.integer(dat[[p]]))), na.rm=TRUE) > 0
    pred.cols.use <- c(top.preds, "PRED_OTHER")
    pred.labels.use <- c(sub("^PREDICTED_","",top.preds), "All other")
    col.pred.base <- setNames(colorRampPalette(c("#A50026","#F46D43","#FEE090",
                                                 "#74ADD1","#4575B4","#762A83",
                                                 "#1B7837","#5AAE61","#F7F7F7"))(max.show.pred), top.preds)
    col.pred <- c(col.pred.base, "PRED_OTHER"=other.pred.col)
  } else {
    pred.cols.use <- nonzero.pred
    pred.labels.use <- sub("^PREDICTED_","",nonzero.pred)
    col.pred <- setNames(colorRampPalette(c("#A50026","#F46D43","#FEE090",
                                            "#74ADD1","#4575B4","#762A83",
                                            "#1B7837","#5AAE61","#F7F7F7"))(length(nonzero.pred)), nonzero.pred)
  }
  pred.cols <- pred.cols.use
  pred.labels <- pred.labels.use
  n.pred <- length(pred.cols)

  has.region <- "REGION" %in% colnames(dat) && any(!is.na(dat$REGION))

  svtype_m <- dat$svtype
  svtype_m[svtype_m %in% c("DEL_SHORT","DEL_SV")] <- "DEL"
  svtype_m[svtype_m %in% c("INS_SHORT","INS_SV")] <- "INS"
  svtype_m[svtype_m %in% c("DUP_SHORT","DUP_SV")] <- "DUP"
  svtype_m[svtype_m %in% c("TR_SNV","TR_INS","TR_DEL")] <- "TR"

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

  panel.w <- max(3600, n.pred * 150)
  png(paste(OUTDIR, "/main_plots/svannotate_distributions.png", sep=""),
      res=300, height=2*1800, width=2*panel.w)
  layout(matrix(1:4, nrow=2, byrow=TRUE))

  # Panel 1: All variants — x=PREDICTED_ category, y=count
  counts.all <- sapply(pred.cols, function(p) sum(as.logical(dat[[p]]), na.rm=T))
  par(bty="n", mar=c(12, 4.5, 4.5, 0.5))
  bp <- barplot(counts.all, names.arg=rep("", n.pred), col=col.pred[pred.cols],
                border=NA, ylim=c(0, max(counts.all)*1.15), yaxt="n")
  axis(2, at=axTicks(2), labels=NA)
  axis(2, at=axTicks(2), tick=F, las=2, cex.axis=0.8, line=-0.4, labels=prettyNum(axTicks(2), big.mark=","))
  mtext(2, text="Count", line=3, cex=0.9)
  mtext(3, text="SVAnnotate Consequences (All Variants)", font=2, line=2.0)
  mtext(3, text=paste("n=", prettyNum(sum(has.pred), big.mark=","), sep=""), line=0.6, cex=0.75)
  text(x=bp, y=par("usr")[3] - diff(par("usr")[3:4])*0.025,
       labels=pred.labels, srt=45, adj=1, xpd=TRUE, cex=0.7)

  # Panel 2: by variant type — x=svtype, stacked by PREDICTED_
  st.labs <- svtypes.merged$svtype
  st.labs <- st.labs[st.labs %in% svtype_m]
  mat.st <- makePredMat(svtype_m, st.labs)
  plotStackedBars(mat=mat.st, colors=col.pred[pred.cols], scaled=F,
                  title="SVAnnotate Consequences by Variant Type")

  # Panel 3: by AF bucket
  af.present <- af.cuts[af.cuts %in% af_grp]
  mat.af <- makePredMat(af_grp[!is.na(af_grp)], af.present)
  if(ncol(mat.af) > 0)
    plotStackedBars(mat=mat.af, colors=col.pred[pred.cols], scaled=F,
                    title="SVAnnotate Consequences by AF")
  else plot.new()

  # Panel 4: by region (if present) or size
  if(has.region && length(unique(dat$REGION[!is.na(dat$REGION)])) > 1){
    reg.labs <- orderRegions(unique(dat$REGION[!is.na(dat$REGION)]))
    mat.reg <- makePredMat(dat$REGION[!is.na(dat$REGION)], reg.labs)
    plotStackedBars(mat=mat.reg, colors=col.pred[pred.cols], scaled=F,
                    title="SVAnnotate Consequences by Region")
  } else {
    sz.labs   <- c("<50bp","50-100bp","100bp-500bp","500bp-5kb","5-50kb",">50kb")
    sz.breaks <- c(-Inf, tiny.max.size, small.max.size, medium.max.size, medlarge.max.size, large.max.size, Inf)
    sz_grp    <- as.character(cut(dat$length, breaks=sz.breaks, labels=sz.labs, include.lowest=TRUE))
    mat.sz    <- makePredMat(sz_grp[!is.na(sz_grp)], sz.labs[sz.labs %in% sz_grp])
    plotStackedBars(mat=mat.sz, colors=col.pred[pred.cols], scaled=F,
                    title="SVAnnotate Consequences by Size")
  }
  dev.off()
}


######################################
#####TRV-specific distribution plots
######################################
# TRV size bins (distinct from main SV size bins, tuned for TR allele length range)
trv.size.breaks <- c(-Inf, 5, 10, 50, 100, Inf)
trv.size.labels <- c("<6bp","6-10bp","11-50bp","51-100bp",">100bp")
trv.motif.breaks <- c(-Inf, 2, 4, 10, Inf)
trv.motif.labels <- c("1-2bp","3-4bp","5-10bp",">10bp")

# Renders 4 panels: All, by TRV size, by REGION, by motif length.
# PNG must be opened before calling and closed afterwards.
plotTrvDistribPanels <- function(trv.dat, vals, xlab, xlim=NULL,
                                  title.all, title.size, title.region, title.motif,
                                  drop.text=NULL, log.y=FALSE){
  has.region <- "REGION" %in% colnames(trv.dat) && any(!is.na(trv.dat$REGION))
  regions    <- if(has.region) orderRegions(unique(trv.dat$REGION[!is.na(trv.dat$REGION)])) else character(0)
  has.motif  <- "max_motif_length" %in% colnames(trv.dat) && any(!is.na(trv.dat$max_motif_length))

  col.size <- colorRampPalette(c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02"))(length(trv.size.labels))
  col.reg  <- colorRampPalette(c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02"))(length(regions))
  col.mot  <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3")[seq_along(trv.motif.labels)]

  if(is.null(xlim)){
    v.f <- vals[!is.na(vals) & is.finite(vals)]
    xlim <- if(length(v.f) > 0) range(v.f) else c(0,1)
  }
  breaks.all <- seq(xlim[1], xlim[2], length.out=51)

  # Panel 1: All TRV
  v.all <- vals[!is.na(vals) & is.finite(vals) & vals >= xlim[1] & vals <= xlim[2]]
  par(bty="n", mar=c(4.5,4,3,1))
  if(length(v.all) >= 2){
    h <- hist(v.all, breaks=breaks.all, plot=F)
    dens.pos <- h$density[h$density > 0]
    ymax <- max(h$density)*1.15
    ymin <- if(log.y && length(dens.pos)>0) max(min(dens.pos)*0.1, 1e-9) else 0
    log.arg <- if(log.y) "y" else ""
    plot(x=xlim, y=c(ymin, if(log.y) ymax*3 else ymax), type="n", xaxt="n", yaxt="n",
         xlab="", ylab="", yaxs="i", log=log.arg)
    abline(v=pretty(xlim), col="gray85", lwd=0.5)
    poly.y <- if(log.y) pmax(c(0,h$density,0), ymin) else c(0,h$density,0)
    polygon(c(h$mids[1],h$mids,h$mids[length(h$mids)]), poly.y,
            col="#4393C3", border="#2166AC")
    axis(1,at=pretty(xlim),labels=NA); axis(1,at=pretty(xlim),tick=F,line=-0.4,cex.axis=0.8)
    if(log.y){
      log.ticks <- 10^seq(floor(log10(ymin)), ceiling(log10(ymax)))
      axis(2,at=log.ticks,labels=NA); axis(2,at=log.ticks,tick=F,las=2,cex.axis=0.8,line=-0.4)
    } else {
      axis(2,at=axTicks(2),labels=NA); axis(2,at=axTicks(2),tick=F,las=2,cex.axis=0.8,line=-0.4)
    }
    mtext(1,text=xlab,line=3,cex=0.9); mtext(2,text="Density",line=2.5,cex=0.9)
    mtext(3,text=title.all,font=2,line=1)
    axis(3,at=mean(xlim),tick=F,line=-0.9,labels=paste("n=",prettyNum(length(v.all),big.mark=","),sep=""))
    if(!is.null(drop.text)){
      axis(3,at=mean(xlim),tick=F,line=-1.9,cex.axis=0.75,labels=drop.text)
    }
  } else {
    plot.new(); mtext(3,text=title.all,font=2,line=1)
  }

  # Panel 2: By TRV size bucket
  size.vals <- lapply(seq_along(trv.size.labels), function(i){
    lo <- trv.size.breaks[i]; hi <- trv.size.breaks[i+1]
    v <- vals[!is.na(vals) & is.finite(vals) & !is.na(trv.dat$length) & trv.dat$length > lo & trv.dat$length <= hi]
    v[v >= xlim[1] & v <= xlim[2]]
  })
  plotDistribOverlaid(sub.list=size.vals, sub.labels=trv.size.labels, sub.colors=col.size,
                      xlab=xlab, title=title.size, xlim=xlim, log.y=log.y)

  # Panel 3: By motif length bucket
  if(has.motif){
    mot.vals <- lapply(seq_along(trv.motif.labels), function(i){
      lo <- trv.motif.breaks[i]; hi <- trv.motif.breaks[i+1]
      ml <- trv.dat$max_motif_length
      v <- vals[!is.na(vals) & is.finite(vals) & !is.na(ml) & ml > lo & ml <= hi]
      v[v >= xlim[1] & v <= xlim[2]]
    })
    plotDistribOverlaid(sub.list=mot.vals, sub.labels=trv.motif.labels, sub.colors=col.mot,
                        xlab=xlab, title=title.motif, xlim=xlim, log.y=log.y)
  } else {
    par(bty="n",mar=c(4.5,4,3,1)); plot.new(); mtext(3,text=title.motif,font=2,line=1)
  }

  # Panel 4: By REGION (far right)
  if(has.region && length(regions) > 0){
    reg.vals <- lapply(regions, function(r){
      v <- vals[!is.na(trv.dat$REGION) & trv.dat$REGION==r & !is.na(vals) & is.finite(vals)]
      v[v >= xlim[1] & v <= xlim[2]]
    })
    plotDistribOverlaid(sub.list=reg.vals, sub.labels=regions, sub.colors=col.reg,
                        xlab=xlab, title=title.region, xlim=xlim, log.y=log.y)
  } else {
    par(bty="n",mar=c(4.5,4,3,1))
    plot(x=c(0,1),y=c(0,1),type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
    text(x=0.5,y=0.5,labels="No region data.")
    mtext(3,text=title.region,font=2,line=1)
  }
}

wrapperPlotTrvAlleleCount <- function(){
  trv.dat <- dat[dat$svtype %in% c("TR_SNV","TR_INS","TR_DEL") & !is.na(dat$AC) & dat$AC > 0, ]
  if(nrow(trv.dat) == 0) return(invisible(NULL))
  vals <- trv.dat$AC
  xlim <- c(1, quantile(vals, 0.99))
  png(paste(OUTDIR,"/main_plots/tr_allele_count_distribution.png",sep=""),
      res=300, height=1800, width=4*1800)
  layout(matrix(1:4, nrow=1))
  plotTrvDistribPanels(trv.dat=trv.dat, vals=vals, xlab="Non-Ref Allele Count", xlim=xlim,
                       title.all="TR Allele Count",
                       title.size="TR Allele Count by Size",
                       title.motif="TR Allele Count by Motif Length",
                       title.region="TR Allele Count by Region")
  dev.off()
}

wrapperPlotTrvSampleCount <- function(){
  trv.dat <- dat[dat$svtype %in% c("TR_SNV","TR_INS","TR_DEL") & !is.na(dat$carriers) & dat$carriers > 0, ]
  if(nrow(trv.dat) == 0) return(invisible(NULL))
  vals <- trv.dat$carriers
  xlim <- c(1, quantile(vals, 0.99))
  png(paste(OUTDIR,"/main_plots/tr_sample_count_distribution.png",sep=""),
      res=300, height=1800, width=4*1800)
  layout(matrix(1:4, nrow=1))
  plotTrvDistribPanels(trv.dat=trv.dat, vals=vals, xlab="Variant Sample Count", xlim=xlim,
                       title.all="TR Sample Count",
                       title.size="TR Sample Count by Size",
                       title.motif="TR Sample Count by Motif Length",
                       title.region="TR Sample Count by Region")
  dev.off()
}

wrapperPlotTrvExpansionRatio <- function(){
  if(!"TRV_EXPANSION_RATIO" %in% colnames(dat)) return(invisible(NULL))
  trv.dat <- dat[dat$svtype %in% c("TR_SNV","TR_INS","TR_DEL") & !is.na(dat$TRV_EXPANSION_RATIO), ]
  if(nrow(trv.dat) == 0) return(invisible(NULL))
  n.trv.total <- sum(dat$svtype %in% c("TR_SNV","TR_INS","TR_DEL"))
  n.trv.dropped <- n.trv.total - nrow(trv.dat)
  trv.drop.text <- if(n.trv.dropped > 0) paste0("(",prettyNum(n.trv.dropped,big.mark=",")," dropped - no alternate allele of distinct size)") else NULL
  vals <- trv.dat$TRV_EXPANSION_RATIO
  xlim <- c(0, 1)
  png(paste(OUTDIR,"/main_plots/tr_expansion_ratio_distribution.png",sep=""),
      res=300, height=1800, width=4*1800)
  layout(matrix(1:4, nrow=1))
  plotTrvDistribPanels(trv.dat=trv.dat, vals=vals, xlab="Expansion Ratio", xlim=xlim,
                       title.all="TR Expansion Ratio",
                       title.size="TR Expansion Ratio by Size",
                       title.motif="TR Expansion Ratio by Motif Length",
                       title.region="TR Expansion Ratio by Region",
                       drop.text=trv.drop.text)
  dev.off()
}


######################################
#####TR loci distribution plot
######################################
wrapperPlotTrLociDistrib <- function(){
  tr.types <- c("TR_SNV","TR_INS","TR_DEL")
  if(!any(dat$svtype %in% tr.types)) return(invisible(NULL))
  if(!"TRID" %in% colnames(dat)) return(invisible(NULL))

  # Build per-locus counts of normalized TR alleles by detailed subtype
  tr.dat <- dat[dat$svtype %in% tr.types,]
  tr.dat$base_vid <- sub("_[0-9]+$", "", tr.dat$VID)
  # Detailed subtype: TR_INS_SHORT, TR_INS_SV, TR_DEL_SHORT, TR_DEL_SV, TR_SNV
  tr.dat$detail_type <- tr.dat$svtype
  tr.dat$detail_type[tr.dat$svtype=="TR_INS" & tr.dat$length < 50] <- "TR_INS_SHORT"
  tr.dat$detail_type[tr.dat$svtype=="TR_INS" & tr.dat$length >= 50] <- "TR_INS_SV"
  tr.dat$detail_type[tr.dat$svtype=="TR_DEL" & tr.dat$length < 50] <- "TR_DEL_SHORT"
  tr.dat$detail_type[tr.dat$svtype=="TR_DEL" & tr.dat$length >= 50] <- "TR_DEL_SV"

  # For each unique TR locus (base_vid), get TRID and count alleles per detail type
  loci <- unique(tr.dat$base_vid)
  trid.map <- setNames(tr.dat$TRID[match(loci, tr.dat$base_vid)], loci)
  detail.types <- c("SNV","INS_SHORT","DEL_SHORT","INS_SV","DEL_SV")
  tr.detail.map <- c("TR_SNV","TR_INS_SHORT","TR_DEL_SHORT","TR_INS_SV","TR_DEL_SV")

  # Count TR alleles per locus per detail type
  tr.allele.counts <- sapply(tr.detail.map, function(dt){
    sapply(loci, function(l) sum(tr.dat$base_vid==l & tr.dat$detail_type==dt))
  })
  if(!is.matrix(tr.allele.counts)) tr.allele.counts <- matrix(tr.allele.counts, nrow=length(loci))
  colnames(tr.allele.counts) <- tr.detail.map

  # Count enveloped non-TR variants per TRID per type
  # Merge DUP into INS for comparison
  non.tr <- dat[!dat$svtype %in% tr.types & !is.na(dat$TRID) & dat$TRID != "" & dat$TRID != ".",]
  non.tr$cmp_type <- non.tr$svtype
  non.tr$cmp_type[non.tr$cmp_type %in% c("DUP_SHORT")] <- "INS_SHORT"
  non.tr$cmp_type[non.tr$cmp_type %in% c("DUP_SV")] <- "INS_SV"
  env.counts <- sapply(detail.types, function(dt){
    sapply(loci, function(l){
      tr.id <- trid.map[l]
      if(is.na(tr.id)) return(0)
      sum(non.tr$TRID==tr.id & non.tr$cmp_type==dt, na.rm=TRUE)
    })
  })
  if(!is.matrix(env.counts)) env.counts <- matrix(env.counts, nrow=length(loci))
  colnames(env.counts) <- detail.types

  # Plot: one scatter per comparison type
  n.panels <- length(detail.types)
  png(paste(OUTDIR,"/main_plots/tr_loci_distribution.png",sep=""),
      res=300, height=1800, width=n.panels*1800)
  par(mfrow=c(1, n.panels))
  for(k in seq_along(detail.types)){
    x <- tr.allele.counts[, tr.detail.map[k]]
    y <- env.counts[, detail.types[k]]
    par(mar=c(4.5,4.5,3,1), bty="n")
    plot(x, y, pch=19, cex=0.4, col=adjustcolor("gray30", alpha=0.5),
         xlab=paste("TR alleles (", tr.detail.map[k], ")", sep=""),
         ylab=paste("Enveloped ", detail.types[k], " variants", sep=""),
         main=detail.types[k])
    abline(0, 1, col="red", lty=2)
    # R-squared
    if(length(x) > 2 && sd(x) > 0 && sd(y) > 0){
      r2 <- cor(x, y)^2
      legend("topleft", bty="n", cex=0.8,
             legend=bquote(R^2 == .(round(r2, 3))))
    }
  }
  dev.off()
}

######################################
#####INS distribution plot
######################################
wrapperPlotInsDistrib <- function(){
  # Filter to variants with 'INS' in their VID
  ins.idx <- grep("INS", dat$VID, ignore.case=FALSE)
  if(length(ins.idx) == 0) return(invisible(NULL))
  ins.dat <- dat[ins.idx,]

  # Get allele_type; derive display labels
  if(!"allele_type" %in% colnames(ins.dat)) return(invisible(NULL))
  atype <- tolower(ins.dat$allele_type)
  # For 'ins' type: split into Unique, TR Enveloped, TR Parsed
  is.ins <- atype == "ins"
  is.tr.env <- is.ins & !is.na(ins.dat$TR_ENVELOPED) & ins.dat$TR_ENVELOPED %in% c(TRUE,"TRUE","1","true")
  is.tr.parsed <- if("TR_PARSED" %in% colnames(ins.dat)){
    is.ins & !is.na(ins.dat$TR_PARSED) & ins.dat$TR_PARSED %in% c(TRUE,"TRUE","1","true")
  } else rep(FALSE, nrow(ins.dat))

  # Build category vector
  cat.vec <- toupper(gsub("_ins", "", atype))
  cat.vec[cat.vec == "INS"] <- "Unique"
  cat.vec[is.tr.env] <- "TR Enveloped"
  cat.vec[is.tr.parsed] <- "TR Parsed"

  # Order: Unique first, then TR Enveloped, TR Parsed, then others alphabetically
  cat.tab <- table(cat.vec)
  fixed.cats <- c("Unique","TR Enveloped","TR Parsed")
  other.cats <- sort(setdiff(names(cat.tab), fixed.cats))
  all.cats <- c(fixed.cats[fixed.cats %in% names(cat.tab)], other.cats)
  cat.counts <- as.numeric(cat.tab[all.cats])
  names(cat.counts) <- all.cats
  total <- sum(cat.counts)

  # Colors
  n.cat <- length(all.cats)
  cat.cols <- setNames(colorRampPalette(c("#440154","#31688E","#35B779","#FDE725","#D95F02"))(n.cat), all.cats)

  # Size buckets for the 6 smaller bars
  sz.labs <- c("<50bp","50-100bp","100bp-500bp","500bp-5kb","5-50kb",">50kb")
  sz.mins <- c(0, tiny.max.size, small.max.size, medium.max.size, medlarge.max.size, large.max.size)
  sz.maxs <- c(tiny.max.size, small.max.size, medium.max.size, medlarge.max.size, large.max.size, huge.max.size)

  # Build stacked bar matrix: rows = categories, cols = ALL + size buckets
  build.col <- function(sub.dat){
    sub.atype <- tolower(sub.dat$allele_type)
    sub.is.ins <- sub.atype == "ins"
    sub.is.env <- sub.is.ins & !is.na(sub.dat$TR_ENVELOPED) & sub.dat$TR_ENVELOPED %in% c(TRUE,"TRUE","1","true")
    sub.is.parsed <- if("TR_PARSED" %in% colnames(sub.dat)){
      sub.is.ins & !is.na(sub.dat$TR_PARSED) & sub.dat$TR_PARSED %in% c(TRUE,"TRUE","1","true")
    } else rep(FALSE, nrow(sub.dat))
    sub.cat <- toupper(gsub("_ins", "", sub.atype))
    sub.cat[sub.cat == "INS"] <- "Unique"
    sub.cat[sub.is.env] <- "TR Enveloped"
    sub.cat[sub.is.parsed] <- "TR Parsed"
    sapply(all.cats, function(c) sum(sub.cat == c))
  }
  mat <- cbind(ALL=cat.counts,
               sapply(seq_along(sz.labs), function(i){
                 sub <- ins.dat[!is.na(ins.dat$length) & ins.dat$length >= sz.mins[i] & ins.dat$length < sz.maxs[i],]
                 if(nrow(sub)==0) return(setNames(rep(0, n.cat), all.cats))
                 build.col(sub)
               }))
  colnames(mat) <- c("ALL", sz.labs)

  # Plot
  png(paste(OUTDIR,"/main_plots/ins_distributions.png",sep=""),
      res=300, height=2400, width=4800)
  layout(matrix(1:7, nrow=1), widths=c(3, rep(1.5, 6)))

  # Main bar (ALL)
  par(mar=c(6,5,3,1), bty="n")
  vals <- mat[,"ALL"]
  starts <- cumsum(c(0, vals[-length(vals)]))
  ends <- cumsum(vals)
  ymax <- total * 1.1
  plot(x=c(0,1), y=c(0, ymax), type="n", xaxt="n", yaxt="n", xlab="", ylab="Count", yaxs="i")
  rect(xleft=0.15, xright=0.85, ybottom=starts, ytop=ends, col=cat.cols, border="white")
  # Labels for each category
  mid.y <- (starts + ends) / 2
  for(j in seq_along(all.cats)){
    if(vals[j] > 0){
      pct <- round(100*vals[j]/total, 1)
      text(0.5, mid.y[j], labels=paste0(all.cats[j], "\n", formatCount(vals[j]), " (", pct, "%)"),
           cex=0.6, col="white", font=2)
    }
  }
  axis(2, at=pretty(c(0, ymax)), labels=NA)
  axis(2, at=pretty(c(0, ymax)), tick=F, las=2, cex.axis=0.8, line=-0.4,
       labels=prettyNum(pretty(c(0, ymax)), big.mark=","))
  mtext(3, text="INS Variant Composition", font=2, line=1)
  axis(3, at=0.5, tick=F, line=-0.9, cex.axis=0.8,
       labels=paste("n=", prettyNum(total, big.mark=","), sep=""))

  # Size-bucketed bars
  for(s in seq_along(sz.labs)){
    par(mar=c(6,3,3,1), bty="n")
    vals.s <- mat[, sz.labs[s]]
    tot.s <- sum(vals.s)
    starts.s <- cumsum(c(0, vals.s[-length(vals.s)]))
    ends.s <- cumsum(vals.s)
    ymax.s <- max(tot.s * 1.1, 1)
    plot(x=c(0,1), y=c(0, ymax.s), type="n", xaxt="n", yaxt="n", xlab="", ylab="", yaxs="i")
    if(tot.s > 0){
      rect(xleft=0.15, xright=0.85, ybottom=starts.s, ytop=ends.s, col=cat.cols, border="white")
    }
    axis(2, at=pretty(c(0, ymax.s)), labels=NA)
    axis(2, at=pretty(c(0, ymax.s)), tick=F, las=2, cex.axis=0.7, line=-0.4,
         labels=prettyNum(pretty(c(0, ymax.s)), big.mark=","))
    mtext(3, text=sz.labs[s], font=2, line=1, cex=0.8)
    axis(3, at=0.5, tick=F, line=-0.9, cex.axis=0.7,
         labels=paste("n=", prettyNum(tot.s, big.mark=","), sep=""))
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

#TRV allele count distribution
wrapperPlotTrvAlleleCount()

#TRV sample count distribution
wrapperPlotTrvSampleCount()

#TRV expansion ratio distribution
wrapperPlotTrvExpansionRatio()

#TR loci distribution
wrapperPlotTrLociDistrib()

#INS distributions
wrapperPlotInsDistrib()
