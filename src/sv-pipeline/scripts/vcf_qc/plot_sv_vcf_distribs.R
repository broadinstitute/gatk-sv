#!/usr/bin/env Rscript

# Copyright (c) 2018 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# Helper script to plot VCF summary stats output from clean_vcf2bed_output.R


###Set master parameters
options(stringsAsFactors=F,scipen=1000)
rare.max.freq <- 0.01
uncommon.max.freq <- 0.1
common.max.freq <- 0.5
major.max.freq <- 1
tiny.max.size <- 100
small.max.size <- 500
medium.max.size <- 2500
medlarge.max.size <- 10000
large.max.size <- 50000
huge.max.size <- 300000000
sex.chroms <- c(1:22, paste("chr", 1:22, sep=""))


###################
###HELPER FUNCTIONS
###################
#General function to plot stacked bars from a matrix
plotStackedBars <- function(mat,colors,scaled=T,title=NULL){
  #Scale columns, if options
  if(scaled==T){
    mat <- apply(mat,2,function(vals){
      vals/sum(vals,na.rm=T)
    })
  }
  
  #Prepare plot area
  ymax <- max(apply(mat,2,sum,na.rm=T),na.rm=T)
  par(mar=c(4,4,2,1),bty="n")
  plot(x=c(0,ncol(mat)),y=c(0,1.02*ymax),type="n",
       xaxt="n",yaxt="n",xaxs="i",yaxs="i",xlab="",ylab="")
  
  #Add axes & title
  axis(2,at=axTicks(2),labels=NA)
  if(scaled==T){
    ylabs <- paste(round(100*axTicks(2),1),"%",sep="")
  }else{
    ylabs <- prettyNum(axTicks(2),big.mark=",")
  }
  axis(2,at=axTicks(2),tick=F,labels=ylabs,las=2,cex.axis=0.8,line=-0.4)
  mtext(3,line=0.5,text=title,font=2)
  
  #Iterate and plot bars
  sapply(1:ncol(mat),function(i){
    #Get plotting values
    vals <- mat[,i]
    starts <- as.numeric(cumsum(c(0,vals[-length(vals)])))
    ends <- as.numeric(cumsum(vals))
    
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
plotSVCountBars <- function(dat,svtypes,title=NULL,ylab="SV Count"){
  #Compute table
  counts <- as.data.frame(t(sapply(svtypes$svtype,function(svtype){
    c(svtype,
      length(which(dat$svtype==svtype)),
      svtypes[which(svtypes$svtype==svtype),2])
  })))
  counts[,2] <- as.numeric(counts[,2])
  
  #Prep plotting area
  par(bty="n",mar=c(3,4.5,2.5,0.5))
  plot(x=c(0,nrow(counts)),y=c(0,1.15*max(counts[,2])),type="n",
       xaxt="n",yaxt="n",xlab="",ylab="",xaxs="i",yaxs="i")
  
  #Add y-axis and title
  axis(2,at=axTicks(2),labels=NA)
  axis(2,at=axTicks(2),las=2,tick=F,line=-0.4,cex.axis=0.7,
       labels=prettyNum(axTicks(2),big.mark=","))
  mtext(2,text=ylab,line=3)
  mtext(3,line=0.5,text=title,font=2)
  
  #Plot per-svtype information
  sapply(1:nrow(counts),function(i){
    #Bars
    rect(xleft=i-0.85,xright=i-0.15,
         ybottom=0,ytop=counts[i,2],
         lwd=0.7,col=counts[i,3])
    #Counts
    text(x=i-0.5,y=counts[i,2],pos=3,col=counts[i,3],
         labels=prettyNum(counts[i,2],big.mark=","),cex=0.7)
    #Labels
    axis(1,at=i-0.5,line=-0.8,tick=F,las=2,cex.axis=0.8,
         labels=counts[i,1],col.axis=counts[i,3])
  })
  
  #Add number of SV to plot
  axis(3,at=mean(par("usr")[1:2]),line=-1.5,tick=F,cex.axis=0.8,
       labels=paste("n=",prettyNum(nrow(dat),big.mark=","),sep=""))
}
#Plot dot for fraction of total SV per chromosome
plotDotsSVperChrom <- function(dat,svtypes,title=NULL,ylab="Fraction of SV Type"){
  #Compute table
  mat <- sapply(svtypes$svtype,function(svtype){
    counts <- sapply(c(1:22,"X","Y"),function(contig){
      length(which(dat$chr==contig & dat$svtype==svtype))
    })
    if(sum(counts,na.rm=T)>0){
      return(counts/sum(counts))
    }else{
      return(rep(0,24))
    }
  })
  
  #Prep plotting area
  par(bty="n",mar=c(3,4.5,2.5,0.5))
  plot(x=c(0,24),y=c(0,1.15*max(mat)),type="n",
       xaxt="n",yaxt="n",xlab="",ylab="",xaxs="i",yaxs="i")
  
  #Add axes and title
  sapply(1:24,function(i){
    axis(1,at=i-0.5,tick=F,labels=c(1:22,"X","Y")[i],line=-0.8,cex.axis=0.7)
  })
  mtext(1,text="Chromosome",line=1.5)
  axis(2,at=axTicks(2),labels=NA)
  axis(2,at=axTicks(2),las=2,tick=F,line=-0.4,cex.axis=0.7,
       labels=paste(round(100*axTicks(2),digits=1),"%",sep=""))
  mtext(2,text=ylab,line=2.5)
  mtext(3,line=0.5,text=title,font=2)
  
  #Plot per-svtype information
  sapply(1:ncol(mat),function(i){
    points(x=(1:24)-seq(0.8,0.2,by=-0.6/(nrow(svtypes)-1))[i],
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
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/sv_count.all_sv.pdf",sep=""),
      height=4,width=2+(nrow(svtypes)/3))
  plotSVCountBars(dat=dat,svtypes=svtypes,
                  title="Variant Count (All SV)")
  dev.off()
  #Singletons
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/sv_count.singletons.pdf",sep=""),
      height=4,width=2+(nrow(svtypes)/3))
  plotSVCountBars(dat=dat[which(dat$AC==1),],svtypes=svtypes,
                  title="Variant Count (AC = 1)")
  dev.off()
  #Rare (>1 & <1%)
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/sv_count.rare_sv.pdf",sep=""),
      height=4,width=2+(nrow(svtypes)/3))
  plotSVCountBars(dat=dat[which(dat$AC>1 & dat$AF<rare.max.freq),],svtypes=svtypes,
                  title="Variant Count (AC > 1, AF < 1%)")
  dev.off()
  #Uncommon (≥1% & <10%)
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/sv_count.uncommon_sv.pdf",sep=""),
      height=4,width=2+(nrow(svtypes)/3))
  plotSVCountBars(dat=dat[which(dat$AF>=rare.max.freq & dat$AF<uncommon.max.freq),],svtypes=svtypes,
                  title="Variant Count (AF 1% - 10%)")
  dev.off()
  #Common (≥10% & <50%)
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/sv_count.common_sv.pdf",sep=""),
      height=4,width=2+(nrow(svtypes)/3))
  plotSVCountBars(dat=dat[which(dat$AF>=uncommon.max.freq & dat$AF<common.max.freq),],svtypes=svtypes,
                  title="Variant Count (AF 10% - 50%)")
  dev.off()
  #Major (≥50%%)
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/sv_count.major_sv.pdf",sep=""),
      height=4,width=2+(nrow(svtypes)/3))
  plotSVCountBars(dat=dat[which(dat$AF>=common.max.freq),],svtypes=svtypes,
                  title="Variant Count (AF > 50%)")
  dev.off()
  #Tiny
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/sv_count.tiny_sv.pdf",sep=""),
      height=4,width=2+(nrow(svtypes)/3))
  plotSVCountBars(dat=dat[which(dat$length<tiny.max.size),],svtypes=svtypes,
                  title="Variant Count (< 100bp)")
  dev.off()
  #Small
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/sv_count.small_sv.pdf",sep=""),
      height=4,width=2+(nrow(svtypes)/3))
  plotSVCountBars(dat=dat[which(dat$length>=tiny.max.size & dat$length<small.max.size),],svtypes=svtypes,
                  title="Variant Count (100 - 500bp)")
  dev.off()
  #Medium
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/sv_count.medium_sv.pdf",sep=""),
      height=4,width=2+(nrow(svtypes)/3))
  plotSVCountBars(dat=dat[which(dat$length>=small.max.size & dat$length<medium.max.size),],svtypes=svtypes,
                  title="Variant Count (500bp - 2.5kb)")
  dev.off()
  #Medium-Large
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/sv_count.medlarge_sv.pdf",sep=""),
      height=4,width=2+(nrow(svtypes)/3))
  plotSVCountBars(dat=dat[which(dat$length>=medium.max.size & dat$length<medlarge.max.size),],svtypes=svtypes,
                  title="Variant Count (2.5 - 10kb)")
  dev.off()
  #Large
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/sv_count.large_sv.pdf",sep=""),
      height=4,width=2+(nrow(svtypes)/3))
  plotSVCountBars(dat=dat[which(dat$length>=medlarge.max.size & dat$length<large.max.size),],svtypes=svtypes,
                  title="Variant Count (10 - 50kb)")
  dev.off()
  #Huge
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/sv_count.huge_sv.pdf",sep=""),
      height=4,width=2+(nrow(svtypes)/3))
  plotSVCountBars(dat=dat[which(dat$length>=large.max.size & dat$length<huge.max.size),],svtypes=svtypes,
                  title="Variant Count (> 50kb)")
  dev.off()
  #Dotplot per chromosome
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/sv_count_per_chromosome.all_sv.pdf",sep=""),
      height=4,width=6)
  plotDotsSVperChrom(dat=dat,svtypes=svtypes,
                     title="Variants per Chromosome (All SV)")
  dev.off()
  #Size x AF Grid
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/sv_counts.size_freq_grid.pdf",sep=""),
      height=6*1.75,width=7*1.75)
  #Iterator data frames
  AF.df <- data.frame("label"=c("AC=1","AF<1%","1-10%",
                                "10-50%",">50%"),
                      "max"=c(1.1/(2*nsamp),rare.max.freq,uncommon.max.freq,
                              common.max.freq,major.max.freq))
  size.df <- data.frame("label"=c("<100bp","100bp-\n500bp","500bp-\n2.5kb",
                                  "2.5-10kb","10-50kb",">50kb"),
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
      plotset <- dat[which(dat$AF>min.AF & dat$AF<=max.AF & 
                             dat$length>min.size & dat$length<=max.size),]
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
        plotSVCountBars(dat=plotset,svtypes=svtypes,
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
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/sv_counts.by_frequency_scaled.pdf",sep=""),
      height=4,width=4)
  plotStackedBars(mat=AF.mat,colors=svtypes$color,scale=T,
                  title="SV Count by Allele Frequency (Scaled)")
  abline(v=1,lty=2,col="gray50")
  axis(1,at=1,labels=NA,col="gray75",tck=-0.22)
  dev.off()  
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/sv_counts.by_frequency_raw.pdf",sep=""),
      height=4,width=4)
  plotStackedBars(mat=AF.mat,colors=svtypes$color,scale=F,
                  title="SV Count by Allele Frequency")
  abline(v=1,lty=2,col="gray50")
  axis(1,at=1,labels=NA,col="gray75",tck=-0.22)
  dev.off()
  #Stacked bars of SV count by size
  size.mat <- as.data.frame(sapply(0:nrow(size.df),function(i){
    #Get cutoffs
    min.size <- max(c(0,size.df[max(c(0,i-1)),2]))
    max.size <- min(c(size.df[i,2],300000000))
    #Get data subset
    plotset <- dat[which(dat$length>min.size & dat$length<=max.size),]
    #Tabulate counts per cutoff
    sapply(svtypes$svtype,function(svtype){
      length(which(plotset$svtype==svtype))
    })
  }))
  colnames(size.mat) <- c("ALL",size.df[,1])
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/sv_counts.by_size_scaled.pdf",sep=""),
      height=4,width=4)
  plotStackedBars(mat=size.mat,colors=svtypes$color,scale=T,
                  title="SV Count by Size (Scaled)")
  abline(v=1,lty=2,col="gray50")
  axis(1,at=1,labels=NA,col="gray75",tck=-0.22)
  dev.off()  
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/sv_counts.by_size_raw.pdf",sep=""),
      height=4,width=4)
  plotStackedBars(mat=size.mat,colors=svtypes$color,scale=F,
                  title="SV Count by Size")
  abline(v=1,lty=2,col="gray50")
  axis(1,at=1,labels=NA,col="gray75",tck=-0.22)
  dev.off()
  pdf(paste(OUTDIR,"/main_plots/VCF_QC.SV_counts.merged.pdf",sep=""),
      height=5,width=9)
  #Merged
  layout(matrix(c(1,2,3,1,4,5),byrow=T,nrow=2),
         widths=c(3,2,2))
  plotSVCountBars(dat=dat,svtypes=svtypes,
                  title="Variant Count (All SV)")
  plotStackedBars(mat=AF.mat,colors=svtypes$color,scale=F,
                  title="SV Count by AF")
  abline(v=1,lty=2,col="gray50")
  axis(1,at=1,labels=NA,col="gray75",tck=-0.22)
  plotStackedBars(mat=AF.mat,colors=svtypes$color,scale=T,
                  title="SV Count by AF (Scaled)")
  abline(v=1,lty=2,col="gray50")
  axis(1,at=1,labels=NA,col="gray75",tck=-0.22)
  plotStackedBars(mat=size.mat,colors=svtypes$color,scale=F,
                  title="SV Count by Size")
  abline(v=1,lty=2,col="gray50")
  axis(1,at=1,labels=NA,col="gray75",tck=-0.22)
  plotStackedBars(mat=size.mat,colors=svtypes$color,scale=T,
                  title="SV Count by Size (Scaled)")
  abline(v=1,lty=2,col="gray50")
  axis(1,at=1,labels=NA,col="gray75",tck=-0.22)
  dev.off()
}


###############
#####Size plots
###############
#Plot single size distribution
plotSizeDistrib <- function(dat, svtypes, n.breaks=150, k=10,
                            min.size=50, max.size=1000000,
                            autosomal=F, biallelic=F,
                            title=NULL, legend=F, lwd.cex=1, text.cex=1){
  #Filter/process sizes & compute range + breaks
  filter.legend <- NULL
  if(autosomal==T){
    dat <- dat[which(dat$chr %in% sex.chroms),]
    filter.legend <- c(filter.legend,"Autosomal SV only")
  }
  if(biallelic==T){
    dat <- dat[which(dat$other_gts==0 & dat$missing_gts<dat$genotyped_samples),]
    filter.legend <- c(filter.legend,"Biallelic SV only")
  }
  sizes <- log10(dat$length)
  if(length(sizes)>0){
    xlims <- range(sizes[which(!is.infinite(sizes))],na.rm=T)
    xlims[1] <- max(c(log10(min.size),xlims[1]))
    xlims[2] <- min(c(log10(max.size),xlims[2]))
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
    mtext(2,text="Fraction of SV",line=2,cex=text.cex)
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
                                  n.breaks=100, min.size=50, max.size=1000000,
                                  autosomal=F, biallelic=T, title=NULL, lwd.cex=1){
  #Process sizes & compute range + breaks
  filter.legend <- NULL
  if(autosomal==T){
    dat <- dat[which(dat$chr %in% sex.chroms),]
    filter.legend <- c(filter.legend,"Autosomal SV only")
  }
  if(biallelic==T){
    dat <- dat[which(dat$other_gts==0 & dat$missing_gts<dat$genotyped_samples),]
    filter.legend <- c(filter.legend,"Biallelic SV only")
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
    xlims <- range(sizes[which(!is.infinite(unlist(sizes)))],na.rm=T)
    xlims[1] <- max(c(log10(min.size),xlims[1]))
    xlims[2] <- min(c(log10(max.size),xlims[2]))
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
    mtext(2,text="Fraction of SV",line=2)
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
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/size_distribution.all_sv.pdf",sep=""),
      height=4,width=6)
  plotSizeDistrib(dat=dat,svtypes=svtypes,
                  title="Size Distribution (All SV)",
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
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/size_distribution.rare_sv.pdf",sep=""),
      height=4,width=6)
  plotSizeDistrib(dat=dat[which(dat$AC>1 & dat$AF<rare.max.freq),],svtypes=svtypes,
                  autosomal=F, biallelic=T,
                  title="Size Distribution (AC > 1, AF < 1%)",
                  legend=T)
  dev.off()
  
  #Uncommon (≥1% & <10%)
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/size_distribution.uncommon_sv.pdf",sep=""),
      height=4,width=6)
  plotSizeDistrib(dat=dat[which(dat$AF>=rare.max.freq & dat$AF<uncommon.max.freq),],svtypes=svtypes,
                  autosomal=F, biallelic=T,
                  title="Size Distribution (AF 1% - 10%)",
                  legend=T)
  dev.off()
  
  #Common (≥10% & <50%)
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/size_distribution.common_sv.pdf",sep=""),
      height=4,width=6)
  plotSizeDistrib(dat=dat[which(dat$AF>=uncommon.max.freq & dat$AF<common.max.freq),],svtypes=svtypes,
                  autosomal=F, biallelic=T,
                  title="Size Distribution (AF 10% - 50%)",
                  legend=T)
  dev.off()
  
  #Major (≥50%%)
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/size_distribution.major_sv.pdf",sep=""),
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
                        title="Size Distributions by Allele Frequency")
  dev.off()
  
  #Merged
  pdf(paste(OUTDIR,"/main_plots/VCF_QC.size_distributions.merged.pdf",sep=""),
      height=6,width=10)
  layout(matrix(c(1,1,1,2,2,3,4,5,6,7),byrow=T,nrow=2),
         heights=c(4,2))
  plotSizeDistrib(dat=dat,svtypes=svtypes,
                  title="Size Distribution (All SV)",
                  legend=T, lwd.cex=1.5)
  plotSizeDistribSeries(dat=dat,svtypes=svtypes,
                        max.AFs=c(1.1/(2*nsamp),rare.max.freq,uncommon.max.freq,
                                  common.max.freq,major.max.freq),
                        legend.labs=c("Singleton","<1%","1-10%","10-50%",">50%"),
                        title="Size Distributions by Allele Frequency",
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
    filter.legend <- c(filter.legend,"Autosomal SV only")
  }
  if(biallelic==T){
    dat <- dat[which(dat$other_gts==0 & dat$missing_gts<dat$genotyped_samples),]
    filter.legend <- c(filter.legend,"Biallelic SV only")
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
    mtext(1,text="Allele Frequency",line=3,cex=lwd.cex)
    axis(2,at=axTicks(2),tck=-0.025,labels=NA)
    axis(2,at=axTicks(2),tick=F,line=-0.4,cex.axis=0.8,las=2,
         labels=paste(round(100*axTicks(2),1),"%",sep=""))
    mtext(2,text="Fraction of SV",line=2.2,cex=lwd.cex)
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
    filter.legend <- c(filter.legend,"Autosomal SV only")
  }
  if(biallelic==T){
    dat <- dat[which(dat$other_gts==0 & dat$missing_gts<dat$genotyped_samples),]
    filter.legend <- c(filter.legend,"Biallelic SV only")
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
    mtext(1,text="Allele Frequency",line=3)
    axis(2,at=axTicks(2),tck=-0.025,labels=NA)
    axis(2,at=axTicks(2),tick=F,line=-0.4,cex.axis=0.8,las=2,
         labels=paste(round(100*axTicks(2),1),"%",sep=""))
    mtext(2,text="Fraction of SV",line=2.2)
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
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/freq_distribution.all_sv.pdf",sep=""),
      height=4,width=4)
  plotFreqDistrib(dat=dat,svtypes=svtypes,
                  title="AF Distribution (All SV)",
                  legend=T)
  dev.off()
  
  #Tiny (<100bp)
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/freq_distribution.tiny_sv.pdf",sep=""),
      height=4,width=4)
  plotFreqDistrib(dat=dat[which(dat$length<tiny.max.size),],svtypes=svtypes,
                  title="AF Distribution (< 100bp)",
                  legend=T)
  dev.off()
  
  #Small (>100bp & <500bp)
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/freq_distribution.small_sv.pdf",sep=""),
      height=4,width=4)
  plotFreqDistrib(dat=dat[which(dat$length>=tiny.max.size & dat$length<small.max.size),],svtypes=svtypes,
                  title="AF Distribution (100bp - 500bp)",
                  legend=T)
  dev.off()
  
  #Medium (>500bp & <2.5kb)
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/freq_distribution.medium_sv.pdf",sep=""),
      height=4,width=4)
  plotFreqDistrib(dat=dat[which(dat$length>=small.max.size & dat$length<medium.max.size),],svtypes=svtypes,
                  title="AF Distribution (500bp - 2.5kb)",
                  legend=T)
  dev.off()
  
  #Med-Large (>2.5kb & <10kb)
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/freq_distribution.medlarge_sv.pdf",sep=""),
      height=4,width=4)
  plotFreqDistrib(dat=dat[which(dat$length>=medium.max.size & dat$length<medlarge.max.size),],svtypes=svtypes,
                  title="AF Distribution (2.5kb - 10kb)",
                  legend=T)
  dev.off()
  
  #Large (>10kb & <50kb)
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/freq_distribution.large_sv.pdf",sep=""),
      height=4,width=4)
  plotFreqDistrib(dat=dat[which(dat$length>=medlarge.max.size & dat$length<large.max.size),],svtypes=svtypes,
                  title="AF Distribution (10kb - 50kb)",
                  legend=T)
  dev.off()
  
  #Huge (>50kb)
  pdf(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/freq_distribution.huge_sv.pdf",sep=""),
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
                        legend.labs=c("<100bp","100-\n500bp","500bp-\n2.5kb",
                                      "2.5-10kb","10kb-50kb",">50kb"),
                        title="AF Distributions by SV Size")
  dev.off()
  
  #Merged
  pdf(paste(OUTDIR,"/main_plots/VCF_QC.freq_distributions.merged.pdf",sep=""),
      height=6,width=10)
  layout(matrix(c(1,1,1,2,2,2,
                  3,4,5,6,7,8),
                byrow=T,nrow=2),
         heights=c(4,2))
  plotFreqDistrib(dat=dat,svtypes=svtypes,
                  title="AF Distribution (All SV)",
                  legend=T)
  plotFreqDistribSeries(dat=dat,svtypes=svtypes,
                        max.sizes=c(tiny.max.size,small.max.size,medium.max.size,
                                    medlarge.max.size,large.max.size,huge.max.size),
                        legend.labs=c("<100bp","100bp-\n500bp","500bp-\n2.5kb",
                                      "2.5-10kb","10kb-50kb",">50kb"),
                        title="AF Distributions by SV Size")
  plotFreqDistrib(dat=dat[which(dat$length<tiny.max.size),],svtypes=svtypes,
                  title="< 100bp",lwd.cex=0.7)
  plotFreqDistrib(dat=dat[which(dat$length>=tiny.max.size & dat$length<small.max.size),],svtypes=svtypes,
                  title="100bp - 500bp",lwd.cex=0.7)
  plotFreqDistrib(dat=dat[which(dat$length>=small.max.size & dat$length<medium.max.size),],svtypes=svtypes,
                  title="500bp - 2.5kb",lwd.cex=0.7)
  plotFreqDistrib(dat=dat[which(dat$length>=medium.max.size & dat$length<medlarge.max.size),],svtypes=svtypes,
                  title="2.5kb - 10kb",lwd.cex=0.7)
  plotFreqDistrib(dat=dat[which(dat$length>=medlarge.max.size & dat$length<large.max.size),],svtypes=svtypes,
                  title="10kb - 50kb",lwd.cex=0.7)
  plotFreqDistrib(dat=dat[which(dat$length>=large.max.size),],svtypes=svtypes,
                  title="> 50kb",lwd.cex=0.7)
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
    
    #Generate HW plot frame
    par(mar=c(1,3.5,3,0.5),bty="n")
    plot(x=1.15*c(-1/sqrt(3),1/sqrt(3)),y=c(-0.15,1.15),type="n",
         xaxt="n",yaxt="n",xlab="",ylab="",xaxs="i",yaxs="i")
    segments(x0=c(-1/sqrt(3),0,1/sqrt(3)),
             x1=c(0,1/sqrt(3),-1/sqrt(3)),
             y0=c(0,1,0),y1=c(1,0,0))
    HWTernaryPlot(X=HW.mat,n=nsamp,newframe=F,
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
    HWTernaryPlot(X=HW.mat,n=nsamp,newframe=F,
                  vbounds=F,mafbounds=F,
                  region=1,vertexlab=NA,
                  alpha=0.05/nrow(HW.mat),
                  curvecols=c("#4DAC26","#AC26A1",NA,NA),
                  pch=21,cex=0.5,signifcolour=F,markercol=HW.cols,
                  markerbgcol=adjustcolor(HW.cols,alpha=0.25))
    
    #Add legend
    n.pass <- length(which(HW.p>=0.05))
    n.nom <- length(which(HW.p<0.05 & HW.p>=0.05/nrow(HW.mat)))
    n.bonf <- length(which(HW.p<0.05/nrow(HW.mat)))
    if(full.legend==T){
      legend("topright",pch=19,col=c("#4DAC26","#81F850","#AC26A1"),pt.cex=2,
             legend=c(paste("SV in H-W equilibrium\n(n=",
                            prettyNum(n.pass,big.mark=","),"; ",
                            round(100*(n.pass/nrow(HW.mat)),2),"%)\n",sep=""),
                      paste("SV not in H-W equilibrium\n(Nominal; n=",
                            prettyNum(n.nom,big.mark=","),"; ",
                            round(100*(n.nom/nrow(HW.mat)),2),"%)\n",sep=""),
                      paste("SV not in H-W equilibrium\n(Bonferroni; n=",
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
           legend=c("Biallelic SV only","Autosomal SV only"))
  }
}
#Correlation of carrier frequency & AF
plotAlleleCarrierCorrelation <- function(dat,autosomal=T,biallelic=T,
                                         title="SV Carrier Freq. vs. Allele Freq."){
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
  mtext(2,text="Allele Frequency",line=2.5)
  mtext(3,text=title,line=0.5,font=2)
  axis(4,at=c(0.5,1),las=2,line=-0.8,tick=F,col.axis="gray50",cex.axis=0.7,
       labels=c("All\nHet.","All\nHom."))
  
  #Plot points & rolling mean
  points(x=CF,y=AF,cex=0.2)
  points(x=rollmean(CF,k=100),
         y=rollmean(AF,k=100),
         type="l",col="red",lwd=1.25)
  
  #Add legend
  legend("bottomright",bg=NA,legend=c("SV Site","Rolling Mean"),cex=0.8,bty="n",
         lwd=c(NA,2),col=c("black","red"),pch=c(1,NA),pt.cex=c(0.4,NA))
  
  #Add filter labels
  if(!is.null(filter.legend)){
    legend("topleft",bg=NA,bty="n",pch=NA,legend=filter.legend,cex=0.8)
  }
}

#Wrapper to plot all HW distributions
wrapperPlotAllHWDistribs <- function(){
  #All SV
  png(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/gt_distribution.all_sv.png",sep=""),
      res=300,height=1800,width=1800)
  plotHWSingle(dat=dat,svtypes=svtypes,
               title="Genotype Distribution (All SV)")
  dev.off()
  #Tiny
  png(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/gt_distribution.tiny_sv.png",sep=""),
      res=300,height=1800,width=1800)
  plotHWSingle(dat=dat[which(dat$length<tiny.max.size),],svtypes=svtypes,
               title="Genotype Distribution (< 100bp)")
  dev.off()
  #Small
  png(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/gt_distribution.small_sv.png",sep=""),
      res=300,height=1800,width=1800)
  plotHWSingle(dat=dat[which(dat$length>=tiny.max.size & dat$length<small.max.size),],svtypes=svtypes,
               title="Genotype Distribution (100 - 500bp)")
  dev.off()
  #Medium
  png(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/gt_distribution.medium_sv.png",sep=""),
      res=300,height=1800,width=1800)
  plotHWSingle(dat=dat[which(dat$length>=small.max.size & dat$length<medium.max.size),],svtypes=svtypes,
               title="Genotype Distribution (500bp - 2.5kb)")
  dev.off()
  #Med-Large
  png(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/gt_distribution.medlarge_sv.png",sep=""),
      res=300,height=1800,width=1800)
  plotHWSingle(dat=dat[which(dat$length>=medium.max.size & dat$length<medlarge.max.size),],svtypes=svtypes,
               title="Genotype Distribution (2.5 - 10kb)")
  dev.off()
  #Large
  png(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/gt_distribution.large_sv.png",sep=""),
      res=300,height=1800,width=1800)
  plotHWSingle(dat=dat[which(dat$length>=medlarge.max.size & dat$length<large.max.size),],svtypes=svtypes,
               title="Genotype Distribution (10 - 50kb)")
  dev.off()
  #Huge
  png(paste(OUTDIR,"/supporting_plots/vcf_summary_plots/gt_distribution.huge_sv.png",sep=""),
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
  png(paste(OUTDIR,"main_plots/VCF_QC.genotype_distributions.png",sep=""),res=300,
      height=1200,width=3.5*1200)
  layout(matrix(c(1,2,3,4,5,
                  1,2,6,7,8),
                byrow=T,nrow=2),
         widths=c(2,2,1,1,1))
  plotAlleleCarrierCorrelation(dat=dat)
  plotHWSingle(dat=dat,svtypes=svtypes,
               title="Genotype Distribution (All SV)")
  plotHWSingle(dat=dat[which(dat$length<tiny.max.size),],svtypes=svtypes,
               title="< 100bp",full.legend=F,lab.cex=0.7)
  plotHWSingle(dat=dat[which(dat$length>=tiny.max.size & dat$length<small.max.size),],svtypes=svtypes,
               title="100 - 500bp",full.legend=F,lab.cex=0.7)
  plotHWSingle(dat=dat[which(dat$length>=small.max.size & dat$length<medium.max.size),],svtypes=svtypes,
               title="500bp - 2.5kb",full.legend=F,lab.cex=0.7)
  plotHWSingle(dat=dat[which(dat$length>=medium.max.size & dat$length<medlarge.max.size),],svtypes=svtypes,
               title="2.5 - 10kb",full.legend=F,lab.cex=0.7)
  plotHWSingle(dat=dat[which(dat$length>=medlarge.max.size & dat$length<large.max.size),],svtypes=svtypes,
               title="10 - 50kb",full.legend=F,lab.cex=0.7)
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

###Plotting block
#SV counts
wrapperPlotAllCountBars()

#SV sizes
wrapperPlotAllSizeDistribs()

#SV frequencies
wrapperPlotAllFreqDistribs()

#Genotype frequencies
wrapperPlotAllHWDistribs()

