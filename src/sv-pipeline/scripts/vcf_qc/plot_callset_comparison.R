#!/usr/bin/env Rscript

# Helper script to plot VCF benchmarking against external dataset


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
ovr.cat.cols <- c("#76E349","#4DAC26","#2C750E",
                  "#FCE744","#EBD321","#D6BD01",
                  "#AC26A1","gray42")


###########################
###GENERAL HELPER FUNCTIONS
###########################
#Convert overlap data into vector of categories
categoryBreakdown <- function(dat,norm=F){
  #Get variant indexes for each group
  a1 <- which(dat$ovr1a != "NO_OVR")
  a2 <- which(dat$ovr2a != "NO_OVR")
  b1 <- which(dat$ovr1b != "NO_OVR")
  b2 <- which(dat$ovr2b != "NO_OVR")
  ab3 <- which(dat$ovr3 != "NO_OVR")
  
  #Get hierarchical counts of variants
  a12.n <- length(intersect(a1,a2))
  a1.n <- length(which(!(a1 %in% a2)))
  a2.n <- length(which(!(a2 %in% a1)))
  b12.n <- length(which(!(intersect(b1,b2) %in% c(a1,a2))))
  b1.n <- length(which(!(b1 %in% c(a1,a2,b2))))
  b2.n <- length(which(!(b2 %in% c(a1,a2,b1))))
  ab3.n <- length(which(!(ab3 %in% c(a1,a2,b1,b2))))
  f.n <- nrow(dat)-sum(c(a12.n,a1.n,a2.n,b12.n,b1.n,b2.n,ab3.n))
  
  #Return vector
  out.vect <- c(a12.n,a1.n,a2.n,b12.n,b1.n,b2.n,ab3.n,f.n)
  if(norm==T){
    out.vect <- out.vect/nrow(dat)
  }
  return(out.vect)
}
#Breakdown data by size
categoryBreakdownBySize <- function(dat,max.sizes,size.labels=NULL,norm=F){
  #Create data frame of size cutoffs & group labels
  if(is.null(size.labels)){
    size.labels <- c("ALL",paste(c(0,max.sizes[-length(max.sizes)]),"-",max.sizes,sep=""))
  }else{
    size.labels <- c("ALL",size.labels)
  }
  sizes.df <- data.frame("min"=c(0,0,max.sizes[-length(max.sizes)]),
                         "max"=c(300000000,max.sizes),
                         "label"=size.labels)
  
  #Create normalized vector of overlap categories per size threshold
  mat <- as.data.frame(sapply(1:nrow(sizes.df),function(i){
    categoryBreakdown(dat[which(dat$length>sizes.df[i,1] & dat$length<=sizes.df[i,2]),],norm=norm)
  }))
  colnames(mat) <- sizes.df$label
  rownames(mat) <- c("Recip. & Bkpt.","Recip.","Bkpt.",
                     "Recip. & Bkpt.\n(Diff. SV types)",
                     "Recip.\n(Diff. SV types)",
                     "Bkpt.\n(Diff. SV types)",
                     "Nearby\n(+/- 250bp)","No Overlap")  
  return(mat)
}
#Breakdown data by frequency
categoryBreakdownByFreq <- function(dat,max.freqs,freq.labels=NULL,norm=F){
  #Create data frame of size cutoffs & group labels
  if(is.null(freq.labels)){
    freq.labels <- c("ALL",paste(c(0,max.freqs[-length(max.freqs)]),"-",max.freqs,sep=""))
  }else{
    freq.labels <- c("ALL",freq.labels)
  }
  freqs.df <- data.frame("min"=c(0,0,max.freqs[-length(max.freqs)]),
                         "max"=c(1,max.freqs),
                         "label"=freq.labels)
  
  #Create normalized vector of overlap categories per freq threshold
  mat <- as.data.frame(sapply(1:nrow(freqs.df),function(i){
    categoryBreakdown(dat[which(dat$AF>freqs.df[i,1] & dat$AF<=freqs.df[i,2]),],norm=norm)
  }))
  colnames(mat) <- freqs.df$label
  rownames(mat) <- c("Recip. & Bkpt.","Recip.","Bkpt.",
                     "Recip. & Bkpt.\n(Diff. SV types)",
                     "Recip.\n(Diff. SV types)",
                     "Bkpt.\n(Diff. SV types)",
                     "Nearby\n(+/- 250bp)","No Overlap")  
  return(mat)
}
#Breakdown data by SV class
categoryBreakdownByClass <- function(dat,norm=F){
  #Create vector of unique SV types
  svtypes <- sort(unique(as.character(dat$svtype)))
  
  #Create normalized vector of overlap categories per sv class
  mat <- as.data.frame(sapply(svtypes,function(svtype){
    categoryBreakdown(dat[which(dat$svtype==svtype),],norm=norm)
  }))
  mat <- cbind(categoryBreakdown(dat,norm=norm),mat)
  colnames(mat) <- c("ALL",svtypes)
  rownames(mat) <- c("Recip. & Bkpt.","Recip.","Bkpt.",
                     "Recip. & Bkpt.\n(Diff. SV types)",
                     "Recip.\n(Diff. SV types)",
                     "Bkpt.\n(Diff. SV types)",
                     "Nearby\n(+/- 250bp)","No Overlap")   
  return(mat)
}
#Extract best-matching AF pairs (after excluding category 3)
getFreqPairs <- function(dat){
  freqPairs <- as.data.frame(t(apply(dat[,7:11],1,function(vals){
    vals <- as.numeric(vals)
    benchmark.AF <- vals[1]
    callset.AFs <- vals[-1]
    deltas <- abs(callset.AFs - benchmark.AF)
    if(all(is.na(deltas))){
      callset.AF <- NA
    }else{
      best.match.AF <- min(deltas, na.rm=T)
      callset.AF <- callset.AFs[head(which(deltas == best.match.AF), 1)]
    }
    return(c(benchmark.AF,callset.AF))
  })))
}


###########
###BARPLOTS
###########
#General function to plot stacked bars from a matrix
plotStackedBars <- function(mat,colors,scaled=T,title=NULL,legend=F){
  #Scale columns, if options
  col.sums <- apply(mat,2,sum)
  col.any.sup <- apply(mat,2,function(vals){
    sum(vals[-length(vals)])
  })
  if(scaled==T){
    mat <- apply(mat,2,function(vals){
      vals/sum(vals,na.rm=T)
    })
  }
  
  #Prepare plot area
  if(legend==T){
    layout(matrix(1:2,nrow=1,byrow=T),
           widths=c(5,1))
  }
  par(mar=c(5,4,3,1),bty="n")
  ymax <- max(apply(mat,2,sum,na.rm=T),na.rm=T)
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
  mtext(3,line=1.5,text=title,font=2)
  
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
    
    #Add category label
    axis(1,at=i-0.5,las=2,line=-0.8,cex.axis=0.8,tick=F,
         labels=paste(colnames(mat)[i],"\n(n=",prettyNum(col.sums[i],big.mark=","),")",sep=""))
    
    #Add top label indicating how many SV had any support
    axis(3,at=i-0.5,line=-1.2,cex.axis=0.7,tick=F,
         labels=paste(prettyNum(col.any.sup[i],big.mark=",")," /\n",
                      prettyNum(col.sums[i],big.mark=","),sep=""))
  })
  
  #Add legend, if optioned
  if(legend==T){
    par(mar=c(4,0.1,2,0.1),bty="n")
    plot(x=c(0,5),y=c(0,nrow(mat)),type="n",
         xaxs="i",yaxs="i",xlab="",ylab="",xaxt="n",yaxt="n")
    points(x=rep(0.4,times=nrow(mat)),
           y=(nrow(mat):1)-0.5,
           pch=15,col=rev(colors),cex=2)
    text(x=rep(0.5,times=nrow(mat)),
         y=(nrow(mat):1)-0.5,
         labels=rev(rownames(mat)),pos=4,cex=0.75)
  }
}
#Master wrapper to generate all barplots
wrapperAllBarplots <- function(){
  #Barplot by size
  size.mat <- categoryBreakdownBySize(dat,norm=F,
                                      max.sizes=c(tiny.max.size,small.max.size,medium.max.size,
                                                  medlarge.max.size,large.max.size,huge.max.size),
                                      size.labels=c("<100bp","100-500bp","500bp-2.5kb",
                                                    "2.5-10kb","10-50kb",">50kb"))
  if(!is.null(prefix)){
    title <- paste("SV Overlap by Size (",prefix,")",sep="")
    pdf.path <- paste(OUTDIR,"/supporting_plots/",prefix,".SV_overlap_barplot.by_size.pdf",sep="")
  }else{
    title <- "SV Overlap by Size"
    pdf.path <- paste(OUTDIR,"/supporting_plots/","SV_overlap_barplot.by_size.pdf",sep="")
  }
  if(carrierFreqs==F){
    title.prefix <- "Allele "
  }else{
    title.prefix <- "Carrier "
  }
  
  pdf(pdf.path,height=5,width=7)
  plotStackedBars(size.mat,colors=ovr.cat.cols,scaled=T,legend=T,title=title)
  dev.off()
  
  #Barplot by frequency
  freq.mat <- categoryBreakdownByFreq(dat,norm=F,
                                      max.freqs=c(rare.max.freq,uncommon.max.freq,
                                                  common.max.freq,major.max.freq),
                                      freq.labels=c("<1%","1-10%","10-50%",">50%"))
  if(!is.null(prefix)){
    title <- paste("SV Overlap by ",title.prefix,"Freq. (",prefix,")",sep="")
    pdf.path <- paste(OUTDIR,"/supporting_plots/",prefix,".SV_overlap_barplot.by_freq.pdf",sep="")
  }else{
    title <- paste("SV Overlap by ",title.prefix,"Freq.",sep="")
    pdf.path <- paste(OUTDIR,"/supporting_plots/","SV_overlap_barplot.by_freq.pdf",sep="")
  }
  pdf(pdf.path,height=5,width=7)
  plotStackedBars(freq.mat,colors=ovr.cat.cols,scaled=T,legend=T,title=title)
  dev.off()
  
  #Barplot by sv type
  svtype.mat <- categoryBreakdownByClass(dat,norm=F)
  if(!is.null(prefix)){
    title <- paste("SV Overlap by Class (",prefix,")",sep="")
    pdf.path <- paste(OUTDIR,"/supporting_plots/",prefix,".SV_overlap_barplot.by_svtype.pdf",sep="")
  }else{
    title <- "SV Overlap by Class"
    pdf.path <- paste(OUTDIR,"/supporting_plots/","SV_overlap_barplot.by_svtype.pdf",sep="")
  }
  pdf(pdf.path,height=5,width=7)
  plotStackedBars(svtype.mat,colors=ovr.cat.cols,scaled=T,legend=T,title=title)
  dev.off()
  
  #Master panel
  if(!is.null(prefix)){
    title.suffix <- paste(" vs. ",prefix,sep="")
    pdf(paste(OUTDIR,"main_plots/VCF_QC.",prefix,".SV_overlap_barplots.pdf",sep=""),
        height=3.5,width=12)
  }else{
    title.suffix <- ""
    pdf(paste(OUTDIR,"main_plots/VCF_QC.SV_overlap_barplots.pdf",sep=""),
        height=3.5,width=12)
  }
  layout(matrix(1:4,nrow=1,byrow=T),widths=c(5,7,5,1.5))
  plotStackedBars(svtype.mat,colors=ovr.cat.cols,scaled=T,legend=F,
                  title=paste("Overlap by Class",title.suffix,sep=""))
  plotStackedBars(size.mat,colors=ovr.cat.cols,scaled=T,legend=F,
                  title=paste("Overlap by Size",title.suffix,sep=""))
  plotStackedBars(freq.mat,colors=ovr.cat.cols,scaled=T,legend=F,
                  title=paste("Overlap by ",title.prefix,"Freq.",title.suffix,sep=""))
  #Legend
  par(mar=c(4,0.1,2,0.1),bty="n")
  plot(x=c(0,5),y=c(0,nrow(size.mat)),type="n",
       xaxs="i",yaxs="i",xlab="",ylab="",xaxt="n",yaxt="n")
  points(x=rep(0.4,times=nrow(size.mat)),
         y=(nrow(size.mat):1)-0.5,
         pch=15,col=rev(ovr.cat.cols),cex=2)
  text(x=rep(0.5,times=nrow(size.mat)),
       y=(nrow(size.mat):1)-0.5,
       labels=rev(rownames(size.mat)),pos=4,cex=0.75)
  dev.off()
}


##################
###AF CORRELATIONS
##################
#General function to make scatterplot from a pair of AF vectors
plotScatterSingle <- function(plot.df,lims=c(-0.01,1.01),colors=NULL,title=NULL,carrierFreqs=F,
                              col.lm=NULL,col.moving=NULL,legend=T,lab.cex=1){
  #Only plot if any values exist
  if(nrow(plot.df)>0 & ncol(plot.df)>0){
    #Order points & restrict to complete observations
    colnames(plot.df) <- c("x","y")
    plot.df <- plot.df[which(!is.na(plot.df[,1]) & !is.na(plot.df[,2])),]
    plot.df <- plot.df[with(plot.df,order(x,y)),]
    x <- plot.df[,1]
    y <- plot.df[,2]
    
    #Get axis limits
    if(is.null(lims)){
      lims <- range(c(0,x,y),na.rm=T)
    }
    #Get point & line colors
    if(is.null(colors)){
      colors <- rep("#062140",times=length(x))
    }
    if(is.null(col.lm)){
      col.lm <- "#FC4523"
    }
    if(is.null(col.moving)){
      col.moving <- "#B53119"
    }
    
    #Prep plot area
    par(mar=c(3.5,4,2.5,1))
    plot(x=lims,y=lims,type="n",
         xaxt="n",yaxt="n",xlab="",ylab="",xaxs="i",yaxs="i")
    
    #Add gridlines
    abline(h=seq(0,1,0.2),v=seq(0,1,0.2),col="gray93")
    abline(0,1,col="gray93")
    
    #Add axes
    axis(1,at=seq(0,1,0.2),labels=NA)
    axis(1,at=seq(0,1,0.2),tick=F,line=-0.4,cex.axis=0.8,
         labels=paste(seq(0,100,20),"%",sep=""))
    if(!is.null(prefix)){
      if(carrierFreqs==F){
        xlab <- paste("Allele Freq. (",prefix,")",sep="")
        ylab <- "Allele Freq. (This Callset)"
      }else{
        xlab <- paste("Carrier Freq. (",prefix,")",sep="")
        ylab <- "Carrier Freq. (This Callset)"
      }
    }else{
      if(carrierFreqs==F){
        xlab <- paste("Allele Freq. (Benchmarking Data)",sep="")
        ylab <- "Allele Freq. (This Callset)"
      }else{
        xlab <- paste("Carrier Freq. (Benchmarking Data)",sep="")
        ylab <- "Carrier Freq. (This Callset)"
      }    }
    mtext(1,text=xlab,line=2,cex=lab.cex)
    axis(2,at=seq(0,1,0.2),labels=NA)
    axis(2,at=seq(0,1,0.2),tick=F,line=-0.4,las=2,cex.axis=0.8,
         labels=paste(seq(0,100,20),"%",sep=""))
    mtext(2,text=ylab,line=2.5,cex=lab.cex)
    
    #Add points, linear fit, and moving average
    points(x,y,pch=21,col=colors,bg=adjustcolor(colors,alpha=0.2),cex=0.25)
    if(nrow(plot.df)>=100 & length(unique(x))>100){
      points(x=rollmean(x,k=ceiling(length(x)/100)),
             y=rollmean(y,k=ceiling(length(x)/100)),
             lwd=1.25,col=col.moving,type="l")
    }
    if(nrow(plot.df)>1){
     #abline(lm(y ~ x),lwd=2,col=col.lm)
      lm_data = lm(y ~ x)
      if(!is.na(is.na(lm_data$coefficients[1])) & !is.na(lm_data$coefficients[2])){
        abline(lm(y ~ x),lwd=2,col=col.lm)
      }
    }
    
    #Add correlation coefficients
    if(nrow(plot.df)>1){
      pval <- suppressWarnings(format(cor.test(x,y,method="spearman")$p.value,
                                      scientific=T,digits=2))
    }else{
      pval <- "NA"
    }
    mtext(3,line=0,cex=0.8*lab.cex,
          text=paste("N=",prettyNum(nrow(plot.df),big.mark=","),
                     "; R=",round(cor(x,y),2),
                     # ";  R^2=",round(cor(x,y)^2,3),
                     "; Rho=",round(cor(x,y,method="spearman"),2),
                     "; P=",pval,sep=""))
    
    
    #Add legend, if optioned
    if(legend==T){
      legend("right",pch=c(21,NA,NA),col=c(colors[1],col.lm,col.moving),
             pt.bg=c(adjustcolor(colors[1],alpha=0.25),NA,NA),
             lwd=c(NA,2,1.25),cex=0.7,
             legend=c("SV","Linear fit","Moving avg."),bg="white")
    }
    
    #Add cleanup box
    box()
  }else{
    par(bty="n",mar=c(3.5,3.5,3,0.5))
    plot(x=c(0,1),y=c(0,1),type="n",
         xaxt="n",yaxt="n",xlab="",ylab="",yaxs="i")
    text(x=0.5,y=0.5,labels="No Data")
  }
  
  #Add title
  mtext(3,text=title,font=2,line=1.25,cex=lab.cex)
}
#Master wrapper to generate all frequency correlations
wrapperAllFreqCors <- function(){
  #Set plotting prefixes & suffixes
  if(!is.null(prefix)){
    title.suffix <- paste("(",prefix,")",sep="")
    png.prefix <- paste(OUTDIR,"/supporting_plots/",prefix,".SV_freq_correlation.",sep="")
  }else{
    title.suffix <- ""
    png.prefix <- paste(OUTDIR,"/supporting_plots/SV_freq_correlation.",sep="")
  }
  if(carrierFreqs==F){
    title.prefix <- "Allele "
  }else{
    title.prefix <- "Carrier "
  }
  
  #All variants
  png(paste(png.prefix,"all_sv.png",sep=""),height=1500,res=300,width=1500)
  plotScatterSingle(plot.df=getFreqPairs(dat),carrierFreqs=carrierFreqs,
                    title=paste(title.prefix,"Freq. Correlation [All SV] ",title.suffix,sep=""))
  dev.off()
  
  #Tiny
  png(paste(png.prefix,"tiny_sv.png",sep=""),height=1500,res=300,width=1500)
  plotScatterSingle(plot.df=getFreqPairs(dat[which(dat$length<=tiny.max.size),]),
                    carrierFreqs=carrierFreqs,
                    title=paste(title.prefix,"Freq. Correlation [<100bp] ",title.suffix,sep=""))
  dev.off()
  #Small
  png(paste(png.prefix,"small_sv.png",sep=""),height=1500,res=300,width=1500)
  plotScatterSingle(plot.df=getFreqPairs(dat[which(dat$length>tiny.max.size & dat$length<=small.max.size),]),
                    carrierFreqs=carrierFreqs,
                    title=paste(title.prefix,"Freq. Correlation [100-500bp] ",title.suffix,sep=""))
  dev.off()
  #Medium
  png(paste(png.prefix,"medium_sv.png",sep=""),height=1500,res=300,width=1500)
  plotScatterSingle(plot.df=getFreqPairs(dat[which(dat$length>small.max.size & dat$length<=medium.max.size),]),
                    carrierFreqs=carrierFreqs,
                    title=paste(title.prefix,"Freq. Correlation [500bp-2.5kb] ",title.suffix,sep=""))
  dev.off()
  #Med-Large
  png(paste(png.prefix,"medlarge_sv.png",sep=""),height=1500,res=300,width=1500)
  plotScatterSingle(plot.df=getFreqPairs(dat[which(dat$length>medium.max.size & dat$length<=medlarge.max.size),]),
                    carrierFreqs=carrierFreqs,
                    title=paste(title.prefix,"Freq. Correlation [2.5-10kb] ",title.suffix,sep=""))
  dev.off()
  #Large
  png(paste(png.prefix,"large_sv.png",sep=""),height=1500,res=300,width=1500)
  plotScatterSingle(plot.df=getFreqPairs(dat[which(dat$length>medlarge.max.size & dat$length<=large.max.size),]),
                    carrierFreqs=carrierFreqs,
                    title=paste(title.prefix,"Freq. Correlation [10-50kb] ",title.suffix,sep=""))
  dev.off()
  #Huge
  png(paste(png.prefix,"huge_sv.png",sep=""),height=1500,res=300,width=1500)
  plotScatterSingle(plot.df=getFreqPairs(dat[which(dat$length>large.max.size & dat$length<=huge.max.size),]),
                    carrierFreqs=carrierFreqs,
                    title=paste(title.prefix,"Freq. Correlation [>50kb] ",title.suffix,sep=""))
  dev.off()
  
  #Iterate over variant classes & one plot per class
  sapply(unique(dat$svtype),function(svtype){
    plot.df <- getFreqPairs(dat[which(dat$svtype==svtype),])
    plot.df <- plot.df[which(!is.na(plot.df[,1]) & !is.na(plot.df[,2])),]
    if(nrow(plot.df)>0 & ncol(plot.df)>0){
      png(paste(png.prefix,svtype,".png",sep=""),height=1500,res=300,width=1500)
      plotScatterSingle(plot.df=plot.df,carrierFreqs=carrierFreqs,
                        title=paste(title.prefix,"Freq. Correlation [",svtype,"] ",title.suffix,sep=""))
      dev.off()
    }
  })
  
  #Master panel
  if(!is.null(prefix)){
    png(paste(OUTDIR,"main_plots/VCF_QC.",prefix,".SV_freq_correlations.png",sep=""),
        height=6*300,width=15*300,res=300)
  }else{
    png(paste(OUTDIR,"main_plots/VCF_QC.SV_freq_correlations.png",sep=""),
        height=6*300,width=15*300,res=300)
  }
  layout(matrix(c(1,2,3,4,
                  1,5,6,7),
                byrow=T,nrow=2),
         widths=c(2,1,1,1))
  plotScatterSingle(plot.df=getFreqPairs(dat),
                    carrierFreqs=carrierFreqs,
                    title=paste(title.prefix,"Freq. Correlation [All SV] ",title.suffix,sep=""))
  plotScatterSingle(plot.df=getFreqPairs(dat[which(dat$length<=tiny.max.size),]),
                    carrierFreqs=carrierFreqs,
                    title="<100bp",legend=F,lab.cex=0.8)
  plotScatterSingle(plot.df=getFreqPairs(dat[which(dat$length>tiny.max.size & dat$length<=small.max.size),]),
                    carrierFreqs=carrierFreqs,
                    title="100-500bp",legend=F,lab.cex=0.8)
  plotScatterSingle(plot.df=getFreqPairs(dat[which(dat$length>small.max.size & dat$length<=medium.max.size),]),
                    carrierFreqs=carrierFreqs,
                    title="500bp-2.5kb",legend=F,lab.cex=0.8)
  plotScatterSingle(plot.df=getFreqPairs(dat[which(dat$length>medium.max.size & dat$length<=medlarge.max.size),]),
                    carrierFreqs=carrierFreqs,
                    title="2.5-10kb",legend=F,lab.cex=0.8)
  plotScatterSingle(plot.df=getFreqPairs(dat[which(dat$length>medlarge.max.size & dat$length<=large.max.size),]),
                    carrierFreqs=carrierFreqs,
                    title="10-50kb",legend=F,lab.cex=0.8)
  plotScatterSingle(plot.df=getFreqPairs(dat[which(dat$length>large.max.size),]),
                    carrierFreqs=carrierFreqs,
                    title=">50kb",legend=F,lab.cex=0.8)
  dev.off()
}


######################
###PAIRWISE BREAKDOWNS
######################
#Breakdown data by size X freq
pairwiseBreakdownBySizeFreq <- function(dat,max.sizes,max.freqs,
                                        size.labels=NULL,freq.labels=NULL){
  #Create data frame of size cutoffs & group labels
  if(is.null(size.labels)){
    size.labels <- c("ALL",paste(c(0,max.sizes[-length(max.sizes)]),"-",max.sizes,sep=""))
  }else{
    size.labels <- c("ALL",size.labels)
  }
  sizes.df <- data.frame("min"=c(0,0,max.sizes[-length(max.sizes)]),
                         "max"=c(300000000,max.sizes),
                         "label"=size.labels)
  
  #Create data frame of size cutoffs & group labels
  if(is.null(freq.labels)){
    freq.labels <- c("ALL",paste(c(0,max.freqs[-length(max.freqs)]),"-",max.freqs,sep=""))
  }else{
    freq.labels <- c("ALL",freq.labels)
  }
  freqs.df <- data.frame("min"=c(0,0,max.freqs[-length(max.freqs)]),
                         "max"=c(1,max.freqs),
                         "label"=freq.labels)
  
  
  #Create matrix of all SV per size/freq threshold pair
  mat.all <- as.data.frame(t(sapply(1:nrow(sizes.df),function(r){
    row <- sapply(1:nrow(freqs.df),function(c){
      vals <- categoryBreakdown(dat[which(dat$length>sizes.df[r,1] & dat$length<=sizes.df[r,2] & 
                                            dat$AF>freqs.df[c,1] & dat$AF<=freqs.df[c,2]),],norm=F)
      return(sum(vals))
    })
    return(row)
  })))
  colnames(mat.all) <- freqs.df$label
  rownames(mat.all) <- sizes.df$label
  
  #Create matrix of overlapping SV per size/freq threshold pair
  mat.hit <- as.data.frame(t(sapply(1:nrow(sizes.df),function(r){
    row <- sapply(1:nrow(freqs.df),function(c){
      vals <- categoryBreakdown(dat[which(dat$length>sizes.df[r,1] & dat$length<=sizes.df[r,2] & 
                                            dat$AF>freqs.df[c,1] & dat$AF<=freqs.df[c,2]),],norm=F)
      return(sum(vals[1:3]))
    })
    return(row)
  })))
  colnames(mat.hit) <- freqs.df$label
  rownames(mat.hit) <- sizes.df$label
  
  #Create matrix of overlap fraction per size/freq threshold pair
  mat.pct <- mat.hit/mat.all
  
  #Return list of matrixes
  return(list(mat.all,mat.hit,mat.pct))
}
#Breakdown data by size X class
pairwiseBreakdownBySizeClass <- function(dat,max.sizes,size.labels=NULL){
  #Create data frame of size cutoffs & group labels
  if(is.null(size.labels)){
    size.labels <- c("ALL",paste(c(0,max.sizes[-length(max.sizes)]),"-",max.sizes,sep=""))
  }else{
    size.labels <- c("ALL",size.labels)
  }
  sizes.df <- data.frame("min"=c(0,0,max.sizes[-length(max.sizes)]),
                         "max"=c(300000000,max.sizes),
                         "label"=size.labels)
  
  #Get vector of SV classes
  svtypes <- sort(unique(dat$svtype))
  
  #Create matrix of all SV per size/class combo
  mat.all <- as.data.frame(t(sapply(1:nrow(sizes.df),function(r){
    row <- sapply(svtypes,function(svtype){
      vals <- categoryBreakdown(dat[which(dat$length>sizes.df[r,1] & dat$length<=sizes.df[r,2] & 
                                            dat$svtype==svtype),],norm=F)
      return(sum(vals))
    })
    return(row)
  })))
  #Add column corresponding to all svtypes
  ALL.all <- sapply(1:nrow(sizes.df),function(r){
    vals <- categoryBreakdown(dat[which(dat$length>sizes.df[r,1] & dat$length<=sizes.df[r,2]),],norm=F)
    return(sum(vals))
  })
  mat.all <- cbind(ALL.all,mat.all)
  rownames(mat.all) <- sizes.df$label
  colnames(mat.all) <- c("ALL",svtypes)
  
  #Create matrix of all SV per size/class combo
  mat.hit <- as.data.frame(t(sapply(1:nrow(sizes.df),function(r){
    row <- sapply(svtypes,function(svtype){
      vals <- categoryBreakdown(dat[which(dat$length>sizes.df[r,1] & dat$length<=sizes.df[r,2] & 
                                            dat$svtype==svtype),],norm=F)
      return(sum(vals[1:3]))
    })
    return(row)
  })))
  #Add column corresponding to all svtypes
  ALL.hit <- sapply(1:nrow(sizes.df),function(r){
    vals <- categoryBreakdown(dat[which(dat$length>sizes.df[r,1] & dat$length<=sizes.df[r,2]),],norm=F)
    return(sum(vals[1:3]))
  })
  mat.hit <- cbind(ALL.hit,mat.hit)
  rownames(mat.hit) <- sizes.df$label
  colnames(mat.hit) <- c("ALL",svtypes)
  
  #Create matrix of fraction of overlapping SV per size/class combo
  mat.pct <- mat.hit/mat.all
  
  return(list(mat.all,mat.hit,mat.pct))
}
#Breakdown data by freq X class
pairwiseBreakdownByFreqClass <- function(dat,max.freqs,freq.labels=NULL){
  #Create data frame of freq cutoffs & group labels
  if(is.null(freq.labels)){
    freq.labels <- c("ALL",paste(c(0,max.freqs[-length(max.freqs)]),"-",max.freqs,sep=""))
  }else{
    freq.labels <- c("ALL",freq.labels)
  }
  freqs.df <- data.frame("min"=c(0,0,max.freqs[-length(max.freqs)]),
                         "max"=c(1,max.freqs),
                         "label"=freq.labels)
  
  #Get vector of SV classes
  svtypes <- sort(unique(dat$svtype))
  
  #Create matrix of overlap categories per freq threshold
  mat.all <- as.data.frame(t(sapply(1:nrow(freqs.df),function(r){
    row <- sapply(svtypes,function(svtype){
      vals <- categoryBreakdown(dat[which(dat$AF>freqs.df[r,1] & dat$AF<=freqs.df[r,2] & 
                                            dat$svtype==svtype),],norm=F)
      return(sum(vals))
    })
    return(row)
  })))
  #Add column corresponding to all svtypes
  ALL.all <- sapply(1:nrow(freqs.df),function(r){
    vals <- categoryBreakdown(dat[which(dat$AF>freqs.df[r,1] & dat$AF<=freqs.df[r,2]),],norm=F)
    return(sum(vals))
  })
  mat.all <- cbind(ALL.all,mat.all)
  rownames(mat.all) <- freqs.df$label
  colnames(mat.all) <- c("ALL",svtypes)
  
  #Create matrix of all SV per freq/class combo
  mat.hit <- as.data.frame(t(sapply(1:nrow(freqs.df),function(r){
    row <- sapply(svtypes,function(svtype){
      vals <- categoryBreakdown(dat[which(dat$AF>freqs.df[r,1] & dat$AF<=freqs.df[r,2] & 
                                            dat$svtype==svtype),],norm=F)
      return(sum(vals[1:3]))
    })
    return(row)
  })))
  #Add column corresponding to all svtypes
  ALL.hit <- sapply(1:nrow(freqs.df),function(r){
    vals <- categoryBreakdown(dat[which(dat$AF>freqs.df[r,1] & dat$AF<=freqs.df[r,2]),],norm=F)
    return(sum(vals[1:3]))
  })
  mat.hit <- cbind(ALL.hit,mat.hit)
  rownames(mat.hit) <- freqs.df$label
  colnames(mat.hit) <- c("ALL",svtypes)
  
  #Create matrix of fraction of overlapping SV per freq/class combo
  mat.pct <- mat.hit/mat.all
  
  return(list(mat.all,mat.hit,mat.pct))
}
#Plot heatmap of fractions
plotHeatmap <- function(mats,col.range=NULL,text.col=NULL,
                        x.labels=NULL,x.title=NULL,
                        y.labels=NULL,y.title=NULL,
                        title=NULL,lab.cex=1){
  #Set values if NULL
  if(is.null(col.range)){
    col.range <- colorRampPalette(c("#440154","#365C8C","#25A584","#FDE725"))(101)
  }
  if(is.null(text.col)){
    text.col <- "gray95"
  }
  if(is.null(x.labels)){
    x.labels <- colnames(mats[[3]])
  }
  if(is.null(y.labels)){
    y.labels <- rownames(mats[[3]])
  }
  
  #Prep plotting area
  par(mar=c(4,4,2,2))
  plot(x=c(0,ncol(mats[[3]])),y=c(0,-nrow(mats[[3]])),type="n",
       xaxt="n",xaxs="i",xlab="",yaxt="n",yaxs="i",ylab="")
  
  #Add axes
  axis(1,at=(1:ncol(mats[[3]]))-0.5,tick=F,line=-0.8,las=2,labels=x.labels,cex.axis=0.7)
  mtext(1,line=2.75,text=x.title,cex=lab.cex)
  axis(2,at=-(1:nrow(mats[[3]]))+0.5,tick=F,line=-0.8,las=2,labels=y.labels,cex.axis=0.7)
  mtext(2,line=2.75,text=y.title,cex=lab.cex)
  mtext(3,line=0.5,text=title,font=2,cex=lab.cex)
  
  #Plot all cells
  sapply(1:ncol(mats[[3]]),function(c){
    sapply(1:nrow(mats[[3]]),function(r){
      #Get & scale value
      val <- round(100*mats[[3]][r,c],0)
      #Get color for shading
      if(mats[[1]][r,c]==0){
        color <- "gray80"
        pct <- "N/A"
        dens <- 12
      }else{
        color <- col.range[val+1]
        pct <- paste(val,"%",sep="")
        dens <- NA
      }
      #Flip text color for high percentages
      if(!is.na(val)){
        if(val>=75){
          text.col <- "gray25"
        }else{
          text.col <- "gray95"
        } 
      }
      #Plot rectangle
      rect(xleft=c-1,xright=c,ybottom=-r,ytop=-(r-1),
           lwd=0.5,border="gray95",col=color,density=dens)
      #Format cell annotation
      n.all <- mats[[1]][r,c]
      n.hit <- mats[[2]][r,c]
      anno.text <- paste(prettyNum(n.hit,big.mark=",")," /\n",
                         prettyNum(n.all,big.mark=","),"\n(",
                         pct,")",sep="")
      
      text(x=c-0.5,y=-(r-0.5),labels=anno.text,
           cex=0.7,col=text.col)
    })
  })
}
#Master wrapper plot for all heatmaps
wrapperAllHeatmaps <- function(){
  #Set plotting prefix
  if(!is.null(prefix)){
    title.suffix <- paste("(",prefix,")",sep="")
    pdf.prefix <- paste(OUTDIR,"/supporting_plots/",prefix,".SV_overlap_heatmap.",sep="")
  }else{
    title.suffix <- ""
    pdf.prefix <- paste(OUTDIR,"/supporting_plots/SV_overlap_heatmap.",sep="")
  }
  if(carrierFreqs==F){
    title.prefix <- "Allele "
  }else{
    title.prefix <- "Carrier "
  }
  
  #Size vs freq
  size_by_freq.mat <- pairwiseBreakdownBySizeFreq(dat=dat,
                                                  max.sizes=c(tiny.max.size,small.max.size,medium.max.size,
                                                              medlarge.max.size,large.max.size,huge.max.size),
                                                  size.labels=c("<100bp","100-\n500bp","500bp-\n2.5kb",
                                                                "2.5-\n10kb","10-\n50kb",">50kb"),
                                                  max.freqs=c(rare.max.freq,uncommon.max.freq,
                                                              common.max.freq,major.max.freq),
                                                  freq.labels=c("<1%","1-10%","10-50%",">50%"))
  pdf(paste(pdf.prefix,"size_vs_frequency.pdf",sep=""),height=5,width=5)
  plotHeatmap(mats=size_by_freq.mat,
              title=paste("SV Overlap, Size vs Freq. [All SV] ",title.suffix,sep=""),
              x.title=paste(title.prefix,"Freq.",sep=""),y.title="Size")
  dev.off()
  
  #Size vs class
  size_by_class.mat <- pairwiseBreakdownBySizeClass(dat=dat,
                                                    max.sizes=c(tiny.max.size,small.max.size,medium.max.size,
                                                                medlarge.max.size,large.max.size,huge.max.size),
                                                    size.labels=c("<100bp","100-\n500bp","500bp-\n2.5kb",
                                                                  "2.5-\n10kb","10-\n50kb",">50kb"))
  pdf(paste(pdf.prefix,"size_vs_svtype.pdf",sep=""),height=5,width=5)
  plotHeatmap(mats=size_by_class.mat,
              title=paste("SV Overlap, Size vs Class [All SV] ",title.suffix,sep=""),
              x.title="SV Class",y.title="Size")
  dev.off()
  
  #Freq vs class
  freq_by_class.mat <- pairwiseBreakdownByFreqClass(dat=dat,
                                                    max.freqs=c(rare.max.freq,uncommon.max.freq,
                                                                common.max.freq,major.max.freq),
                                                    freq.labels=c("<1%","1-10%","10-50%",">50%"))
  pdf(paste(pdf.prefix,"frequency_vs_svtype.pdf",sep=""),height=5,width=5)
  plotHeatmap(mats=freq_by_class.mat,
              title=paste("SV Overlap, Freq. vs Class [All SV] ",title.suffix,sep=""),
              x.title="SV Class",y.title=paste(title.prefix,"Freq.",sep=""))
  dev.off()
  
  #Iterate over variant classes & one size X freq heatmap per class
  sapply(unique(dat$svtype),function(svtype){
    plot.df <- pairwiseBreakdownBySizeFreq(dat=dat[which(dat$svtype==svtype),],
                                           max.sizes=c(tiny.max.size,small.max.size,medium.max.size,
                                                       medlarge.max.size,large.max.size,huge.max.size),
                                           size.labels=c("<100bp","100-\n500bp","500bp-\n2.5kb",
                                                         "2.5-\n10kb","10-\n50kb",">50kb"),
                                           max.freqs=c(rare.max.freq,uncommon.max.freq,
                                                       common.max.freq,major.max.freq),
                                           freq.labels=c("<1%","1-10%","10-50%",">50%"))
    pdf(paste(pdf.prefix,"size_vs_frequency.",svtype,".pdf",sep=""),height=5,width=5)
    plotHeatmap(mat=plot.df,
                title=paste("SV Overlap, Freq. vs Class [",svtype,"] ",title.suffix,sep=""),
                x.title="SV Class",y.title=paste(title.prefix,"Freq.",sep=""))
    dev.off()
  })
  
  #Master panel
  if(!is.null(prefix)){
    pdf(paste(OUTDIR,"main_plots/VCF_QC.",prefix,".SV_overlap_heatmaps.pdf",sep=""),
        height=5,width=15)
  }else{
    pdf(paste(OUTDIR,"main_plots/VCF_QC.SV_freq_correlations.pdf",sep=""),
        height=5,width=15)
  }
  par(mfrow=c(1,3))
  plotHeatmap(mats=size_by_freq.mat,
              title=paste("SV Overlap, Size vs Freq. [All SV] ",title.suffix,sep=""),
              x.title=paste(title.prefix,"Freq.",sep=""),y.title="Size",lab.cex=0.8)
  plotHeatmap(mats=size_by_class.mat,
              title=paste("SV Overlap, Size vs Class [All SV] ",title.suffix,sep=""),
              x.title="SV Class",y.title="Size",lab.cex=0.8)
  plotHeatmap(mats=freq_by_class.mat,
              title=paste("SV Overlap, Freq. vs Class [All SV] ",title.suffix,sep=""),
              x.title="SV Class",y.title=paste(title.prefix,"Freq.",sep=""),lab.cex=0.8)
  dev.off()
}


######################
###MASTER SUMMARY PLOT
######################
#Master wrapper to generate summary plot for whole analysis
masterWrapper <- function(){
  #Gather data required for barplots
  size.mat.bar <- categoryBreakdownBySize(dat,norm=F,
                                          max.sizes=c(tiny.max.size,small.max.size,medium.max.size,
                                                      medlarge.max.size,large.max.size,huge.max.size),
                                          size.labels=c("<100bp","100-500bp","500bp-2.5kb",
                                                        "2.5-10kb","10-50kb",">50kb"))
  freq.mat.bar <- categoryBreakdownByFreq(dat,norm=F,
                                          max.freqs=c(rare.max.freq,uncommon.max.freq,
                                                      common.max.freq,major.max.freq),
                                          freq.labels=c("<1%","1-10%","10-50%",">50%"))
  svtype.mat.bar <- categoryBreakdownByClass(dat,norm=F)
  
  #Gather data required for heatmaps
  size_by_freq.mat <- pairwiseBreakdownBySizeFreq(dat=dat,
                                                  max.sizes=c(tiny.max.size,small.max.size,medium.max.size,
                                                              medlarge.max.size,large.max.size,huge.max.size),
                                                  size.labels=c("<100bp","100-\n500bp","500bp-\n2.5kb",
                                                                "2.5-\n10kb","10-\n50kb",">50kb"),
                                                  max.freqs=c(rare.max.freq,uncommon.max.freq,
                                                              common.max.freq,major.max.freq),
                                                  freq.labels=c("<1%","1-10%","10-50%",">50%"))
  size_by_class.mat <- pairwiseBreakdownBySizeClass(dat=dat,
                                                    max.sizes=c(tiny.max.size,small.max.size,medium.max.size,
                                                                medlarge.max.size,large.max.size,huge.max.size),
                                                    size.labels=c("<100bp","100-\n500bp","500bp-\n2.5kb",
                                                                  "2.5-\n10kb","10-\n50kb",">50kb"))
  freq_by_class.mat <- pairwiseBreakdownByFreqClass(dat=dat,
                                                    max.freqs=c(rare.max.freq,uncommon.max.freq,
                                                                common.max.freq,major.max.freq),
                                                    freq.labels=c("<1%","1-10%","10-50%",">50%"))
  
  #Master panel
  if(!is.null(prefix)){
    title.suffix <- paste("(",prefix,")",sep="")
    png(paste(OUTDIR,"main_plots/VCF_QC.",prefix,".callset_benchmarking.png",sep=""),
        height=6*300,width=14*300,res=300)
  }else{
    title.suffix <- ""
    png(paste(OUTDIR,"main_plots/VCF_QC.callset_benchmarking.png",sep=""),
        height=6*300,width=14*300,res=300)
  }
  if(carrierFreqs==F){
    title.prefix <- "Allele "
  }else{
    title.prefix <- "Carrier "
  }
  layout(matrix(c(1,2,9,5,6,
                  1,2,4,5,6,
                  1,3,4,7,8,
                  1,3,9,7,8),
                byrow=T,nrow=4),
         widths=c(4,3,1,3,3))
  
  #Barplot by class
  plotStackedBars(svtype.mat.bar,colors=ovr.cat.cols,scaled=T,legend=F,
                  title=paste("SV Overlap by Class ",title.suffix,sep=""))
  
  #Barplot by size
  plotStackedBars(size.mat.bar,colors=ovr.cat.cols,scaled=T,legend=F,
                  title="Overlap by Size")
  
  #Barplot by frequency
  plotStackedBars(freq.mat.bar,colors=ovr.cat.cols,scaled=T,legend=F,
                  title=paste("Overlap by ",title.prefix,"Freq.",sep=""))
  
  #Barplot legend
  par(mar=c(4,0.1,2,0.1),bty="n")
  plot(x=c(0,5),y=c(0,nrow(size.mat.bar)),type="n",
       xaxs="i",yaxs="i",xlab="",ylab="",xaxt="n",yaxt="n")
  points(x=rep(0.4,times=nrow(size.mat.bar)),
         y=(nrow(size.mat.bar):1)-0.5,
         pch=15,col=rev(ovr.cat.cols),cex=2)
  text(x=rep(0.5,times=nrow(size.mat.bar)),
       y=(nrow(size.mat.bar):1)-0.5,
       labels=rev(rownames(size.mat.bar)),pos=4,cex=0.75)
  
  #Frequency correlation scatterplot
  plotScatterSingle(plot.df=getFreqPairs(dat),lab.cex=0.8,carrierFreqs=carrierFreqs,
                    title=paste(title.prefix,"Freq. Correlation",sep=""))
  
  #Heatmap size vs freq
  plotHeatmap(mats=size_by_freq.mat,
              title="Overlap [Size x Freq.]",
              x.title=paste(title.prefix,"Freq.",sep=""),y.title="Size",lab.cex=0.8)
  
  #Heatmap size vs class
  plotHeatmap(mats=size_by_class.mat,
              title="Overlap [Size x Class]",
              x.title="SV Class",y.title="Size",lab.cex=0.8)
  
  #Heatmap freq vs class
  plotHeatmap(mats=freq_by_class.mat,
              title="Overlap [Freq. x Class]",
              x.title="SV Class",y.title=paste(title.prefix,"Freq.",sep=""),lab.cex=0.8)
  
  #Close device
  dev.off()
}


########################
###RSCRIPT FUNCTIONALITY
########################
###Load libraries as needed
require(optparse)
require(RColorBrewer)
require(zoo)

###List of command-line options
option_list <- list(
  make_option(c("-p", "--prefix"), type="character", default=NULL,
              help="prefix to be appended to all outputs [default %default]",
              metavar="character"),
  make_option(c("-C", "--carrierFreqs"), type="logical", action="store_true", default=FALSE,
              help="flag to indicate use of carrier frequencies [default %default]",
              metavar="logical")
)

###Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog INFILE OUTDIR",
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
prefix <- opts$prefix
carrierFreqs <- opts$carrierFreqs

###Prepares I/O files
#Read & clean data
dat <- read.table(INFILE,comment.char="",sep="\t",header=T,check.names=F)
colnames(dat)[1] <- "chr"
#Create output directory structure, if necessary
if(!dir.exists(OUTDIR)){
  dir.create(OUTDIR)
}
if(!dir.exists(paste(OUTDIR,"main_plots/",sep=""))){
  dir.create(paste(OUTDIR,"main_plots/",sep=""))
}
if(!dir.exists(paste(OUTDIR,"supporting_plots/",sep=""))){
  dir.create(paste(OUTDIR,"supporting_plots/",sep=""))
}


###Plotting block
#Summary plot
masterWrapper()
#Barplots
wrapperAllBarplots()
#Frequency correlations
wrapperAllFreqCors()
#Pairwise heatmaps
wrapperAllHeatmaps()
