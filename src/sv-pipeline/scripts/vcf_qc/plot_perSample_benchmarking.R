#!/usr/bin/env Rscript

# Helper script to plot per-sample SV benchmarking vs. an external dataset


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
#Read overlap data for a list of samples
readMultiSampDat <- function(samples, measurement="sensitivity"){
  #Iterate over samples
  dat <- lapply(samples,function(ID){
    #Set path
    path <- paste(perSampDir,"/",ID,".",measurement,".bed.gz",sep="")
    #Read & process data if file exists
    if(file.exists(path)){
      #Read & clean data
      x <- read.table(path,header=T,comment.char="",check.names=F)
      colnames(x)[1] <- "chr"
      
      #Return data
      return(x)
    }else{
      warning(paste("VID file not found for sample ",ID,sep=""))
      return(NULL)
    }
  })
  
  #Clean list of overlaps & exclude NULL samples, and returns
  names(dat) <- samples
  dat <- dat[which(!unlist(lapply(dat,is.null)))]
  return(dat)
}
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
  if(sum(c(a12.n,a1.n,a2.n,b12.n,b1.n,b2.n,ab3.n,f.n))>0){
    ovr <- (a12.n+a1.n+a2.n)/sum(c(a12.n,a1.n,a2.n,b12.n,b1.n,b2.n,ab3.n,f.n))
  }else{
    ovr <- NA
  }
  
  #Return vector
  out.vect <- c(a12.n,a1.n,a2.n,b12.n,b1.n,b2.n,ab3.n,f.n)
  if(norm==T){
    out.vect <- out.vect/nrow(dat)
  }
  out.vect <- c(out.vect,ovr)
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
                     "Nearby\n(+/- 250bp)","No Overlap","Frac")  
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
                     "Nearby\n(+/- 250bp)","No Overlap","Frac")  
  return(mat)
}
#Breakdown data by SV class
categoryBreakdownByClass <- function(dat,norm=F){
  #Create normalized vector of overlap categories per sv class
  mat <- as.data.frame(sapply(svtypes$svtype,function(svtype){
    categoryBreakdown(dat[which(dat$svtype==svtype),],norm=norm)
  }))
  mat <- cbind(categoryBreakdown(dat,norm=norm),mat)
  colnames(mat) <- c("ALL",svtypes$svtype)
  rownames(mat) <- c("Recip. & Bkpt.","Recip.","Bkpt.",
                     "Recip. & Bkpt.\n(Diff. SV types)",
                     "Recip.\n(Diff. SV types)",
                     "Bkpt.\n(Diff. SV types)",
                     "Nearby\n(+/- 250bp)","No Overlap","Frac")
  return(mat)
}
#Extract best matching AF pairs
getFreqPairs <- function(dat){
  freqPairs <- as.data.frame(t(apply(dat[,7:ncol(dat)],1,function(vals){
    benchmark.AF <- as.numeric(vals[1])
    callset.AF <- NULL
    # for(i in c(2,4,3,5,6)){
    for(i in c(4,2)){
      if(is.null(callset.AF)){
        if(!(vals[i]) %in% c("NS","NO_OVR")){
          callset.AF <- as.numeric(vals[i])
        }
      }
    }
    if(is.null(callset.AF)){
      callset.AF <- NA
    }
    return(c(benchmark.AF,callset.AF))
  })))
}
#Create matrix of Frac values per sample given list of per-sample overlap summaries
makeFracMat <- function(ovrlist){
  #Get matrix of fractions per sample
  fracMat <- as.data.frame(matrix(unlist(lapply(ovrlist,function(df){
    fracs <- df[which(rownames(df)=="Frac"),]
  })),nrow=length(ovrlist),byrow=T))
  
  #Clean & return matrix
  colnames(fracMat) <- colnames(ovrlist[[1]])
  rownames(fracMat) <- names(ovrlist)
  return(fracMat)
}
#Compute variance of weighted mean
#Taken from: https://stats.stackexchange.com/questions/25895/computing-standard-error-in-weighted-mean-estimation
weighted.var.se <- function(x,w,na.rm=F){
  if (na.rm) { w <- w[i <- !is.na(x)]; x <- x[i] }
  n = length(w)
  xWbar = weighted.mean(x,w,na.rm=na.rm)
  wbar = mean(w)
  out = n/((n-1)*sum(w)^2)*(sum((w*x-wbar*xWbar)^2)-2*xWbar*sum((w-wbar)*(w*x-wbar*xWbar))+xWbar^2*sum((w-wbar)^2))
  return(out)
}
#Summarize per-sample fracs into two smaller matrixes of mean and 95% CI
#Note: weight by number of SV per genome
makeMeanStderrMats <- function(ovrlist){
  #Get values per sample
  means.noweight <- makeFracMat(ovrlist)
  
  #Get weights per sample
  weights <- as.data.frame(matrix(unlist(lapply(ovrlist,function(df){
    return(apply(df[-which(rownames(df)=="Frac"),],2,sum,na.rm=T))
  })),nrow=length(ovrlist),byrow=T))
  
  #Get weighted means
  means.weighted <- sapply(1:ncol(means.noweight),function(i){
    wm <- weighted.mean(x=means.noweight[,i],
                        w=weights[,i],
                        na.rm=T)
    if(is.nan(wm)){
      wm <- NA
    }
    return(wm)
  })
  
  #Get weighted 95% CI adjustment
  ci_adjustment.weighted <- sapply(1:ncol(means.noweight),function(i){
    wsd <- sqrt(weighted.var.se(x=means.noweight[,i],
                                w=weights[,i],
                                na.rm=T))
    if(is.nan(wsd)){
      wsd <- NA
    }else{
      1.96*(wsd/sqrt(nrow(means.noweight)))
    }
    return(wsd)
  })
  
  #Return vectors as list
  return(list("mean"=means.weighted,
              "ci_adj"=ci_adjustment.weighted))
}


######################
###PAIRWISE BREAKDOWNS
######################
#Breakdown data by size X freq
pairwiseBreakdownBySizeFreq <- function(dat,max.sizes,size.labels=NULL,
                                        max.freqs,freq.labels=NULL){
  #Create data frame of size cutoffs & group labels
  if(is.null(size.labels)){
    size.labels <- c("ALL",paste(c(0,max.sizes[-length(max.sizes)]),"-",max.sizes,sep=""))
  }else{
    size.labels <- c("ALL",size.labels)
  }
  sizes.df <- data.frame("min"=c(0,0,max.sizes[-length(max.sizes)]),
                         "max"=c(300000000,max.sizes),
                         "label"=size.labels)
  
  #Create data frame of freq cutoffs & group labels
  if(is.null(freq.labels)){
    freq.labels <- c("ALL",paste(c(0,max.freqs[-length(max.freqs)]),"-",max.freqs,sep=""))
  }else{
    freq.labels <- c("ALL",freq.labels)
  }
  freqs.df <- data.frame("min"=c(0,0,max.freqs[-length(max.freqs)]),
                         "max"=c(1,max.freqs),
                         "label"=freq.labels)
  
  ###Gather mean overlap stats per freq bin
  ovrstats <- lapply(1:nrow(freqs.df),function(f){
    #Iterate per sample & compute overlap matrix list after restricting on size
    ovrlist <- lapply(dat,function(df){
      samp.mat <- sapply(1:nrow(sizes.df),function(r){
        categoryBreakdown(df[which(df$length>sizes.df[r,1] & df$length<=sizes.df[r,2] &
                                     df$AF>freqs.df[f,1] & df$AF<=freqs.df[f,2]),],norm=F)
      })
      colnames(samp.mat) <- sizes.df$label
      rownames(samp.mat) <- c("Recip. & Bkpt.","Recip.","Bkpt.",
                              "Recip. & Bkpt.\n(Diff. SV types)",
                              "Recip.\n(Diff. SV types)",
                              "Bkpt.\n(Diff. SV types)",
                              "Nearby\n(+/- 250bp)","No Overlap","Frac")
      return(samp.mat)
    })
    
    #Get weighted mean for this sv type
    makeMeanStderrMats(ovrlist)
  })
  
  #Create matrix of weighted means (col: sizes, row: freqs)
  means <- matrix(unlist(lapply(ovrstats,function(vects){return(vects$mean)})),
                  nrow=nrow(freqs.df),byrow=T)
  colnames(means) <- sizes.df$label
  colnames(means)[3:4] <- c("100-\n500bp","500bp\n-2.5kb")
  rownames(means) <- freqs.df$label
  
  #Create matrix of CI adjustments
  ci_adj <- matrix(unlist(lapply(ovrstats,function(vects){return(vects$ci_adj)})),
                   nrow=nrow(freqs.df),byrow=T)
  colnames(ci_adj) <- sizes.df$label
  colnames(ci_adj)[3:4] <- c("100-\n500bp","500bp\n-2.5kb")
  rownames(ci_adj) <- freqs.df$label
  
  #Format & return output list of two matrixes
  return(list("mean"=means,"ci_adj"=ci_adj))
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
  
  ###Gather mean overlap stats per svtype
  #First for ALL SV (simplest case)
  allsv.stats <- makeMeanStderrMats(lapply(dat,categoryBreakdownBySize,
                                           max.sizes=sizes.df$max[-1],
                                           size.labels=sizes.df$label[-1],
                                           norm=F))
  #Next for each svtype individually (more complicated)
  ovrstats <- lapply(svtypes$svtype,function(svtype){
    #Iterate per sample & compute overlap matrix list after restricting on size
    ovrlist <- lapply(dat,function(df){
      samp.mat <- sapply(1:nrow(sizes.df),function(r){
        categoryBreakdown(df[which(df$length>sizes.df[r,1] & df$length<=sizes.df[r,2] &
                                     df$svtype==svtype),],norm=F)
      })
      colnames(samp.mat) <- sizes.df$label
      rownames(samp.mat) <- c("Recip. & Bkpt.","Recip.","Bkpt.",
                              "Recip. & Bkpt.\n(Diff. SV types)",
                              "Recip.\n(Diff. SV types)",
                              "Bkpt.\n(Diff. SV types)",
                              "Nearby\n(+/- 250bp)","No Overlap","Frac")
      return(samp.mat)
    })
    
    #Get weighted mean for this sv type
    makeMeanStderrMats(ovrlist)
  })
  
  #Create matrix of weighted means (col: sizes, row: svtypes)
  means <- matrix(unlist(lapply(ovrstats,function(vects){return(vects$mean)})),
                  nrow=nrow(svtypes),byrow=T)
  means <- rbind(allsv.stats$mean,means)
  colnames(means) <- sizes.df$label
  colnames(means)[3:4] <- c("100-\n500bp","500bp\n-2.5kb")
  rownames(means) <- c("ALL",svtypes$svtype)
  
  #Create matrix of CI adjustments
  ci_adj <- matrix(unlist(lapply(ovrstats,function(vects){return(vects$ci_adj)})),
                   nrow=nrow(svtypes),byrow=T)
  ci_adj <- rbind(allsv.stats$ci_adj,ci_adj)
  colnames(ci_adj) <- sizes.df$label
  colnames(ci_adj)[3:4] <- c("100-\n500bp","500bp\n-2.5kb")
  rownames(ci_adj) <- c("ALL",svtypes$svtype)
  
  #Format & return output list of two matrixes
  return(list("mean"=means,"ci_adj"=ci_adj))
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
  
  ###Gather mean overlap stats per svtype
  #First for ALL SV (simplest case)
  allsv.stats <- makeMeanStderrMats(lapply(dat,categoryBreakdownByFreq,
                                           max.freqs=freqs.df$max[-1],
                                           freq.labels=freqs.df$label[-1],
                                           norm=F))
  #Next for each svtype individually (more complicated)
  ovrstats <- lapply(svtypes$svtype,function(svtype){
    #Iterate per sample & compute overlap matrix list after restricting on freq
    ovrlist <- lapply(dat,function(df){
      samp.mat <- sapply(1:nrow(freqs.df),function(r){
        categoryBreakdown(df[which(df$AF>freqs.df[r,1] & df$AF<=freqs.df[r,2] &
                                     df$svtype==svtype),],norm=F)
      })
      colnames(samp.mat) <- freqs.df$label
      rownames(samp.mat) <- c("Recip. & Bkpt.","Recip.","Bkpt.",
                              "Recip. & Bkpt.\n(Diff. SV types)",
                              "Recip.\n(Diff. SV types)",
                              "Bkpt.\n(Diff. SV types)",
                              "Nearby\n(+/- 250bp)","No Overlap","Frac")
      return(samp.mat)
    })
    
    #Get weighted mean for this sv type
    makeMeanStderrMats(ovrlist)
  })
  
  #Create matrix of weighted means (col: freqs, row: svtypes)
  means <- matrix(unlist(lapply(ovrstats,function(vects){return(vects$mean)})),
                  nrow=nrow(svtypes),byrow=T)
  means <- rbind(allsv.stats$mean,means)
  colnames(means) <- freqs.df$label
  rownames(means) <- c("ALL",svtypes$svtype)
  
  #Create matrix of CI adjustments
  ci_adj <- matrix(unlist(lapply(ovrstats,function(vects){return(vects$ci_adj)})),
                   nrow=nrow(svtypes),byrow=T)
  ci_adj <- rbind(allsv.stats$ci_adj,ci_adj)
  colnames(ci_adj) <- freqs.df$label
  rownames(ci_adj) <- c("ALL",svtypes$svtype)
  
  #Format & return output list of two matrixes
  return(list("mean"=means,"ci_adj"=ci_adj))
}



############################
###PLOTTING HELPER FUNCTIONS
############################
#Violin plot of fractions on y axis
plotViolins <- function(mat,colors=NULL,log=F,ylims=c(0,1),
                        xlab=NULL,ylab=NULL,title=NULL,
                        gridlines=T,violin=NULL,lab.cex=1){
  
  #Log-transform matrix, if optioned
  if(log==T){
    mat <- apply(mat,2,log10)
    y.suffix <- " (log10-scaled)"
  }else{
    y.suffix <- NULL
  }
  
  #Set arguments as needed
  if(is.null(colors)){
    colors <- rep("gray60",times=ncol(mat))
  }else{
    colors <- c("gray35",colors)
  }
  if(is.null(ylims)){
    ylims <- range(mat)
  }
  #If more than 1000 samples, use violin plots
  if(is.null(violin)){
    if(nrow(mat)>1000){
      violin <- T
    }else{
      violin <- F
    }
  }
  
  #Prep plotting area
  par(mar=c(4.5,4,2,1),bty="n")
  plot(x=c(0,ncol(mat)),y=ylims,type="n",
       xaxt="n",xaxs="i",yaxt="n",xlab="",ylab="")
  
  #Add gridlines, if optioned
  if(gridlines==T){
    abline(h=seq(0,1,0.2),col="gray92")
  }
  
  #Add axes & labels
  mtext(1,line=3-max(c(0,(1-lab.cex))),text=xlab,cex=lab.cex)
  axis(2,at=seq(0,1,0.2),labels=NA)
  axis(2,at=seq(0,1,0.2),tick=F,line=-0.4,cex.axis=0.8,las=2,
       labels=paste(seq(0,100,20),"%",sep=""))
  mtext(2,line=2.75-max(c(0,(1-lab.cex))),text=ylab,cex=lab.cex)
  mtext(3,line=0,cex=0.7*lab.cex,
        text=paste("N=",prettyNum(nrow(mat),big.mark=",")," Samples",sep=""))
  mtext(3,line=0.8,text=paste(title,y.suffix,sep=""),font=2,cex=lab.cex)
  
  #Iterate over columns per class and plot distributions
  sapply(1:ncol(mat),function(i){
    #Get values
    vals <- mat[,i]
    
    #Only plot distribution if any values are non-zero
    if(any(vals>0 & !is.infinite(as.vector(vals)) & !is.na(as.vector(vals)))){
      if(violin==T){
        #Get inliers & outliers
        IQR <- IQR(vals,na.rm=T)
        median <- median(vals,na.rm=T)
        inliers <- vals[which(vals>=median-(3*IQR) & vals<=median+(3*IQR))]
        outliers <- vals[which(vals<=median-(3*IQR) | vals>=median+(3*IQR))]
        
        #Plot violin of inliers
        if(length(unique(inliers))>1){
          vioplot(inliers,add=T,at=i-0.5,wex=0.35,
                  col=colors[i],drawRect=F,border=NA)
        }else{
          points(x=i-0.5,y=unique(inliers),pch=18,cex=1.5,col=colors[i])
        }
        
        #Plot points for outliers
        points(x=jitter(rep(i-0.5,length(outliers)),amount=0.01),
               y=outliers,pch=21,cex=0.175,lwd=0.7,
               col=colors[i],bg=adjustcolor(colors[i],alpha=0.3))
        #Otherwise, use swarmplots
      }else{
        beeswarm(vals,add=T,at=i-0.5,pch=21,cex=0.25,
                 col=colors[i],bg=adjustcolor(colors[i],alpha=0.3),
                 corral="wrap",corralWidth=0.35)
      }
      
      #Add bar for median & vertical rule for IQR
      segments(x0=i-0.6,x1=i-0.4,lwd=2,
               y0=median(vals,na.rm=T))
      segments(x0=i-0.5,x1=i-0.5,
               y0=quantile(vals,probs=0.25,na.rm=T),
               y1=quantile(vals,probs=0.75,na.rm=T))
      
      #Add mean on x-axis
      axis(1,at=i-0.3,
           labels=bquote(mu == .(paste(round(100*mean(vals,na.rm=T),1),"%",sep=""))),
           tick=F,las=2,line=-0.8,cex.axis=0.6) 
    }
    
    #Add category label on x-axis
    axis(1,at=i-0.5,
         labels=paste(colnames(mat)[i],"\n",sep=""),
         tick=F,las=2,line=-0.8,cex.axis=0.8,
         col.axis=colors[i],font=2)
  })
  
  #Iterate over columns and add text labels for median values
  sapply(1:ncol(mat),function(i){
    vals <- mat[,i]
    
    #Only plot distribution if any values are non-zero
    if(any(vals>0 & !is.infinite(as.vector(vals)) & !is.na(as.vector(vals)))){
      vals <- mat[,i]
      text(x=i-0.48,y=median(vals,na.rm=T),pos=4,cex=0.5,font=4,
           col=colors[i],labels=paste(round(100*median(vals,na.rm=T),1),"%",sep=""))
    }
  })
  
  #Clean up boxes
  rect(xleft=par("usr")[1],xright=par("usr")[2],
       ybottom=par("usr")[3],ytop=par("usr")[4],
       col=NA)
}
#Lineplots of fractions on y-axis split by SV class
plotLinesByClass <- function(means,ci_adj,nsamp,xlab=NULL,ylab=NULL,title=NULL,
                             legend=T,gridlines=T,lab.cex=1){
  #Prep plot area
  par(mar=c(4.5,4,2,1),bty="n")
  plot(x=c(0,ncol(means)),y=c(0,1),type="n",
       xaxt="n",yaxt="n",xlab="",ylab="",yaxs="i")
  
  #Set type colors & midpoints
  colors <- c("gray15",svtypes$color)
  lwds <- c(3,rep(2,times=nrow(svtypes)))
  mids <- (1:ncol(means))-0.5
  
  #Add gridlines, if optioned
  if(gridlines==T){
    abline(h=seq(0,1,0.2),col="gray92")
  }
  
  #Add axes & labels
  mtext(1,line=3-max(c(0,(1-lab.cex))),text=xlab,cex=lab.cex)
  axis(1,at=mids,line=-0.8,tick=F,labels=colnames(means),las=2,cex.axis=0.8)
  axis(2,at=seq(0,1,0.2),labels=NA)
  axis(2,at=seq(0,1,0.2),tick=F,line=-0.4,cex.axis=0.8,las=2,
       labels=paste(seq(0,100,20),"%",sep=""))
  mtext(2,line=2.75-max(c(0,(1-lab.cex))),text=ylab,cex=lab.cex)
  mtext(3,line=0,cex=0.7*lab.cex,
        text=paste("Weighted Mean of N=",prettyNum(nsamp,big.mark=",")," Samples",sep=""))
  mtext(3,line=0.8,text=title,font=2,cex=lab.cex)
  
  # # Iterate over classes and plot background shading for 95% CI
  #Add background shading for height of mean 
  sapply(nrow(means):1,function(i){
    #Get values
    vals.upper <- as.numeric(means[i,-1])
    vals.lower <- rep(0,times=length(vals.upper))
    # vals.upper <- as.numeric(means[i,-1])+as.numeric(ci_adj[i,-1])
    
    #Plot shading
    polygon(x=c(mids[-1][which(!is.na(means[i,-1]))],
                rev(mids[-1][which(!is.na(means[i,-1]))])),
            y=c(vals.lower[which(!is.na(means[i,-1]))],
                rev(vals.upper[which(!is.na(means[i,-1]))])),
            bty="n",border=NA,col=adjustcolor(colors[i],alpha=0.12))
  })
  
  #Iterate over classes and plot solid lines for means
  sapply(nrow(means):1,function(i){
    #Get values
    vals <- as.numeric(means[i,])
    
    #Plot lines of all but the first column
    points(x=mids[-1],y=vals[-1],type="l",lwd=lwds[i],col=colors[i])
    
    #Plot solid points of all
    points(x=mids,y=vals,pch=19,cex=0.5*lwds[i],col=colors[i])
  })
  
  #Add legend
  if(legend==T){
    idx.for.legend <- which(apply(means,1,function(vals){any(!is.na(vals))}))
    if(length(idx.for.legend) > 0){
      legend("bottomleft",bg="white",pch=19,cex=0.7*lab.cex,lwd=2,
             legend=rownames(means)[idx.for.legend],
             col=colors[idx.for.legend])
    }
  }
  
  #Add cleanup boxes
  rect(xleft=c(par("usr")[1],0.8,ncol(means)-0.2),
       xright=c(0.2,1.2,par("usr")[2]),
       ybottom=par("usr")[3],ytop=par("usr")[4],border="white",col="white")
  rect(xleft=0.2,xright=0.8,ybottom=par("usr")[3],ytop=par("usr")[4],col=NA)
  rect(xleft=1.2,xright=ncol(means)-0.2,ybottom=par("usr")[3],ytop=par("usr")[4],col=NA)
  axis(2,at=c(0,1),labels=NA,tck=0)
}
#Plot heatmap of means
plotHeatmap <- function(mat,nsamp,col.range=NULL,text.col=NULL,
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
    x.labels <- colnames(mat)
  }
  if(is.null(y.labels)){
    y.labels <- rownames(mat)
  }
  
  #Prep plotting area
  par(mar=c(4,4,2,2))
  plot(x=c(0,ncol(mat)),y=c(0,-nrow(mat)),type="n",
       xaxt="n",xaxs="i",xlab="",yaxt="n",yaxs="i",ylab="")
  
  #Add axes
  axis(1,at=(1:ncol(mat))-0.5,tick=F,line=-0.8,las=2,labels=x.labels,cex.axis=0.7)
  mtext(1,line=2.75,text=x.title,cex=lab.cex)
  axis(2,at=-(1:nrow(mat))+0.5,tick=F,line=-0.8,las=2,labels=y.labels,cex.axis=0.7)
  mtext(2,line=2.75,text=y.title,cex=lab.cex)
  mtext(3,line=0,cex=0.7*lab.cex,
        text=paste("Weighted Mean of N=",prettyNum(nsamp,big.mark=",")," Samples",sep=""))
  mtext(3,line=0.8,text=title,font=2,cex=lab.cex)
  
  #Plot all cells
  sapply(1:ncol(mat),function(c){
    sapply(1:nrow(mat),function(r){
      #Get & scale value
      val <- round(100*mat[r,c],0)
      #Get color for shading
      if(is.na(mat[r,c]) | mat[r,c]==0){
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
      text(x=c-0.5,y=-(r-0.5),labels=pct,
           cex=0.75,col=text.col)
    })
  })
}


############################
###MASTER PLOTTING FUNCTIONS
############################
#Compute master plotting data list from list of per-sample overlaps
gatherPlottingData <- function(dat){
  ###One-variable breakdowns
  #Breakdown by size
  datBySize <- lapply(dat,categoryBreakdownBySize,norm=F,
                      max.sizes=c(tiny.max.size,small.max.size,medium.max.size,
                                  medlarge.max.size,large.max.size,huge.max.size),
                      size.labels=c("<100bp","100-500bp","500bp-2.5kb",
                                    "2.5-10kb","10-50kb",">50kb"))
  datBySize.mat <- makeFracMat(datBySize)
  #Breakdown by freq
  datByFreq <- lapply(dat,categoryBreakdownByFreq,norm=F,
                      max.freqs=c(rare.max.freq,uncommon.max.freq,
                                  common.max.freq,major.max.freq),
                      freq.labels=c("<1%","1-10%","10-50%",">50%"))
  datByFreq.mat <- makeFracMat(datByFreq)
  #Breakdown by class
  datByClass <- lapply(dat,categoryBreakdownByClass,norm=F)
  datByClass.mat <- makeFracMat(datByClass)
  
  ###Two-variable breakdowns
  datBySizeClass.mat <- pairwiseBreakdownBySizeClass(dat=dat,
                                                     max.sizes=c(tiny.max.size,small.max.size,medium.max.size,
                                                                 medlarge.max.size,large.max.size,huge.max.size),
                                                     size.labels=c("<100bp","100-500bp","500bp-2.5kb",
                                                                   "2.5-10kb","10-50kb",">50kb"))
  datByFreqClass.mat <- pairwiseBreakdownByFreqClass(dat=dat,
                                                     max.freqs=c(rare.max.freq,uncommon.max.freq,
                                                                 common.max.freq,major.max.freq),
                                                     freq.labels=c("<1%","1-10%","10-50%",">50%"))
  datBySizeFreq.mat <- pairwiseBreakdownBySizeFreq(dat=dat,
                                                   max.sizes=c(tiny.max.size,small.max.size,medium.max.size,
                                                               medlarge.max.size,large.max.size,huge.max.size),
                                                   size.labels=c("<100bp","100-500bp","500bp-2.5kb",
                                                                 "2.5-10kb","10-50kb",">50kb"),
                                                   max.freqs=c(rare.max.freq,uncommon.max.freq,
                                                               common.max.freq,major.max.freq),
                                                   freq.labels=c("<1%","1-10%","10-50%",">50%"))
  
  
  ###Format & return list
  plotDat <- list("datBySize"=datBySize,
                  "datBySize.mat"=datBySize.mat,
                  "datByFreq"=datByFreq,
                  "datByFreq.mat"=datByFreq.mat,
                  "datByClass"=datByClass,
                  "datByClass.mat"=datByClass.mat,
                  "datBySizeClass.mat"=datBySizeClass.mat,
                  "datByFreqClass.mat"=datByFreqClass.mat,
                  "datBySizeFreq.mat"=datBySizeFreq.mat)
  return(plotDat)
}
#Wrapper for all violin plots
wrapperViolinPlots <- function(plotDat,compset.prefix,ylab=NULL){
  #Violin by class
  png(paste(OUTDIR,"/supporting_plots/per_sample_benchmarking_",compset.prefix,"/",
            ylab,"_per_genome_by_class.png",sep=""),
      height=1200,width=1500,res=300)
  plotViolins(mat=plotDat$datByClass.mat,
              colors=svtypes$color,
              xlab="SV Class",ylab=ylab,
              title=paste(ylab," by Class vs. ",compset.prefix,sep=""))
  dev.off()
  
  #Violin by size
  png(paste(OUTDIR,"/supporting_plots/per_sample_benchmarking_",compset.prefix,"/",
            ylab,"_per_genome_by_size.png",sep=""),
      height=1200,width=1500,res=300)
  colnames(plotDat$datBySize.mat)[3:4] <- c("100-\n500bp","500bp\n-2.5kb")
  plotViolins(mat=plotDat$datBySize.mat,
              colors=colorRampPalette(c("#440154","#365C8C","#25A584","#FDE725"))(ncol(plotDat$datBySize.mat)-1),
              xlab="SV Size",ylab=ylab,
              title=paste(ylab," by Size vs. ",compset.prefix,sep=""))
  dev.off()
  
  #Violin by freq
  png(paste(OUTDIR,"/supporting_plots/per_sample_benchmarking_",compset.prefix,"/",
            ylab,"_per_genome_by_freq.png",sep=""),
      height=1200,width=1500,res=300)
  # colnames(plotDat$datByFreq.mat)[3:4] <- c("100-\n500bp","500bp\n-2.5kb")
  plotViolins(mat=plotDat$datByFreq.mat,
              colors=rev(colorRampPalette(c("#440154","#365C8C","#25A584","#FDE725"))(ncol(plotDat$datByFreq.mat)-1)),
              xlab="SV Frequency",ylab=ylab,
              title=paste(ylab," by Freq. vs. ",compset.prefix,sep=""))
  dev.off()
}
#Wrapper for all lineplots
wrapperLinePlots <- function(plotDat,compset.prefix,ylab=NULL,nsamp){
  #Size x class
  pdf(paste(OUTDIR,"/supporting_plots/per_sample_benchmarking_",compset.prefix,"/",
            "mean_",ylab,"_per_genome_by_size_class.pdf",sep=""),
      height=4,width=6)
  plotLinesByClass(means=plotDat$datBySizeClass.mat$mean,
                   ci_adj=plotDat$datBySizeClass.mat$ci_adj,
                   nsamp=nsamp,xlab="SV Size",ylab=ylab,
                   title=paste("Mean ",ylab," by Size vs. ",compset.prefix,sep=""))
  dev.off()
  
  #Size x class
  pdf(paste(OUTDIR,"/supporting_plots/per_sample_benchmarking_",compset.prefix,"/",
            "mean_",ylab,"_per_genome_by_freq_class.pdf",sep=""),
      height=4,width=6)
  plotLinesByClass(means=plotDat$datByFreqClass.mat$mean,
                   ci_adj=plotDat$datByFreqClass.mat$ci_adj,
                   nsamp=nsamp,xlab="SV Frequency",ylab=ylab,
                   title=paste("Mean ",ylab," by Freq. vs. ",compset.prefix,sep=""))
  dev.off()
  
  #HEATMAP of Size x class
  pdf(paste(OUTDIR,"/supporting_plots/per_sample_benchmarking_",compset.prefix,"/",
            "mean_",ylab,"_per_genome_by_size_freq.pdf",sep=""),
      height=4,width=5)
  plotHeatmap(mat=plotDat$datBySizeFreq.mat$mean,x.title="SV Size",y.title="SV Freq.",
              title=paste("Mean ",ylab," by Size x Freq. vs. ",compset.prefix,sep=""),
              nsamp=nsamp)
  dev.off()
}
#Master plot wrapper
masterWrapper <- function(plotDat.all,compset.prefix){
  #Set shared params
  nsamp <- length(plotDat.all$sensitivity$datBySize)
  lab.cex <- 0.8
  
  #Prepare plot layout
  png(paste(OUTDIR,"main_plots/VCF_QC.",compset.prefix,".perSample_benchmarking.png",sep=""),
      height=5*300,width=10.5*300,res=300)
  layout(matrix(c(1,2,3,4,
                  5,6,7,8),
                nrow=2,byrow=T),
         widths=c(5,4,4,4))
  
  ###Top row: sensitivity
  #Set params
  ylab <- "Sensitivity"
  plotDat <- plotDat.all$sensitivity
  #Violins by class
  plotViolins(mat=plotDat$datByClass.mat,
              colors=svtypes$color,
              xlab="SV Class",ylab=ylab,lab.cex=lab.cex,
              title=paste(ylab," vs. ",compset.prefix,sep=""))
  #Lineplot by size
  plotLinesByClass(means=plotDat$datBySizeClass.mat$mean,
                   ci_adj=plotDat$datBySizeClass.mat$ci_adj,
                   nsamp=nsamp,xlab="SV Size",ylab=ylab,lab.cex=lab.cex,
                   title=paste("Mean ",ylab," by Size",sep=""),legend=F)
  #Lineplot by freq
  plotLinesByClass(means=plotDat$datByFreqClass.mat$mean,
                   ci_adj=plotDat$datByFreqClass.mat$ci_adj,
                   nsamp=nsamp,xlab="SV Frequency",ylab=ylab,lab.cex=lab.cex,
                   title=paste("Mean ",ylab," by Freq. ",sep=""),legend=F)
  #Heatmap of size v freq
  plotHeatmap(mat=plotDat$datBySizeFreq.mat$mean,x.title="SV Size",y.title="SV Freq.",
              title=paste("Mean ",ylab," by Size x Freq.",sep=""),lab.cex=lab.cex,nsamp=nsamp)
  
  #Set params
  ylab <- "Specificity"
  plotDat <- plotDat.all$specificity
  #Violins by class
  plotViolins(mat=plotDat$datByClass.mat,
              colors=svtypes$color,
              xlab="SV Class",ylab=ylab,lab.cex=lab.cex,
              title=paste(ylab," vs. ",compset.prefix,sep=""))
  #Lineplot by size
  plotLinesByClass(means=plotDat$datBySizeClass.mat$mean,
                   ci_adj=plotDat$datBySizeClass.mat$ci_adj,
                   nsamp=nsamp,xlab="SV Size",ylab=ylab,lab.cex=lab.cex,
                   title=paste("Mean ",ylab," by Size",sep=""),legend=F)
  #Lineplot by freq
  plotLinesByClass(means=plotDat$datByFreqClass.mat$mean,
                   ci_adj=plotDat$datByFreqClass.mat$ci_adj,
                   nsamp=nsamp,xlab="SV Frequency",ylab=ylab,lab.cex=lab.cex,
                   title=paste("Mean ",ylab," by Freq. ",sep=""),legend=F)
  #Heatmap of size v freq
  plotHeatmap(mat=plotDat$datBySizeFreq.mat$mean,x.title="SV Size",y.title="SV Freq.",
              title=paste("Mean ",ylab," by Size x Freq.",sep=""),lab.cex=lab.cex,nsamp=nsamp)
  
  #Close device
  dev.off()
}


########################
###RSCRIPT FUNCTIONALITY
########################
###Load libraries as needed
require(optparse, quietly=T)
require(vioplot, quietly=T)
require(beeswarm, quietly=T)

###List of command-line options
option_list <- list(
  make_option(c("-c", "--comparisonSetName"), type="character", default="external",
              help="name of comparison set used for prefixing output [default %default]",
              metavar="character")
)

###Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog perSampleDir samples.list svtypes OUTDIR",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

###Checks for appropriate positional arguments
if(length(args$args) != 4){
  stop("Incorrect number of required positional arguments\n")
}

###Writes args & opts to vars
perSampDir <- args$args[1]
samples.in <- args$args[2]
svtypes.file <- args$args[3]
OUTDIR <- args$args[4]
compset.prefix <- opts$comparisonSetName

# #Dev parameters
# perSampDir <- "/Users/collins/scratch/gnomAD-SV_v3.chr19_to_22.v1/"
# samples.in <- "/Users/collins/scratch/gnomAD-SV_v3.chr19_to_22.v1.chr19.shard.shard_.analysis_samples.list"
# OUTDIR <- "~/scratch/perSample_benchmarking_plots_test/"
# svtypes.file <- "/Users/collins/Desktop/Collins/Talkowski/NGS/SV_Projects/gnomAD_v3/gnomad-sv-v3-qc//src/sv-pipeline/scripts/vcf_qc/SV_colors.txt"
# compset.prefix <- "HGSV_Ebert_perSample"

###Read & process input data
#Read list of samples
samples <- as.character(read.table(samples.in,check.names=F)[,1])
#Sets sv types & colors
svtypes <- read.table(svtypes.file,sep="\t",header=F,comment.char="",check.names=F)
svtypes <- as.data.frame(apply(svtypes,2,as.character))
colnames(svtypes) <- c("svtype","color")
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
if(!dir.exists(paste(OUTDIR,"/supporting_plots/per_sample_benchmarking_",compset.prefix,"/",sep=""))){
  dir.create(paste(OUTDIR,"/supporting_plots/per_sample_benchmarking_",compset.prefix,"/",sep=""))
}


###Iterates over sensitivity and specificity and generates plots
plotDat.all <- lapply(list("sensitivity","specificity"),function(measurement){
  #Formats variable
  measurement.lab <- paste(toupper(substr(measurement,1,1)),substr(measurement,2,nchar(measurement)),sep="")
  
  #Reads list of per-sample overlap data
  dat <- readMultiSampDat(samples, measurement=measurement)
  
  #Get number of samples
  nsamp <- length(dat)
  
  #Only plot if any samples exist
  if(nsamp>0){
    #Computes plot data
    plotDat <- gatherPlottingData(dat)
    
    #Master plotting block
    wrapperViolinPlots(plotDat=plotDat,compset.prefix=compset.prefix,ylab=measurement.lab)
    wrapperLinePlots(plotDat=plotDat,compset.prefix=compset.prefix,ylab=measurement.lab,nsamp=nsamp)
    
    #Return plotting data
    return(plotDat)
  }else{
    return(NULL)
  }
})
names(plotDat.all) <- c("sensitivity","specificity")

###Generates final master panel
if(!is.null(plotDat.all$sensitivity)){
  masterWrapper(plotDat.all=plotDat.all,compset.prefix=compset.prefix)
}

