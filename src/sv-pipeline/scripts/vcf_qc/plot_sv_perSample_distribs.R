#!/usr/bin/env Rscript

# Helper script to plot per-sample VCF stats output from collectQC.sh


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


###########################
###GENERIC HELPER FUNCTIONS
###########################
#Read & clean list of variant IDs & genotypes per sample
readDatPerSample <- function(ID){
  #Set path
  path <- paste(perSampDir,"/",ID,".VIDs_genotypes.txt.gz",sep="")
  #Read & process data if file exists
  if(file.exists(path)){
    #Check number of lines in file
    nlines <- as.integer(system(paste("gunzip -c ",path," | wc -l",sep=""),wait=T,intern=T,ignore.stderr=T))
    
    #Only run if number of lines > 0
    if(nlines>0){
      #Read data
      x <- read.table(path,header=F,check.names=F)
      
      #Check for correct number of columns
      if(ncol(x)==3){
        #Convert genotypes to number of alleles
        x[,2] <- sapply(x[,2],function(gt){
          #return NA for number of alleles for no-calls
          if(gt=="./."){
            return(NA)
          }else{
            #Check for multiallelic ./N notation, and compute # of alleles as divergence from diploid
            if(length(grep(".",as.character(gt),fixed=T))>0){
              abs(2-sum(as.numeric(gsub(".","",unlist(strsplit(as.character(gt),split="/")),fixed=T)),na.rm=T))
            #Otherwise, sum number of alleles
            }else{
              sum(as.numeric(gsub(".","",unlist(strsplit(as.character(gt),split="/")),fixed=T)))
            }
          }
        })
        
        #Format output data & remove variants with no observed alleles
        x[,2:3] <- apply(x[,2:3],2,as.numeric)
        colnames(x) <- c("VID","alleles","GQ")
        x <- x[which(x[,2]>0 & !is.na(x[,2])),]
        
        #Return data
        return(x)
      }else{
        warning(paste("VID file for sample ",ID," does not have correct number of columns",sep=""))
        return(NULL)
      }
    }else{
      warning(paste("VID file for sample ",ID," has no lines",sep=""))
      return(NULL)
    }
  }else{
    warning(paste("VID file not found for sample ",ID,sep=""))
    return(NULL)
  }
}
#Filters VIDs for a single sample based on GQ
GQfilterVIDs <- function(ID,min.GQ=0,max.GQ=99){
  #Check input variant list
  if(ID %in% names(VID.lists)){
    vlist <- VID.lists[[which(names(VID.lists)==ID)]]
  }else{
    vlist <- NULL
  }
  if(!is.null(vlist)){
    #Filter variant list on GQ
    vlist <- vlist[which(vlist$GQ>=min.GQ & vlist$GQ<=max.GQ),]
  }
  #Return result
  return(vlist)
}
#Count variants or alleles per class given a list of VIDs
countVarsSingle <- function(dat,vlist,count="variants"){
  #Sanity check data & args
  if(!(count %in% c("variants","alleles"))){
    stop("Must specify counting of 'variants' or 'alleles'")
  }
  if(!is.null(vlist)){
    if(count=="alleles" & !any(colnames(vlist)=="alleles")){
      stop("Counting 'alleles' specified, but no 'alleles' column present in input data")
    }
    
    #Iterate over svtypes and count number of entries
    if(count=="variants"){
      res <- sapply(svtypes$svtype,function(svtype){
        length(which(dat$svtype==svtype & dat$VID %in% vlist$VID))
      })
    }else{
      res <- sapply(svtypes$svtype,function(svtype){
        sum(vlist[which(vlist$VID %in% dat[which(dat$svtype==svtype),]$VID),]$alleles)
      })
    }
    #Return data
    return(res)
  }
}
#Create matrix of variant/allele counts for a list of samples
countVarsMulti <- function(dat,samples,count="variants",
                           min.GQ=0,max.GQ=99,biallelic=F){
  #Exclude non-biallelic sites, if optioned
  if(biallelic==T){
    dat <- dat[which(dat$carriers>0 & dat$other_gts==0),]
  }
  
  #Iterate over samples and count variants/alleles
  mat <- as.data.frame(t(sapply(samples,function(ID){
    vlist <- GQfilterVIDs(ID=ID,min.GQ=min.GQ,max.GQ=max.GQ)
    countVarsSingle(dat=dat,vlist=vlist,count=count)
  })))

  #Format & return data
  colnames(mat) <- svtypes$svtype
  rownames(mat) <- NULL
  mat <- cbind("ALL"=apply(mat,1,sum),mat)
  rownames(mat) <- samples
  return(mat)
}
#Gather median count of variants per sample per size bin
mediansBySize <- function(dat,samples,svtypes,max.sizes,size.labs,count="variants",biallelic=F){
  #Add buffer to max sizes
  max.sizes <- c(0,as.numeric(max.sizes))
  
  #Iterate over sizes and compose matrix of medians
  median.mat <- t(sapply(2:length(max.sizes),function(i){
    #Set size range
    min.size <- max.sizes[i-1]
    max.size <- max.sizes[i]
    
    #Get matrix of variants per sample
    mat.sub <- countVarsMulti(dat=dat[which(dat$length>min.size & dat$length<=max.size),],
                              samples=samples,count=count,biallelic=biallelic)
    
    #Get median per sample
    medians <- apply(mat.sub,2,median,na.rm=T)
    
    #Return medians
    return(medians)
  }))
  
  #Add size thresholds as row names
  rownames(median.mat) <- size.labs
  return(median.mat)
}
#Gather median count of variants per sample per freq bin
mediansByFreq <- function(dat,samples,svtypes,max.freqs,freq.labs,count="variants",biallelic=F){
  #Add buffer to max freqs
  max.freqs <- c(0,as.numeric(max.freqs))
  
  #Get index for carrier vs allele freq
  if(count=="variants"){
    freq.idx <- which(colnames(dat)=="carrierFreq")
  }else{
    freq.idx <- which(colnames(dat)=="AF")
  }
  
  #Iterate over freqs and compose matrix of medians
  median.mat <- t(sapply(2:length(max.freqs),function(i){
    #Set freq range
    min.freq <- max.freqs[i-1]
    max.freq <- max.freqs[i]
    
    #Get matrix of variants per sample
    mat.sub <- countVarsMulti(dat=dat[which(dat[,freq.idx]>min.freq & dat[,freq.idx]<=max.freq),],
                              samples=samples,count=count,biallelic=biallelic)
    
    #Get median per sample
    medians <- apply(mat.sub,2,median,na.rm=T)
    
    #Return medians
    return(medians)
  }))
  
  #Add freq thresholds as row names
  rownames(median.mat) <- freq.labs
  return(median.mat)
}
#Gather median count of variants per sample by minimum GQ
mediansByGQ <- function(dat,samples,svtypes,min.GQs=seq(0,90,10),count="variants",biallelic=F,max.GQ=99){
  #Iterate over min GQs and compose matrix of medians
  median.mat <- t(sapply(min.GQs,function(min.GQ){
    #Get matrix of variants per sample
    mat.sub <- countVarsMulti(dat=dat,samples=samples,count=count,min.GQ=min.GQ,biallelic=biallelic,max.GQ=max.GQ)
    
    #Get median per sample
    medians <- apply(mat.sub,2,median,na.rm=T)
    
    #Return medians
    return(medians)
  }))
  
  #Add size thresholds as row names
  rownames(median.mat) <- paste(">",min.GQs,sep="")
  return(median.mat)
}


############################
###PLOTTING HELPER FUNCTIONS
############################
#Violin plot of a matrix of counts
plotViolins <- function(mat,colors=NULL,log=F,ylims=NULL,
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
  
  #Split matrix into all vals and per-class vals for plotting
  mat.all <- mat[,1]
  mat.classes <- mat[,-1]
  #Scale to match mat.all
  if(any(mat.classes>0)){
    all.scalar <- (max(mat.all,na.rm=T)/max(mat.classes,na.rm=T))
  }else{
    all.scalar <- 1
  }
  mat.classes.scaled <- mat.classes*all.scalar
  
  #Prep plotting area
  right.sep <- 3
  text.buffer <- 0.4
  par(mar=c(4.5,4,2,1),bty="n")
  plot(x=c(-0.15-right.sep,ncol(mat)-1+text.buffer),y=ylims,type="n",
       xaxt="n",xaxs="i",yaxt="n",xlab="",ylab="")
  
  #Add gridlines, if optioned
  if(gridlines==T){
    abline(h=axTicks(2),col="gray92")
  }
  rect(xleft=-right.sep+1+text.buffer,xright=0,
       ybottom=par("usr")[3],ytop=par("usr")[4],
       border=NA,col="white")
  
  #Add shared axes & labels
  mtext(1,line=3-max(c(0,(1-lab.cex))),text=xlab,cex=lab.cex)
  mtext(3,line=0,cex=0.7*lab.cex,
        text=paste("N=",prettyNum(nrow(mat),big.mark=",")," Samples",sep=""))
  mtext(3,line=0.8,text=paste(title,y.suffix,sep=""),font=2,cex=lab.cex)
  
  #Add left-panel y-axis
  axis(2,at=axTicks(2),labels=NA)
  axis(2,at=axTicks(2),tick=F,line=-0.4,cex.axis=0.8,las=2,
       labels=prettyNum(round(axTicks(2),2),big.mark=","))
  mtext(2,line=2.75-max(c(0,(1-lab.cex))),text=ylab,cex=lab.cex)
  
  #Add right panel axes
  segments(x0=0,x1=0,y0=min(axTicks(2)),y1=max(axTicks(2)))
  segments(x0=-0.1,x1=0,y0=axTicks(2),y1=axTicks(2))
  text(x=0,y=axTicks(2),labels=prettyNum(round(axTicks(2)/all.scalar,0),big.mark=","),
       pos=2,cex=0.8)
  
  #Plot distribution of all variants on left panel
  for(i in 1){
    #Get values
    vals <- mat.all
    
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
          vioplot(inliers,add=T,at=i-right.sep-0.5,wex=0.35,
                  col=colors[i],drawRect=F,border=NA)
        }else{
          points(x=i-right.sep-0.5,y=unique(inliers),pch=18,cex=1.5,col=colors[i])
        }
        
        #Plot points for outliers
        points(x=jitter(rep(i-0.5,length(outliers)),amount=0.01),
               y=outliers,pch=21,cex=0.175,lwd=0.7,
               col=colors[i],bg=adjustcolor(colors[i],alpha=0.3))
        #Otherwise, use swarmplots
      }else{
        beeswarm(vals,add=T,at=i-right.sep-0.5,pch=21,cex=0.25,
                 col=colors[i],bg=adjustcolor(colors[i],alpha=0.3),
                 corral="wrap",corralWidth=0.35)
      }
      
      #Add bar for median & vertical rule for IQR
      segments(x0=i-right.sep-0.6,x1=i-right.sep-0.4,lwd=2,
               y0=median(vals,na.rm=T))
      segments(x0=i-right.sep-0.5,x1=i-right.sep-0.5,
               y0=quantile(vals,probs=0.25,na.rm=T),
               y1=quantile(vals,probs=0.75,na.rm=T))
      
      #Add mean on x-axis
      axis(1,at=i-right.sep-0.3,
           labels=bquote(mu == .(prettyNum(round(mean(vals,na.rm=T),1)))),
           tick=F,las=2,line=-0.8,cex.axis=0.6)
      
      #Add text label for median value
      text(x=i-right.sep-0.48,y=median(vals),pos=4,cex=0.5,font=4,
           col=colors[i],labels=prettyNum(round(median(vals),0),big.mark=","))
    }
    
    #Add category label on x-axis
    axis(1,at=i-right.sep-0.5,
         labels=paste(colnames(mat)[i],"\n",sep=""),
         tick=F,las=2,line=-0.8,cex.axis=0.8,
         col.axis=colors[i],font=2)
  }
  
  #Iterate over columns per class and plot distributions
  sapply(1:ncol(mat.classes),function(i){
    #Get values
    vals <- mat.classes[,i]
    
    #Only plot distribution if any values are non-zero
    if(any(vals>0 & !is.infinite(as.vector(vals)) & !is.na(as.vector(vals)))){
      vals.scaled <- mat.classes.scaled[,i]
      
      if(violin==T){
        #Get inliers & outliers
        IQR <- IQR(vals.scaled,na.rm=T)
        median <- median(vals.scaled,na.rm=T)
        inliers <- vals.scaled[which(vals.scaled>=median-(3*IQR) & vals.scaled<=median+(3*IQR))]
        outliers <- vals.scaled[which(vals.scaled<=median-(3*IQR) | vals.scaled>=median+(3*IQR))]
        
        #Plot violin of inliers
        if(length(unique(inliers))>1){
          vioplot(inliers,add=T,at=i-0.5,wex=0.35,
                  col=colors[i+1],drawRect=F,border=NA)
        }else{
          points(x=i-0.5,y=unique(inliers),pch=18,cex=1.5,col=colors[i])
        }
        
        #Plot points for outliers
        points(x=jitter(rep(i-0.5,length(outliers)),amount=0.01),
               y=outliers,pch=21,cex=0.175,lwd=0.7,
               col=colors[i+1],bg=adjustcolor(colors[i],alpha=0.3))
        #Otherwise, use swarmplots
      }else{
        beeswarm(vals.scaled,add=T,at=i-0.5,pch=21,cex=0.25,
                 col=colors[i+1],bg=adjustcolor(colors[i],alpha=0.3),
                 corral="wrap",corralWidth=0.35)
      }
      
      #Add bar for median & vertical rule for IQR
      segments(x0=i-0.6,x1=i-0.4,lwd=2,
               y0=median(vals.scaled,na.rm=T))
      segments(x0=i-0.5,x1=i-0.5,
               y0=quantile(vals.scaled,probs=0.25,na.rm=T),
               y1=quantile(vals.scaled,probs=0.75,na.rm=T))
      
      #Add mean on x-axis
      axis(1,at=i-0.3,
           labels=bquote(mu == .(prettyNum(round(mean(vals,na.rm=T),1)))),
           tick=F,las=2,line=-0.8,cex.axis=0.6) 
    }
    
    #Add category label on x-axis
    axis(1,at=i-0.5,
         labels=paste(colnames(mat)[i+1],"\n",sep=""),
         tick=F,las=2,line=-0.8,cex.axis=0.8,
         col.axis=colors[i+1],font=2)
  })
  
  #Iterate over columns and add text labels for median values
  sapply(1:ncol(mat.classes),function(i){
    vals <- mat.classes[,i]
    
    #Only plot distribution if any values are non-zero
    if(any(vals>0 & !is.infinite(as.vector(vals)) & !is.na(as.vector(vals)))){
      vals.scaled <- mat.classes.scaled[,i]
      text(x=i-0.48,y=median(vals.scaled),pos=4,cex=0.5,font=4,
           col=colors[i+1],labels=prettyNum(round(median(vals),0),big.mark=","))
    }
  })
  
  #Clean up boxes
  rect(xleft=par("usr")[1],xright=-right.sep+1+text.buffer,
       ybottom=par("usr")[3],ytop=par("usr")[4],col=NA)
  rect(xleft=0,xright=par("usr")[2],ybottom=par("usr")[3],ytop=par("usr")[4],col=NA)
}
#Generic heatmap plotting function
plotHeatmap <- function(mat,base.cols,
                        x.labels=NULL,x.title=NULL,
                        y.labels=NULL,y.title=NULL,
                        title=NULL,lab.cex=1){
  #Set values if NULL
  if(is.null(x.labels)){
    x.labels <- colnames(mat)
  }
  if(is.null(y.labels)){
    y.labels <- rownames(mat)
  }
  
  #Convert medians to fractions
  mat.frac <- apply(mat,2,function(vals){
    if(any(vals>0)){
      newvals <- vals/sum(vals,na.rm=T)
    }else{
      newvals <- rep(NA,times=length(vals))
    }
    return(newvals)
  })
  
  #Prep plotting area
  par(mar=c(4,4,2,2))
  plot(x=c(0,ncol(mat)),y=c(0,-nrow(mat)),type="n",
       xaxt="n",xaxs="i",xlab="",yaxt="n",yaxs="i",ylab="")
  
  #Add axes
  sapply(1:ncol(mat),function(i){
    axis(1,at=i-0.5,tick=F,line=-0.8,las=2,labels=x.labels[i],font=2,
         cex.axis=0.7,col.axis=base.cols[i])
  })
  mtext(1,line=2.75-max(c(0,(1-lab.cex))),text=x.title,cex=lab.cex)
  axis(2,at=-(1:nrow(mat))+0.5,tick=F,line=-0.8,las=2,labels=y.labels,cex.axis=0.7)
  mtext(2,line=2.75-max(c(0,(1-lab.cex))),text=y.title,cex=lab.cex)
  mtext(3,line=0,cex=0.7*lab.cex,
        text=paste("Median of N=",prettyNum(length(samples),big.mark=",")," Samples",sep=""))
  mtext(3,line=0.8,text=title,font=2,cex=lab.cex)
  
  #Plot all cells
  sapply(1:nrow(mat),function(r){
    sapply(1:ncol(mat),function(c){
      #Get color range
      col.range <- colorRampPalette(c("white",base.cols[c]))(101)
      #Get & scale value
      val <- ceiling(mat[r,c])
      #Get color for shading
      if(mat[r,c]==0 | is.na(mat[r,c])){
        color <- "gray80"
        pct <- "N/A"
        dens <- 12
      }else{
        frac <- as.numeric(round(100*mat.frac[r,c],0))
        color <- col.range[min((2*frac)+1,101)]
        pct <- paste(frac,"%",sep="")
        dens <- NA
      }
      #Plot rectangle
      rect(xleft=c-1,xright=c,ybottom=-r,ytop=-(r-1),
           lwd=0.5,border="gray95",col=color,density=dens)
      #Format cell annotation
      anno.text <- paste(prettyNum(val,big.mark=","),"\n(",
                         pct,")",sep="")
      text(x=c-0.5,y=-(r-0.5),labels=anno.text,
           cex=0.6)
    })
  })
  
  #Clean up box
  box()
}
#General function to plot stacked bars from a matrix
plotStackedBars <- function(mat,colors,scaled=T,lab.cex=1,
                            xlabel=NULL,ylabel=NULL,title=NULL,subtitle=NULL){
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
  axis(2,at=axTicks(2),tick=F,labels=ylabs,las=2,cex.axis=0.8*lab.cex,line=-0.4)
  mtext(3,line=1,text=title,font=2,cex=lab.cex)
  mtext(3,line=0.1,text=subtitle,cex=0.7*lab.cex)
  mtext(1,line=2.5-max(c(0,(1-lab.cex))),text=xlabel,cex=lab.cex)
  mtext(2,line=2.75-max(c(0,(1-lab.cex))),text=ylabel,cex=lab.cex)
  
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
    axis(1,at=i-0.5,las=2,line=-0.8,labels=colnames(mat)[i],cex.axis=0.8*lab.cex,tick=F)
  })
}

# #Get count of het/hom/other genotypes per sample
# countGTs <- function(sample){
#   #Restrict to gts from biallelic variants
#   biallelic.vars <- dat$VID[which(dat$other_gts==0)]
#   VID.list <- VID.lists[[which(names(VID.lists)==sample)]]
#   gts <- VID.list[which(VID.list[,1] %in% biallelic.vars),2]
#   return(gts)
# }


#############################
###VARIANT/ALLELE COUNT PLOTS
#############################
#Master wrapper for all variant/allele count violins
wrapperVariantCountViolins <- function(count="variants"){
  #Subset plotting data
  plot.data.sub <- plot.data[[which(names(plot.data)==count)]]
  
  #Set plotting label & prefix
  if(count=="variants"){
    label.prefix <- "SV Sites"
    freq.idx <- which(colnames(dat)=="carrierFreq")
  }else{
    label.prefix <- "SV Alleles"
    freq.idx <- which(colnames(dat)=="AF")
  }
  
  #All variants
  png(paste(OUTDIR,"/supporting_plots/per_sample_plots/sv_counts.",count,"_per_genome.all_sv.png",sep=""),
      height=1200,width=1500,res=300)
  plotViolins(mat=plot.data.sub$count.all,
              colors=svtypes$color,log=F,
              xlab="SV Classes",ylab=paste(label.prefix," per Genome",sep=""),
              title=paste(label.prefix," per Genome [All SV]",sep=""))
  dev.off()
  
  #Tiny
  png(paste(OUTDIR,"/supporting_plots/per_sample_plots/sv_counts.",count,"_per_genome.tiny_sv.png",sep=""),
      height=1200,width=1500,res=300)
  plotViolins(mat=plot.data.sub$count.tiny,
              colors=svtypes$color,log=F,
              xlab="SV Classes",ylab=paste(label.prefix," per Genome",sep=""),
              title=paste(label.prefix," per Genome [< 100bp]",sep=""))
  dev.off()
  #Small
  png(paste(OUTDIR,"/supporting_plots/per_sample_plots/sv_counts.",count,"_per_genome.small_sv.png",sep=""),
      height=1200,width=1500,res=300)
  plotViolins(mat=plot.data.sub$count.small,
              colors=svtypes$color,log=F,
              xlab="SV Classes",ylab=paste(label.prefix," per Genome",sep=""),
              title=paste(label.prefix," per Genome [100-500bp]",sep=""))
  dev.off()
  #Medium
  png(paste(OUTDIR,"/supporting_plots/per_sample_plots/sv_counts.",count,"_per_genome.medium_sv.png",sep=""),
      height=1200,width=1500,res=300)
  plotViolins(mat=plot.data.sub$count.medium,
              colors=svtypes$color,log=F,
              xlab="SV Classes",ylab=paste(label.prefix," per Genome",sep=""),
              title=paste(label.prefix," per Genome [500bp-2.5kb]",sep=""))
  dev.off()
  #Medium-large
  png(paste(OUTDIR,"/supporting_plots/per_sample_plots/sv_counts.",count,"_per_genome.medlarge_sv.png",sep=""),
      height=1200,width=1500,res=300)
  plotViolins(mat=plot.data.sub$count.medlarge,
              colors=svtypes$color,log=F,
              xlab="SV Classes",ylab=paste(label.prefix," per Genome",sep=""),
              title=paste(label.prefix," per Genome [2.5-10kb]",sep=""))
  dev.off()
  #Large
  png(paste(OUTDIR,"/supporting_plots/per_sample_plots/sv_counts.",count,"_per_genome.large_sv.png",sep=""),
      height=1200,width=1500,res=300)
  plotViolins(mat=plot.data.sub$count.large,
              colors=svtypes$color,log=F,
              xlab="SV Classes",ylab=paste(label.prefix," per Genome",sep=""),
              title=paste(label.prefix," per Genome [10-50kb]",sep=""))
  dev.off()
  #Huge
  png(paste(OUTDIR,"/supporting_plots/per_sample_plots/sv_counts.",count,"_per_genome.huge_sv.png",sep=""),
      height=1200,width=1500,res=300)
  plotViolins(mat=plot.data.sub$count.huge,
              colors=svtypes$color,log=F,
              xlab="SV Classes",ylab=paste(label.prefix," per Genome",sep=""),
              title=paste(label.prefix," per Genome [>50kb]",sep=""))
  dev.off()
  
  #Singleton
  png(paste(OUTDIR,"/supporting_plots/per_sample_plots/sv_counts.",count,"_per_genome.singleton_sv.png",sep=""),
      height=1200,width=1500,res=300)
  plotViolins(mat=plot.data.sub$count.singleton,
              colors=svtypes$color,log=F,
              xlab="SV Classes",ylab=paste(label.prefix," per Genome",sep=""),
              title=paste(label.prefix," per Genome [AC=1]",sep=""))
  dev.off()
  #Rare
  png(paste(OUTDIR,"/supporting_plots/per_sample_plots/sv_counts.",count,"_per_genome.rare_sv.png",sep=""),
      height=1200,width=1500,res=300)
  plotViolins(mat=plot.data.sub$count.rare,
              colors=svtypes$color,log=F,
              xlab="SV Classes",ylab=paste(label.prefix," per Genome",sep=""),
              title=paste(label.prefix," per Genome [AC>1 & Freq. <1%]",sep=""))
  dev.off()
  #Uncommon
  png(paste(OUTDIR,"/supporting_plots/per_sample_plots/sv_counts.",count,"_per_genome.uncommon_sv.png",sep=""),
      height=1200,width=1500,res=300)
  plotViolins(mat=plot.data.sub$count.uncommon,
              colors=svtypes$color,log=F,
              xlab="SV Classes",ylab=paste(label.prefix," per Genome",sep=""),
              title=paste(label.prefix," per Genome [Freq. 1-10%]",sep=""))
  dev.off()
  #Common
  png(paste(OUTDIR,"/supporting_plots/per_sample_plots/sv_counts.",count,"_per_genome.common_sv.png",sep=""),
      height=1200,width=1500,res=300)
  plotViolins(mat=plot.data.sub$count.common,
              colors=svtypes$color,log=F,
              xlab="SV Classes",ylab=paste(label.prefix," per Genome",sep=""),
              title=paste(label.prefix," per Genome [Freq. 10-50%]",sep=""))
  dev.off()
  #Major
  png(paste(OUTDIR,"/supporting_plots/per_sample_plots/sv_counts.",count,"_per_genome.major_sv.png",sep=""),
      height=1200,width=1500,res=300)
  plotViolins(mat=plot.data.sub$count.major,
              colors=svtypes$color,log=F,
              xlab="SV Classes",ylab=paste(label.prefix," per Genome",sep=""),
              title=paste(label.prefix," per Genome [Freq. >50%]",sep=""))
  dev.off()
}
#Master wrapper for all variant/allele count barplots
wrapperVariantCountBarplots <- function(count="variants"){
  #Subset plotting data
  plot.data.sub <- plot.data[[which(names(plot.data)==count)]]
  
  #Set plotting label & prefix
  if(count=="variants"){
    label.prefix <- "SV Sites per Genome"
    freq.lab <- "Carrier Frequency"
  }else{
    label.prefix <- "SV Alleles per Genome"
    freq.lab <- "Allele Frequency"
  }
  
  #Raw counts vs size
  png(paste(OUTDIR,"/supporting_plots/per_sample_plots/sv_counts.",count,"_per_genome.raw_barplots_by_size.png",sep=""),
      height=1200,width=1200,res=300)
  plotStackedBars(mat=t(plot.data.sub$median.size)[-1,],
                  colors=svtypes$color,scaled=F,
                  xlabel="SV Size",ylabel="Count per Genome",
                  title=paste(label.prefix," vs. Size",sep=""),
                  subtitle=paste("Median of ",prettyNum(length(samples),big.mark=",")," Samples",sep=""))
  dev.off()
  
  #Scaled counts vs size
  png(paste(OUTDIR,"/supporting_plots/per_sample_plots/sv_counts.",count,"_per_genome.scaled_barplots_by_size.png",sep=""),
      height=1200,width=1200,res=300)
  plotStackedBars(mat=t(plot.data.sub$median.size)[-1,],
                  colors=svtypes$color,scaled=T,
                  xlabel="SV Size",ylabel="Pct. per Genome",
                  title=paste("Pct. of ",label.prefix," vs. Size",sep=""),
                  subtitle=paste("Median of ",prettyNum(length(samples),big.mark=",")," Samples",sep=""))
  dev.off()
  
  #Raw counts vs freq
  png(paste(OUTDIR,"/supporting_plots/per_sample_plots/sv_counts.",count,"_per_genome.raw_barplots_by_freq.png",sep=""),
      height=1200,width=1200,res=300)
  plotStackedBars(mat=t(plot.data.sub$median.freq)[-1,],
                  colors=svtypes$color,scaled=F,
                  xlabel=freq.lab,ylabel="Count per Genome",
                  title=paste(label.prefix," vs. Freq.",sep=""),
                  subtitle=paste("Median of ",prettyNum(length(samples),big.mark=",")," Samples",sep=""))
  dev.off()
  
  #Scaled counts vs freq
  png(paste(OUTDIR,"/supporting_plots/per_sample_plots/sv_counts.",count,"_per_genome.scaled_barplots_by_freq.png",sep=""),
      height=1200,width=1200,res=300)
  plotStackedBars(mat=t(plot.data.sub$median.freq)[-1,],
                  colors=svtypes$color,scaled=T,
                  xlabel=freq.lab,ylabel="Pct. per Genome",
                  title=paste("Pct. of ",label.prefix," vs. Freq.",sep=""),
                  subtitle=paste("Median of ",prettyNum(length(samples),big.mark=",")," Samples",sep=""))
  dev.off()
  
  #Raw counts vs GQ
  png(paste(OUTDIR,"/supporting_plots/per_sample_plots/sv_counts.",count,"_per_genome.raw_barplots_by_GQ.png",sep=""),
      height=1200,width=1200,res=300)
  plotStackedBars(mat=plot.data.sub$median.GQ[-1,],
                  colors=svtypes$color,scaled=F,
                  xlabel="Min GQ",ylabel="Count per Genome",
                  title=paste(label.prefix," vs. GQ",sep=""),
                  subtitle=paste("Median of ",prettyNum(length(samples),big.mark=",")," Samples",sep=""))
  dev.off()
  
  #Scaled counts vs GQ
  png(paste(OUTDIR,"/supporting_plots/per_sample_plots/sv_counts.",count,"_per_genome.scaled_barplots_by_GQ.png",sep=""),
      height=1200,width=1200,res=300)
  plotStackedBars(mat=plot.data.sub$median.GQ[-1,],
                  colors=svtypes$color,scaled=T,
                  xlabel="Min GQ",ylabel="Pct. per Genome",
                  title=paste("Pct. of ",label.prefix," vs. GQ",sep=""),
                  subtitle=paste("Median of ",prettyNum(length(samples),big.mark=",")," Samples",sep=""))
  dev.off()
}
#Master wrapper for all variant/allele count heatmaps
wrapperVariantCountHeats <- function(count="variants"){
  #Subset plotting data
  plot.data.sub <- plot.data[[which(names(plot.data)==count)]]
  
  #Set plotting label & prefix
  if(count=="variants"){
    label.prefix <- "Median SV Sites"
    ylab.prefix <- "Carrier"
  }else{
    label.prefix <- "Median SV Alleles"
    ylab.prefix <- "Allele"
  }
  
  #SV count by size
  pdf(paste(OUTDIR,"/supporting_plots/per_sample_plots/sv_counts.",count,"_per_genome.heatmap_by_size.pdf",sep=""),
      height=5,width=5)
  plotHeatmap(mat=plot.data.sub$median.size,
              base.cols=c("gray15",svtypes$color),
              x.title="SV Classes",y.title="SV Size",
              title=paste(label.prefix," per Genome, by Size",sep=""))
  dev.off()
  
  #SV count by frequency
  pdf(paste(OUTDIR,"/supporting_plots/per_sample_plots/sv_counts.",count,"_per_genome.heatmap_by_freq.pdf",sep=""),
      height=5,width=5)
  plotHeatmap(mat=plot.data.sub$median.freq,
              base.cols=c("gray15",svtypes$color),
              x.title="SV Classes",y.title=paste(ylab.prefix," Freq.",sep=""),
              title=paste(label.prefix," per Genome, by Freq.",sep=""))
  dev.off()
}


#######################
###MASTER PANEL WRAPPER
#######################
#Master wrapper for final top-level plot
masterWrapperSummaryPlot <- function(){
  #Prep plot area
  png(paste(OUTDIR,"/main_plots/VCF_QC.SV_per_genome.png",sep=""),
      height=5*300,width=11*300,res=300)
  layout(matrix(c(1,2,3,4,5,
                  6,7,8,9,10),nrow=2,byrow=T),
         widths=c(4,1.5,1.5,3,3))
  
  #Violin plot of SV sites per sample
  plotViolins(mat=plot.data$variants$count.all,
              colors=svtypes$color,log=F,
              xlab="SV Classes",ylab=paste("Sites",sep=""),
              title="SV Sites per Genome",lab.cex=0.75)
  
  #Raw counts of sites vs. GQ
  plotStackedBars(mat=plot.data$variants$median.GQ[-1,],
                  colors=svtypes$color,scaled=F,lab.cex=0.75,
                  xlabel="Min GQ",ylabel="Sites per Genome",
                  title="Sites vs. GQ",
                  subtitle=paste("Median of ",prettyNum(length(samples),big.mark=",")," Samples",sep=""))
  
  #Scaled counts of sites vs. GQ
  plotStackedBars(mat=plot.data$variants$median.GQ[-1,],
                  colors=svtypes$color,scaled=T,lab.cex=0.75,
                  xlabel="Min GQ",ylabel="Pct. of Sites per Genome",
                  title="Pct. Sites vs. GQ",
                  subtitle=paste("Median of ",prettyNum(length(samples),big.mark=",")," Samples",sep=""))
  
  #Heatmap of SV sites per sample by size
  plotHeatmap(mat=plot.data$variants$median.size,
              base.cols=c("gray15",svtypes$color),
              x.title="SV Classes",y.title="SV Size",
              title="Sites per Genome, by Size",lab.cex=0.75)
  
  #Heatmap of SV sites per sample by frequency
  plotHeatmap(mat=plot.data$variants$median.freq,
              base.cols=c("gray15",svtypes$color),
              x.title="SV Classes",y.title="Carrier Frequency",
              title="Sites per Genome, by Freq.",lab.cex=0.75)
  
  #Violin plot of SV alleles per sample
  plotViolins(mat=plot.data$alleles$count.all,
              colors=svtypes$color,log=F,
              xlab="SV Classes",ylab=paste("Alleles",sep=""),
              title="SV Alleles per Genome",lab.cex=0.75)
  
  #Raw counts of sites vs. GQ
  plotStackedBars(mat=plot.data$alleles$median.GQ[-1,],
                  colors=svtypes$color,scaled=F,lab.cex=0.75,
                  xlabel="Min GQ",ylabel="Alleles per Genome",
                  title="Alleles vs. GQ",
                  subtitle=paste("Median of ",prettyNum(length(samples),big.mark=",")," Samples",sep=""))
  
  #Scaled counts of sites vs. GQ
  plotStackedBars(mat=plot.data$alleles$median.GQ[-1,],
                  colors=svtypes$color,scaled=T,lab.cex=0.75,
                  xlabel="Min GQ",ylabel="Pct. of Alleles per Genome",
                  title="Pct. Alleles vs. GQ",
                  subtitle=paste("Median of ",prettyNum(length(samples),big.mark=",")," Samples",sep=""))
  
  #Heatmap of SV alleles per sample by size
  plotHeatmap(mat=plot.data$alleles$median.size,
              base.cols=c("gray15",svtypes$color),
              x.title="SV Classes",y.title="SV Size",
              title="Alleles per Genome by Size",lab.cex=0.75)
  
  #Heatmap of SV alleles per sample by frequency
  plotHeatmap(mat=plot.data$alleles$median.freq,
              base.cols=c("gray15",svtypes$color),
              x.title="SV Classes",y.title="Allele Frequency",
              title="Alleles per Genome, by Freq.",lab.cex=0.75)
  
  #Close device
  dev.off()
}


########################
###RSCRIPT FUNCTIONALITY
########################
###Load libraries as needed
require(optparse)
require(beeswarm)
require(vioplot)

###List of command-line options
option_list <- list(
  make_option(c("-S", "--svtypes"), type="character", default=NULL,
              help="tab-delimited file specifying SV types and HEX colors [default %default]",
              metavar="character"),
  make_option(c("-D", "--downsample"), type="integer", default=1000,
              help="restrict analyses to a random subset of N samples (for speed) [default %default]",
              metavar="integer"),
  make_option(c("-G", "--maxgq"), type="integer", default=99,
              help="Max GQ value, ie. 99 for GQ on a scale of [0,99]",
              metavar="integer")
)

###Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog svstats.bed samples.list perSampleDir OUTDIR",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options
downsamp <- opts$downsample
maxgq <- opts$maxgq

###Checks for appropriate positional arguments
if(length(args$args) != 4){
  stop("Incorrect number of required positional arguments\n")
}

###Writes args & opts to vars
dat.in <- args$args[1]
samples.in <- args$args[2]
perSampDir <- args$args[3]
OUTDIR <- args$args[4]
svtypes.file <- opts$svtypes

# #Dev parameters
# dat.in <- "~/scratch/xfer/gnomAD_v2_SV_test_cohort_v2_postCPX_VCF.VCF_sites.stats.bed.gz"
# samples.in <- "~/scratch/xfer/samples.list"
# perSampDir <- "~/scratch/xfer/gnomAD_v2_SV_test_cohort_v2_perSample_VIDs/"
# OUTDIR <- "~/scratch/perSample_plots_test/"
# # OUTDIR <- "~/scratch/melt_debug_plots/"
# svtypes.file <- "~/Desktop/Collins/Talkowski/code/sv-pipeline/ref/vcf_qc_refs/SV_colors.txt"
# downsamp <- 100


###Prepares I/O files
#Read & clean SV stats data
dat <- read.table(dat.in,comment.char="",sep="\t",header=T,check.names=F)
colnames(dat)[1] <- "chr"
#Read list of samples
samples <- unique(as.character(read.table(samples.in,check.names=F)[,1]))
#Downsample list of all samples, if optioned
if(downsamp<length(samples)){
  cat(paste("NOTE: downsampling to ",prettyNum(downsamp,big.mark=",")," samples, as specified. ",
            prettyNum(length(samples)-downsamp,big.mark=",")," samples excluded.\n",sep=""))
  samples <- sample(samples,size=downsamp,replace=F)
}
#Get number of samples
nsamp <- length(samples)
#Read list of variant IDs and allele counts per sample
VID.lists <- lapply(samples,readDatPerSample)
names(VID.lists) <- samples
#Exclude samples with NULL VIDs
exclude <- as.numeric(which(unlist(lapply(VID.lists,is.null))))
if(length(exclude)>0){
  cat(paste("WARNING: ",prettyNum(length(exclude),big.mark=","),"/",
            prettyNum(length(VID.lists),big.mark=",")," (",
            round(100*length(exclude)/length(VID.lists),1),"%) of samples have no observed variants.",
            " Excluding these samples from all per-sample analyses.\n",sep=""))
  VID.lists <- VID.lists[-exclude]
  samples <- samples[-exclude]
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
if(!dir.exists(paste(OUTDIR,"/supporting_plots/per_sample_plots/",sep=""))){
  dir.create(paste(OUTDIR,"/supporting_plots/per_sample_plots/",sep=""))
}


###Compute all plotting data
plot.data <- lapply(list("variants","alleles"),function(count){
  #Get freq.idx
  if(count=="variants"){
    freq.idx <- which(colnames(dat)=="carrierFreq")
  }else{
    freq.idx <- which(colnames(dat)=="AF")
  }
  
  #Counts per sample
  count.all <- countVarsMulti(dat=dat,samples=samples,count=count)
  count.tiny <- countVarsMulti(dat=dat[which(dat$length<=tiny.max.size),],
                               samples=samples,count=count)
  count.small <- countVarsMulti(dat=dat[which(dat$length>tiny.max.size & dat$length<=small.max.size),],
                                samples=samples,count=count)
  count.medium <- countVarsMulti(dat=dat[which(dat$length>small.max.size & dat$length<=medium.max.size),],
                                 samples=samples,count=count)
  count.medlarge <- countVarsMulti(dat=dat[which(dat$length>medium.max.size & dat$length<=medlarge.max.size),],
                                   samples=samples,count=count)
  count.large <- countVarsMulti(dat=dat[which(dat$length>medlarge.max.size & dat$length<=large.max.size),],
                                samples=samples,count=count)
  count.huge <- countVarsMulti(dat=dat[which(dat$length>large.max.size),],
                               samples=samples,count=count)
  count.singleton <- countVarsMulti(dat=dat[which(dat$AC==1),],
                                    samples=samples,count=count)
  count.rare <- countVarsMulti(dat=dat[which(dat$AC>1 & dat[,freq.idx]<=rare.max.freq),],
                               samples=samples,count=count)
  count.uncommon <- countVarsMulti(dat=dat[which(dat[,freq.idx]>rare.max.freq & dat[,freq.idx]<=uncommon.max.freq),],
                                   samples=samples,count=count)
  count.common <- countVarsMulti(dat=dat[which(dat[,freq.idx]>uncommon.max.freq & dat[,freq.idx]<=common.max.freq),],
                                 samples=samples,count=count)
  count.major <- countVarsMulti(dat=dat[which(dat[,freq.idx]>common.max.freq & dat[,freq.idx]<=major.max.freq),],
                                samples=samples,count=count)
  
  #Medians per sample
  stepgq <- ceiling(maxgq / 10)
  median.GQ <- t(mediansByGQ(dat=dat,samples=samples,count=count,max.GQ=maxgq,min.GQs=seq(0,stepgq*9,stepgq)))
  median.size <- mediansBySize(dat=dat,samples=samples,count=count,
                               max.sizes=c(tiny.max.size,small.max.size,medium.max.size,
                                           medlarge.max.size,large.max.size,huge.max.size),
                               size.labs=c("<100bp","100-\n500bp","500bp-\n2.5kb",
                                           "2.5-10kb","10kb-50kb",">50kb"))
  median.freq <- mediansByFreq(dat=dat,samples=samples,count=count,
                               max.freqs=c(1.1/nsamp,rare.max.freq,uncommon.max.freq,
                                           common.max.freq,major.max.freq),
                               freq.labs=c("AC=1","AC>1 &\nCF<1%","1-10%","10-50%",">50%"))
  
  #Format master list output
  out.list <- list("count.all"=count.all,
                   "count.tiny"=count.tiny,
                   "count.small"=count.small,
                   "count.medium"=count.medium,
                   "count.medlarge"=count.medlarge,
                   "count.large"=count.large,
                   "count.huge"=count.huge,
                   "count.singleton"=count.singleton,
                   "count.rare"=count.rare,
                   "count.uncommon"=count.uncommon,
                   "count.common"=count.common,
                   "count.major"=count.major,
                   "median.GQ"=median.GQ,
                   "median.size"=median.size,
                   "median.freq"=median.freq)
  return(out.list)
})
names(plot.data) <- c("variants","alleles")


###Master plotting block
#Summary panels
masterWrapperSummaryPlot()
#Variant sites per sample
wrapperVariantCountViolins(count="variants")
wrapperVariantCountBarplots(count="variants")
wrapperVariantCountHeats(count="variants")
#Variant alleles per sample
wrapperVariantCountViolins(count="alleles")
wrapperVariantCountBarplots(count="alleles")
wrapperVariantCountHeats(count="alleles")

