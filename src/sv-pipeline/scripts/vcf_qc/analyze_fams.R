#!/usr/bin/env Rscript

# Copyright (c) 2018 Talkowski Laboratory
# Contact: Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

# Helper script to perform family-based VCF QC & plot results


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
nocall.placeholder <- 9999


###################
###HELPER FUNCTIONS
###################
#Read & clean list of variant IDs & genotypes per sample
readDatPerSample <- function(ID,nocall.placeholder=9999){
  #Set path
  path <- paste(perSampDir,"/",ID,".VIDs_genotypes.txt.gz",sep="")
  #Read & process data if file exists
  if(file.exists(path)){
    #Read data
    x <- read.table(path,header=F,check.names=F)
    
    #Convert genotypes to number of alleles
    x[,2] <- sapply(x[,2],function(gt){
      #Return nocall.placeholder for no-calls
      if(gt=="./."){
        return(nocall.placeholder)
      }else{
        sum(as.numeric(gsub(".","",unlist(strsplit(as.character(gt),split="/")),fixed=T)))
      }
    })
    
    #Format output data
    x[,2:3] <- apply(x[,2:3],2,as.numeric)
    colnames(x) <- c("VID","alleles","GQ")
    
    #Return data
    return(x)
  }else{
    warning(paste("VID file not found for sample ",ID,sep=""))
    return(NULL)
  }
}
#Subset SV stats data
subsetDat <- function(dat,vlist,biallelic=T,
                      min.GQ.pro=0,max.GQ.pro=999,
                      min.GQ.par=0,max.GQ.par=999){
  #Check input variant list
  #Subset dat to variants found in sample & append number of alleles in sample
  x <- merge(dat,vlist,by="VID",sort=F)
  
  #Reorder columns
  x <- x[,c(2:4,1,5:ncol(x))]
  
  #Exclude non-biallelic sites, if optioned
  if(biallelic==T){
    x <- x[which(x$carriers>0 & x$other_gts==0),]
  }
  
  #Filter sites to specified GQ ranges
  x <- x[which(is.na(x$pro.GQ) | (x$pro.GQ>=min.GQ.pro & x$pro.GQ<=max.GQ.pro)),]
  # parent.GQ.range <- as.data.frame(t(apply(x[,which(colnames(x) %in% c("fa.GQ","mo.GQ")),],1,function(vals){
  #   if(any(!is.na(vals))){
  #     return(range(vals,na.rm=T))
  #   }else{
  #     return(c(NA,NA))
  #   }
  # })))
  # colnames(parent.GQ.range) <- c("min","max")
  # x <- x[which((is.na(parent.GQ.range$min) | parent.GQ.range$min>=min.GQ.par) & 
  #                (is.na(parent.GQ.range$max) | parent.GQ.range$max<=max.GQ.par)),]
  
  #Return result
  return(x)
}
#Gather matrix of SVs in any member of a family with information on allele counts in child & parent(s)
getFamDat <- function(dat,proband,father=NA,mother=NA,biallelic=T,nocall.placeholder=9999){
  #Clean parent IDs
  if(is.na(father)){
    father <- NULL
  }
  if(is.na(mother)){
    mother <- NULL
  }
  #Read VID lists for family members
  VID.lists <- lapply(c(proband,father,mother),readDatPerSample)
  names(VID.lists) <- c(proband,father,mother)
  
  #Get master VID list of all family members
  if(!is.null(father)){
    vlist <- merge(VID.lists[[which(names(VID.lists)==proband)]],
                   VID.lists[[which(names(VID.lists)==father)]],
                   sort=F,by="VID",all=T,suffixes=c(".pro",".fa"))
    if(!is.null(mother)){
      vlist <- merge(vlist,
                     VID.lists[[which(names(VID.lists)==mother)]],
                     sort=F,by="VID",all=T)
      colnames(vlist)[6:7] <- c("alleles.mo","GQ.mo")
    }else{
      vlist$alleles.mo <- NA
      vlist$GQ.mo <- NA
    }
  }else{
    vlist <- VID.lists[[which(names(VID.lists)==proband)]]
    colnames(vlist[which(colnames(vlist)=="alleles")]) <- "alleles.pro"
    colnames(vlist[which(colnames(vlist)=="GQ")]) <- "GQ.pro"
    vlist$alleles.fa <- 0
    vlist$GQ.fa <- NA
    vlist <- merge(vlist,
                   VID.lists[[which(names(VID.lists)==mother)]],
                   sort=F,by="VID",all=T,suffixes=c(".pro",".mo"))
  }
  #Only retain sites where all three samples are not null genotype (no-call, ./., nocall.placeholder)
  exclude <- sapply(vlist[,c(2,4,6)],
                          function(vals){
    which(as.numeric(vals)==nocall.placeholder)
  })

  exclude_out = c()
  for(exclude_x in exclude){  exclude_out=c(exclude_out,exclude_x) }
  exclude_out=unique(exclude_out)

  if(length(exclude_out) > 0){
    vlist <- vlist[-exclude_out,]
  }
  #Convert remaining NA allele counts to 0s
  vlist[,c(2,4,6)] <- apply(vlist[,c(2,4,6)],2,function(vals){
    vals[which(is.na(vals))] <- 0
    return(vals)
  })
  
  #Add transmission information to vlist
  trans <- t(apply(vlist[,c(2,4,6)],1,function(alleles){
    #Convert allele counts to numeric
    sapply(alleles,function(vals){
      vals[which(is.na(vals))] <- 0
    })
    #Get allele counts
    pro <- as.numeric(alleles[1])
    fa <- as.numeric(alleles[2])
    mo <- as.numeric(alleles[3])
    
    #Infer child inheritance status
    pro.denovo <- max(c(pro-(fa+mo),0))
    pro.inherited <- pro-pro.denovo
    
    #Divide credit for inherited variants between parents based on ratio of parent allele counts
    parental.alleles <- sum(c(fa,mo),na.rm=T)
    if(fa>0){
      p.fa <- fa/parental.alleles
    }else{
      p.fa <- 0
    }
    if(mo>0){
      p.mo <- mo/parental.alleles
    }else{
      p.mo <- 0
    }
    fa.transmitted <- pro.inherited*p.fa
    fa.untransmitted <- max(c(0,fa-fa.transmitted))
    mo.transmitted <- pro.inherited*p.mo
    mo.untransmitted <- max(c(0,mo-mo.transmitted))
    
    #Return vector of relevant transmission allele counts
    return(c(pro,pro.inherited,pro.denovo,
             fa,fa.transmitted,fa.untransmitted,
             mo,mo.transmitted,mo.untransmitted))
  }))
  trans <- as.data.frame(cbind(vlist[,1],trans,vlist[,c(3,5,7)]))
  colnames(trans) <- c("VID","pro.all","pro.inherited","pro.denovo",
                       "fa.all","fa.trans","fa.untrans",
                       "mo.all","mo.trans","mo.untrans",
                       "pro.GQ","fa.GQ","mo.GQ")
  
  #Subset data & add transmission data
  dat.fam <- subsetDat(dat=dat,vlist=trans,biallelic=biallelic)
  dat.fam[,(ncol(dat.fam)-9):ncol(dat.fam)] <- apply(dat.fam[,(ncol(dat.fam)-9):ncol(dat.fam)],2,as.numeric)
  return(dat.fam)
}
#Compute inheritance stats from a dat.fam dataframe
computeInheritance <- function(dat.fam,VIDs=NULL){
  #Clean dat.fam
  dat.fam <- as.data.frame(dat.fam)
  
  #Subset data frame based on list of VIDs
  if(!is.null(VIDs)){
    dat.fam <- dat.fam[which(as.character(dat.fam$VID) %in% as.character(VIDs)),]
  }
  #Compute allele-based inheritance rates
  pro.a.all <- sum(dat.fam$pro.all)
  pro.a.denovo <- sum(dat.fam$pro.denovo)
  pro.a.denovorate <- pro.a.denovo/pro.a.all
  pro.a.inh <- pro.a.all-pro.a.denovo
  pro.a.inhrate <- pro.a.inh/pro.a.all
  fa.a.all <- sum(dat.fam$fa.all)
  fa.a.trans <- sum(dat.fam$fa.trans)
  fa.a.transrate <- fa.a.trans/fa.a.all
  fa.a.untrans <- fa.a.all-fa.a.trans
  fa.a.untransrate <- fa.a.untrans/fa.a.all
  pro.a.patfrac <- fa.a.trans/pro.a.inh
  mo.a.all <- sum(dat.fam$mo.all)
  mo.a.trans <- sum(dat.fam$mo.trans)
  mo.a.transrate <- mo.a.trans/mo.a.all
  mo.a.untrans <- mo.a.all-mo.a.trans
  mo.a.untransrate <- mo.a.untrans/mo.a.all
  pro.a.matfrac <- mo.a.trans/pro.a.inh
  pro.a.patmatratio <- pro.a.patfrac/(pro.a.patfrac+pro.a.matfrac)
  
  #Compute variant-based inheritance rates
  pro.v.all <- length(which(dat.fam$pro.all>0))
  pro.v.denovo <- length(which(dat.fam$pro.all>0 & dat.fam$pro.inherited==0))
  pro.v.denovorate <- pro.v.denovo/pro.v.all
  pro.v.inh <- pro.v.all-pro.v.denovo
  pro.v.inhrate <- pro.v.inh/pro.v.all
  fa.v.all <- length(which(dat.fam$fa.all>0))
  fa.v.trans <- length(which(dat.fam$pro.all>0 & dat.fam$fa.all>0))
  fa.v.transrate <- fa.v.trans/fa.v.all
  fa.v.untrans <- fa.v.all-fa.v.trans
  fa.v.untransrate <- fa.v.untrans/fa.v.all
  pro.v.patfrac <- fa.v.trans/pro.v.inh
  mo.v.all <- length(which(dat.fam$mo.all>0))
  mo.v.trans <- length(which(dat.fam$pro.all>0 & dat.fam$mo.all>0))
  mo.v.transrate <- mo.v.trans/mo.v.all
  mo.v.untrans <- mo.v.all-mo.v.trans
  mo.v.untransrate <- mo.v.untrans/mo.v.all
  pro.v.matfrac <- mo.v.trans/pro.v.inh
  pro.v.patmatratio <- pro.a.patfrac/(pro.a.patfrac+pro.a.matfrac)
  
  #Format & return vector of rates
  return(c("pro.allele.all"=pro.a.all,
           "pro.allele.inh"=pro.a.inh,
           "pro.allele.inhrate"=pro.a.inhrate,
           "pro.allele.patfrac"=pro.a.patfrac,
           "pro.allele.matfrac"=pro.a.matfrac,
           "pro.allele.patmatratio"=pro.a.patmatratio,
           "pro.allele.denovo"=pro.a.denovo,
           "pro.allele.denovorate"=pro.a.denovorate,
           "fa.allele.all"=fa.a.all,
           "fa.allele.trans"=fa.a.trans,
           "fa.allele.transrate"=fa.a.transrate,
           "fa.allele.untrans"=fa.a.untrans,
           "fa.allele.untransrate"=fa.a.untransrate,
           "mo.allele.all"=mo.a.all,
           "mo.allele.trans"=mo.a.trans,
           "mo.allele.transrate"=mo.a.transrate,
           "mo.allele.untrans"=mo.a.untrans,
           "mo.allele.untransrate"=mo.a.untransrate,
           "pro.site.all"=pro.v.all,
           "pro.site.inh"=pro.v.inh,
           "pro.site.inhrate"=pro.v.inhrate,
           "pro.site.patfrac"=pro.v.patfrac,
           "pro.site.matfrac"=pro.v.matfrac,
           "pro.site.patmatratio"=pro.v.patmatratio,
           "pro.site.denovo"=pro.v.denovo,
           "pro.site.denovorate"=pro.v.denovorate,
           "fa.site.all"=fa.v.all,
           "fa.site.trans"=fa.v.trans,
           "fa.site.transrate"=fa.v.transrate,
           "fa.site.untrans"=fa.v.untrans,
           "fa.site.untransrate"=fa.v.untransrate,
           "mo.site.all"=mo.v.all,
           "mo.site.trans"=mo.v.trans,
           "mo.site.transrate"=mo.v.transrate,
           "mo.site.untrans"=mo.v.untrans,
           "mo.site.untransrate"=mo.v.untransrate))
}
#Compute inheritance for a list of trios and return as a data frame
computeInheritanceMulti <- function(trio.dat.list,VIDs=NULL){
  #Iterate over trios and compute inheritance
  res <- as.data.frame(t(sapply(trio.dat.list,computeInheritance,VIDs=VIDs)))
  return(res)
}
#Collect de novo rate per SV class
deNovoRateByClass <- function(trio.dat.list,VIDs=NULL){
  #Collect median DNR across all classes
  all.dat <- computeInheritanceMulti(trio.dat.list=trio.dat.list,VIDs=VIDs)
  all.dnrs <- c(median(all.dat$pro.site.denovorate,na.rm=T),
                median(all.dat$pro.allele.denovorate,na.rm=T))
  
  #Iterate over classes and return DNRs
  res <- sapply(svtypes$svtype,function(svtype){
    if(is.null(VIDs)){
      VIDs <- dat$VID
    }
    sub.dat <- computeInheritanceMulti(trio.dat.list=trio.dat.list,
                                       VIDs=dat$VID[which(dat$VID %in% VIDs & dat$svtype==svtype)])
    sub.dnrs <- c(median(sub.dat$pro.site.denovorate,na.rm=T),
                  median(sub.dat$pro.allele.denovorate,na.rm=T))
  })
  
  #Format & return results
  res <- as.data.frame(cbind(all.dnrs,res))
  colnames(res) <- c("ALL",svtypes$svtype)
  rownames(res) <- c("variants","alleles")
  return(res)
}
#Collect matrix of de novo rates by class by freq
deNovoRateByFreq <- function(trio.dat.list,freq.bins=40,count="variants"){
  #Get frequency index
  if(count=="variants"){
    freq.idx <- which(colnames(dat)=="carrierFreq")
  }else{
    freq.idx <- which(colnames(dat)=="AF")
  }
  
  #Create evenly spaced freq bins on log10-scale
  logfreq.min <- log10(min(dat[,freq.idx]))
  logfreq.max <- log10(1)
  logfreq.steps <- seq(logfreq.min,logfreq.max,by=(logfreq.max-logfreq.min)/(freq.bins-1))
  freq.df <- data.frame("min.freq"=c(0,10^logfreq.steps[-length(logfreq.steps)]),
                        "max.freq"=10^logfreq.steps)
  rownames(freq.df) <- paste(round(100*freq.df[,1],4),"-",
                             round(100*freq.df[,2],4),"%",sep="")
  
  #Iterate over frequency bins and gather de novo rates
  DNRs <- apply(freq.df,1,function(bounds){
    dnrs <- deNovoRateByClass(trio.dat.list=trio.dat.list,
                              VIDs=dat$VID[which(dat[,freq.idx]>bounds[1] & dat[,freq.idx]<=bounds[2])])
    return(as.numeric(dnrs[which(rownames(dnrs)==count),]))
  })
  
  #Format & return DNRs & freq.df
  DNRs <- as.data.frame(DNRs)
  DNRs <- apply(DNRs,2,as.numeric)
  rownames(DNRs) <- c("ALL",svtypes$svtype)
  return(list("DNRs"=DNRs,"bins"=freq.df))
}
#Collect matrix of de novo rates by class by size
deNovoRateBySize <- function(trio.dat.list,size.bins=40,count="variants"){
  #Create evenly spaced size bins on log10-scale
  logsize.min <- log10(50)
  logsize.max <- log10(1000000)
  logsize.steps <- seq(logsize.min,logsize.max,by=(logsize.max-logsize.min)/(size.bins-2))
  size.df <- data.frame("min.size"=c(0,10^logsize.steps),
                        "max.size"=c(10^logsize.steps,300000000))
  rownames(size.df) <- paste("10^",round(log10(size.df[,1]),1),
                             "-",round(log10(size.df[,2]),1),sep="")
  
  #Iterate over sizeuency bins and gather de novo rates
  DNRs <- apply(size.df,1,function(bounds){
    dnrs <- deNovoRateByClass(trio.dat.list=trio.dat.list,
                              VIDs=dat$VID[which(dat$length>bounds[1] & dat$length<=bounds[2])])
    return(as.numeric(dnrs[which(rownames(dnrs)==count),]))
  })
  
  #Format & return DNRs & size.df
  DNRs <- as.data.frame(DNRs)
  DNRs <- apply(DNRs,2,as.numeric)
  rownames(DNRs) <- c("ALL",svtypes$svtype)
  return(list("DNRs"=DNRs,"bins"=size.df))
}
#Collect matrix of de novo rates by size & freq combination
deNovoRateBySizeFreq <- function(trio.dat.list,VIDs=NULL,count="variants",
                                 max.sizes,size.labs,max.freqs,freq.labs){
  #Get frequency index
  if(count=="variants"){
    freq.idx <- which(colnames(dat)=="carrierFreq")
  }else{
    freq.idx <- which(colnames(dat)=="AF")
  }
  
  #Create size & freq bins
  size.df <- data.frame("min.size"=c(0,0,max.sizes),
                        "max.size"=c(300000000,max.sizes,300000000))
  rownames(size.df) <- c("ALL",size.labs)
  freq.df <- data.frame("min.freq"=c(0,0,max.freqs),
                        "max.freq"=c(1,max.freqs,1))
  rownames(freq.df) <- c("ALL",freq.labs)
  
  #Instantiate VIDs if necessary
  if(is.null(VIDs)){
    VIDs <- dat$VID
  }
  
  #Iterate over size bins & create DNR df for all SV
  DNRs <- as.data.frame(t(sapply(1:nrow(size.df),function(s){
    #Iterate over frequency bins
    sapply(1:nrow(freq.df),function(f){
      #Get de novo rate
      DNR <- deNovoRateByClass(trio.dat.list,
                               VIDs=dat$VID[which(dat$VID %in% VIDs & 
                                                    dat$length>size.df[s,1] & dat$length<=size.df[s,2] & 
                                                    dat[,freq.idx]>freq.df[f,1] & dat[,freq.idx]<=freq.df[f,2])])
      return(DNR$ALL[which(rownames(DNR)==count)])
    })
  })))
  colnames(DNRs) <- rownames(freq.df)
  rownames(DNRs) <- rownames(size.df)
  
  #Iterate over SV classes and create DNR df for each class
  DNRs.byClass <- lapply(svtypes$svtype,function(svtype){
    DNRs <- as.data.frame(t(sapply(1:nrow(size.df),function(s){
      #Iterate over frequency bins
      sapply(1:nrow(freq.df),function(f){
        #Get de novo rate
        DNR <- deNovoRateByClass(trio.dat.list,
                                 VIDs=dat$VID[which(dat$VID %in% VIDs & dat$svtype==svtype &
                                                      dat$length>size.df[s,1] & dat$length<=size.df[s,2] & 
                                                      dat[,freq.idx]>freq.df[f,1] & dat[,freq.idx]<=freq.df[f,2])])
        return(DNR$ALL[which(rownames(DNR)==count)])
      })
    })))
    colnames(DNRs) <- rownames(freq.df)
    rownames(DNRs) <- rownames(size.df)
    return(DNRs)
  })
  names(DNRs.byClass) <- svtypes$svtype
  
  #Combine all DNR dfs & return
  DNRs.all <- c(list(DNRs),DNRs.byClass)
  names(DNRs.all)[1] <- "ALL"
  return(DNRs.all)
}
#Collect matrix of de novo rates by class by size
deNovoRateBySize <- function(trio.dat.list,size.bins=40,count="variants"){
  #Create evenly spaced size bins on log10-scale
  logsize.min <- log10(50)
  logsize.max <- log10(1000000)
  logsize.steps <- seq(logsize.min,logsize.max,by=(logsize.max-logsize.min)/(size.bins-2))
  size.df <- data.frame("min.size"=c(0,10^logsize.steps),
                        "max.size"=c(10^logsize.steps,300000000))
  rownames(size.df) <- paste("10^",round(log10(size.df[,1]),1),
                             "-",round(log10(size.df[,2]),1),sep="")
  
  #Iterate over sizeuency bins and gather de novo rates
  DNRs <- apply(size.df,1,function(bounds){
    dnrs <- deNovoRateByClass(trio.dat.list=trio.dat.list,
                              VIDs=dat$VID[which(dat$length>bounds[1] & dat$length<=bounds[2])])
    return(as.numeric(dnrs[which(rownames(dnrs)==count),]))
  })
  
  #Format & return DNRs & size.df
  DNRs <- as.data.frame(DNRs)
  DNRs <- apply(DNRs,2,as.numeric)
  rownames(DNRs) <- c("ALL",svtypes$svtype)
  return(list("DNRs"=DNRs,"bins"=size.df))
}
#Collect matrix of de novo rates by class by minimum proband GQ
deNovoRateByProGQ <- function(trio.dat.list,GQ.bins=40,count="variants"){
  #Create evenly spaced GQ bins
  GQ.steps <- seq(0,1000,by=1000/GQ.bins)
  
  #Iterate over min GQs and gather de novo rates
  DNRs <- sapply(GQ.steps,function(min.GQ){
    tdl.tmp <- lapply(trio.dat.list,function(df){
      return(df[which(df$pro.GQ>=min.GQ),])
    })
    dnrs <- deNovoRateByClass(trio.dat.list=tdl.tmp)
    return(as.numeric(dnrs[which(rownames(dnrs)==count),]))
  })
  
  #Format & return DNRs & GQ.df
  DNRs <- as.data.frame(DNRs)
  DNRs <- apply(DNRs,2,as.numeric)
  rownames(DNRs) <- c("ALL",svtypes$svtype)
  colnames(DNRs) <- paste("gt",GQ.steps,sep="")
  return(list("DNRs"=DNRs,"bins"=GQ.steps))
}


############################
###PLOTTING HELPER FUNCTIONS
############################
#Generate main inheritance plot
plotInhStats <- function(inh.stats,count="variants",title=NULL,cex.lab=1){
  #Subset inh.stats to relevant columns
  if(count=="variants"){
    plot.df <- data.frame(inh.stats$pro.site.inhrate,
                          inh.stats$pro.site.patfrac,
                          inh.stats$pro.site.matfrac,
                          inh.stats$pro.site.patmatratio,
                          inh.stats$pro.site.denovorate,
                          "total.site.transrate"=(inh.stats$fa.site.trans+inh.stats$mo.site.trans)/(inh.stats$fa.site.all+inh.stats$mo.site.all),
                          inh.stats$fa.site.transrate,
                          inh.stats$mo.site.transrate)
  }else{
    plot.df <- data.frame(inh.stats$pro.allele.inhrate,
                          inh.stats$pro.allele.patfrac,
                          inh.stats$pro.allele.matfrac,
                          inh.stats$pro.allele.patmatratio,
                          inh.stats$pro.allele.denovorate,
                          "total.allele.transrate"=(inh.stats$fa.allele.trans+inh.stats$mo.allele.trans)/(inh.stats$fa.allele.all+inh.stats$mo.allele.all),
                          inh.stats$fa.allele.transrate,
                          inh.stats$mo.allele.transrate)
  }
  
  
  #Create vector of median fractions as representative for each category
  if(count=="variants"){
    median.counts <- c(paste(prettyNum(round(median(inh.stats$pro.site.inh),0),big.mark=",")," / ",
                             prettyNum(round(median(inh.stats$pro.site.all),0),big.mark=","),sep=""),
                       paste(prettyNum(round(median(inh.stats$fa.site.trans),0),big.mark=",")," / ",
                             prettyNum(round(median(inh.stats$pro.site.inh),0),big.mark=","),sep=""),
                       paste(prettyNum(round(median(inh.stats$mo.site.trans),0),big.mark=",")," / ",
                             prettyNum(round(median(inh.stats$pro.site.inh),0),big.mark=","),sep=""),
                       paste(round(100*median(inh.stats$pro.site.patfrac),0),"% : ",
                             round(100*median(inh.stats$pro.site.matfrac),0),"%",sep=""),
                       paste(prettyNum(round(median(inh.stats$pro.site.denovo),0),big.mark=",")," / ",
                             prettyNum(round(median(inh.stats$pro.site.all),0),big.mark=","),sep=""),
                       paste(prettyNum(round(median(inh.stats$fa.site.trans+inh.stats$mo.site.trans),0),big.mark=",")," / ",
                             prettyNum(round(median(inh.stats$fa.site.all+inh.stats$mo.site.all),0),big.mark=","),sep=""),
                       paste(prettyNum(round(median(inh.stats$fa.site.trans),0),big.mark=",")," / ",
                             prettyNum(round(median(inh.stats$fa.site.all),0),big.mark=","),sep=""),
                       paste(prettyNum(round(median(inh.stats$mo.site.trans),0),big.mark=",")," / ",
                             prettyNum(round(median(inh.stats$mo.site.all),0),big.mark=","),sep=""))
  }else{
    median.counts <- c(paste(prettyNum(round(median(inh.stats$pro.allele.inh),0),big.mark=",")," / ",
                             prettyNum(round(median(inh.stats$pro.allele.all),0),big.mark=","),sep=""),
                       paste(prettyNum(round(median(inh.stats$fa.allele.trans),0),big.mark=",")," / ",
                             prettyNum(round(median(inh.stats$pro.allele.inh),0),big.mark=","),sep=""),
                       paste(prettyNum(round(median(inh.stats$mo.allele.trans),0),big.mark=",")," / ",
                             prettyNum(round(median(inh.stats$pro.allele.inh),0),big.mark=","),sep=""),
                       paste(round(100*median(inh.stats$pro.allele.patfrac),0),"% : ",
                             round(100*median(inh.stats$pro.allele.matfrac),0),"%",sep=""),
                       paste(prettyNum(round(median(inh.stats$pro.allele.denovo),0),big.mark=",")," / ",
                             prettyNum(round(median(inh.stats$pro.allele.all),0),big.mark=","),sep=""),
                       paste(prettyNum(round(median(inh.stats$fa.allele.trans+inh.stats$mo.allele.trans),0),big.mark=",")," / ",
                             prettyNum(round(median(inh.stats$fa.allele.all+inh.stats$mo.allele.all),0),big.mark=","),sep=""),
                       paste(prettyNum(round(median(inh.stats$fa.allele.trans),0),big.mark=",")," / ",
                             prettyNum(round(median(inh.stats$fa.allele.all),0),big.mark=","),sep=""),
                       paste(prettyNum(round(median(inh.stats$mo.allele.trans),0),big.mark=",")," / ",
                             prettyNum(round(median(inh.stats$mo.allele.all),0),big.mark=","),sep=""))
  }
  
  #Set plot colors
  col.pat <- "#00B0CF"
  col.mat <- "#F064A5"
  col.dn <- "#F13D15"
  col.other <- "gray35"
  plot.cols <- c(col.other,col.pat,col.mat,col.other,
                 col.dn,col.other,col.pat,col.mat)
  lab.cols <- plot.cols
  lab.cols[which(lab.cols==col.other)] <- "black"
  lab.fonts <- c(2,3,3,3,2,2,3,3)
  
  #Prep plot area
  par(mar=c(1,6.5,3,4.5))
  plot(x=c(-0.05,1.05),y=c(0,-8.25),type="n",
       xaxt="n",yaxt="n",xlab="",ylab="",yaxs="i")
  
  #Dress up plot
  abline(v=seq(0,1,0.2),lty=2,col="gray85")
  # abline(h=c(0,-9),lwd=2,col="gray80")
  abline(v=c(0,1),col="gray75")
  # text(x=0.5,y=-0.2,pos=3,labels="SV Sites",font=4,cex=0.9)
  # text(x=0.5,y=-9.2,pos=3,labels="SV Alleles",font=4,cex=0.9)
  
  #Add category axes
  cat.labels <- c("Proband\nInheritance Rate","Inh. Rate   \n(Paternal)   ","Inh. Rate   \n(Maternal)   ",
                  "Pat:Mat Ratio   ","Proband\nDe Novo Rate","Parental\nTransmission Rate",
                  "Trans. Rate   \n(Paternal)   ","Trans. Rate   \n(Maternal)   ")
  sapply(1:8,function(i){
    axis(2,at=-i+0.3,line=-0.8,tick=F,las=2,cex.axis=0.7,font=lab.fonts[i],
         labels=cat.labels[i],col.axis=lab.cols[i])
  })
  # sapply(1:8,function(i){
  #   axis(2,at=-i-8.7,line=-0.8,tick=F,las=2,cex.axis=0.7,font=lab.fonts[i],
  #        labels=cat.labels[i],col.axis=lab.cols[i])
  # })
  axis(2,at=c(-1,-5,-6,-9,-10,-14,-15,-18)+0.75,labels=NA,tck=-0.1,col="gray80")
  
  #Add other axes & title
  # axis(2,at=0.3,tick=F,las=2,line=-0.5,cex.axis=0.8,font=2,
  #      label=paste("n=",prettyNum(nrow(plot.df),big.mark=","),
  #                  " ",fam.type,"s",sep=""))
  axis(3,at=seq(0,1,0.2),labels=NA)
  axis(3,at=seq(0,1,0.2),tick=F,line=-0.4,cex.axis=0.7,
       labels=paste(seq(0,100,20),"%",sep=""))
  mtext(3,line=1.5,text=title,font=2,cex=cex.lab)
  axis(4,at=-0.1,tick=F,las=2,line=-0.5,cex.axis=0.8,font=2,
       label="Median")
  
  #Plot points & category info
  sapply(1:ncol(plot.df),function(i){
    if(any(!is.na(plot.df[,i]))){
      #Add shading rect
      rect(xleft=par("usr")[1],xright=par("usr")[2],
           ybottom=-i+0.05,ytop=-i+0.45,bty="n",border=NA,
           col=adjustcolor(plot.cols[i],alpha=0.15))
      
      #Add points
      beeswarm(plot.df[,i],add=T,horizontal=T,at=-i+0.25,
               pch=21,cex=0.5,bg=plot.cols[i],pt.lwd=0.1,
               corral="wrap",corralWidth=0.2)
      
      #Add thick line & label for mean
      cat.mean <- mean(plot.df[,i],na.rm=T)
      segments(x0=cat.mean,x1=cat.mean,y0=-i,y1=-i+0.5,
               lend="round",lwd=4,col=plot.cols[i])
      text(x=cat.mean,y=-i+0.35,pos=3,cex=0.6,
           labels=paste(round(100*cat.mean,1),"%",sep=""))
      
      #Add median to right margin
      axis(4,at=-i+0.3,line=-0.8,las=2,tick=F,cex.axis=0.6,
           col.axis=lab.cols[i],labels=median.counts[i])
    }
  })
  
  #Add clean-up box
  box()
}
#Generate size distribution frame with log10 scaling
prepSizePlot <- function(xlims=c(50,1000000),cex.lab=1){
  #Prep plot area
  plot(x=log10(xlims),y=c(0,1),type="n",
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
  axis(1,at=logscale.minor,tick=F,cex.axis=0.8,line=-0.3,las=2,
       labels=logscale.minor.labs)
  mtext(1,text="Size",line=2.25,cex=cex.lab)
  axis(2,at=seq(0,1,0.2),tck=-0.025,labels=NA)
  axis(2,at=seq(0,1,0.2),tick=F,line=-0.4,cex.axis=0.8,las=2,
       labels=paste(seq(0,100,20),"%",sep=""))
  mtext(2,text="De Novo Rate",line=2.2,cex=cex.lab)
  
  #Add cleanup box
  box()
}
#Plot DNRs vs size for all classes
plotDNRvsSize <- function(DNRs,bins,k=4,title=NULL,legend=T,fam.type="families",nfams,cex.lab=1){
  #Prep plot area
  par(mar=c(3.5,3.5,2.5,1))
  prepSizePlot(cex.lab=cex.lab)
  
  #Get midpoints for lines
  mids <- log10(c(bins[1,2],
                  (bins[-c(1,nrow(bins)),1]+bins[-c(1,nrow(bins)),2])/2,
                  bins[nrow(bins),1]))
  
  #Set type colors
  colors <- c("gray15",svtypes$color)
  lwds <- c(3,rep(2,times=nrow(svtypes)))
  
  #Iterate over DNRs and plot per class
  sapply(nrow(DNRs):1,function(i){
    #Get values
    vals <- as.numeric(DNRs[i,])
    
    #Plot line & points
    points(x=mids,y=vals,pch=19,cex=0.4,col=colors[i])
    points(x=mids,
           y=rollapply(vals,k,mean,partial=T,na.rm=T),
           type="l",lwd=lwds[i],col=colors[i])
  })
  
  #Add legend
  if(legend==T){
    idx.for.legend <- which(apply(DNRs,1,function(vals){any(!is.na(vals))}))
    legend("topright",bg="white",pch=19,cex=0.8*cex.lab,lwd=2,
           legend=rownames(DNRs)[idx.for.legend],
           col=colors[idx.for.legend])
  }
  
  #Add title & number of families
  mtext(3,line=0.2,cex=0.8*cex.lab,text=paste("n=",prettyNum(nfams,big.mark=","),
                                              " ",fam.type,"s",sep=""))
  mtext(3,line=1,text=title,font=2,cex=cex.lab)
  
  #Add cleanup box
  box()
}
#Generate freq distribution frame with log10 scaling
prepFreqPlot <- function(xlims=c(1/10000,1),xlabel="Frequency",cex.lab=1){
  #Prep plot area
  plot(x=floor(log10(xlims)),y=c(0,1),type="n",
       xaxt="n",yaxt="n",xlab="",ylab="",yaxs="i")
  
  #Add vertical gridlines
  logscale.all <- log10(as.numeric(sapply(0:8,function(i){(1:9)/10^i})))
  logscale.minor <- log10(as.numeric(sapply(0:8,function(i){c(5,10)/10^i})))
  logscale.major <- log10(as.numeric(1/10^(0:8)))
  abline(v=logscale.all,col="gray97")
  abline(v=logscale.minor,col="gray92")
  abline(v=logscale.major,col="gray85")
  
  #Add axes, title, and Alu/SVA/L1 ticks
  axis(1,at=logscale.all,tck=-0.015,col="gray50",labels=NA)
  axis(1,at=logscale.minor,tck=-0.0225,col="gray20",labels=NA)
  axis(1,at=logscale.major,tck=-0.03,labels=NA)
  for(i in -8:0){
    axis(1,at=i,tick=F,cex.axis=0.8,line=-0.2,
         labels=bquote(10^{.(i)}))
  }
  mtext(1,text=xlabel,line=2.25,cex=cex.lab)
  axis(2,at=seq(0,1,0.2),tck=-0.025,labels=NA)
  axis(2,at=seq(0,1,0.2),tick=F,line=-0.4,cex.axis=0.8,las=2,
       labels=paste(seq(0,100,20),"%",sep=""))
  mtext(2,text="De Novo Rate",line=2.2,cex=cex.lab)
  
  #Add cleanup box
  box()
}
#Plot DNRs vs freq for all classes
plotDNRvsFreq <- function(DNRs,bins,k=4,title=NULL,legend=T,fam.type="familie",nfams,count="variants",cex.lab=1){
  #Get frequency index & x axis title
  if(count=="variants"){
    freq.idx <- which(colnames(dat)=="carrierFreq")
    x.label <- "Carrier Frequency"
  }else{
    freq.idx <- which(colnames(dat)=="AF")
    x.label <- "Allele Frequency"
  }
  
  #Prep plot area
  par(mar=c(3.5,3.5,2.5,1))
  prepFreqPlot(xlims=c(min(dat[,freq.idx],na.rm=T),1),
               xlabel=x.label,cex.lab=cex.lab)
  
  #Get midpoints for lines
  mids <- log10(c(bins[1,2],(bins[-nrow(bins),1]+bins[-1,2])/2))
  
  #Set type colors
  colors <- c("gray15",svtypes$color)
  lwds <- c(3,rep(2,times=nrow(svtypes)))
  
  #Iterate over DNRs and plot per class
  sapply(nrow(DNRs):1,function(i){
    #Get values
    vals <- as.numeric(DNRs[i,])
    
    #Plot line & points
    points(x=mids,y=vals,pch=19,cex=0.4,col=colors[i])
    points(x=mids,
           y=rollapply(vals,k,mean,partial=T,na.rm=T),
           type="l",lwd=lwds[i],col=colors[i])
  })
  
  #Add legend
  if(legend==T){
    idx.for.legend <- which(apply(DNRs,1,function(vals){any(!is.na(vals))}))
    legend("topright",bg="white",pch=19,cex=0.7,lwd=3,
           legend=rownames(DNRs)[idx.for.legend],
           col=colors[idx.for.legend])
  }
  
  #Add title & number of families
  mtext(3,line=0.2,cex=0.8*cex.lab,text=paste("n=",prettyNum(nfams,big.mark=","),
                                              " ",fam.type,"s",sep=""))
  mtext(3,line=1,text=title,font=2,cex=cex.lab)
  
  #Add cleanup box
  box()
}
#Plot DNRs vs GQ for all classes
plotDNRvsGQ <- function(DNRs,bins,k=4,title=NULL,xlabel="Mininum GQ",
                        legend=T,fam.type="familie",nfams,count="variants",cex.lab=1){
  #Get x axis title
  if(count=="variants"){
    x.label <- "Carrier Frequency"
  }else{
    x.label <- "Allele Frequency"
  }
  
  #Prep plot area
  par(mar=c(3.5,3.5,2.5,1))
  plot(x=range(bins),y=c(0,1),type="n",
       xaxt="n",yaxt="n",xlab="",ylab="",yaxs="i")
  
  #Add vertical gridlines
  abline(v=seq(0,1000,50),col="gray92")
  abline(v=seq(0,1000,100),col="gray85")
  
  #Add axes & title
  axis(1,at=seq(0,1000,100),tck=-0.03,labels=NA)
  axis(1,at=seq(0,1000,100),tick=F,cex.axis=0.7*cex.lab,line=-0.4,
       las=2,labels=paste(">",seq(0,1000,100),sep=""))
  mtext(1,text=xlabel,line=2.25,cex=cex.lab)
  axis(2,at=seq(0,1,0.2),tck=-0.025,labels=NA)
  axis(2,at=seq(0,1,0.2),tick=F,line=-0.4,cex.axis=0.8,las=2,
       labels=paste(seq(0,100,20),"%",sep=""))
  mtext(2,text="De Novo Rate",line=2.2,cex=cex.lab)
  
  #Set type colors
  colors <- c("gray15",svtypes$color)
  lwds <- c(3,rep(2,times=nrow(svtypes)))
  
  #Iterate over DNRs and plot per class
  sapply(nrow(DNRs):1,function(i){
    #Get values
    vals <- as.numeric(DNRs[i,])
    
    #Plot line & points
    points(x=bins,y=vals,pch=19,cex=0.4,col=colors[i])
    points(x=bins,
           y=rollapply(vals,k,mean,partial=T,na.rm=T),
           type="l",lwd=lwds[i],col=colors[i])
  })
  
  #Add legend
  if(legend==T){
    idx.for.legend <- which(apply(DNRs,1,function(vals){any(!is.na(vals))}))
    legend("topright",bg="white",pch=19,cex=0.7,lwd=3,
           legend=rownames(DNRs)[idx.for.legend],
           col=colors[idx.for.legend])
  }
  
  #Add title & number of families
  mtext(3,line=0.2,cex=0.8*cex.lab,text=paste("n=",prettyNum(nfams,big.mark=","),
                                              " ",fam.type,"s",sep=""))
  mtext(3,line=1,text=title,font=2,cex=cex.lab)
  
  #Add cleanup box
  box()
}
#Generic heatmap function
plotHeatmap <- function(mat,nfams,fam.type,
                        x.labels=NULL,x.title=NULL,
                        y.labels=NULL,y.title=NULL,
                        title=NULL,cex.lab=1){
  #Set values if NULL
  if(is.null(x.labels)){
    x.labels <- colnames(mat)
  }
  if(is.null(y.labels)){
    y.labels <- rownames(mat)
  }
  
  #Prep plotting area
  par(mar=c(2,4,4,2))
  plot(x=c(0,ncol(mat)),y=c(0,-nrow(mat)),type="n",
       xaxt="n",xaxs="i",xlab="",yaxt="n",yaxs="i",ylab="")
  
  #Add axes
  sapply(1:ncol(mat),function(i){
    axis(3,at=i-0.5,tick=F,line=-0.8,las=2,labels=x.labels[i],cex.axis=0.7)
  })
  # mtext(1,line=2.75,text=x.title,cex=cex.lab)
  axis(2,at=-(1:nrow(mat))+0.5,tick=F,line=-0.8,las=2,labels=y.labels,cex.axis=0.7)
  # mtext(2,line=2.75,text=y.title,cex=cex.lab)
  mtext(1,line=0,cex=0.7*cex.lab,
        text=paste("Median of N=",prettyNum(nfams,big.mark=",")," ",fam.type,"s",sep=""))
  mtext(3,line=2.25,text=title,font=2,cex=cex.lab)
  
  #Plot all cells
  sapply(1:nrow(mat),function(r){
    sapply(1:ncol(mat),function(c){
      #Get color range
      col.range <- colorRampPalette(c("#FFFFFF","#FBDB69","#EF9C4B",
                                      "#E45F30","#8B412B","#000000"))(101)
      
      #Get & scale value
      val <- mat[r,c]
      pct <- round(100*val,0)
      
      #Get color for shading
      if(is.na(val)){
        color <- "gray80"
        label <- "N/A"
        dens <- 12
      }else{
        color <- col.range[pct+1]
        label <- paste(pct,"%",sep="")
        dens <- NA
      }
      #Get text color
      if(is.na(val)){
        text.col <- "gray60"
      }else{
        if(val>0.5){
          text.col <- "white"
        }else{
          text.col <- "black"
        }
      }
      
      #Plot rectangle
      rect(xleft=c-1,xright=c,ybottom=-r,ytop=-(r-1),
           lwd=0.5,border="gray95",col=color,density=dens)
      #Format cell annotation
      text(x=c-0.5,y=-(r-0.5),labels=label,
           cex=0.8,col=text.col)
    })
  })
  
  #Clean up box
  box()
}


############################
###INHERITANCE PLOT WRAPPERS
############################
#Wrapper for all standard inheritance plots
wrapperInheritancePlots <- function(fam.dat.list,fam.type,count="variants"){
  #Set title prefix & suffix and freq filter index
  if(count=="variants"){
    title.prefix <- "SV Site "
    freq.idx <- which(colnames(dat)=="carrierFreq")
    freq.lab <- "CF"
  }else{
    title.prefix <- "SV Allele "
    freq.idx <- which(colnames(dat)=="AF")
    freq.lab <- "AF"
  }
  if(fam.type=="trio"){
    title.suffix <- paste("(Trios; n=",
                          prettyNum(length(fam.dat.list),big.mark=","),
                          ")",sep="")
  }else{
    if(fam.type=="duo"){
      title.suffix <- paste("(Duos; n=",
                            prettyNum(length(fam.dat.list),big.mark=","),
                            ")",sep="")
    }else{
      title.suffix <- paste("(Families; n=",
                            prettyNum(length(fam.dat.list),big.mark=","),
                            ")",sep="")
    }
  }
  
  #All variants
  pdf(paste(OUTDIR,"/supporting_plots/sv_inheritance_plots/sv_inheritance.",fam.type,"s.",count,".all_sv.pdf",sep=""),
      height=3.75,width=4.5)
  plotInhStats(inh.stats=computeInheritanceMulti(trio.dat.list=fam.dat.list,
                                                 VIDs=NULL),
               title=paste(title.prefix,"Inheritance [All SV] ",title.suffix,sep=""),
               count=count)
  dev.off()  
  
  #Variants by class
  sapply(svtypes$svtype,function(svtype){
    pdf(paste(OUTDIR,"/supporting_plots/sv_inheritance_plots/sv_inheritance.",fam.type,"s.",count,".",svtype,".pdf",sep=""),
        height=3.75,width=4.5)
    plotInhStats(inh.stats=computeInheritanceMulti(trio.dat.list=fam.dat.list,
                                                   VIDs=dat$VID[which(dat$svtype==svtype)]),
                 title=paste(title.prefix,"Inheritance [",svtype,"] ",title.suffix,sep=""),
                 count=count)
    dev.off()  
  })
  
  #Tiny
  pdf(paste(OUTDIR,"/supporting_plots/sv_inheritance_plots/sv_inheritance.",fam.type,"s.",count,".tiny_sv.pdf",sep=""),
      height=3.75,width=4.5)
  plotInhStats(inh.stats=computeInheritanceMulti(trio.dat.list=fam.dat.list,
                                                 VIDs=dat$VID[which(dat$length<=tiny.max.size)]),
               title=paste(title.prefix,"Inheritance [<100bp] ",title.suffix,sep=""),
               count=count)
  dev.off()  
  #Small
  pdf(paste(OUTDIR,"/supporting_plots/sv_inheritance_plots/sv_inheritance.",fam.type,"s.",count,".small_sv.pdf",sep=""),
      height=3.75,width=4.5)
  plotInhStats(inh.stats=computeInheritanceMulti(trio.dat.list=fam.dat.list,
                                                 VIDs=dat$VID[which(dat$length>tiny.max.size & dat$length<=small.max.size)]),
               title=paste(title.prefix,"Inheritance [100-500bp] ",title.suffix,sep=""),
               count=count)
  dev.off()  
  #Medium
  pdf(paste(OUTDIR,"/supporting_plots/sv_inheritance_plots/sv_inheritance.",fam.type,"s.",count,".medium_sv.pdf",sep=""),
      height=3.75,width=4.5)
  plotInhStats(inh.stats=computeInheritanceMulti(trio.dat.list=fam.dat.list,
                                                 VIDs=dat$VID[which(dat$length>small.max.size & dat$length<=medium.max.size)]),
               title=paste(title.prefix,"Inheritance [500bp-2.5kb] ",title.suffix,sep=""),
               count=count)
  dev.off()  
  #Med-large
  pdf(paste(OUTDIR,"/supporting_plots/sv_inheritance_plots/sv_inheritance.",fam.type,"s.",count,".medlarge_sv.pdf",sep=""),
      height=3.75,width=4.5)
  plotInhStats(inh.stats=computeInheritanceMulti(trio.dat.list=fam.dat.list,
                                                 VIDs=dat$VID[which(dat$length>medium.max.size & dat$length<=medlarge.max.size)]),
               title=paste(title.prefix,"Inheritance [2.5-10kb] ",title.suffix,sep=""),
               count=count)
  dev.off()  
  #Large
  pdf(paste(OUTDIR,"/supporting_plots/sv_inheritance_plots/sv_inheritance.",fam.type,"s.",count,".large_sv.pdf",sep=""),
      height=3.75,width=4.5)
  plotInhStats(inh.stats=computeInheritanceMulti(trio.dat.list=fam.dat.list,
                                                 VIDs=dat$VID[which(dat$length>medlarge.max.size & dat$length<=large.max.size)]),
               title=paste(title.prefix,"Inheritance [10-50kb] ",title.suffix,sep=""),
               count=count)
  dev.off()  
  #Huge
  pdf(paste(OUTDIR,"/supporting_plots/sv_inheritance_plots/sv_inheritance.",fam.type,"s.",count,".huge_sv.pdf",sep=""),
      height=3.75,width=4.5)
  plotInhStats(inh.stats=computeInheritanceMulti(trio.dat.list=fam.dat.list,
                                                 VIDs=dat$VID[which(dat$length>large.max.size)]),
               title=paste(title.prefix,"Inheritance [>50kb] ",title.suffix,sep=""),
               count=count)
  dev.off()  
  
  #Rare
  pdf(paste(OUTDIR,"/supporting_plots/sv_inheritance_plots/sv_inheritance.",fam.type,"s.",count,".rare_sv.pdf",sep=""),
      height=3.75,width=4.5)
  plotInhStats(inh.stats=computeInheritanceMulti(trio.dat.list=fam.dat.list,
                                                 VIDs=dat$VID[which(dat[,freq.idx]<=rare.max.freq)]),
               title=paste(title.prefix,"Inheritance [",freq.lab,"<1%] ",title.suffix,sep=""),
               count=count)
  dev.off()  
  #Uncommon
  pdf(paste(OUTDIR,"/supporting_plots/sv_inheritance_plots/sv_inheritance.",fam.type,"s.",count,".uncommon_sv.pdf",sep=""),
      height=3.75,width=4.5)
  plotInhStats(inh.stats=computeInheritanceMulti(trio.dat.list=fam.dat.list,
                                                 VIDs=dat$VID[which(dat[,freq.idx]>rare.max.freq & dat[,freq.idx]<=uncommon.max.freq)]),
               title=paste(title.prefix,"Inheritance [",freq.lab," 1-10%] ",title.suffix,sep=""),
               count=count)
  dev.off()  
  #Common
  pdf(paste(OUTDIR,"/supporting_plots/sv_inheritance_plots/sv_inheritance.",fam.type,"s.",count,".common_sv.pdf",sep=""),
      height=3.75,width=4.5)
  plotInhStats(inh.stats=computeInheritanceMulti(trio.dat.list=fam.dat.list,
                                                 VIDs=dat$VID[which(dat[,freq.idx]>uncommon.max.freq & dat[,freq.idx]<=common.max.freq)]),
               title=paste(title.prefix,"Inheritance [",freq.lab," 10-50%] ",title.suffix,sep=""),
               count=count)
  dev.off()  
  #Major
  pdf(paste(OUTDIR,"/supporting_plots/sv_inheritance_plots/sv_inheritance.",fam.type,"s.",count,".major_sv.pdf",sep=""),
      height=3.75,width=4.5)
  plotInhStats(inh.stats=computeInheritanceMulti(trio.dat.list=fam.dat.list,
                                                 VIDs=dat$VID[which(dat[,freq.idx]>common.max.freq)]),
               title=paste(title.prefix,"Inheritance [",freq.lab,">50%] ",title.suffix,sep=""),
               count=count)
  dev.off()
}
#Wrapper for de novo rate lineplots
wrapperDeNovoRateLines <- function(fam.dat.list,fam.type,count="variants",gq=F){
  #Set title prefix
  if(count=="variants"){
    title.prefix <- "SV Site "
  }else{
    title.prefix <- "SV Allele "
  }
  
  #DNR by Size
  size.dat <- deNovoRateBySize(trio.dat.list=fam.dat.list,size.bins=40,count=count)
  pdf(paste(OUTDIR,"/supporting_plots/sv_inheritance_plots/sv_de_novo_rate.",fam.type,"s.",count,".by_size.pdf",sep=""),
      height=4,width=5)
  plotDNRvsSize(DNRs=size.dat$DNRs,bins=size.dat$bins,k=4,nfams=length(fam.dat.list),
                title=paste(title.prefix,"De Novo Rate by Size",sep=""),
                fam.type=fam.type,legend=T)
  dev.off()  
  
  #DNR by Freq
  freq.dat <- deNovoRateByFreq(trio.dat.list=fam.dat.list,freq.bins=40,count=count)
  pdf(paste(OUTDIR,"/supporting_plots/sv_inheritance_plots/sv_de_novo_rate.",fam.type,"s.",count,".by_frequency.pdf",sep=""),
      height=4,width=5)
  plotDNRvsFreq(DNRs=freq.dat$DNRs,bins=freq.dat$bins,k=4,nfams=length(fam.dat.list),
                title=paste(title.prefix,"De Novo Rate by Freq.",sep=""),
                count=count,fam.type=fam.type,legend=T)
  dev.off()
  
  #DNR by Proband GQ
  if(gq) {
    GQ.dat <- deNovoRateByProGQ(trio.dat.list=fam.dat.list,GQ.bins=40,count=count)
    pdf(paste(OUTDIR,"/supporting_plots/sv_inheritance_plots/sv_de_novo_rate.",fam.type,"s.",count,".by_proband_GQ.pdf",sep=""),
        height=4,width=5)
    plotDNRvsGQ(DNRs=GQ.dat$DNRs,bins=GQ.dat$bins,k=4,nfams=length(fam.dat.list),
                title=paste(title.prefix,"De Novo Rate by Min. Proband GQ",sep=""),
                count=count,fam.type=fam.type,legend=T,xlab="Min. Proband GQ")
    dev.off()
  }
}
#Wrapper for de novo rate heatmaps
wrapperDeNovoRateHeats <- function(fam.dat.list,fam.type,count="variants"){
  #Set title prefix
  if(count=="variants"){
    title.prefix <- "SV Site "
  }else{
    title.prefix <- "SV Allele "
  }
  
  #Gather DNR data
  DNR.dat <- deNovoRateBySizeFreq(trio.dat.list=fam.dat.list,count=count,
                                  max.sizes=c(tiny.max.size,small.max.size,medium.max.size,
                                              medlarge.max.size,large.max.size),
                                  size.labs=c("<100bp","100-\n500bp","500bp-\n2.5kb",
                                              "2.5-10kb","10kb-50kb",">50kb"),
                                  max.freqs=c(0.01,0.05,0.10,0.50),
                                  freq.labs=c("<1%","1-5%","5-10%","10-50%",">50%"))
  
  #Plot one heatmap for all variants
  pdf(paste(OUTDIR,"/supporting_plots/sv_inheritance_plots/sv_de_novo_rate.",fam.type,"s.",count,".size_vs_freq.all_sv.pdf",sep=""),
      height=5,width=5)
  plotHeatmap(mat=DNR.dat$ALL,nfams=length(fam.dat.list),fam.type=fam.type,
              title=paste(title.prefix,"De Novo Rate, Size vs. Freq. [All SV]",sep=""))
  dev.off()
  
  #Plot one heatmap per variant class
  sapply(svtypes$svtype,function(svtype){
    pdf(paste(OUTDIR,"/supporting_plots/sv_inheritance_plots/sv_de_novo_rate.",fam.type,"s.",count,".size_vs_freq.",svtype,".pdf",sep=""),
        height=5,width=5)
    plotHeatmap(mat=DNR.dat[[which(names(DNR.dat)==svtype)]],nfams=length(fam.dat.list),fam.type=fam.type,
                title=paste(title.prefix,"De Novo Rate, Size vs. Freq. [",svtype,"]",sep=""))
    dev.off()
  })
}
#Wrapper for master summary panel
masterInhWrapper <- function(fam.dat.list,fam.type, gq=T){
  #Set title suffix
  if(fam.type=="trio"){
    title.suffix <- paste("(Trios; n=",
                          prettyNum(length(fam.dat.list),big.mark=","),
                          ")",sep="")
  }else{
    if(fam.type=="duo"){
      title.suffix <- paste("(Duos; n=",
                            prettyNum(length(fam.dat.list),big.mark=","),
                            ")",sep="")
    }else{
      title.suffix <- paste("(Families; n=",
                            prettyNum(length(fam.dat.list),big.mark=","),
                            ")",sep="")
    }
  }
  
  #Prepare plot area
  width <- ifelse(gq, 12, 10)
  pdf(paste(OUTDIR,"/main_plots/VCF_QC.SV_",fam.type,"_inheritance.pdf",sep=""),
      height=5,width=width)
  if(gq) {
    layout(matrix(c(1,2,3,4,5,
                    6,7,8,9,10),
                  byrow=T,nrow=2))
  } else {
    layout(matrix(c(1,2,3,4,
                    5,6,7,8),
                  byrow=T,nrow=2))
  }
  
  
  #Set global cex.lab
  cex.lab <- 0.75
  
  ###Top row: SV sites
  #Master inheritance plot
  plotInhStats(inh.stats=computeInheritanceMulti(trio.dat.list=fam.dat.list),
               title=paste("SV Site Inheritance (n=",
                           prettyNum(length(fam.dat.list),big.mark=","),
                           " ",fam.type,"s)",sep=""),
               count="variants",cex.lab=cex.lab)
  #DNR vs size
  size.dat.v <- deNovoRateBySize(trio.dat.list=fam.dat.list,size.bins=40,count="variants")
  plotDNRvsSize(DNRs=size.dat.v$DNRs,bins=size.dat.v$bins,k=4,nfams=length(fam.dat.list),
                title=paste("Site De Novo Rate by Size",sep=""),
                fam.type=fam.type,legend=T,cex.lab=cex.lab)
  #DNR vs frequency
  freq.dat.v <- deNovoRateByFreq(trio.dat.list=fam.dat.list,freq.bins=40,count="variants")
  plotDNRvsFreq(DNRs=freq.dat.v$DNRs,bins=freq.dat.v$bins,k=4,nfams=length(fam.dat.list),
                title=paste("Site De Novo Rate by Freq.",sep=""),
                count="variants",fam.type=fam.type,legend=F,cex.lab=cex.lab)
  #DNR vs min proband GQ
  if(gq) {
    GQ.dat.v <- deNovoRateByProGQ(trio.dat.list=fam.dat.list,GQ.bins=40,count="variants")
    plotDNRvsGQ(DNRs=GQ.dat.v$DNRs,bins=GQ.dat.v$bins,k=4,nfams=length(fam.dat.list),
                title=paste("Site De Novo Rate by GQ",sep=""),
                count="variants",fam.type=fam.type,legend=F,cex.lab=cex.lab,
                xlab="Min. Proband GQ")
  }
  #DNR heatmap (size vs freq.)
  DNR.dat.v <- deNovoRateBySizeFreq(trio.dat.list=fam.dat.list,count="variants",
                                    max.sizes=c(tiny.max.size,small.max.size,medium.max.size,
                                                medlarge.max.size,large.max.size),
                                    size.labs=c("<100bp","100-\n500bp","500bp-\n2.5kb",
                                                "2.5-10kb","10kb-50kb",">50kb"),
                                    max.freqs=c(0.01,0.05,0.10,0.50),
                                    freq.labs=c("<1%","1-5%","5-\n10%","10-\n50%",">50%"))
  plotHeatmap(mat=DNR.dat.v$ALL,nfams=length(fam.dat.list),fam.type=fam.type,
              title=paste("Site De Novo Rate, Size vs. Freq.",sep=""),cex.lab=cex.lab)
  
  ###Bottom row: SV alleles
  #Master inheritance plot
  plotInhStats(inh.stats=computeInheritanceMulti(trio.dat.list=fam.dat.list),
               title=paste("SV Allele Inheritance (n=",
                           prettyNum(length(fam.dat.list),big.mark=","),
                           " ",fam.type,"s)",sep=""),
               count="alleles",cex.lab=cex.lab)
  #DNR vs size
  size.dat.a <- deNovoRateBySize(trio.dat.list=fam.dat.list,size.bins=40,count="alleles")
  plotDNRvsSize(DNRs=size.dat.a$DNRs,bins=size.dat.a$bins,k=4,nfams=length(fam.dat.list),
                title=paste("Allele De Novo Rate by Size",sep=""),
                fam.type=fam.type,legend=F,cex.lab=cex.lab)
  #DNR vs frequency
  freq.dat.a <- deNovoRateByFreq(trio.dat.list=fam.dat.list,freq.bins=40,count="alleles")
  plotDNRvsFreq(DNRs=freq.dat.a$DNRs,bins=freq.dat.a$bins,k=4,nfams=length(fam.dat.list),
                title=paste("Allele De Novo Rate by Freq.",sep=""),
                count="alleles",fam.type=fam.type,legend=F,cex.lab=cex.lab)
  #DNR vs min proband GQ
  if(gq){
    GQ.dat.a <- deNovoRateByProGQ(trio.dat.list=fam.dat.list,GQ.bins=40,count="alleles")
    plotDNRvsGQ(DNRs=GQ.dat.a$DNRs,bins=GQ.dat.a$bins,k=4,nfams=length(fam.dat.list),
                title=paste("Allele De Novo Rate by GQ",sep=""),
                count="alleles",fam.type=fam.type,legend=F,cex.lab=cex.lab,
                xlabel="Min. Proband GQ")
  }
  #DNR heatmap (size vs freq.)
  DNR.dat.a <- deNovoRateBySizeFreq(trio.dat.list=fam.dat.list,count="alleles",
                                    max.sizes=c(tiny.max.size,small.max.size,medium.max.size,
                                                medlarge.max.size,large.max.size),
                                    size.labs=c("<100bp","100-\n500bp","500bp-\n2.5kb",
                                                "2.5-10kb","10kb-50kb",">50kb"),
                                    max.freqs=c(0.01,0.05,0.10,0.50),
                                    freq.labs=c("<1%","1-5%","5-\n10%","10-\n50%",">50%"))
  plotHeatmap(mat=DNR.dat.a$ALL,nfams=length(fam.dat.list),fam.type=fam.type,
              title=paste("Allele De Novo Rate, Size vs. Freq.",sep=""),cex.lab=cex.lab)
  
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
require(zoo)

###List of command-line options
option_list <- list(
  make_option(c("-S", "--svtypes"), type="character", default=NULL,
              help="tab-delimited file specifying SV types and HEX colors [default %default]",
              metavar="character"),
  make_option(c("-M", "--multiallelics"), type="logical", default=FALSE,
              help="include multiallelic sites in inheritance calculations [default %default]",
              metavar="logical")
)

###Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog svstats.bed famfile perSampleDir OUTDIR",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

###Checks for appropriate positional arguments
if(length(args$args) != 4){
  stop("Incorrect number of required positional arguments\n")
}

###Writes args & opts to vars
dat.in <- args$args[1]
famfile.in <- args$args[2]
perSampDir <- args$args[3]
OUTDIR <- args$args[4]
svtypes.file <- opts$svtypes
multiallelics <- opts$multiallelics

# #Dev parameters
# dat.in <- "~/scratch/xfer/gnomAD_v2_SV_MASTER_resolved_VCF.VCF_sites.stats.bed.gz"
# famfile.in <- "~/scratch/xfer/cleaned.fam"
# perSampDir <- "~/scratch/xfer/gnomAD_v2_SV_MASTER_resolved_VCF_perSample_VIDs_merged/"
# OUTDIR <- "~/scratch/famQC_plots_test/"
# # OUTDIR <- "~/scratch/VCF_plots_test/"
# svtypes.file <- "~/Desktop/Collins/Talkowski/code/sv-pipeline/ref/vcf_qc_refs/SV_colors.txt"
# multiallelics <- F

###Prepares I/O files
#Read & clean SV stats data
dat <- read.table(dat.in,comment.char="",sep="\t",header=T,check.names=F)
colnames(dat)[1] <- "chr"
#Restrict data to autosomes only, and exclude multiallelics (if optioned)
allosome.exclude.idx <- which(!(dat$chr %in% c(1:22,paste("chr",1:22,sep=""))))
multi.exclude.idx <- which(dat$other_gts>0)
cat(paste("NOTE: only autosomes considered during transmission analyses. Excluded ",
          prettyNum(length(allosome.exclude.idx),big.mark=","),"/",
          prettyNum(nrow(dat),big.mark=",")," (",
          round(100*length(allosome.exclude.idx)/nrow(dat),1),
          "%) of all variants as non-autosomal.\n",sep=""))
if(multiallelics==F){
  cat(paste("NOTE: only biallelic variants considered during transmission analyses. Excluded ",
            prettyNum(length(multi.exclude.idx),big.mark=","),"/",
            prettyNum(nrow(dat),big.mark=",")," (",
            round(100*length(multi.exclude.idx)/nrow(dat),1),
            "%) of all variants as multiallelic.\n",sep=""))
  all.exclude.idx <- unique(c(allosome.exclude.idx,multi.exclude.idx))
  cat(paste("NOTE: excluded a nonredundant total of ",
            prettyNum(length(all.exclude.idx),big.mark=","),"/",
            prettyNum(nrow(dat),big.mark=",")," (",
            round(100*length(all.exclude.idx)/nrow(dat),1),
            "%) of all variants due to autosomal and/or multiallelic filters.\n",sep=""))
  
  dat <- dat[-all.exclude.idx,]
}else{
  dat <- dat[-allosome.exclude.idx,]
}
cat(paste("NOTE: retained ",
          prettyNum(nrow(dat),big.mark=","),
          " variants for transmission analyses.\n",sep=""))
#Read fam file and splits into duos and trios
fams <- read.table(famfile.in,comment.char="",header=T,check.names=F)
colnames(fams)[1] <- "FAM_ID"
# duos <- fams[grep("DUO_",fams$FAM_ID,fixed=T),]
# duos$FATHER[which(duos$FATHER==".")] <- NA
# duos$MOTHER[which(duos$MOTHER==".")] <- NA
trios <- fams[grep("TRIO_",fams$FAM_ID,fixed=T),]
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
if(!dir.exists(paste(OUTDIR,"/supporting_plots/sv_inheritance_plots/",sep=""))){
  dir.create(paste(OUTDIR,"/supporting_plots/sv_inheritance_plots/",sep=""))
}

###Performs trio analyses, if any trios exist
if(nrow(trios)>0){
  #Read data
  trio.dat <- apply(trios[,2:4],1,function(IDs){
    IDs <- as.character(IDs)
    return(getFamDat(dat=dat,proband=IDs[1],father=IDs[2],mother=IDs[3],biallelic=!multiallelics))
  })
  names(trio.dat) <- trios[,1]
  
  # if there are no GQ values in any of the trios, do not make GQ plots
  gq <- any(unlist(lapply(trio.dat, function(trio){ sum(!is.na(trio$pro.GQ)) > 0 })))
  
  #Master wrapper
  masterInhWrapper(fam.dat.list=trio.dat,fam.type="trio", gq=gq)
  #Standard inheritance panels
  sapply(c("variants","alleles"),function(count){
    wrapperInheritancePlots(fam.dat.list=trio.dat,
                            fam.type="trio",
                            count=count)  
  })
  #De novo rate panels
  sapply(c("variants","alleles"),function(count){
    wrapperDeNovoRateLines(fam.dat.list=trio.dat,
                           fam.type="trio",
                           count=count,
                           gq=gq)  
  })
  #De novo rate heatmaps
  sapply(c("variants","alleles"),function(count){
    wrapperDeNovoRateHeats(fam.dat.list=trio.dat,
                           fam.type="trio",
                           count=count)  
  })
}


