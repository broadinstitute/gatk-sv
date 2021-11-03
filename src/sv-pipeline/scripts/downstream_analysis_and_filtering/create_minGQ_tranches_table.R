#!/usr/bin/env Rscript

# Helper script to generate table of minGQ filtering conditions


###Set master parameters
options(stringsAsFactors=F,scipen=1000)


###################
###HELPER FUNCTIONS
###################
###Format size range table
formatSizes <- function(min.sizes,max.sizes){
  #Vectorize inputs
  min.vect <- as.integer(unlist(strsplit(min.sizes,split=";")))
  # min.vect[which(min.vect==0)] <- -2
  max.vect <- as.integer(unlist(strsplit(max.sizes,split=";")))
  
  #Sanity check correct formatting of vectorized inputs
  if(length(min.vect) != length(max.vect)){
    stop(paste("Number of minimum and maximum SVLEN cutoffs must match.\n",
               "   --min.sizes as parsed:  ",paste(prettyNum(min.vect,big.mark=","),collapse="  "),"\n",
               "   --max.sizes as parsed:  ",paste(prettyNum(max.vect,big.mark=","),collapse="  "),
               sep=""))
  }
  if(any(is.na(min.vect)) | any(is.nan(min.vect)) | any(is.infinite(min.vect))){
    stop(paste("Minimum SVLEN cutoffs must all be integers.\n",
               "   --min.sizes as parsed:  ",paste(prettyNum(min.vect,big.mark=","),collapse="  "),"\n",
               sep=""))
  }
  if(any(is.na(max.vect)) | any(is.nan(max.vect)) | any(is.infinite(max.vect))){
    stop(paste("maximum SVLEN cutoffs must all be integers.\n",
               "   --max.sizes as parsed:  ",paste(prettyNum(max.vect,big.mark=","),collapse="  "),
               sep=""))
  }
  
  #Create size range table
  sizes.df <- data.frame("minSVLEN"=min.vect,"maxSVLEN"=max.vect)
  size.ranges <- paste("\n   ",paste(prettyNum(min.vect,big.mark=","),
                                     prettyNum(max.vect,big.mark=","),
                                     sep="-"),sep="")
  
  #Sanity check size range table
  apply(sizes.df,1,function(vals){
    if(vals[1]>vals[2]){
      stop(paste("Minimum SVLEN cutoffs must be <= their corresponding maximum SVLEN cutoffs.\n",
                 "   SVLEN ranges as parsed [n=",prettyNum(length(size.ranges),big.mark=","),
                 "]:  ",paste(size.ranges,collapse=""),sep=""))
    }
  })
  
  #Return results
  cat(paste("\nSVLEN ranges as parsed [n=",prettyNum(length(size.ranges),big.mark=","),
            "]:  ",paste(size.ranges,collapse=""),"\n",sep=""))
  return(sizes.df)
}

###Format freq range table
formatFreqs <- function(min.freqs,max.freqs){
  #Vectorize inputs
  min.vect <- as.numeric((unlist(strsplit(min.freqs,split=";"))))
  max.vect <- as.numeric((unlist(strsplit(max.freqs,split=";"))))
  
  #Sanity check correct formatting of vectorized inputs
  if(length(min.vect) != length(max.vect)){
    stop(paste("Number of minimum and maximum AF cutoffs must match.\n",
               "   --min.freqs as parsed:  ",paste(paste(round(100*min.vect,2),"%",sep=""),collapse="  "),"\n",
               "   --max.freqs as parsed:  ",paste(paste(round(100*max.vect,2),"%",sep=""),collapse="  "),
               sep=""))
  }
  if(any(is.na(min.vect)) | any(is.nan(min.vect))){
    stop(paste("Minimum AF cutoffs must all be numeric\n",
               "   --min.freqs as parsed:  ",paste(paste(round(100*min.vect,2),"%",sep=""),collapse="  "),"\n",
               sep=""))
  }
  if(any(is.na(max.vect)) | any(is.nan(max.vect))){
    stop(paste("maximum AF cutoffs must all be numeric\n",
               "   --max.freqs as parsed:  ",paste(paste(round(100*max.vect,2),"%",sep=""),collapse="  "),
               sep=""))
  }
  
  #Create freq range table
  freqs.df <- data.frame("minAF"=min.vect,"maxAF"=max.vect)
  freq.ranges <- paste("\n   ",paste(round(100*min.vect,2),
                                     paste(round(100*max.vect,2),"%",sep=""),
                                     sep="-"),sep="")
  
  #Sanity check freq range table
  apply(freqs.df,1,function(vals){
    if(vals[1]>vals[2]){
      stop(paste("Minimum AF cutoffs must be <= their corresponding maximum AF cutoffs.\n",
                 "   AF ranges as parsed [n=",prettyNum(length(freq.ranges),big.mark=","),
                 "]:  ",paste(freq.ranges,collapse=""),sep=""))
    }
  })
  
  #Return results
  cat(paste("\nAF ranges as parsed [n=",prettyNum(length(freq.ranges),big.mark=","),
            "]:  ",paste(freq.ranges,collapse=""),"\n",sep=""))
  return(freqs.df)
}

###Format SVTYPE table
formatSVTYPEs <- function(include,exclude){
  #Vectorize inputs
  if(!is.null(include)){
    include.vect <- unlist(strsplit(include,split=";"))
    include.vect[which(include.vect=="")] <- "None"
  }
  if(!is.null(exclude)){
    exclude.vect <- unlist(strsplit(exclude,split=";"))
    exclude.vect[which(exclude.vect=="")] <- "None"
  }
  
  if(is.null(include)){
    include.vect <- rep("None",times=length(exclude.vect))
  }
  if(is.null(exclude)){
    exclude.vect <- rep("None",times=length(include.vect))
  }
  
  #Sanity check correct formatting of vectorized inputs
  if(length(include.vect) != length(exclude.vect)){
    stop(paste("Number of inclusion and exclusion SVTYPEs must match.\n",
               "   --svtype.include as parsed:  ",paste(include.vect,collapse="  "),"\n",
               "   --svtype.exclude as parsed:  ",paste(exclude.vect,collapse="  "),
               sep=""))
  }
  if(any(is.na(include.vect))){
    stop(paste("At least one SVTYPE to include is malformed.\n",
               "   --svtype.include as parsed:  ",paste(include.vect,collapse="  "),"\n",
               sep=""))
  }
  if(any(is.na(exclude.vect))){
    stop(paste("At least one SVTYPE to include is malformed.\n",
               "   --svtype.exclude as parsed:  ",paste(exclude.vect,collapse="  "),
               sep=""))
  }
  
  #Create SVTYPE table
  svtypes.df <- data.frame("includeSVTYPE"=include.vect,
                           "excludeSVTYPE"=exclude.vect)
  svtype.pairings <- apply(svtypes.df,1,function(vals){
    return(paste("\n  -Include",vals[1],"& Exclude",vals[2],sep=" "))
  })
  
  #Return results
  cat(paste("\nSVTYPE pairings as parsed [n=",prettyNum(length(svtype.pairings),big.mark=","),
            "]:  ",paste(svtype.pairings,collapse=""),"\n",sep=""))
  return(svtypes.df)
}

###Format FILTER table
formatFILTERs <- function(include,exclude){
  #Vectorize inputs
  if(!is.null(include)){
    include.vect <- unlist(strsplit(include,split=";"))
    include.vect[which(include.vect=="")] <- "None"
  }
  if(!is.null(exclude)){
    exclude.vect <- unlist(strsplit(exclude,split=";"))
    exclude.vect[which(exclude.vect=="")] <- "None"
  }
  
  if(is.null(include)){
    include.vect <- rep("None",times=length(exclude.vect))
  }
  if(is.null(exclude)){
    exclude.vect <- rep("None",times=length(include.vect))
  }
  
  #Sanity check correct formatting of vectorized inputs
  if(length(include.vect) != length(exclude.vect)){
    stop(paste("Number of inclusion and exclusion VCF FILTER statuses must match.\n",
               "   --filter.include as parsed:  ",paste(include.vect,collapse="  "),"\n",
               "   --filter.exclude as parsed:  ",paste(exclude.vect,collapse="  "),
               sep=""))
  }
  if(any(is.na(include.vect))){
    stop(paste("At least one VCF FILTER status to include is malformed.\n",
               "   --filter.include as parsed:  ",paste(include.vect,collapse="  "),"\n",
               sep=""))
  }
  if(any(is.na(exclude.vect))){
    stop(paste("At least one VCF FILTER status to include is malformed.\n",
               "   --filter.exclude as parsed:  ",paste(exclude.vect,collapse="  "),
               sep=""))
  }
  
  #Create FILTER table
  filters.df <- data.frame("includeFILTER"=include.vect,
                           "excludeFILTER"=exclude.vect)
  filter.pairings <- apply(filters.df,1,function(vals){
    return(paste("\n  -Include",vals[1],"& Exclude",vals[2],sep=" "))
  })
  
  #Return results
  cat(paste("\nVCF FILTER status pairings as parsed [n=",prettyNum(length(filter.pairings),big.mark=","),
            "]:  ",paste(filter.pairings,collapse=""),"\n",sep=""))
  return(filters.df)
}

###Format EV table
formatEVs <- function(include,exclude){
  #Vectorize inputs
  if(!is.null(include)){
    include.vect <- unlist(strsplit(include,split=";"))
    include.vect[which(include.vect=="")] <- "None"
  }
  if(!is.null(exclude)){
    exclude.vect <- unlist(strsplit(exclude,split=";"))
    exclude.vect[which(exclude.vect=="")] <- "None"
  }
  
  if(is.null(include)){
    include.vect <- rep("None",times=length(exclude.vect))
  }
  if(is.null(exclude)){
    exclude.vect <- rep("None",times=length(include.vect))
  }
  
  #Sanity check correct formatting of vectorized inputs
  if(length(include.vect) != length(exclude.vect)){
    stop(paste("Number of inclusion and exclusion per-sample genotype EVs must match.\n",
               "   --ev.include as parsed:  ",paste(include.vect,collapse="  "),"\n",
               "   --ev.exclude as parsed:  ",paste(exclude.vect,collapse="  "),
               sep=""))
  }
  if(any(is.na(include.vect))){
    stop(paste("At least one per-sample genotype EV to include is malformed.\n",
               "   --ev.include as parsed:  ",paste(include.vect,collapse="  "),"\n",
               sep=""))
  }
  if(any(is.na(exclude.vect))){
    stop(paste("At least one per-sample genotype EV to include is malformed.\n",
               "   --ev.exclude as parsed:  ",paste(exclude.vect,collapse="  "),
               sep=""))
  }
  
  #Create EV table
  evs.df <- data.frame("includeEV"=include.vect,
                       "excludeEV"=exclude.vect)
  ev.pairings <- apply(evs.df,1,function(vals){
    return(paste("\n  -Include",vals[1],"& Exclude",vals[2],sep=" "))
  })
  
  #Return results
  cat(paste("\nPer-sample genotype EV pairings as parsed [n=",prettyNum(length(ev.pairings),big.mark=","),
            "]:  ",paste(ev.pairings,collapse=""),"\n",sep=""))
  return(evs.df)
}


################
###RSCRIPT BLOCK
################
require(optparse,quietly=T)
###List of command-line options
option_list <- list(
  make_option(c("--min.sizes"), type="character", default=NULL,
              help="minimum SVLEN cutoffs [default %default]"),
  make_option(c("--max.sizes"), type="character", default=NULL,
              help="maximum SVLEN cutoffs [default %default]"),
  make_option(c("--min.freqs"), type="character", default=NULL,
              help="minimum AF cutoffs [default %default]"),
  make_option(c("--max.freqs"), type="character", default=NULL,
              help="maximum AF cutoffs [default %default]"),
  make_option(c("--svtype.include"), type="character", default=NULL,
              help="SVTYPE to include [default %default]"),
  make_option(c("--svtype.exclude"), type="character", default=NULL,
              help="SVTYPE to exclude [default %default]"),
  make_option(c("--filter.include"), type="character", default=NULL,
              help="VCF FILTER status to include [default %default]"),
  make_option(c("--filter.exclude"), type="character", default=NULL,
              help="VCF FILTER status to exclude [default %default]"),
  make_option(c("--ev.include"), type="character", default=NULL,
              help="Per-sample genotype EV to include [default %default]"),
  make_option(c("--ev.exclude"), type="character", default=NULL,
              help="Per-sample genotype EV to exclude [default %default]")
)

###Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog OUTFILE",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

###Checks for appropriate positional arguments
if(length(args$args) != 1){
  stop("Incorrect number of required positional arguments\n")
}

###Writes args & opts to vars
OUTFILE <- args$args[1]

min.sizes <- opts$min.sizes
max.sizes <- opts$max.sizes
if(is.null(min.sizes) & is.null(max.sizes)){
  sizes.df <- data.frame("minSVLEN"=-2,"maxSVLEN"=400000000)
  cat("\nNo SVLEN restrictions detected.\n")
}else{
  sizes.df <- suppressWarnings(formatSizes(min.sizes,max.sizes))
}

min.freqs <- opts$min.freqs
max.freqs <- opts$max.freqs
if(is.null(min.freqs) & is.null(max.freqs)){
  freqs.df <- data.frame("minAF"=0,"maxAF"=1)
  cat("\nNo AF restrictions detected.\n")
}else{
  freqs.df <- suppressWarnings(formatFreqs(min.freqs,max.freqs))
}

svtype.include <- opts$svtype.include
svtype.exclude <- opts$svtype.exclude
if(is.null(svtype.include) & is.null(svtype.exclude)){
  svtypes.df <- data.frame("includeSVTYPE"="None","excludeSVTYPE"="None")
  cat("\nNo SVTYPE restrictions detected.\n")
}else{
  svtypes.df <- suppressWarnings(formatSVTYPEs(svtype.include,svtype.exclude))
}

filter.include <- opts$filter.include
filter.exclude <- opts$filter.exclude
if(is.null(filter.include) & is.null(filter.exclude)){
  filters.df <- data.frame("includeFILTER"="None","excludeFILTER"="None")
  cat("\nNo VCF FILTER status restrictions detected.\n")
}else{
  filters.df <- suppressWarnings(formatFILTERs(filter.include,filter.exclude))
}

ev.include <- opts$ev.include
ev.exclude <- opts$ev.exclude
if(is.null(ev.include) & is.null(ev.exclude)){
  evs.df <- data.frame("includeEV"="None","excludeEV"="None")
  cat("\nNo per-sample genotype EV restrictions detected.\n")
}else{
  evs.df <- suppressWarnings(formatEVs(ev.include,ev.exclude))
}

###Builds table of conditions
#Enumerate all possible conditions
out.vect <- as.vector(unlist(apply(sizes.df,1,function(sizes){
  unlist(apply(freqs.df,1,function(freqs){
    unlist(apply(svtypes.df,1,function(svtypes){
      unlist(apply(filters.df,1,function(filters){
        unlist(apply(evs.df,1,function(evs){
          return(paste(c(sizes,freqs,svtypes,filters,evs),collapse="\t"))
        }))
      }))
    }))
  }))
})))
out.table <- t(as.data.frame(sapply(out.vect,strsplit,split="\t")))
rownames(out.table) <- 1:nrow(out.table)
colnames(out.table) <- unlist(lapply(list(sizes.df,freqs.df,svtypes.df,filters.df,evs.df),colnames))
out.table <- as.data.frame(out.table)
cat(paste("\nIdentified a total of ",prettyNum(nrow(out.table),big.mark=","),
          " distinct minGQ evaluation conditions.\n",sep=""))

#Exclude various impossible conditions
# - RD-only evidence for non-CNVs
# - RD-only evidence for CNVs with max size <= 5kb
# - SR-only evidence for CNVs with min size >= 5kb
# - PESR_GT_OVERDISPERSION or HIGH_SR_BACKGROUND for RD-only
cond.excl <- sort(unique(c(which(out.table$includeEV=="RD" & !(out.table$includeSVTYPE %in% c("DEL","DUP"))),
                      which(out.table$includeEV=="RD" & out.table$includeSVTYPE %in% c("DEL","DUP") & as.numeric(out.table$maxSVLEN)<=5000),
                      which(out.table$includeEV=="SR" & out.table$includeSVTYPE %in% c("DEL","DUP") & as.numeric(out.table$minSVLEN)>=5000),
                      which(out.table$includeEV=="RD" & out.table$includeFILTER %in% c("HIGH_SR_BACKGROUND","PESR_GT_OVERDISPERSION")))))
cat(paste("\nExcluded ",prettyNum(length(cond.excl),big.mark=","),"/",
          prettyNum(nrow(out.table),big.mark=","),
          " minGQ evaluation conditions due to impossible combinations of parameters (e.g. RD-only inversions).\n",sep=""))
out.table <- out.table[-cond.excl,]

#Write table out
out.table <- cbind(data.frame(paste("minGQ_condition",1:nrow(out.table),sep="_")),
                   out.table)
colnames(out.table)[1] <- "#condition"
write.table(out.table,OUTFILE,col.names=T,row.names=F,sep="\t",quote=F)
cat(paste("\nFinished. Wrote a total of ",prettyNum(nrow(out.table),big.mark=","),
          " minGQ evaluation conditions to the following file:\n  ",
          OUTFILE,"\n\n",sep=""))


