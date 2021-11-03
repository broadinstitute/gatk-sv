#!/usr/bin/env Rscript

# Helper script to filter trio-based SV data to meet
# a set of conditions for minGQ parameterization


###Set master parameters
options(stringsAsFactors=F,scipen=1000)


###################
###HELPER FUNCTIONS
###################
#Check & clean optional arguments
checkOpt <- function(value,numeric=F){
  if(!is.null(value) & value!="None"){
    if(numeric==T){
      cleaned <- as.numeric(value)
    }else{
      cleaned <- unlist(strsplit(value,split=","))
    }
  }else{
    cleaned <- NA
  }
  return(cleaned)
}
#Filter a numeric vector
filterNumeric <- function(vals,cutoff,greater=T){
  if(length(vals)>0){
    if(!is.na(cutoff)){
      if(greater==T){
        return(which(as.numeric(vals)>=as.numeric(cutoff)))
      }else{
        return(which(as.numeric(vals)<as.numeric(cutoff)))
      }
    }else{
      return(1:length(vals))
    }
  }else{
    return(NULL)
  }
}
#Filter a categorical vector
filterCategorical <- function(vals,criteria,include=T){
  if(length(vals)>1){
    if(any(sapply(criteria,is.na))){
      return(1:length(vals))
    }else{
      vals <- strsplit(vals,split=",")
      matches <- unlist(lapply(vals,function(val){
        if(any(unlist(val) %in% criteria)){
          return(TRUE)
        }else{
          return(FALSE)
        }
      }))
      if(include==T){
        return(which(matches))
      }else{
        return(which(!matches))
      }
    }
  }else{
    return(NULL)
  }
}
#Filter variant table
filterVariants <- function(dat,min.size,max.size,min.freq,max.freq,
                           svtype.include,svtype.exclude,
                           filter.include,filter.exclude,
                           ev.include,ev.exclude,max.variants){
  #Filter on SVLEN
  dat <- dat[filterNumeric(vals=dat$SVLEN,cutoff=min.size,greater=T),]
  dat <- dat[filterNumeric(vals=dat$SVLEN,cutoff=max.size,greater=F),]
  
  #Filter on AF
  dat <- dat[filterNumeric(vals=dat$AF,cutoff=min.freq,greater=T),]
  dat <- dat[filterNumeric(vals=dat$AF,cutoff=max.freq,greater=F),]
  
  #Filter on SVTYPE
  dat <- dat[filterCategorical(vals=dat$SVTYPE,criteria=svtype.include,include=T),]
  dat <- dat[filterCategorical(vals=dat$SVTYPE,criteria=svtype.exclude,include=F),]
  
  #Filter on FILTER
  dat <- dat[filterCategorical(vals=dat$FILTER,criteria=filter.include,include=T),]
  dat <- dat[filterCategorical(vals=dat$FILTER,criteria=filter.exclude,include=F),]
  
  #Filter on pro_EV
  dat <- dat[filterCategorical(vals=dat$pro_EV,criteria=ev.include,include=T),]
  dat <- dat[filterCategorical(vals=dat$pro_EV,criteria=ev.exclude,include=F),]
  
  #Downsample to max.variants if > max.variants
  if(nrow(dat) > max.variants){
    set.seed(123456789)
    dat <- dat[sample(1:nrow(dat),size=max.variants,replace=F),]
  }
  
  #Return filtered data
  return(dat)
}


################
###RSCRIPT BLOCK
################
require(optparse,quietly=T)
###List of command-line options
option_list <- list(
  make_option(c("--min.size"), type="character", default=NULL,
              help="minimum SVLEN [default %default]"),
  make_option(c("--max.size"), type="character", default=NULL,
              help="maximum SVLEN [default %default]"),
  make_option(c("--min.freq"), type="character", default=NULL,
              help="minimum AF [default %default]"),
  make_option(c("--max.freq"), type="character", default=NULL,
              help="maximum AF [default %default]"),
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
              help="Per-sample genotype EV to exclude [default %default]"),
  make_option(c("--max.variants"), type="integer", default=1000000,
              help="Maximum number of variants to retain per trio [default %default]")
)

###Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog INFILE OUTFILE",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options


###Checks for appropriate positional arguments
if(length(args$args) != 2){
  stop("Incorrect number of required positional arguments\n")
}

###Writes args & opts to vars
INFILE <- args$args[1]
OUTFILE <- args$args[2]
min.size <- checkOpt(opts$min.size,numeric=T)
if(!is.na(min.size)){
  if(min.size==0){
    min.size <- -2
  }
}
max.size <- checkOpt(opts$max.size,numeric=T)
min.freq <- checkOpt(opts$min.freq,numeric=T)
max.freq <- checkOpt(opts$max.freq,numeric=T)
svtype.include <- checkOpt(opts$svtype.include)
svtype.exclude <- checkOpt(opts$svtype.exclude)
filter.include <- checkOpt(opts$filter.include)
filter.exclude <- checkOpt(opts$filter.exclude)
ev.include <- checkOpt(opts$ev.include)
ev.exclude <- checkOpt(opts$ev.exclude)
max.variants <- as.numeric(opts$max.variants)

###Read input data
dat <- read.table(INFILE,header=T,comment.char="",check.names=F)
colnames(dat)[1] <- "famID"

###Filter input data
if(nrow(dat)>0){
  dat <- filterVariants(dat=dat,
                        min.size=min.size,
                        max.size=max.size,
                        min.freq=min.freq,
                        max.freq=max.freq,
                        svtype.include=svtype.include,
                        svtype.exclude=svtype.exclude,
                        filter.include=filter.include,
                        filter.exclude=filter.exclude,
                        ev.include=ev.include,
                        ev.exclude=ev.exclude,
                        max.variants=max.variants)
}
dat <- data.frame(dat$famID,dat$VID,dat$pro_AC,dat$fa_AC,dat$mo_AC,dat$pro_GQ,dat$fa_GQ,dat$mo_GQ)

###Return filtered data
write.table(dat,OUTFILE,col.names=F,row.names=F,quote=F,sep="\t")
