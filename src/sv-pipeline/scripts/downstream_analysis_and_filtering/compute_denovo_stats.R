#!/usr/bin/env Rscript

# Helper script to titrate proband variants & de novo rate given a minimum GQ


###Set master parameters
options(stringsAsFactors=F,scipen=1000)
svtypes <- c("DEL","DUP","INS","INV","CPX","BND")


###########################
###GENERAL HELPER FUNCTIONS
###########################
#Read & clean trio data
#Restrict to heterozygous calls in proband only for this analysis
readTrioDat <- function(trio_dat, famID){
  dat <- read.table(trio_dat,header=T,comment.char="")
  dat[,(ncol(dat)-5):ncol(dat)] <- apply(dat[,(ncol(dat)-5):ncol(dat)],2,as.integer)
  dat <- dat[which(dat$pro_AC==1),]
  colnames(dat)[1] <- "VID"
  dat$dn <- dat$pro_AC - (dat$fa_AC + dat$mo_AC)
  return(dat)
}
#Get counts of all proband variants & de novo proband variants given min GQ
filterGQStats <- function(dat, minGQ=0, svtype=NULL){
  if(!is.null(svtype)){
    max.pro <- nrow(dat[which(dat$svtype==svtype),])
  }else{
    max.pro <- nrow(dat)
  }
  sub.pro <- dat[which(dat$pro_GQ>=minGQ),]
  sub.trio <- dat[which(dat$pro_GQ>=minGQ & dat$fa_GQ>=minGQ & dat$mo_GQ>=minGQ),]
  if(!is.null(svtype)){
    sub.pro <- sub.pro[which(sub.pro$svtype==svtype),]
    sub.trio <- sub.trio[which(sub.trio$svtype==svtype),]
  }
  pro.all <- nrow(sub.pro)
  pro.dn <- nrow(sub.trio[which(sub.trio$dn>0),])
  if(max.pro>0){
    pro.all.ret <- pro.all/max.pro
  }else{
    pro.all.ret <- NA
  }
  if(pro.all>0){
    pro.dnr <- pro.dn/pro.all
  }else{
    pro.dnr <- NA
  }
  return(c(pro.all.ret,pro.dnr))
}
#Master wrapper to iterate over all SV types across a range of GQs and return a table of counts
getAllGQStats <- function(dat, rangeGQ, svtypes){
  iter.res <- as.data.frame(t(sapply(rangeGQ,function(minGQ){
    ALL <- filterGQStats(dat,minGQ)
    bySVTYPE <- as.vector(sapply(svtypes,function(svtype){
      filterGQStats(dat,minGQ,svtype)
    }))
    out.v <- c(ALL,bySVTYPE)
    return(out.v)
  })))
  colnames(iter.res) <- as.vector(sapply(c("ALL",svtypes),function(svtype){return(paste(svtype,c(".het",".dn"),sep=""))}))
  iter.res <- cbind("minGQ"=rangeGQ,iter.res)
  return(iter.res)
}


################
###RSCRIPT BLOCK
################
require(optparse)
###List of command-line options
option_list <- list(
  make_option(c("--famID"), type="character", default="trio",
              help="name of trio used for prefixing output [default %default]",
              metavar="character"),
  make_option(c("--minGQ"), type="integer", default=0,
              help="minimum GQ in window to test [default %default]",
              metavar="integer"),
  make_option(c("--maxGQ"), type="integer", default=200,
              help="maximum GQ in window to test [default %default]",
              metavar="integer"),
  make_option(c("--step"), type="integer", default=1,
              help="GQ increment step size [default %default]",
              metavar="integer")
)

###Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog trio_dat outfile",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

###Checks for appropriate positional arguments
if(length(args$args) != 2){
  stop("Incorrect number of required positional arguments\n")
}

###Writes args & opts to vars
trio_dat <- args$args[1]
outfile <- args$args[2]
famID <- opts$famID
minGQ <- opts$minGQ
maxGQ <- opts$maxGQ
step <- opts$step

###Read & clean data
dat <- readTrioDat(trio_dat, famID)

###Get master GQ table
GQrange <- seq(minGQ,maxGQ,by=step)
res <- getAllGQStats(dat,GQrange,svtypes)

###Write table to outfile
write.table(res,outfile,col.names=T,row.names=F,sep="\t",quote=F)
