#!/usr/bin/env Rscript


# evenSplitter: helper script for splitting tab-delimited files


####################################
#####Set parameters & load libraries
####################################
options(scipen=1000,stringsAsFactors=F)
require(optparse)


##########################################
#####Read command-line arguments & options
##########################################
#List of command-line options
option_list <- list(
  make_option(c("-L", "--targetLines"),type="integer",default=NULL,
              help="target lines per split [default: %default]",
              metavar="integer"),
  make_option(c("-S", "--targetSplits"),type="integer",default=NULL,
              help="target splits [default: %default]",
              metavar="integer"),
  make_option(c("--shuffle"),action="store_true",default=FALSE,
              help="randomly shuffle lines in input file before splitting [default: %default]"),
  make_option(c("-q","--quiet"),action="store_true",default=FALSE,
              help="suppress output messages [default: %default]")
)

#Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog [options] INFILE PREFIX",
                                option_list=option_list),
                   positional_arguments=TRUE)

#Sanity-check arguments
if(length(args$args) != 2){
  stop("Must supply an input file and output prefix\n")
}

#Clean arguments
INFILE <- args$args[1]
PREFIX <- args$args[2]

#Clean options
shuf <- args$options$shuffle
quiet <- args$options$quiet
target.lines <- args$options$targetLines
target.splits <- args$options$targetSplits


#######################
#####Sanity check input
#######################
# n.modes.set <- length(which(!sapply(c(min.lines,max.lines,target.lines,
#                                        min.splits,max.splits,target.splits),
#                                     is.null)))
n.modes.set <- length(which(!sapply(c(target.lines,target.splits),is.null)))

if(n.modes.set < 1 ){
  # stop("Must specify a splitting mode (min/max/target lines or splits)")
  stop("Must specify a splitting mode (lines or splits)\n")
}
if(n.modes.set > 1 ){
  # stop("Cannot specify more than one splitting mode (min/max/target lines or splits)")
  stop("Cannot specify more than one splitting mode (lines or splits)\n")
}


#############################
#####Read & format input file
#############################
#Read file as data frame
dat <- as.data.frame(read.table(INFILE,sep="\t",header=F))

#Count number of lines
nlines <- nrow(dat)

#Randomly shuffle lines (if optioned)
if(shuf == T){
  if(nlines > 0){
    dat <- as.data.frame(dat[sample(1:nlines,size=nlines,replace=F),])
  }
}


##############################################################
#####Determine number of splits, and number of lines per split
##############################################################
#Perform different operations based on input mode
if(!is.null(target.lines)){
  #Determine optimum number of splits
  target.splits <- max(c(1,round(nlines/target.lines)))
}

#Determine minimum number of lines per file
min.lines <- floor(nlines/target.splits)

#Create split data frame
splits.df <- data.frame("split.no"=1:target.splits,
                        "lines"=min.lines)

#Determine excess number of lines to be distributed
excess.lines <- nlines-sum(splits.df$lines)

#Increment random splits by 1
excess.sinks <- sample(1:nrow(splits.df),size=excess.lines,replace=F)
splits.df[excess.sinks,2] <- splits.df[excess.sinks,2]+1


##############################
#####Report plan to split file
##############################
if(quiet == F){
  schema <- table(splits.df[,2])
  cat("SPLITTING SCHEMA:\n")
  for(i in 1:length(schema)){
    cat(paste(schema[i]," splits @ ",names(schema)[i]," lines\n",sep=""))
  }
}


###############
#####Split file
###############
#Set start and stop lines for splits
splits.df$start <- NA
splits.df$end <- cumsum(splits.df$lines)
if(nrow(splits.df)==1){
  splits.df$start <- 1
}else{
  splits.df$start <- c(1,splits.df$end[1:(nrow(splits.df)-1)]+1)
}


#Iterate over splits and write to file
apply(splits.df,1,function(vals){
  vals <- as.integer(vals)
  chunk <- as.data.frame(dat[vals[3]:vals[4],])
  write.table(chunk,paste(PREFIX,vals[1],sep=""),sep="\t",quote=F,
              row.names=F,col.names=F)
})
