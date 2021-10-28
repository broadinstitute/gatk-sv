#! /usr/bin/env Rscript

# Wrapper for WGD R companion library (WGDR) function
# WGD.sample.summaryplots.R to generate WGD summary plots
# for a set of libraries in a given coverage matrix

#installs optparse if not already installed
if("optparse" %in% rownames(installed.packages()) == FALSE)
{install.packages("optparse",repos="http://cran.rstudio.com")}
suppressWarnings(suppressPackageStartupMessages(library(optparse)))

#installs WGDR if not already installed
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
if("WGDR" %in% rownames(installed.packages()) == FALSE)
{install.packages(paste(script.basename,"/WGDR.tar.gz",sep=""),repos=NULL,type="source")}
suppressWarnings(suppressPackageStartupMessages(library(WGDR)))

#Disables factor default
options(stringsAsFactors=F)

#list of Rscript options
option_list <- list(
  make_option(c("-s", "--samples"), type="character", default=NULL,
              help="list of samples to evaluate [default: all samples]",
              metavar="character"),
  make_option(c("-t","--tranche"), action="store_false", default=TRUE,
              help="exclude most extreme 0.5% of bins genome-wide [default %default]"),
  make_option(c("-n","--normalize"), action="store_true", default=FALSE,
              help="normalize raw coverage values; only necessary when passing raw binCov matrix of counts [default %default]"),
  make_option(c("-q","--quiet"), action="store_true", default=FALSE,
              help="disable verbose output [default %default]")
)

#Get command-line arguments & options
parser <- OptionParser(usage="%prog [options] matrix outdir",
                       option_list=option_list,add_help_option=T)
args <- parse_args(parser,positional_arguments=TRUE)
opts <- args$options

#checks for appropriate positional arguments
if(length(args$args) != 2) {
  print_help(parser)
  stop("Incorrect number of required positional arguments")
}

#Prints startup message
if(opts$quiet==F){cat(paste("+------------------+\n",
                            "|       WGDR       |\n",
                            "| Sample Summaries |\n",
                            "|     (c) 2016     |\n",
                            "+------------------+\n",
                            sep=""))}

#Processes arguments & cleans options
if(!(file.exists(args$args[1]))){
  cat(paste("WGDR::ERROR [",
            strsplit(as.character(Sys.time()),split=" ")[[1]][2],
            "]: specified input matrix does not exist (",args$args[1],")\n",sep=""))

  stop()
}
if(!(file.exists(args$args[2]))){
  cat(paste("WGDR::STATUS [",
            strsplit(as.character(Sys.time()),split=" ")[[1]][2],
            "]: specified output directory does not exist; attempting to create\n",sep=""))
  dir.create(args$args[2])
}
if(is.null(opts$samples)==T){
  if(opts$quiet==F){
    cat(paste("WGDR::STATUS [",
              strsplit(as.character(Sys.time()),split=" ")[[1]][2],
              "]: processing all sample IDs in ",args$args[1],"\n",sep=""))
  }
  samples <- NULL
}else{
  if(opts$quiet==F){
    cat(paste("WGDR::STATUS [",
              strsplit(as.character(Sys.time()),split=" ")[[1]][2],
              "]: reading sample IDs from ",opts$samples,"\n",sep=""))
  }
  samples <- read.table(opts$samples)[,1]
}

#Reads matrix
mat <- WGD.matrix.read(args$args[1],
                       norm=opts$normalize)

#Tranches matrix (if optioned)
if(opts$tranche==T){
  mat <- WGD.matrix.tranche(mat)
}

#Assigns vector of sample IDs
if(is.null(samples)){
  samples <- names(mat$mat)[-c(1:3)]
}

#Plots sample summary per library
for(ID in samples){
  cat(paste("WGDR::STATUS [",
            strsplit(as.character(Sys.time()),split=" ")[[1]][2],
            "]: plotting summary for sample ",ID,"\n",sep=""))
  WGD.sample.plot(mat,ID,args$args[2])
}
