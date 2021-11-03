#!/usr/bin/env Rscript

# Script to run cn.MOPS on a coverage matrix as output by makeMatrix.sh

####Load optparse
require(optparse)

####List of command-line options
option_list <- list(
  make_option(c("-I", "--ID"), type="character", default="CNMOPS",
              help="sample group ID (used for output filenames) [default %default]",
              metavar="character")
)

####Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog [options] covMatrix.bed OUTDIR",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

####Checks for appropriate positional arguments
if(length(args$args) != 2){
  stop("Incorrect number of required positional arguments\n")
}

####Writes args & opts to vars
path.to.matrix <- args$args[1]
OUTDIR <- args$args[2]
ID <- opts$ID

####Loads remiaining packages
require(cn.mops)
require(rtracklayer)

####Helper function for coercing data frames to GRanges
dataframe2GRanges <- function (df, keepColumns = TRUE, ignoreStrand = TRUE){
  stopifnot(class(df) == "data.frame")
  stopifnot(all(c("Start","Stop") %in% names(df)))
  stopifnot(any(c("Chr","Seqnames") %in% names(df)))
  if("Seqnames" %in% names(df)){
    names(df)[names(df)=="Seqnames"] <- "Chr"
  }
  if(!ignoreStrand && "Strand" %in% names(df)){
    if (is.numeric(df$Strand)) {
      Strand <- ifelse(df$Strand==1,"+","*")
      Strand[df$Strand==-1] <- "-"
      df$Strand <- Strand
    }
    gr <- GRanges(seqnames=df$Chr,
                  ranges=IRanges(start=df$Start,
                                 end=df$Stop),
                  strand=df$Strand)
  }else{
    gr <- GRanges(seqnames=df$Chr,
                  ranges=IRanges(start=df$Start,
                                 end=df$Stop))
  }
  if (keepColumns){
    dt <- as(df[,setdiff(names(df),c("Chr","Start","Stop","Strand"))],"DataFrame")
    elementMetadata(gr) <- dt
  }
  names(gr) <- rownames(df)
  return(gr)
}

####Loads coverage matrix & coerces to GRange
cov <- read.table(path.to.matrix,header=T,sep="\t",stringsAsFactors=F,comment.char="",check.names=F)
cov[2:ncol(cov)] <- apply(cov[2:ncol(cov)],2,function(vals){
  return(as.numeric(as.character(vals)))
})
colnames(cov)[1:3] <- c("Chr","Start","Stop")
cov[is.na(cov)] <- 0
#Exclude samples with median coverage = 0
median.per.sample <- apply(cov[,-c(1:3)],2,median,na.rm=T)
if(any(median.per.sample==0)){
  samples.to.exclude <- which(median.per.sample==0)
  if(length(samples.to.exclude)>0){
    samples.to.exclude <- samples.to.exclude+3
    cat(paste("\nEXCLUDING SAMPLES DUE TO NO COVERAGE:\n",
              paste(colnames(cov)[samples.to.exclude],collapse="\n"),sep=""))
    cov <- cov[,-samples.to.exclude]
    #Check to make sure at least three samples don't have zero coverage
    if(ncol(cov)<6){
      stop("runcnMOPS.R: Less than three samples retained in coverage matrix after filtering on median coverage")
    }
  }
}
#Coerce to GRange
cov <- dataframe2GRanges(cov)

#Runs cn.mops
res <- cn.mops(cov)
res <- calcIntegerCopyNumbers(res)

#Exports CNVs
export(cnvs(res),paste(OUTDIR,"/",ID,".cnMOPS.gff",sep=""),"GFF3")
