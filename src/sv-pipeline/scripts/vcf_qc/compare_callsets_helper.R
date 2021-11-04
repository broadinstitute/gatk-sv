#!/usr/bin/env Rscript

# Helper script to clean final table of vcf vs vcf comparisons

###Set parameters
options(stringsAsFactors=F)

###Helper function to clean an input overlap file
cleanData <- function(dat){
  #Get list of repeated VIDs
  repeat.IDs <- names(which(table(dat[,1])>1))
  #Iterate over repeated VIDs and resolve
  resolved.vals <- as.data.frame(t(sapply(repeat.IDs,function(VID){
    #Get values for repeated variant
    vals <- dat[which(dat[,1]==VID),2]
    #Test if any values remain after excluding NS & NO_OVR
    numeric.vals <- as.numeric(vals[which(!(vals %in% c("NO_OVR","NS")))])
    #If any exist, take mean as final value
    if(length(numeric.vals)>0){
      final <- mean(numeric.vals)
    }else{
      #Otherwise, report "NS" if any vals are NS
      if(any(vals) %in% "NS"){
        final <- "NS"
      }else{
        #Otherwise, report no overlap
        final <- "NO_OVR"
      }
    }
    #Return new row
    return(c(VID,final))
  })))
  #Make new data frame of single vals
  cleaned.dat <- rbind(dat[which(!(dat[,1] %in% repeat.IDs)),],
        resolved.vals)
  #Clean & return cleaned data
  rownames(cleaned.dat) <- NULL
  colnames(cleaned.dat) <- c("VID","OVR")
  return(cleaned.dat)
}

###Read arguments
args <- commandArgs(trailingOnly=T)
#Checks for appropriate number of positional arguments
if(length(args) != 7){
  stop("Incorrect number of required positional arguments\n")
}
#Assigns args
set2.in <- args[1]
ovr1a.in <- args[2]
ovr1b.in <- args[3]
ovr2a.in <- args[4]
ovr2b.in <- args[5]
ovr3.in <- args[6]
OUTFILE <- args[7]

###Reads & cleans all input files
#Read original set2 bed file
set2 <- read.table(set2.in,header=T,comment.char="",sep="\t",check.names=F)
#Read & clean overlap data
ovr.dat <- lapply(c(ovr1a.in,ovr1b.in,ovr2a.in,ovr2b.in,ovr3.in),function(path){
  #Read data
  dat <- read.table(path,header=F,sep="\t",check.names=F)
  #Clean data
  dat <- cleanData(dat)
  #Return data
  return(dat)
})
#Check to make sure all variants represented once and only once in each overlap dataset
match.checks <- unlist(lapply(ovr.dat,function(df){
  matches <- length(which(names(table(df[,1])==1) %in% as.character(set2[,4])))
  if(matches!=nrow(set2)){
    return(1)
  }else{
    return(0)
  }
}))
if(any(match.checks==1)){
  stop("ERROR: some benchmarking SV do not appear in overlap files")
}


###Formats final table & writes out
#Iteratively merge files
merged <- set2
for(i in 1:length(ovr.dat)){
  merged <- suppressWarnings(merge(merged,ovr.dat[[i]],by="VID"))
}
#Reformat merged files
merged <- merged[,c(2:4,1,5:ncol(merged))]
colnames(merged)[-(1:7)] <- c("ovr1a","ovr1b","ovr2a","ovr2b","ovr3")
colnames(merged)[1] <- "chr"
#Check to make sure no rows were gained or lost during processing
if(nrow(merged)!=nrow(set2)){
  stop("ERROR: final table does not match input benchmarking SV set")
}
#Coordinate-based sorting
merged <- merged[with(merged,order(chr,start,end,VID)),]
colnames(merged)[1] <- "#chr"
#Write out merged table
write.table(merged,OUTFILE,col.names=T,row.names=F,quote=F,sep="\t")

