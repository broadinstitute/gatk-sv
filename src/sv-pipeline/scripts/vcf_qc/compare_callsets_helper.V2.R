#!/usr/bin/env Rscript

# Helper script to clean final table of vcf vs vcf comparisons

###Set parameters
options(stringsAsFactors=F)

###Helper function to clean an input overlap file
cleanData <- function(dat){
  uniq = dat[dat[,2]=='NO_OVR',]
  dat = dat[dat[,2]!='NO_OVR',]

  dat[,2]=as.numeric(as.character(dat[,2]))
  #Get list of repeated VIDs
  repeat.IDs <- names(which(table(dat[,1])>1))
  out = dat[which(!(dat[,1] %in% repeat.IDs)),]
  for(i in sort(unique(table(dat[,1])))){
  	if (i>1){
	  	print(i)
	  	tmp_IDs <- names(which(table(dat[,1])==i))
	  	if (length(tmp_IDs)>0){
		  	tmp_table <- dat[dat[,1]%in%tmp_IDs,]
		  	tmp_tmp <- tmp_table[c(1:length(tmp_IDs))*(i)-i+1,]
		  	for(j in c(2:i)){
		  		tmp_tmp = cbind(tmp_tmp, tmp_table[c(1:length(tmp_IDs))*(i)-i+j,2])
		  	}
		  	tmp_tmp[,1]=as.character(tmp_tmp[,1])
			tmp_table <-cbind(tmp_tmp[,1],rowMeans(tmp_tmp[,c(2:ncol(tmp_tmp))]))
		  	colnames(tmp_table)=c('V1','V2')
		  	out = rbind(out, tmp_table)
	  	}
  	}
  }

  uniq2 = uniq[!uniq[,1]%in%out[,1],]
  cleaned.dat=rbind(out, uniq2)
  #Iterate over repeated VIDs and resolve

  #Make new data frame of single vals
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
    test = set2[!set2[,4]%in%df[,1],]
    if(nrow(test)>0){
        print(test) 
        return(1)
  }
    if (nrow(test)==0){
    	return(0)
    }
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

