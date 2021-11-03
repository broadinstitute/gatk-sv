#!/usr/bin/env Rscript

# Helper script to write list of variant IDs, genotypes, and GQs per sample

###Set parameters
options(stringsAsFactors=F)

###Read arguments
args <- commandArgs(trailingOnly=T)
#Checks for appropriate number of positional arguments
if(length(args) != 3){
  stop("Incorrect number of required positional arguments\n")
}
#Assigns args
vcf.in <- args[1]
samples.in <- args[2]
OUTDIR <- args[3]

###Ensure OUTDIR is created
if(!dir.exists(OUTDIR)){
  dir.create(OUTDIR)
}

###Read data
if(vcf.in=="/dev/stdin"){
  vcf <- read.table(file("stdin"),comment.char="",header=T,sep="\t",check.names=F)
}else{
  vcf <- read.table(vcf.in,comment.char="",header=T,sep="\t",check.names=F)
}
samples <- read.table(samples.in,header=F,check.names=F)[,1]

###Iterate over samples, parse genotypes, and write genotypes, GQs, and VIDs to file
sapply(samples,function(ID){
  #Get column index
  idx <- which(colnames(vcf) %in% ID)
  #Only run if sample ID present in VCF header
  if(!is.null(idx)){
    #Parse genotypes & GQs for SVs other than MCNV
	tmp = vcf[vcf[,5]!='<CNV>',c(3,idx)]
	tmp[,3]=apply(tmp,1,function(x){strsplit(as.character(x[2]),':')[[1]][1]})
	tmp[,4]=apply(tmp,1,function(x){strsplit(as.character(x[2]),':')[[1]][2]})
	keep_1 = tmp[!tmp[,3] %in% c('0/0'),c(1,3,4)]
	colnames(keep_1)=c('SVID','GT','GQ')

	keep_2 <- NULL
	if (any(vcf[,5]=='<CNV>')) {
        	tmp_mcnv = vcf[vcf[,5]=='<CNV>',c(3,9,idx)]
        	tmp_mcnv[,4]=apply(tmp_mcnv,1,function(x){match('CN',strsplit(as.character(x[2]),':')[[1]])})
        	tmp_mcnv[,5]=apply(tmp_mcnv,1,function(x){strsplit(as.character(x[3]),':')[[1]][as.integer(x[4])]})
        	tmp_mcnv[,6]='0/1'
        	tmp_mcnv[,7]=apply(tmp_mcnv,1,function(x){match('CNQ',strsplit(as.character(x[2]),':')[[1]])})
        	tmp_mcnv[,8]=apply(tmp_mcnv,1,function(x){strsplit(as.character(x[3]),':')[[1]][as.integer(x[4])]})
        	keep_2 = tmp_mcnv[tmp_mcnv[,5]!=2,c(1,6,8)]
		colnames(keep_2)=c('SVID','GT','GQ')
	}
    write.table(rbind(keep_1, keep_2),
                paste(OUTDIR,"/",ID,".VIDs_genotypes.txt",sep=""),
                col.names=F,row.names=F,sep="\t",quote=F)
    }
})

