#!/usr/bin/env Rscript
library("optparse")

option_list = list(
        make_option(c( "--anno"), type="character", default=NULL, help="bed file annotated with all annotations", metavar="character"),
        make_option(c( "--rd"), type="character", default=NULL, help="bed file annotated with bincov rd", metavar="character"),
        make_option(c( "--rd_le"), type="character", default=NULL, help="left flank of SVs in bed file annotated with bincov rd", metavar="character"),
        make_option(c( "--rd_ri"), type="character", default=NULL, help="right flank of SVs in bed file annotated with bincov rd", metavar="character"),
        make_option(c( "--pesr"), type="character", default=NULL, help="SVs in bed file annotated with pe sr counts", metavar="character"),
        #make_option(c( "--info"), type="character", default=NULL, help="SVID with annotations such as SVTYPE SVLEN ALGORITHMS EVIDENCE FILTER", metavar="character"),
        make_option(c( "--gt"), type="character", default=NULL, help="SVID with annotations such as GT and GQ", metavar="character"),
        make_option(c( "--raw_manta"), type="character", default=NULL, help="comparison results of SV vs. raw manta SVs", metavar="character"),
        make_option(c( "--raw_wham"), type="character", default=NULL, help="comparison results of SV vs. raw wham SVs", metavar="character"),
        make_option(c( "--raw_melt"), type="character", default=NULL, help="comparison results of SV vs. raw melt SVs", metavar="character"),
        make_option(c( "--denovo"), type="character", default=NULL, help="two column file with SVID and de novo rate", metavar="character"),
        make_option(c( "--output"), type="character", default=NULL, help="output file", metavar="character")
 );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

extract_duphold_info<-function(duphold_name){
 duphold=read.table(duphold_name)
 out=duphold[,c(3,9,10)]
 for(i in c('DHFC','DHFFC','DHBFC')){
  out[,ncol(out)+1]=apply(out,1,function(x){strsplit(as.character(x[3]),':')[[1]][match(i,strsplit(as.character(x[2]),':')[[1]])]})
  colnames(out)[ncol(out)]=i
 }
 colnames(out)[1]='SVID'
 return(out[,c(1,4:ncol(out))])
}

anno=read.table(opt$anno, header=T)

rd=read.table(opt$rd)
colnames(rd)[c(1:4,ncol(rd)-2,ncol(rd)-1,ncol(rd))]=c('#chr','pos','end','SVID','rd_median','rd_mean','rd_std')
dat=merge(rd[,c(1:4,ncol(rd)-2,ncol(rd)-1,ncol(rd))],anno,  by='SVID')

rd_le=read.table(opt$rd_le)
colnames(rd_le)[c(4,ncol(rd_le)-2,ncol(rd_le)-1,ncol(rd_le))]=c('SVID','rd_median_le','rd_mean_le','rd_std_le')
dat=merge(dat, rd_le[,c(4,ncol(rd_le)-2,ncol(rd_le)-1,ncol(rd_le))], by='SVID')

rd_ri=read.table(opt$rd_ri)
colnames(rd_ri)[c(4,ncol(rd_ri)-2,ncol(rd_ri)-1,ncol(rd_ri))]=c('SVID','rd_median_ri','rd_mean_ri','rd_std_ri')
dat=merge(dat, rd_ri[,c(4,ncol(rd_ri)-2,ncol(rd_ri)-1,ncol(rd_ri))], by='SVID')

pesr=read.table(opt$pesr)
colnames(pesr)[c(4,ncol(pesr)-5,ncol(pesr)-4,ncol(pesr)-3,ncol(pesr)-2,ncol(pesr)-1,ncol(pesr))]=c('SVID','PE_le','PE_ri','SR_le','SR_ri','SR_le_V2','SR_ri_V2')
dat=merge(dat, pesr[,c(4,ncol(pesr)-5,ncol(pesr)-4,ncol(pesr)-3,ncol(pesr)-2,ncol(pesr)-1,ncol(pesr))], by='SVID')

dat[,ncol(dat)+1]=apply(dat[,c('PE_le','PE_ri')],1,max)
colnames(dat)[ncol(dat)]='PE_max'
dat[,ncol(dat)+1]=apply(dat[,c('PE_le','PE_ri')],1,min)
colnames(dat)[ncol(dat)]='PE_min'
dat[,ncol(dat)+1]=apply(dat[,c('SR_le','SR_ri')],1,max)
colnames(dat)[ncol(dat)]='SR_max'
dat[,ncol(dat)+1]=apply(dat[,c('SR_le','SR_ri')],1,min)
colnames(dat)[ncol(dat)]='SR_min'


gtgq=read.table(opt$gt, header =T)
dat=merge(dat, gtgq, by='SVID')

dnv = read.table(opt$denovo, header =T)
dnv[,7]=rowSums(dnv[,c(2:5)])
colnames(dnv)[7]='inheri_trios'
dat=merge(dat, dnv[,c('SVID','denovo_rate','inheri_trios')], by='SVID')


vs_manta=read.table(opt$raw_manta, comment.char="", header=T, sep='\t')
colnames(vs_manta)[c(4,8,9)]=c('SVID','vs_raw_manta_ovr1a','vs_raw_manta_ovr1b')
dat=merge(dat, vs_manta[,c(4,8,9)], by='SVID')
vs_wham=read.table(opt$raw_wham, comment.char="", header=T, sep='\t')
colnames(vs_wham)[c(4,8,9)]=c('SVID','vs_raw_wham_ovr1a','vs_raw_wham_ovr1b')
dat=merge(dat, vs_wham[,c(4,8,9)], by='SVID')
vs_melt=read.table(opt$raw_melt, comment.char="", header=T, sep='\t')
colnames(vs_melt)[c(4,8,9)]=c('SVID','vs_raw_melt_ovr1a','vs_raw_melt_ovr1b')
dat=merge(dat, vs_melt[,c(4,8,9)], by='SVID')

write.table(dat[,c(2:4,1,5:ncol(dat))],opt$output, quote =F, sep='\t', col.names=T, row.names=F)
