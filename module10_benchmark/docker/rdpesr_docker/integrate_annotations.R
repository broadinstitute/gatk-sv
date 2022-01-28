#!/usr/bin/env Rscript
library("optparse")

option_list = list(
        make_option(c( "--bed"), type="character", default=NULL, help="bed file with basic information", metavar="character"),
        make_option(c( "--gc_anno"), type="character", default=NULL, help="bed file annotated with genomic content", metavar="character"),
        make_option(c( "--duphold_il"), type="character", default=NULL, help="vcf file annotated with duphold", metavar="character"),
        make_option(c( "--duphold_il_le"), type="character", default=NULL, help="left flank of SVs in vcf annotated with duphold", metavar="character"),
        make_option(c( "--duphold_il_ri"), type="character", default=NULL, help="right flank of SVs in vcf file annotated with duphold", metavar="character"),
        make_option(c( "--rd"), type="character", default=NULL, help="bed file annotated with bincov rd", metavar="character"),
        make_option(c( "--rd_le"), type="character", default=NULL, help="left flank of SVs in bed file annotated with bincov rd", metavar="character"),
        make_option(c( "--rd_ri"), type="character", default=NULL, help="right flank of SVs in bed file annotated with bincov rd", metavar="character"),
        make_option(c( "--pesr"), type="character", default=NULL, help="SVs in bed file annotated with pe sr counts", metavar="character"),
        make_option(c( "--info"), type="character", default=NULL, help="SVID with annotations such as SVTYPE SVLEN ALGORITHMS EVIDENCE FILTER", metavar="character"),
        make_option(c( "--gt"), type="character", default=NULL, help="SVID with annotations such as GT and GQ", metavar="character"),
        make_option(c( "--raw_manta"), type="character", default=NULL, help="comparison results of SV vs. raw manta SVs", metavar="character"),
        make_option(c( "--raw_wham"), type="character", default=NULL, help="comparison results of SV vs. raw wham SVs", metavar="character"),
        make_option(c( "--raw_melt"), type="character", default=NULL, help="comparison results of SV vs. raw melt SVs", metavar="character"),
        make_option(c( "--raw_depth"), type="character", default=NULL, help="comparison results of SV vs. raw depth SVs", metavar="character"),
        make_option(c( "--vs_hgsv"), type="character", default=NULL, help="comparison results of SV vs. hgsv SVs", metavar="character"),
        make_option(c( "--vs_pacbio"), type="character", default=NULL, help="comparison results of SV vs. pacbio SVs", metavar="character"),
        make_option(c( "--vs_bionano"), type="character", default=NULL, help="comparison results of SV vs. bionano SVs", metavar="character"),
        make_option(c( "--vs_array"), type="character", default=NULL, help="comparison results of SV vs. array SVs", metavar="character"),
        make_option(c( "--vapor"), type="character", default=NULL, help="vapor evaluation results", metavar="character"),
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

dat=read.table(opt$bed, sep='\t')
colnames(dat)=c('#CHR','POS','END', 'SVID','SVTYPE','svtype','length')

#categorize SV sizes:
dat[,ncol(dat)+1] = 's4_over5Kb'
colnames(dat)[ncol(dat)]='size_cate'
if(nrow(dat[dat$length<5000,])>0){
    dat[dat$length<5000,][,ncol(dat)]='s3_1to5Kb'
}
if(nrow(dat[dat$length<1000,])>0){
    dat[dat$length<1000,][,ncol(dat)]='s2_250bpto1Kb'
}
if(nrow(dat[dat$length<250,])){
    dat[dat$length<250,][,ncol(dat)]='s1_under250bp'
}

if(!is.null(opt$gc_anno)){
    gc_anno=read.table(opt$gc_anno,sep='\t')
    gc_anno=gc_anno[gc_anno[,6]!="",]
    colnames(gc_anno)=c('#CHR','POS','END','SVID','SVTYPE','sample','svtype','length','GC')
    dat = merge(dat, gc_anno[,c('SVID','GC')], by='SVID', all=T)
}

if(!is.null(opt$duphold_il)){
    duphold = extract_duphold_info(opt$duphold_il)
    colnames(duphold)=c('SVID','DHFC_IL','DHFFC_IL','DHBFC_IL')
    dat=merge(dat, duphold, by='SVID')
}

if(!is.null(opt$duphold_il_le)){
    duphold_le=extract_duphold_info(opt$duphold_il_le)
    colnames(duphold_le)=c('SVID','DHFC_IL_le','DHFFC_IL_le','DHBFC_IL_le')
    dat=merge(dat, duphold_le, by='SVID')
}

if(!is.null(opt$duphold_il_ri)){
    duphold_ri=extract_duphold_info(opt$duphold_il_ri)
    colnames(duphold_ri)=c('SVID','DHFC_IL_ri','DHFFC_IL_ri','DHBFC_IL_ri')
    dat=merge(dat, duphold_ri, by='SVID')
}

if(!is.null(opt$rd)){
    rd=read.table(opt$rd)
    colnames(rd)[c(4,ncol(rd)-2,ncol(rd)-1,ncol(rd))]=c('SVID','rd_median','rd_mean','rd_std')
    dat=merge(dat, rd[,c(4,ncol(rd)-2,ncol(rd)-1,ncol(rd))], by='SVID')
}

if(!is.null(opt$rd_le)){
    rd_le=read.table(opt$rd_le)
    colnames(rd_le)[c(4,ncol(rd_le)-2,ncol(rd_le)-1,ncol(rd_le))]=c('SVID','rd_median_le','rd_mean_le','rd_std_le')
    dat=merge(dat, rd_le[,c(4,ncol(rd_le)-2,ncol(rd_le)-1,ncol(rd_le))], by='SVID')
}

if(!is.null(opt$rd_ri)){
    rd_ri=read.table(opt$rd_ri)
    colnames(rd_ri)[c(4,ncol(rd_ri)-2,ncol(rd_ri)-1,ncol(rd_ri))]=c('SVID','rd_median_ri','rd_mean_ri','rd_std_ri')
    dat=merge(dat, rd_ri[,c(4,ncol(rd_ri)-2,ncol(rd_ri)-1,ncol(rd_ri))], by='SVID')
}

if(!is.null(opt$pesr)){
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
}

if(!is.null(opt$info)){
    info=read.table(opt$info,sep='\t', header =T)
    colnames(info)[1]='SVID'
    info[,ncol(info)+1] = apply(info,1,function(x){grepl('BOTHSIDES_SUPPORT',as.character(x[ncol(info)]))})
    colnames(info)[ncol(info)]='BothSideSupp'
    dat=merge(dat, info[,c('SVID','ALGORITHMS','EVIDENCE','FILTER','BothSideSupp')], by='SVID')
}

if(!is.null(opt$gt)){
    gtgq=read.table(opt$gt, header =T)
    dat=merge(dat, gtgq, by='SVID')
}

if(!is.null(opt$denovo)){
    dnv = read.table(opt$denovo, header =T)
    dat=merge(dat, dnv, by='SVID')
}

if(!is.null(opt$raw_manta)){
    vs_manta=read.table(opt$raw_manta, comment.char="", header=T, sep='\t')
    colnames(vs_manta)[c(4,8,9)]=c('SVID','vs_raw_manta_ovr1a','vs_raw_manta_ovr1b')
    dat=merge(dat, vs_manta[,c(4,8,9)], by='SVID')
}

if(!is.null(opt$raw_wham)){
    vs_wham=read.table(opt$raw_wham, comment.char="", header=T, sep='\t')
    colnames(vs_wham)[c(4,8,9)]=c('SVID','vs_raw_wham_ovr1a','vs_raw_wham_ovr1b')
    dat=merge(dat, vs_wham[,c(4,8,9)], by='SVID')
}

if(!is.null(opt$raw_melt)){
    vs_melt=read.table(opt$raw_melt, comment.char="", header=T, sep='\t')
    colnames(vs_melt)[c(4,8,9)]=c('SVID','vs_raw_melt_ovr1a','vs_raw_melt_ovr1b')
    dat=merge(dat, vs_melt[,c(4,8,9)], by='SVID')
}

if(!is.null(opt$raw_depth)){
    vs_depth=read.table(opt$raw_depth, comment.char="", header=T, sep='\t')
    colnames(vs_depth)[c(4,8,9)]=c('SVID','vs_raw_depth_ovr1a','vs_raw_depth_ovr1b')
    dat=merge(dat, vs_depth[,c(4,8,9)], by='SVID')
}

if(!is.null(opt$vs_hgsv)){
    vs_hgsv = read.table(opt$vs_hgsv, comment.char="", header=T, sep='\t')
    colnames(vs_hgsv)[c(4,8,9)]=c('SVID','vs_hgsv_ovr1a','vs_hgsv_ovr1b')
    dat=merge(dat, vs_hgsv[,c(4,8,9)], by='SVID')
}

if(!is.null(opt$vs_pacbio)){
    vs_pacbio = read.table(opt$vs_pacbio, comment.char="", header=T, sep='\t')
    colnames(vs_pacbio)[c(4,8,9)]=c('SVID','vs_pacbio_ovr1a','vs_pacbio_ovr1b')
    dat=merge(dat, vs_pacbio[,c(4,8,9)], by='SVID')
}

if(!is.null(opt$vs_bionano)){
    vs_bionano = read.table(opt$vs_bionano, comment.char="", header=T, sep='\t')
    colnames(vs_bionano)[c(4,8,9)]=c('SVID','vs_bionano_ovr1a','vs_bionano_ovr1b')
    dat=merge(dat, vs_bionano[,c(4,8,9)], by='SVID')
}

if(!is.null(opt$vs_array)){
    vs_bionano = read.table(opt$vs_array, comment.char="", header=T, sep='\t')
    colnames(vs_bionano)[c(4,8,9)]=c('SVID','vs_array_ovr1a','vs_array_ovr1b')
    dat=merge(dat, vs_bionano[,c(4,8,9)], by='SVID')
}

if(!is.null(opt$vapor)){
    vapor = read.table(opt$vapor, comment.char="", header=T, sep='\t')
    dat=merge(dat, vapor[,c(5:ncol(vapor))], by='SVID', all=T)
}

write.table(dat[,c(2:4,1,5:ncol(dat))],opt$output, quote =F, sep='\t', col.names=T, row.names=F)








