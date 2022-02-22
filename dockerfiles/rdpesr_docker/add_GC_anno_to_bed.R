#!/usr/bin/env Rscript
library("optparse")

option_list = list(
        make_option(c("-b", "--bedfile"), type="character", default=NULL,
              help="name of input bed file", metavar="character"),
        make_option(c("-o", "--output"), type="character", default=NULL,
              help="name of output file", metavar="character"),

        make_option(c( "--left_vs_SR"), type="character", default=NULL,
              help="", metavar="character"),
        make_option(c( "--left_vs_SD"), type="character", default=NULL,
              help="", metavar="character"),
        make_option(c( "--left_vs_RM"), type="character", default=NULL,
              help="", metavar="character"),

        make_option(c( "--right_vs_SR"), type="character", default=NULL,
              help="", metavar="character"),
        make_option(c( "--right_vs_SD"), type="character", default=NULL,
              help="", metavar="character"),
        make_option(c( "--right_vs_RM"), type="character", default=NULL,
              help="", metavar="character")
 );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

filename=opt$input

dat=read.table(opt$bedfile,sep='\t')
sr_le=read.table(opt$left_vs_SR,sep='\t')
sd_le=read.table(opt$left_vs_SD,sep='\t')
rm_le=read.table(opt$left_vs_RM,sep='\t')

sr_ri=read.table(opt$right_vs_SR,sep='\t')
sd_ri=read.table(opt$right_vs_SD,sep='\t')
rm_ri=read.table(opt$right_vs_RM,sep='\t')

dat[,ncol(dat)+1] = 'US'
dat[dat[,4]%in%rm_le[,4],][,ncol(dat)]='RM'
dat[dat[,4]%in%rm_ri[,4] & !dat[,5]%in%c('INS','ALU','LINE1','MEI','SVA'),][,ncol(dat)]='RM'

dat[dat[,4]%in%sd_le[,4],][,ncol(dat)]='SD'
dat[dat[,4]%in%sd_ri[,4] & !dat[,5]%in%c('INS','ALU','LINE1','MEI','SVA'),][,ncol(dat)]='SD'

dat[dat[,4]%in%sr_le[,4],][,ncol(dat)]='SR'
dat[dat[,4]%in%sr_ri[,4] & !dat[,5]%in%c('INS','ALU','LINE1','MEI','SVA'),][,ncol(dat)]='SR'

write.table(dat,opt$output, quote=F, sep='\t', col.names=F, row.names=F)


