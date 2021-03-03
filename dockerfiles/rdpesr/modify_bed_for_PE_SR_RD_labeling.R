#!/usr/bin/env Rscript
library("optparse")

option_list = list(
        make_option(c("-i", "--input"), type="character", default=NULL,
              help="name of input bed file", metavar="character"),
        make_option(c("--le_bp"), type="character", default=NULL,
              help="name of output left bp file", metavar="character"),
        make_option(c("--ri_bp"), type="character", default=NULL,
              help="name of output right bp file", metavar="character"),
        make_option(c("--le_flank"), type="character", default=NULL,
              help="name of output left flank file", metavar="character"),
        make_option(c("--ri_flank"), type="character", default=NULL,
              help="name of output right flank file", metavar="character")

 );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

filename=opt$input

flank_length=1000
dat=read.table(filename,sep='\t')
dat2=dat[,c(1,2,2,4)]
dat3=dat[,c(1,3,3,4)]
dat4=dat2
dat4[,2]=dat4[,2]-flank_length
dat5=dat3
dat5[,3]=dat5[,3]+flank_length
write.table(dat2, opt$le_bp, quote=F, sep='\t', col.names=F, row.names=F)
write.table(dat3, opt$ri_bp, quote=F, sep='\t', col.names=F, row.names=F)
write.table(dat4, opt$le_flank, quote=F, sep='\t', col.names=F, row.names=F)
write.table(dat5, opt$ri_flank, quote=F, sep='\t', col.names=F, row.names=F)

