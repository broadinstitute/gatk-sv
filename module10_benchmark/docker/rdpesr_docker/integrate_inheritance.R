#!/usr/bin/env Rscript
library("optparse")

option_list = list(
        make_option(c( "--bed"), type="character", default=NULL, help="bed file with basic information", metavar="character"),
        make_option(c( "--pb_vs_fa"), type="character", default=NULL, help="comparison results between proband and father", metavar="character"),
        make_option(c( "--pb_vs_mo"), type="character", default=NULL, help="comparison results between proband and mother", metavar="character"),
        make_option(c( "--output"), type="character", default=NULL, help="output file", metavar="character")
 );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

dat=read.table(opt$bed, header=T, sep='\t', comment.char="")
colnames(dat)[1]='#CHROM'
pb_vs_fa = read.table(opt$pb_vs_fa, sep='\t')
pb_vs_mo = read.table(opt$pb_vs_mo, sep='\t')
dat[,ncol(dat)+1] = 'denovo'
dat[dat[,4]%in%pb_vs_fa[,4],][,ncol(dat)]='pb_fa'
dat[dat[,4]%in%pb_vs_mo[,4],][,ncol(dat)]='pb_mo'
dat[dat[,4]%in%pb_vs_fa[,4] & dat[,4]%in%pb_vs_mo[,4],][,ncol(dat)]='pb_fa_mo'
colnames(dat)[ncol(dat)] = 'inheritance'
write.table(dat, opt$output, quote=F, sep='\t', col.names=T, row.names=F)
