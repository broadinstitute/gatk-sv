#!/usr/bin/env Rscript
library("optparse")

option_list = list(
        make_option(c("-i", "--input"), type="character", default=NULL,
              help="name of input bincov tsv file", metavar="character")
 );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

filein=opt$input
dat=read.table(filein)
median=median(dat[,4])
dat[,4]=dat[,4]/median
write.table(dat, gsub('.gz','',gsub('bincov','normCov',filein)), quote=F, sep='\t', col.names=F, row.names=F)


