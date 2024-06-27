#!R
library("optparse")

option_list = list(
        make_option(c("--sv_vs_DHS_mamm"), type="character", default=NULL, help="permutation_round", metavar="character"),
        make_option(c("--sv_vs_DHS_prim"), type="character", default=NULL, help="permutation_round", metavar="character"),
        make_option(c("--sv_vs_footprint_mamm"), type="character", default=NULL, help="permutation_round", metavar="character"),
        make_option(c("--sv_vs_footprint_prim"), type="character", default=NULL, help="permutation_round", metavar="character"),
        make_option(c("--sv_vs_UCE"), type="character", default=NULL, help="permutation_round", metavar="character"),
        make_option(c("--sv_vs_phastCons100way"), type="character", default=NULL, help="permutation_round", metavar="character"),
        make_option(c("--sv_vs_phyloP100way"), type="character", default=NULL, help="permutation_round", metavar="character"),
        make_option(c("--sv_vs_z_over_2"), type="character", default=NULL, help="permutation_round", metavar="character"),
        make_option(c("--sv_vs_z_over_4"), type="character", default=NULL, help="permutation_round", metavar="character"),
        make_option(c("--output"), type="character", default=NULL, help="permutation_round", metavar="character")
 );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


d1=read.table(opt$sv_vs_DHS_mamm)
d1=d1[,c(4,ncol(d1))]
colnames(d1) = c('name','X239prim_DHS_mamm')

d2=read.table(opt$sv_vs_DHS_prim)
d2=d2[,c(4,ncol(d2))]
colnames(d2) = c('name','X239prim_DHS_prim')

d3=read.table(opt$sv_vs_footprint_mamm)
d3=d3[,c(4,ncol(d3))]
colnames(d3) = c('name','X239prim_footprint_mamm')

d4=read.table(opt$sv_vs_footprint_prim)
d4=d4[,c(4,ncol(d4))]
colnames(d4) = c('name','X239prim_footprint_prim')

d5=read.table(opt$sv_vs_UCE)
d5=d5[,c(4,ncol(d5))]
colnames(d5) = c('name','X239prim_uce')

d6=read.table(opt$sv_vs_phastCons100way)
d6=d6[,c(4,ncol(d6))]
colnames(d6) = c('name','phastCons100way')

d7=read.table(opt$sv_vs_phyloP100way)
d7=d7[,c(4,ncol(d7))]
colnames(d7) = c('name','phyloP100way')

dat=merge(d1,  d2, by='name')
dat=merge(dat, d3, by='name')
dat=merge(dat, d4, by='name')
dat=merge(dat, d5, by='name')
dat=merge(dat, d6, by='name')
dat=merge(dat, d7, by='name')

d8=read.table(opt$sv_vs_z_over_2)
d9=read.table(opt$sv_vs_z_over_4)

dat[,ncol(dat)+1] = 0
colnames(dat)[ncol(dat)] = 'constraint_z'
dat[dat$name%in%d8[d8[,ncol(d8)]>0, 4],]$constraint_z = 'z_over_2'
dat[dat$name%in%d9[d9[,ncol(d9)]>0, 4],]$constraint_z = 'z_over_4'
write.table(dat, opt$output, quote=F, sep='\t', col.names=T, row.names=F)

