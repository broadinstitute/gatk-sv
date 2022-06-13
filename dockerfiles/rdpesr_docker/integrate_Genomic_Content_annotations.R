#!R 
#integrate GC of each breakpoints:
add_genomic_content_anno<-function(dat, le_bp_vs_sr, le_bp_vs_sd, le_bp_vs_rm, ri_bp_vs_sr, ri_bp_vs_sd, ri_bp_vs_rm, lg_cnv_vs_sr, lg_cnv_vs_sd, lg_cnv_vs_rm){
        dat=read.table(dat)
        sr_le=read.table(le_bp_vs_sr)
        sd_le=read.table(le_bp_vs_sd)
        rm_le=read.table(le_bp_vs_rm)

        sr_ri=read.table(ri_bp_vs_sr)
        sd_ri=read.table(ri_bp_vs_sd)
        rm_ri=read.table(ri_bp_vs_rm)

          
        dat[,ncol(dat)+1] = 'US'
        dat[dat[,4]%in%rm_le[,4],][,ncol(dat)]='RM'
        dat[dat[,4]%in%rm_ri[,4],][,ncol(dat)]='RM'
        dat[dat[,4]%in%sd_le[,4],][,ncol(dat)]='SD'
        dat[dat[,4]%in%sd_ri[,4],][,ncol(dat)]='SD'
        dat[dat[,4]%in%sr_le[,4],][,ncol(dat)]='SR'
        dat[dat[,4]%in%sr_ri[,4],][,ncol(dat)]='SR'

        lg_cnv_sr=read.table(lg_cnv_vs_sr)
        lg_cnv_sd=read.table(lg_cnv_vs_sd)
        lg_cnv_rm=read.table(lg_cnv_vs_rm)

        lg_cnv = cbind(lg_cnv_sr[,c(4,ncol(lg_cnv_sr))], lg_cnv_sd[,ncol(lg_cnv_sd)], lg_cnv_rm[,ncol(lg_cnv_rm)])
        lg_cnv[,ncol(lg_cnv)+1] = 1-rowSums(lg_cnv[,c(2:4)])
        colnames(lg_cnv)=c('SVID','SR','SD','RM','US')
        lg_cnv[,ncol(lg_cnv)+1] = 'US'
        lg_cnv[lg_cnv$RM > .5,][,ncol(lg_cnv)]='RM'
        lg_cnv[lg_cnv$SD > .5,][,ncol(lg_cnv)]='SD'
        lg_cnv[lg_cnv$SR > .5,][,ncol(lg_cnv)]='SR'

        colnames(dat)[4]='SVID'
        dat=merge(dat, lg_cnv[,c(1,6)], by='SVID',all=T)
        dat[!is.na(dat[,ncol(dat)]),][,ncol(dat)-1] = dat[!is.na(dat[,ncol(dat)]),][,ncol(dat)] 
        colnames(dat)[ncol(dat)-1]='GC'

        return(dat[,c("SVID","GC")])
}


#!/usr/bin/env Rscript
library("optparse")

option_list = list(
        make_option(c( "--bed"), type="character", default=NULL, help="bed file", metavar="character"),
        make_option(c( "--out"), type="character", default=NULL, help="name of output file", metavar="character"),
        make_option(c( "--le_bp_vs_sr"), type="character", default=NULL, help="overlap between left breakpoint and simple repeats", metavar="character"),
        make_option(c( "--le_bp_vs_sd"), type="character", default=NULL, help="overlap between left breakpoint and segmetnal duplicates", metavar="character"),
        make_option(c( "--le_bp_vs_rm"), type="character", default=NULL, help="overlap between left breakpoint and repeat masks", metavar="character"),
        make_option(c( "--ri_bp_vs_sr"), type="character", default=NULL, help="overlap between right breakpoint and simple repeats", metavar="character"),
        make_option(c( "--ri_bp_vs_sd"), type="character", default=NULL, help="overlap between right breakpoint and segmetnal duplicates", metavar="character"),
        make_option(c( "--ri_bp_vs_rm"), type="character", default=NULL, help="overlap between right breakpoint and repeat masks", metavar="character"),
        make_option(c( "--lg_cnv_vs_sr"), type="character", default=NULL, help="overlap between large CNVs and simple repeats", metavar="character"),
        make_option(c( "--lg_cnv_vs_sd"), type="character", default=NULL, help="overlap between large CNVs and segmetnal duplicates", metavar="character"),
        make_option(c( "--lg_cnv_vs_rm"), type="character", default=NULL, help="overlap between large CNVs and repeat masks", metavar="character")
 );


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

annotated_SV = add_genomic_content_anno(opt$bed, opt$le_bp_vs_sr, opt$le_bp_vs_sd, opt$le_bp_vs_rm, opt$ri_bp_vs_sr, opt$ri_bp_vs_sd, opt$ri_bp_vs_rm, opt$lg_cnv_vs_sr, opt$lg_cnv_vs_sd, opt$lg_cnv_vs_rm)

write.table(annotated_SV, opt$out, quote=F, sep='\t', col.names=T, row.names=F)





