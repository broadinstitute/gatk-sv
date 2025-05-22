#!R 
library("optparse")

option_list = list(
        make_option(c('-i', "--input"), type="character", default=NULL, help="input file", metavar="character"),
        make_option(c('-o', "--output"), type="character", default=NULL, help="output file", metavar="character"),
        make_option(c('-p', "--path"), type="character", default=NULL, help="input file", metavar="character") );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


#integrate GC of each breakpoints:
	dat=read.table(opt$input,sep='\t')
	sr_le=read.table(paste(opt$path,'il_inte.le_bp.vs.SR', sep='/'))
	sd_le=read.table(paste(opt$path,'il_inte.le_bp.vs.SD', sep='/'))
	rm_le=read.table(paste(opt$path,'il_inte.le_bp.vs.RM', sep='/'))

	sr_ri=read.table(paste(opt$path,'il_inte.ri_bp.vs.SR', sep='/'))
	sd_ri=read.table(paste(opt$path,'il_inte.ri_bp.vs.SD', sep='/'))
	rm_ri=read.table(paste(opt$path,'il_inte.ri_bp.vs.RM', sep='/'))
	dat[,ncol(dat)+1] = 'US'
	dat[dat[,4]%in%rm_le[,4],][,ncol(dat)]='RM'
	dat[dat[,4]%in%rm_ri[,4],][,ncol(dat)]='RM'
	dat[dat[,4]%in%sd_le[,4],][,ncol(dat)]='SD'
	dat[dat[,4]%in%sd_ri[,4],][,ncol(dat)]='SD'
	dat[dat[,4]%in%sr_le[,4],][,ncol(dat)]='SR'
	dat[dat[,4]%in%sr_ri[,4],][,ncol(dat)]='SR'

	lg_cnv_sr=read.table(paste(opt$path,'il_inte.lg_cnv.vs.SR', sep='/'))
	lg_cnv_sd=read.table(paste(opt$path,'il_inte.lg_cnv.vs.SD', sep='/'))
	lg_cnv_rm=read.table(paste(opt$path,'il_inte.lg_cnv.vs.RM', sep='/'))
	lg_cnv = cbind(lg_cnv_sr[,c(4,ncol(lg_cnv_sr))], lg_cnv_sd[,ncol(lg_cnv_sd)], lg_cnv_rm[,ncol(lg_cnv_rm)])
	lg_cnv[,ncol(lg_cnv)+1] = 1-rowSums(lg_cnv[,c(2:4)])
	colnames(lg_cnv)=c('SVID','SR','SD','RM','US')
	lg_cnv[,ncol(lg_cnv)+1] = 'US'
	lg_cnv[lg_cnv$RM>.5,][,ncol(lg_cnv)]='RM'
	lg_cnv[lg_cnv$SD>.5,][,ncol(lg_cnv)]='SD'
	lg_cnv[lg_cnv$SR>.5,][,ncol(lg_cnv)]='SR'

	colnames(dat)[4]='SVID'
	dat=merge(dat, lg_cnv[,c(1,6)], by='SVID',all=T)
	dat[!is.na(dat[,ncol(dat)]),][,ncol(dat)-1] = dat[!is.na(dat[,ncol(dat)]),][,ncol(dat)] 
	colnames(dat)[ncol(dat)-1]='GC'
	write.table(dat[,c(1,ncol(dat)-1)],opt$output, quote=F, sep='\t', col.names=F, row.names=F)

