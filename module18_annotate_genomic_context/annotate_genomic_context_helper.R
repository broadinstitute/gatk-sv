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
	if(nrow(rm_le)>0){
			dat[dat[,4]%in%rm_le[,4],][,ncol(dat)]='RM'
	}
	if(nrow(rm_ri)>0){
		dat[dat[,4]%in%rm_ri[,4],][,ncol(dat)]='RM'
	}
	if(nrow(sd_le)>0){
		dat[dat[,4]%in%sd_le[,4],][,ncol(dat)]='SD'
	}
	if(nrow(sd_ri)>0){
		dat[dat[,4]%in%sd_ri[,4],][,ncol(dat)]='SD'
	}
	if(nrow(sr_le)>0){
		dat[dat[,4]%in%sr_le[,4],][,ncol(dat)]='SR'
	}
	if(nrow(sr_ri)>0){
		dat[dat[,4]%in%sr_ri[,4],][,ncol(dat)]='SR'
	}


	colnames(dat)[ncol(dat)] = 'GC'
        colnames(dat)[4]='SVID'
 
	if (file.exists(paste(opt$path,'il_inte.lg_cnv.vs.SR', sep='/'))) {
		lg_cnv_sr=read.table(paste(opt$path,'il_inte.lg_cnv.vs.SR', sep='/'))
		lg_cnv_sd=read.table(paste(opt$path,'il_inte.lg_cnv.vs.SD', sep='/'))
		lg_cnv_rm=read.table(paste(opt$path,'il_inte.lg_cnv.vs.RM', sep='/'))
		lg_cnv = cbind(lg_cnv_sr[,c(4,ncol(lg_cnv_sr))], lg_cnv_sd[,ncol(lg_cnv_sd)], lg_cnv_rm[,ncol(lg_cnv_rm)])
		lg_cnv[,ncol(lg_cnv)+1] = 1-rowSums(lg_cnv[,c(2:4)])
		colnames(lg_cnv)=c('SVID','SR','SD','RM','US')
		lg_cnv[,ncol(lg_cnv)+1] = 'US'
		if(nrow(lg_cnv[lg_cnv$RM>.5,])>0){
			lg_cnv[lg_cnv$RM>.5,][,ncol(lg_cnv)]='RM'
		}
		if(nrow(lg_cnv[lg_cnv$SD>.5,])>0){
			lg_cnv[lg_cnv$SD>.5,][,ncol(lg_cnv)]='SD'
		}
		if(nrow(lg_cnv[lg_cnv$SR>.5,])>0){
			lg_cnv[lg_cnv$SR>.5,][,ncol(lg_cnv)]='SR'
		}

 
		dat=merge(dat, lg_cnv[,c(1,6)], by='SVID',all=T)
		dat[!is.na(dat[,ncol(dat)]),][,ncol(dat)-1] = dat[!is.na(dat[,ncol(dat)]),][,ncol(dat)] 
	}
	
	write.table(dat[,c('SVID','GC')],opt$output, quote=F, sep='\t', col.names=T, row.names=F)