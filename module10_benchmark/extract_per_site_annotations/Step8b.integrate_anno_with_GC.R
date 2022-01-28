#!R 
#integrate GC of each breakpoints:
add_genomic_content_anno<-function(chr){
        dat=read.table(paste('gnomad-sv-v3.',chr,'.final_cleanup.info3', sep=''),sep='\t', header=T, comment.char="")
        sr_le=read.table(paste('US_RM_SD_SR/gnomad-sv-v3.',chr,'.le_bp.vs.SR', sep=''))
        sd_le=read.table(paste('US_RM_SD_SR/gnomad-sv-v3.',chr,'.le_bp.vs.SD', sep=''))
        rm_le=read.table(paste('US_RM_SD_SR/gnomad-sv-v3.',chr,'.le_bp.vs.RM', sep=''))

        sr_ri=read.table(paste('US_RM_SD_SR/gnomad-sv-v3.',chr,'.ri_bp.vs.SR', sep=''))
        sd_ri=read.table(paste('US_RM_SD_SR/gnomad-sv-v3.',chr,'.ri_bp.vs.SD', sep=''))
        rm_ri=read.table(paste('US_RM_SD_SR/gnomad-sv-v3.',chr,'.ri_bp.vs.RM', sep=''))
        dat[,ncol(dat)+1] = 'US'
        dat[dat[,4]%in%rm_le[,4],][,ncol(dat)]='RM'
        dat[dat[,4]%in%rm_ri[,4],][,ncol(dat)]='RM'
        dat[dat[,4]%in%sd_le[,4],][,ncol(dat)]='SD'
        dat[dat[,4]%in%sd_ri[,4],][,ncol(dat)]='SD'
        dat[dat[,4]%in%sr_le[,4],][,ncol(dat)]='SR'
        dat[dat[,4]%in%sr_ri[,4],][,ncol(dat)]='SR'

        lg_cnv_sr=read.table(paste('US_RM_SD_SR/gnomad-sv-v3.',chr,'.lg_cnv.vs.SR', sep=''))
        lg_cnv_sd=read.table(paste('US_RM_SD_SR/gnomad-sv-v3.',chr,'.lg_cnv.vs.SD', sep=''))
        lg_cnv_rm=read.table(paste('US_RM_SD_SR/gnomad-sv-v3.',chr,'.lg_cnv.vs.RM', sep=''))
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
        write.table(dat[,c(1:(ncol(dat)-1))],paste('gnomad-sv-v3.',chr,'.final_cleanup.info4',sep=''), quote=F, sep='\t', col.names=T, row.names=F)
}

for (i in c('chr2','chr3','chr4','chr5','chr6','chr7', 'chr8', 'chr9', 'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17', 'chr18', 'chr19', 'chr20','chr21','chr22','chrY')){
        print(i)
        add_genomic_content_anno(i)
}

