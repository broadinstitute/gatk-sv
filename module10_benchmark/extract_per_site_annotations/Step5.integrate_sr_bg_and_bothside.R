add_bg_fail_and_sr_bg_info<-function(chr){
        d1=read.table(paste('gnomad-sv-v3.',chr,'.final_cleanup.info', sep=''), header=T, sep='\t', comment.char="")
        d2=read.table(paste('drop_redun/gnomad-sv-v3.',chr,'.drop_redundant_cnvs.sorted.reheadered.sites.gz', sep=''), header=T, sep='\t', comment.char="")
        d3=read.table(paste('stitch_fix/gnomad-sv-v3.',chr,'.stitch_fragmented_cnvs.sites.gz', sep=''),header=T, sep='\t', comment.char="")
        dat=merge(d1, d3, by=c('X.chrom','start','end','CHR2'))
        bg_fail=read.table(paste('bg_fail/gnomad-sv-v3.',chr,'.',chr,'.sr_background_fail.updated3.txt', sep=''))
        bs_supp=read.table(paste('sr_bothside/gnomad-sv-v3.',chr,'.',chr,'.sr_bothside_pass.updated3.txt', sep=''))

        #bg_fail_svs = dat[dat$MEMBERS%in%bg_fail[,ncol(bg_fail)-1],]
        #bs_supp_svs = dat[dat$MEMBERS%in%bs_supp[,ncol(bs_supp)-1],]

        bg_fail_svs = dat[dat$name.y%in%bg_fail[,ncol(bg_fail)],]
        bs_supp_svs = dat[dat$name.y%in%bs_supp[,ncol(bs_supp)],]


        d1[,ncol(d1)+1] = 0
        d1[d1[,4]%in%bg_fail_svs[,5],][,ncol(d1)]=1
        colnames(d1)[ncol(d1)]='sr_background_fail'

        d1[,ncol(d1)+1] = 0
        d1[d1[,4]%in%bs_supp_svs[,5],][,ncol(d1)]=1
        colnames(d1)[ncol(d1)]='sr_bothside_support'

        colnames(d1)[1]='#chrom'
        write.table(d1, paste('gnomad-sv-v3.',chr,'.final_cleanup.info2', sep=''), quote=F, sep='\t', col.names=T, row.names=F)
}

for(chr in c("chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11",'chr13',"chr14","chr15","chr17","chr18","chr19","chr20","chr21","chr22","chrY")){
        print(chr)
        add_bg_fail_and_sr_bg_info(chr)
}


