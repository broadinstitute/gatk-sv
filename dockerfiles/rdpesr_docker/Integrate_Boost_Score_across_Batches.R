#!R
#!/usr/bin/env Rscript
library("optparse")

option_list = list(
        make_option(c( "--chr"), type="character", default=NULL, help="chromosome name", metavar="character"),
        make_option(c( "--prefix"), type="character", default=NULL, help="prefix of output file name", metavar="character")
 );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

chr_name = opt$chr
stat_list=c()
for(i in list.files('./')){
	if(length(strsplit(as.character(i),'[.]')[[1]]) > 1){
		if(strsplit(as.character(i),'[.]')[[1]][length(strsplit(as.character(i),'[.]')[[1]])]=='stat' & grepl(chr_name,as.character(i))) {
			stat_list = c(stat_list, i)
		}
	}
}

stat = read.table(stat_list[1], header=F)
stat[,ncol(stat)+1 ] = paste(stat[,1],stat[,2],stat[,3],sep='.')
colnames(stat)[ncol(stat)] = 'SVID'
stat=stat[,c(4:ncol(stat))]
for(i in stat_list[2:length(stat_list)]){
	print(i)
	stat2 = read.table(i, header=F)
	stat2[,ncol(stat2)+1 ] = paste(stat2[,1],stat2[,2],stat2[,3],sep='.')
	colnames(stat2)[ncol(stat2)] = 'SVID'
	stat=merge(stat, stat2[,c(4:ncol(stat2))], by='SVID', all=T)
	stat[is.na(stat)]=0
	stat[,2] = stat[,2]+stat[,4]
	stat[,3] = stat[,3]+stat[,5]
	stat=stat[,c(1:3)]
}
colnames(stat)[c(2:3)] = c('count_samp_fail','count_samp_pass')
write.table(stat,paste(opt$prefix, opt$chr, 'stat', sep='.'), quote=F, sep='\t', col.names=T, row.names=F)

