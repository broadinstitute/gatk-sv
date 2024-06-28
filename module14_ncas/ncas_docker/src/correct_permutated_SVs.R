#!R
#script to correct permutated SV

library("optparse")

option_list = list(
        make_option(c('-i', "--input"), type="character", default=NULL, help="input file : permutated SVs", metavar="character"),
        make_option(c('-r', "--real_data"), type="character", default=NULL, help="input file : permutated SVs", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

sv=read.table(opt$real_data, header=T, comment.char="")
sv[,ncol(sv)+1] = sv[,3] - sv[,2]
colnames(sv)[ncol(sv)] = 'dis_bp'
permu=read.table(opt$input, header=T, comment.char="")
permu=merge(permu, sv[,c('name','dis_bp')])
permu$end = permu$start + permu$dis_bp
out = permu[,c(2:4,1,5:(ncol(permu)-1))]
colnames(out)[1] = '#chr'
out = out[order(out[,3]),]
out = out[order(out[,2]),]
out = out[order(out[,1]),]
write.table(out, gsub('.gz','.corrected',opt$input), quote=F, sep='\t', col.names=T, row.names=F)

