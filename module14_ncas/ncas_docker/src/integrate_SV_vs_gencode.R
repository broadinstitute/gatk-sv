#reorganize SV vs. gencode
#!R
library("optparse")

option_list = list(
        make_option(c('-i', "--input"), type="character", default=NULL, help="input file", metavar="character"),
        make_option(c('-o', "--output"), type="character", default=NULL, help="output file", metavar="character")
 );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

sv_vs_gencode=read.table(opt$input, header=T, comment.char="")
gencode_list=unique(sv_vs_gencode$gencode_cate)
out = data.frame(table(sv_vs_gencode[sv_vs_gencode$gencode_cate==gencode_list[1],]$name))
colnames(out)=c('name',paste('gencode', gencode_list[1],sep='.'))
for(i in c(2:length(gencode_list))){
	print(i)
	tmp = data.frame(table(sv_vs_gencode[sv_vs_gencode$gencode_cate==gencode_list[i],]$name))
	colnames(tmp) = c('name',paste('gencode', gencode_list[i],sep='.'))
	out = merge(out, tmp, by='name', all=T)
}
write.table(out, opt$output, quote=F, sep='\t', col.names=T, row.names=F)

