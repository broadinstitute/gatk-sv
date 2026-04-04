#integrate SVs vs. gene overlaps
#!R
library("optparse")


unify_SVID_gene_list<-function(d1, gene){
    d1[,1] = paste(d1[,1],d1[,2],sep=',')
    d1 = d1[,-2]
    stat=data.frame(table(d1[,1]))
    uni = d1[d1[,1]%in%stat[stat[,2]==1,][,1],]
    uni[,2] = as.character(uni[,2])
    uni[,3] = as.character(uni[,3])
    print(table(stat[,2]))
    colnames(uni) = c('key', 'gene_id', 'gene_symbol')
    for(i in unique(stat[,2])){
            print(i)
            if(i>1){
				# Instead of growing inside the loop:
                svid_list = stat[stat[,2]==i,][,1]
				new_rows <- lapply(svid_list, function(svid) {
				  tmp <- d1[d1[,1] == svid, c(2,3)]
				  # ... merge and order ...
				  data.frame(key = svid,
				             gene_id     = paste(tmp[,1], collapse = ","),
				             gene_symbol = paste(tmp[,2], collapse = ","))
				})
				uni <- rbind(uni, do.call(rbind, new_rows))
            }
    }
    uni[,4] = apply(uni,1,function(x){strsplit(as.character(x[1]),',')[[1]][1]})
    uni[,5] = apply(uni,1,function(x){strsplit(as.character(x[1]),',')[[1]][2]})
    uni = uni[,c(4,5,2,3)]
    return(uni)
}



option_list = list(
        make_option(c('-i', "--input"), type="character", default=NULL, help="input file", metavar="character"),
        make_option(c('-g', "--gtf"), type="character", default=NULL, help="gtf file", metavar="character"),
        make_option(c('-o', "--output"), type="character", default=NULL, help="output file", metavar="character")
 );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


svid_als_file = opt$aps

d1=read.table(opt$input)
gene=read.table(opt$gtf)
colnames(gene) = c('gene_chr','gene_pos','gene_end','gene_strands','gene_id','gene_symbol')
d1_uniq = unify_SVID_gene_list(d1, gene)
colnames(d1_uniq) = c('SVID','svtype','gene_id','gene_name')
write.table(d1_uniq, opt$output, quote=F, sep='\t', col.names=T, row.names=F)


