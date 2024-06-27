#!R
#Rscript to generate permutation of SVs

library("optparse")

option_list = list(
  make_option(c('-p', "--permutate"), type="character", default=NULL, help="permutation round", metavar="character"),
  make_option(c('-g', "--genome"), type="character", default=NULL, help="genome file", metavar="character"),
  make_option(c('-i', "--input"), type="character", default=NULL, help="name of input file", metavar="character"),
  make_option(c('-o', "--output"), type="character", default=NULL, help="name of output file", metavar="character")
 );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
permu = paste('permu',opt$permutate, sep='_')

reorganize_chromosome_arms<-function(genome_length){
  out = genome_length[1,]
  for(i in unique(genome_length[,1])){
    tmp = genome_length[genome_length[,1]==i,]
    tmp = tmp[tmp[,4]%in%c('telomere','centromere'),]
    tmp = tmp[order(tmp[,2]),]
    out[nrow(out)+1,1] = i
    out[nrow(out),2] = tmp[1,3]
    out[nrow(out),3] = tmp[2,2]
    out[nrow(out),4] = 'p'
    out[nrow(out)+1,1] = i
    out[nrow(out),2] = tmp[2,3]
    out[nrow(out),3] = tmp[3,2]
    out[nrow(out),4] = 'q'
  }
  out = out[-1,]
  return(out)
}


#readin transcripts
print("reading in SV coordinates ... ")
trans=read.table(opt$input, header=T,comment.char="")
#readin telomere and centeromere regions:
genome_length = read.table(opt$genome, sep = '\t')
chr_arms = reorganize_chromosome_arms(genome_length)


permutate_round = as.integer(opt$permutate)
trans_permu = trans[1,]
set.seed(permutate_round)
permu_set = 1000
rec_flag = 0
while(TRUE){
        print(nrow(trans_permu))
        if(nrow(trans)> permu_set){
                out = trans[sample(1:nrow(trans), permu_set, replace=F),]
                chrom = chr_arms[sample(c(1:nrow(chr_arms)),1),]
                out[,1] = chrom[1,1]
                pos_list = sample(c(chrom[1,2]:chrom[1,3]), nrow(out), replace=F)
                out[,2] = pos_list
                out[,3] = out[,2] +out$SVLEN     
                out = out[out$end< chrom[1,3],] 
                trans_permu = rbind(trans_permu, out)
                trans = trans[!trans$name%in%trans_permu$name,]    
        }
        if(nrow(trans)<=permu_set){
                out = trans
                chrom = chr_arms[sample(c(1:nrow(chr_arms)),1),]
                out[,1] = chrom[1,1]
                pos_list = sample(c(chrom[1,2]:chrom[1,3]), nrow(out), replace=F)
                out[,2] = pos_list
                out[,3] = out[,2] +out$SVLEN     
                out = out[out$end< chrom[1,3],] 
                trans_permu = rbind(trans_permu, out)
                trans = trans[!trans$name%in%trans_permu$name,]    
        }
        if(nrow(trans)==0){ break }

}
trans_permu = trans_permu[-1,]
trans_permu = trans_permu[order(trans_permu[,3]),]
trans_permu = trans_permu[order(trans_permu[,2]),]
trans_permu = trans_permu[order(trans_permu[,1]),]
colnames(trans_permu)[1] = '#chr'
write.table(trans_permu, opt$output, quote=F, sep='\t', col.names=T, row.names=F)

