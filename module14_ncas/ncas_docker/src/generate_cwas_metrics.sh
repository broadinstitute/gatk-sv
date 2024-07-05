#!R
#Rscript to integrate SVs, APS, conservative, gencode and non coding elements to form the cwas data frame

library("optparse")

option_list = list(
        make_option(c("--permu"), type="character", default=NULL, help="permutation_round", metavar="character"),
        make_option(c("--sv_file_real"), type="character", default=NULL, help="permutation_round", metavar="character"),
        make_option(c("--sv_file_permu"), type="character", default=NULL, help="permutation_round", metavar="character"),
        make_option(c("--sv_vs_gencode"), type="character", default=NULL, help="permutation_round", metavar="character"),
        make_option(c("--sv_vs_conserve"), type="character", default=NULL, help="permutation_round",metavar="character"),
        make_option(c("--sv_vs_noncoding"), type="character", default=NULL, help="permutation_round", metavar="character"),
        make_option(c("--sv_vs_gene"), type="character", default=NULL, help="permutation_round", metavar="character"),
        make_option(c("--sv_vs_coding"), type="character", default=NULL, help="permutation_round", metavar="character"),
        make_option(c("--SVID_aps"), type="character", default=NULL, help="permutation_round", metavar="character"),
        make_option(c("--SVID_genomic_context"), type="character", default=NULL, help="permutation_round", metavar="character"),
        make_option(c("--output"), type="character", default=NULL, help="permutation_round", metavar="character")

 );


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

permu = opt$permu
sv_file_real = opt$sv_file_real
sv_file_permu = opt$sv_file_permu
sv_vs_gene_file = opt$sv_vs_gene
sv_vs_coding_file = opt$sv_vs_coding
sv_vs_gencode_file = opt$sv_vs_gencode
sv_vs_conserve_file = opt$sv_vs_conserve
sv_vs_noncoding_file = opt$sv_vs_noncoding
SVID_aps_file = opt$SVID_aps
SVID_genomic_context = opt$SVID_genomic_context

print('read in permutated SVs ...')
dat=read.table(sv_file_permu, header=T, comment.char="")
print('read in SV vs. gencode ...')
sv_vs_gencode=read.table(sv_vs_gencode_file, header=T, comment.char="")
print('read in SV vs. conserved elements ...')
sv_vs_conserve=read.table(sv_vs_conserve_file, header=T, comment.char="")
print('read in SV vs. nonciding elements ... ')
sv_vs_noncoding = read.table(sv_vs_noncoding_file, header=T, comment.char="")
colnames(sv_vs_noncoding)[c(2:ncol(sv_vs_noncoding))] = paste('nc',colnames(sv_vs_noncoding)[c(2:ncol(sv_vs_noncoding))] , sep='.')
print('read in SVID aps ...')
aps=read.table(SVID_aps_file, header=T)
print('read in SV vs. genes ... ')
vs_gene=read.table(sv_vs_gene_file, sep='')
print('read in SV vs. coding sequences ... ')
vs_coding=read.table(sv_vs_coding_file, sep='')
print('read in real SV calls ... ')
d1=read.table(sv_file_real, header=T, comment.char="")

if(length(colnames(d1)[colnames(d1)=='PREDICTED_INTERGENIC'])==0){
        d1[,ncol(d1)+1] = 'intergenic'
        colnames(d1)[ncol(d1)] = 'PREDICTED_INTERGENIC'
}
if(length(colnames(d1)[colnames(d1)=='PREDICTED_INTERGENIC'])==1){
        d1$PREDICTED_INTERGENIC = as.character(d1$PREDICTED_INTERGENIC)
        d1$PREDICTED_INTERGENIC = 'intergenic'
}
d1[d1$name%in%vs_gene[,4],]$PREDICTED_INTERGENIC = 'genic'
d1[d1$name%in%vs_coding[,4],]$PREDICTED_INTERGENIC = 'coding'

dat2 = merge(dat[,c('X.chr','start','end','name','SVLEN','SVTYPE')], d1[,c('name','svtype','AC','AF','PREDICTED_INTERGENIC')], by='name')
sv_key = dat2[,c('X.chr','start','end','name','svtype','AC','AF','PREDICTED_INTERGENIC','SVLEN','SVTYPE')]
sv_key = merge(sv_key, sv_vs_gencode, by='name',all=T)
sv_key = merge(sv_key, sv_vs_conserve, by='name',all=T)
sv_key = merge(sv_key, sv_vs_noncoding, by='name',all=T)
sv_key[is.na(sv_key)] = 0
sv_key=merge(sv_key, aps, by='name', all=T)
sv_key = sv_key[!is.na(sv_key$start),]

#define list of conservative, gencode and non coding elements
gencode_list = colnames(sv_vs_gencode)[2:ncol(sv_vs_gencode)]
nc_list = colnames(sv_vs_noncoding)[2:ncol(sv_vs_noncoding)]

conserv_list = c("X239prim_DHS_mamm", "X239prim_DHS_prim", "X239prim_footprint_mamm", "X239prim_footprint_prim", "X239prim_uce", "UCE_481", "UCNE", "phastCons100way", "phyloP100way", 'zoonomia_highly_conserved', 'zoonomia_TFBSs','constraint_z')
anti_conserv_list = c('HAR', 'zooHAR', 'zooCHAR', 'zoonomia_actively_evolving','zoonomia_primate_spec')


#generate gencode.and and gencode.non columns
print('calculating parental stats for gencode ... ')
sv_key[,ncol(sv_key)+1] = apply(sv_key[,gencode_list], 1, sum)
colnames(sv_key)[ncol(sv_key)] = 'gencode.any'
sv_key[,ncol(sv_key)+1] = 0
colnames(sv_key)[ncol(sv_key)] = 'gencode.none'
sv_key[sv_key$gencode.any==0,]$gencode.none = 1
sv_key[,ncol(sv_key)+1] = 1
colnames(sv_key)[ncol(sv_key)] = 'gencode.NA'

#generate noncoding.and and noncoding.non columns
print('calculating parental stats for noncoding ... ')
sv_key[,ncol(sv_key)+1] = apply(sv_key[,nc_list], 1, sum)
colnames(sv_key)[ncol(sv_key)] = 'noncoding.any'
sv_key[,ncol(sv_key)+1] = 0
colnames(sv_key)[ncol(sv_key)] = 'noncoding.none'
sv_key[sv_key$noncoding.any==0,]$noncoding.none = 1
sv_key[,ncol(sv_key)+1] = 1
colnames(sv_key)[ncol(sv_key)] = 'noncoding.NA'

#generate additional conservaite columns
print('calculating parental stats for conserved elements ... ')
sv_key[,ncol(sv_key)+1] = 0
colnames(sv_key)[ncol(sv_key)] = 'constrained_z_cate'
sv_key[sv_key$constraint_z=='z_over_2',]$constrained_z_cate = 1
sv_key[sv_key$constraint_z=='z_over_4',]$constrained_z_cate = 2

sv_key[,ncol(sv_key)+1]= apply(sv_key[,conserv_list], 1, sum)
colnames(sv_key)[ncol(sv_key)] = 'conserve.any'
sv_key[,ncol(sv_key)+1] = 0
colnames(sv_key)[ncol(sv_key)] = 'conserve.none'
sv_key[sv_key$conserve.any==0,]$conserve.none = 1

sv_key[,ncol(sv_key)+1]= apply(sv_key[,anti_conserv_list], 1, sum)
colnames(sv_key)[ncol(sv_key)] = 'anti_conserve.any'
sv_key[,ncol(sv_key)+1] = 0
colnames(sv_key)[ncol(sv_key)] = 'anti_conserve.none'
sv_key[sv_key$anti_conserve.any==0,]$anti_conserve.none = 1

sv_key[,ncol(sv_key)+1] = 1
colnames(sv_key)[ncol(sv_key)] = 'conserve.NA'

#add genomic context
print('read in genomic context ... ')
sv_gc=read.table(SVID_genomic_context)
colnames(sv_gc) = c('name','genomic_context')
sv_key=merge(sv_key, sv_gc, by='name')

#save the sv_key metrics
print('write output file ...')
save(sv_key, file = opt$output)
