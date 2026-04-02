#integrate SVs across gene overlaps:
#!R

library("optparse")

option_list = list(
        make_option(c('-p', "--prefix"), type="character", default=NULL, help="prefix of input and output files", metavar="character"),
        make_option(c('-o', "--output"), type="character", default=NULL, help="output file name", metavar="character")
 );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

d1a=read.table(paste(opt$prefix, '.inside_exons.reorganized', sep=''), header=T)
colnames(d1a)[4] = 'inside_exons'
d1b=read.table(paste(opt$prefix, '.whole_transcript_overlap.reorganized', sep=''), header=T)
colnames(d1b)[4] = 'whole_transcript_overlap'
d1=merge(d1a[,c(1,2,4)], d1b[,c(1,2,4)], by=c('SVID','svtype'), all=T)
print(dim(d1))

d2=read.table(paste(opt$prefix, '.intact_exon_overlap.reorganized', sep=''), header=T)
colnames(d2)[4] = 'intact_exon_overlap'
dat=merge(d1, d2[,c(1,2,4)], by=c('SVID','svtype'), all=T)

d3=read.table(paste(opt$prefix, '.partial_exon_overlap.reorganized', sep=''), header=T)
colnames(d3)[4] = 'partial_exon_overlap'
dat=merge(dat, d3[,c(1,2,4)], by=c('SVID','svtype'), all=T)

d4b=read.table(paste(opt$prefix, '.tss_transcripts_overlap.reorganized', sep=''), header=T)
colnames(d4b)[4] = 'tss_transcripts_overlap'
dat=merge(dat, d4b[,c(1,2,4)], by=c('SVID','svtype'), all=T)

d5b=read.table(paste(opt$prefix, '.partial_transcripts_overlap.reorganized', sep=''), header=T)
colnames(d5b)[4] = 'partial_transcripts_overlap'
dat=merge(dat, d5b[,c(1,2,4)], by=c('SVID','svtype'), all=T)

d6 = read.table(paste(opt$prefix, '.5_prime_utr.reorganized', sep=''), header=T)
colnames(d6)[4] = 'X5_prime_utr'
dat=merge(dat, d6[,c(1,2,4)], by=c('SVID','svtype'), all=T)

d7 = read.table(paste(opt$prefix, '.3_prime_utr.reorganized', sep=''), header=T)
colnames(d7)[4] = 'X3_prime_utr'
dat=merge(dat, d7[,c(1,2,4)], by=c('SVID','svtype'), all=T)

d8=read.table(paste(opt$prefix, '.inside_introns.reorganized', sep=''), header=T)
colnames(d8)[4] = 'inside_introns'
dat=merge(dat, d8[,c(1,2,4)], by=c('SVID','svtype'), all=T)

d9=read.table(paste(opt$prefix, '.promoter.reorganized', sep=''), header=T)
colnames(d9)[4] = 'promoter'
dat=merge(dat, d9[,c(1,2,4)], by=c('SVID','svtype'), all=T)

for(i in c(3:ncol(dat))){dat[,i] = as.character(dat[,i])}

dat[!is.na(dat$partial_exon_overlap),]$partial_exon_overlap = apply(dat[!is.na(dat$partial_exon_overlap),c('partial_exon_overlap','inside_exons')],1,function(x){paste(strsplit(as.character(x[1]),',')[[1]][!strsplit(as.character(x[1]),',')[[1]]%in%strsplit(as.character(x[2]),',')[[1]]], collapse=',')})
dat[!is.na(dat$partial_transcripts_overlap),]$partial_transcripts_overlap = apply(dat[!is.na(dat$partial_transcripts_overlap),c('partial_transcripts_overlap','X3_prime_utr')],1,function(x){paste(strsplit(as.character(x[1]),',')[[1]][!strsplit(as.character(x[1]),',')[[1]]%in%strsplit(as.character(x[2]),',')[[1]]], collapse=',')})
dat[!is.na(dat$tss_transcripts_overlap),]$tss_transcripts_overlap = apply(dat[!is.na(dat$tss_transcripts_overlap),c('tss_transcripts_overlap','X5_prime_utr')],1,function(x){paste(strsplit(as.character(x[1]),',')[[1]][!strsplit(as.character(x[1]),',')[[1]]%in%strsplit(as.character(x[2]),',')[[1]]], collapse=',')})
dat[!is.na(dat$promoter_overlap),]$promoter_overlap = apply(dat[!is.na(dat$promoter_overlap),c('promoter_overlap','X5_prime_utr')],1,function(x){paste(strsplit(as.character(x[1]),',')[[1]][!strsplit(as.character(x[1]),',')[[1]]%in%strsplit(as.character(x[2]),',')[[1]]], collapse=',')})
dat[!is.na(dat$promoter_overlap),]$promoter_overlap = apply(dat[!is.na(dat$promoter_overlap),c('promoter_overlap','tss_transcripts_overlap')],1,function(x){paste(strsplit(as.character(x[1]),',')[[1]][!strsplit(as.character(x[1]),',')[[1]]%in%strsplit(as.character(x[2]),',')[[1]]], collapse=',')})
dat[!is.na(dat$promoter_overlap),]$promoter_overlap = apply(dat[!is.na(dat$promoter_overlap),c('promoter_overlap','whole_transcript_overlap')],1,function(x){paste(strsplit(as.character(x[1]),',')[[1]][!strsplit(as.character(x[1]),',')[[1]]%in%strsplit(as.character(x[2]),',')[[1]]], collapse=',')})


dat[dat==""] = NA

write.table(dat, opt$output, sep='\t', quote=F, col.names=T, row.names=F)







