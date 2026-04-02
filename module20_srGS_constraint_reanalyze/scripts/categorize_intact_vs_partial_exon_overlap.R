
#!extract these overlaps from exon overlap:
#!R

library("optparse")

option_list = list(
        make_option(c('-c', "--CDS"), type="character", default=NULL, help="CDS input: comparison results between SVs and CDSs in gtf", metavar="character"),
        make_option(c('-g', "--sv_gene"), type="character", default=NULL, help="3 column input linking genes with the SVs that entirely fell within the transcripts", metavar="character"),
        make_option(c('-p', "--prefix"), type="character", default=NULL, help="prefix of output file", metavar="character")
 );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

d1=read.table(opt$CDS)
gene_id=read.table(opt$sv_gene)
colnames(d1)[c(4,5,10,11)] = c('SVID','SVTYPE','gene_id','gene_name')
colnames(gene_id) = c('SVID','SVTYPE','gene_id','gene_name')

dat = merge(d1, gene_id, by=c('SVID','SVTYPE','gene_id','gene_name'))
dat[,ncol(dat)+1] = 'intact_exon_overlap'
dat[dat[,6]>dat[,9] & dat[,6]<dat[,10],][,ncol(dat)] = 'partial_exon_overlap'
dat[dat[,7]>dat[,9] & dat[,7]<dat[,10],][,ncol(dat)] = 'partial_exon_overlap'
colnames(dat)[ncol(dat)] = 'overlap_category'
dat[,ncol(dat)+1] = paste(dat$SVID, dat$gene_id, sep=',')
colnames(dat)[ncol(dat)] = 'svid_geneid'
partial_overlap_svid_gene = dat[dat$overlap_category=='partial_exon_overlap',]
dat[dat$svid_geneid%in%partial_overlap_svid_gene$svid_geneid,]$overlap_category = 'partial_exon_overlap'

sv_exon_stat = data.frame(table(dat$svid_geneid))
colnames(sv_exon_stat) = c('svid_geneid','count_exons')
out = unique(dat[,c('SVID','SVTYPE','gene_id','gene_name','overlap_category','svid_geneid')])
out = merge(out, sv_exon_stat, by='svid_geneid')

partial_overlap = out[out$overlap_category =='partial_exon_overlap',c('SVID','SVTYPE','gene_id','gene_name','count_exons')]
intact_overlap  = out[out$overlap_category =='intact_exon_overlap',c('SVID','SVTYPE','gene_id','gene_name','count_exons')]
write.table(partial_overlap, paste(opt$prefix, '.partial_exon_overlap', sep=''), quote=F, sep='\t', col.names=F, row.names=F)
write.table(intact_overlap,  paste(opt$prefix, '.intact_exon_overlap',  sep=''), quote=F, sep='\t', col.names=F, row.names=F)

