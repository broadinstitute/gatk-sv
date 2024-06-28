#!R
#Rscript to integrate SVs, APS, conservative, gencode and non coding elements to form the cwas data frame

library("optparse")

option_list = list(
        make_option(c("--sv_file_real"), type="character", default=NULL, help="input SV file", metavar="character"),
        make_option(c("--sv_file_permu"), type="character", default=NULL, help="input SV file", metavar="character"),
        make_option(c("--sv_vs_gencode"), type="character", default=NULL, help="input", metavar="character"),
        make_option(c("--sv_vs_conserve"), type="character", default=NULL, help="input", metavar="character"),
        make_option(c("--sv_vs_noncoding"), type="character", default=NULL, help="input", metavar="character"),
        make_option(c("--sv_vs_gene"), type="character", default=NULL, help="input", metavar="character"),
        make_option(c("--sv_vs_coding"), type="character", default=NULL, help="input", metavar="character"),
        make_option(c("--aps"), type="character", default=NULL, help="input", metavar="character"),
        make_option(c("--output"), type="character", default=NULL, help="output rData", metavar="character")

 );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

permu = opt$permu
print('read in permutated SVs ...')
dat=read.table(opt$sv_file_permu, header=T, comment.char="")
print('read in SV vs. conservative elements ...')
sv_vs_conserve  = read.table(opt$sv_vs_conserve,  header=T, comment.char="")
print('read in SV vs. gencode elements ...')
sv_vs_gencode   = read.table(opt$sv_vs_gencode,   header=T, comment.char="")
print('read in SV vs. noncoding elements ...')
sv_vs_noncoding = read.table(opt$sv_vs_noncoding, header=T, comment.char="")
colnames(sv_vs_noncoding)[c(2:ncol(sv_vs_noncoding))] = paste('nc',colnames(sv_vs_noncoding)[c(2:ncol(sv_vs_noncoding))] , sep='.')
print('read in SVID vs. APS ...')
aps=read.table(opt$aps, header=T)

print('read in SVs ...')
d1=read.table(opt$sv_file_permu, header=T)
print('read in SV vs. genes ...')
vs_gene=read.table(opt$sv_vs_gene)
print('read in SV vs. coding')
vs_coding=read.table(opt$sv_vs_coding)
d1$PREDICTED_INTERGENIC = 'intergenic'
d1[d1$name%in%vs_gene[,4],]$PREDICTED_INTERGENIC = 'genic'
d1[d1$name%in%vs_coding[,4],]$PREDICTED_INTERGENIC = 'coding'

dat2 = merge(dat, d1[,c('name','svtype','AC','AF','PREDICTED_INTERGENIC')], by='name')
sv_key = dat2[,c('X.chr','start','end','name','svtype','AC','AF','PREDICTED_INTERGENIC','SVLEN','SVTYPE')]
sv_key = merge(sv_key, sv_vs_conserve, by='name',all=T)
sv_key = merge(sv_key, sv_vs_gencode, by='name',all=T)
sv_key = merge(sv_key, sv_vs_noncoding, by='name',all=T)
sv_key[is.na(sv_key)] = 0
sv_key=merge(sv_key, aps, by='name', all=T)
colnames(sv_key) = gsub('.over_50perc_ovr','', gsub(paste('gnomAD_SV_v3.',permu,'.vs.', sep=''),'', colnames(sv_key)))
#define list of conservative, gencode and non coding elements
conserv_list = c("X239prim_DHS_mamm","X239prim_DHS_prim","X239prim_footprint_mamm","X239prim_footprint_prim","X239prim_uce","phastCons100way","phyloP100way","constraint_z")
gencode_list = c("gencode.pseudogene_gene","gencode.pseudogene_exon","gencode.lncRNA.gene","gencode.miRNA.gene","gencode.lncRNA.exon","gencode.miRNA.exon","gencode.intron","gencode.protein_coding.exon","gencode.UTR","gencode.protein_coding.gene","gencode.snRNA.gene","gencode.snRNA.exon","gencode.misc_RNA.gene","gencode.misc_RNA.exon","gencode.tRNA","gencode.snoRNA.gene","gencode.snoRNA.exon")
nc_list = c("nc.abc.cognate.loeuf_constraint","nc.abc.cognate.loeuf_middle","nc.abc.cognate.loeuf.phi.pts.tall_v0.1","nc.abc.cognate.loeuf_unconstraint","nc.abc.cognate.phaplo_constraint","nc.abc.cognate.phaplo_middle","nc.abc.cognate.phaplo_unconstraint","nc.abc.cognate.ptriplo_constraint","nc.abc.cognate.ptriplo_middle","nc.abc.cognate.ptriplo_unconstraint","nc.atac.seq_OCRs_2020.autosomes","nc.atac.seq_pREs_enhancers_2020.autosomes","nc.chromHMM.Enh.txEnh","nc.chromHMM.insulator","nc.chromHMM.polycomb","nc.chromHMM.promoter.TSS","nc.DCR2.adrenal_glands","nc.DCR2.brain","nc.DCR2.eye","nc.DCR2.face","nc.DCR2.general","nc.DCR2.heart","nc.DCR2.intestine","nc.DCR2.kidney","nc.DCR2.limb","nc.DCR2.liver","nc.DCR2.lung","nc.DCR2.muscle","nc.DCR2.neural_tube","nc.DCR2.placenta","nc.DCR2.skin","nc.DCR2.stomach","nc.DNaseHS.UCSC","nc.Encode3_ccREs.CTCF.only.CTCF.bound","nc.Encode3_ccREs.dELS.CTCF.bound","nc.Encode3_ccREs.dELS","nc.Encode3_ccREs.DNase.H3K4me3.CTCF.bound","nc.Encode3_ccREs.DNase.H3K4me3","nc.Encode3_ccREs.pELS.CTCF.bound","nc.Encode3_ccREs.pELS","nc.Encode3_ccREs.PLS.CTCF.bound","nc.Encode3_ccREs.PLS","nc.encode4.Bivalent","nc.encode4.ConstitutiveHet","nc.encode4.CTCF","nc.encode4.EnhancerLow","nc.encode4.Enhancer","nc.encode4.FacultativeHet","nc.encode4.K9K36","nc.encode4.PromoterFlanking","nc.encode4.Promoter","nc.encode4.Transcribed","nc.enhancerAtlasV2","nc.enhancer_DCR2v3.autosomes","nc.fantom5_enhancers.all_cells","nc.fantom5_enhancers.all_organs","nc.fantom5_enhancers.enhancer_clusters","nc.fantom5_enhancers.permissive_enhancers","nc.fantom5_enhancers.robust_enhancers","nc.fantom5_enhancers.specific_enhancers_cells","nc.fantom5_enhancers.specific_enhancers_organs","nc.fantom5_enhancers.ubiquitous_enhancers_cells","nc.fantom5_enhancers.ubiquitous_enhancers_organs","nc.fantom_enhancers.enhancers.siwei_paper","nc.FetalBrain_2685.autosomes","nc.HAR.2023","nc.human_fetal_brain_TF_footprints.autosomes","nc.masked.regions.LCRs.gaps.autosomes","nc.NDD_664.autosomes","nc.Neuron_PROcap_enhancers.autosomes","nc.SE_package.lft38","nc.Shi_Lab.TAD_Boundry","nc.TFChip.UCSC","nc.UCNE_coord","nc.vista_cores.neg","nc.vista_cores.pos")

sv_key = sv_key[!is.na(sv_key$start),]
#generate gencode.and and gencode.non columns
sv_key[,ncol(sv_key)+1] = apply(sv_key[,gencode_list], 1, sum)
colnames(sv_key)[ncol(sv_key)] = 'gencode.any'
sv_key[,ncol(sv_key)+1] = 0
colnames(sv_key)[ncol(sv_key)] = 'gencode.none'
sv_key[sv_key$gencode.any==0,]$gencode.none = 1
sv_key[,ncol(sv_key)+1] = 1
colnames(sv_key)[ncol(sv_key)] = 'gencode.NA'
#generate noncoding.and and noncoding.non columns
sv_key[,ncol(sv_key)+1] = apply(sv_key[,nc_list], 1, sum)
colnames(sv_key)[ncol(sv_key)] = 'noncoding.any'
sv_key[,ncol(sv_key)+1] = 0
colnames(sv_key)[ncol(sv_key)] = 'noncoding.none'
sv_key[sv_key$noncoding.any==0,]$noncoding.none = 1
sv_key[,ncol(sv_key)+1] = 1
colnames(sv_key)[ncol(sv_key)] = 'noncoding.NA'
#generate additional conservaite columns
sv_key[,ncol(sv_key)+1] = 0
colnames(sv_key)[ncol(sv_key)] = 'constrained_z_cate'
sv_key[sv_key$constraint_z=='z_over_2',]$constrained_z_cate = 1
sv_key[sv_key$constraint_z=='z_over_4',]$constrained_z_cate = 2
conserv_list_numeric = c("X239prim_DHS_mamm","X239prim_DHS_prim","X239prim_footprint_mamm","X239prim_footprint_prim","X239prim_uce","phastCons100way","phyloP100way","constrained_z_cate")
sv_key[,ncol(sv_key)+1]= apply(sv_key[,conserv_list_numeric], 1, sum)
colnames(sv_key)[ncol(sv_key)] = 'conserve.any'
sv_key[,ncol(sv_key)+1] = 0
colnames(sv_key)[ncol(sv_key)] = 'conserve.none'
sv_key[sv_key$conserve.any==0,]$conserve.none = 1
sv_key[,ncol(sv_key)+1] = 1
colnames(sv_key)[ncol(sv_key)] = 'conserve.NA'
#add genomic context
sv_gc=read.table('../../../gnomAD_SV_v3.sites.Genomic_Context.gz')
colnames(sv_gc) = c('name','genomic_context')
sv_key=merge(sv_key, sv_gc, by='name')
#save the sv_key metricsi
print('generating final output ... ')
save(sv_key, file = opt$output)


