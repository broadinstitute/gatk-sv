#!R
#!Rscript to integrate batch effect comparison statistics
library("optparse")

option_list = list(
        make_option(c( "--stat_plus_vs_minus"), type="character", default=NULL, help="list of all samples", metavar="character"),
        make_option(c( "--stat_pre_vs_post"), type="character", default=NULL, help="list of all samples", metavar="character"),
        make_option(c( "--stat_pairwise_minus"), type="character", default=NULL, help="loose union of all filters", metavar="character"),
        make_option(c( "--stat_pairwise_plus"), type="character", default=NULL, help="loose union of all filters", metavar="character"),
        make_option(c( "--stat_one_vs_all_minus"), type="character", default=NULL, help="loose union of all filters", metavar="character"),
        make_option(c( "--stat_one_vs_all_plus"), type="character", default=NULL, help="loose union of all filters", metavar="character"),
        make_option(c( "--out_fail_info"), type="character", default=NULL, help="loose union of all filters", metavar="character"),
        make_option(c( "--out_fail_filter"), type="character", default=NULL, help="loose union of all filters", metavar="character"),
        make_option(c( "--out_fail_format"), type="character", default=NULL, help="loose union of all filters", metavar="character")
 );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

add_bonferroni_correction<-function(stat1, pv_cff = .05){
	correction_fac = nrow(stat1[!is.na(stat1[,2]),])
	colnames(stat1)[1] = 'VID'
	fail_info = stat1[!is.na(stat1[,2]) & stat1[,2]< pv_cff/correction_fac , ]
	fail_filter = stat1[!is.na(stat1[,2]) & stat1[,2]< pv_cff/(correction_fac*5) , ]
	
	stat1[,3] = 'pass'
	stat1[stat1[,1]%in%fail_info[,1],][,3] = 'fail'
	stat1[,4] = 'pass'
	stat1[stat1[,1]%in%fail_filter[,1],][,4] = 'fail'
	colnames(stat1)[3] = 'info'
	colnames(stat1)[4] = 'filter'

	return(stat1)
}

add_bonferroni_per_batch<-function(stat4){
	out_metrics = data.frame(stat4[,c(1,2)])
	out_metrics = add_bonferroni_correction(out_metrics)
	colnames(out_metrics)[3] = gsub('chisq_p_','',colnames(out_metrics)[2])
	out_metrics = out_metrics[,c(1,3)]
	out_metrics[out_metrics[,2]=='fail',][,2] = colnames(out_metrics)[2]
	for(i in c(3:ncol(stat4))){
		tmp_metrics = data.frame(stat4[,c(1,i)])
		tmp_metrics = add_bonferroni_correction(tmp_metrics)
		colnames(tmp_metrics)[3] = gsub('chisq_p_','',colnames(tmp_metrics)[2])
		tmp_metrics = tmp_metrics[,c(1,3)]
		tmp_metrics[tmp_metrics[,2]=='fail',][,2] = colnames(tmp_metrics)[2]

		out_metrics = cbind(out_metrics,tmp_metrics[,2])
		colnames(out_metrics)[ncol(out_metrics)] = colnames(tmp_metrics)[ncol(tmp_metrics)]
	}
	return(out_metrics)
}

#read in batch effect comparison results as a 2-column file with SVID and p-value
print("readin PCR+ vs. PCR- comparison statistics ... ")
stat1_pre=read.table(opt$stat_plus_vs_minus, header=T, comment.char="")

print("readin pre- vs. post- filtering comparison statistics ... ")
stat2=read.table(opt$stat_pre_vs_post, header=T, comment.char="")
stat2a_pre = stat2[,c('X.VID','p.PCRMINUS')]
stat2b_pre = stat2[,c('X.VID','p.PCRPLUS')]

print("readin pairwise comparison statistics ... ")
stat3a_pre = read.table(opt$stat_pairwise_minus, header=T, comment.char="")
stat3b_pre = read.table(opt$stat_pairwise_plus, header=T, comment.char="")

print("readin one vs. all comparison statistics ... ")
stat4a_pre = read.table(opt$stat_one_vs_all_minus, header=T, comment.char="")
stat4b_pre = read.table(opt$stat_one_vs_all_plus, header=T, comment.char="")

#seek for SVs with significant P value after bonferroni correction
print("identify SVs with significant bonferroni corrected p-values ... ")
stat1 = add_bonferroni_correction(stat1_pre)
stat2a = add_bonferroni_correction(stat2a_pre)
stat2b = add_bonferroni_correction(stat2b_pre)
stat3a = add_bonferroni_correction(stat3a_pre)
stat3b = add_bonferroni_correction(stat3b_pre)
stat4a = add_bonferroni_per_batch(stat4a_pre)
stat4b = add_bonferroni_per_batch(stat4b_pre)

#collect failed SVID
fail_SVID_plus_vs_minus_info = data.frame(stat1[stat1$info=='fail',][,c(1,2)])
fail_SVID_pre_vs_post_minus_info = stat2a[stat2a$info=='fail',][,c(1,2)]
fail_SVID_pre_vs_post_plus_info = stat2b[stat2b$info=='fail',][,c(1,2)]
fail_SVID_pairwise_minus_info = stat3a[stat3a$info=='fail',][,c(1,2)]
fail_SVID_pairwise_plus_info = stat3b[stat3b$info=='fail',][,c(1,2)]

#label SVID with the reason for failure
fail_SVID_plus_vs_minus_info[,2] = 'PCR_AF_BIAS'
fail_SVID_pre_vs_post_minus_info[,2] = 'UNSTABLE_AF_ESTIMATE_PCRMINUS_pre_post'
fail_SVID_pre_vs_post_plus_info[,2]  = 'UNSTABLE_AF_ESTIMATE_PCRPLUS_pre_post'
fail_SVID_pairwise_minus_info[,2] = 'UNSTABLE_AF_ESTIMATE_PCRMINUS_pairwise'
fail_SVID_pairwise_plus_info[,2] = 'UNSTABLE_AF_ESTIMATE_PCRPLUS_pairwise'

#integrate and print significant results to output file
print("write SVs with significant bonferroni corrected p-values ... ")
info_out = merge(fail_SVID_plus_vs_minus_info,fail_SVID_pre_vs_post_minus_info,by='VID', all=T)
info_out = merge(info_out,fail_SVID_pre_vs_post_plus_info,by='VID', all=T)
info_out = merge(info_out,fail_SVID_pairwise_minus_info,by='VID', all=T)
info_out = merge(info_out,fail_SVID_pairwise_plus_info,by='VID', all=T)
info_out[,ncol(info_out)+1] = apply(info_out[,c(2,5,6)],1,function(x){paste(x[!is.na(x)],collapse = ',')})
colnames(info_out)[ncol(info_out)] = 'info_wo_pre_post'
info_out[,ncol(info_out)+1] = apply(info_out[,c(2:6)],1,function(x){paste(x[!is.na(x)],collapse = ',')})
colnames(info_out)[ncol(info_out)] = 'info_with_pre_post'
info_out_a = info_out[,c('VID','info_wo_pre_post')]
info_out_b = info_out[,c('VID','info_with_pre_post')]
info_out_a = info_out_a[info_out_a[,2]!="",]
info_out_b = info_out_b[info_out_b[,2]!="",]
write.table(info_out_a, opt$out_fail_info, quote=F, sep='\t',col.names=F, row.names=F)

#integrate SVID that failed batch effect and need filter column to be revised
print("write SVs with significant bias between batches ... ")
fail_SVID_plus_vs_minus_filter = stat1[stat1$filter=='fail',][,c(1,2)]
fail_SVID_pre_vs_post_minus_filter = stat2a[stat2a$filter=='fail',][,c(1,2)]
fail_SVID_pre_vs_post_plus_filter = stat2b[stat2b$filter=='fail',][,c(1,2)]
fail_SVID_pairwise_minus_filter = stat3a[stat3a$filter=='fail',][,c(1,2)]
fail_SVID_pairwise_plus_filter = stat3b[stat3b$filter=='fail',][,c(1,2)]
fail_filter = unique(c(fail_SVID_plus_vs_minus_filter, fail_SVID_pre_vs_post_minus_filter,fail_SVID_pairwise_minus_filter))
fail_SVID_plus_vs_minus_filter[,2] = 'PCR_AF_BIAS'
fail_SVID_pairwise_minus_filter[,2] = 'UNSTABLE_AF_ESTIMATE_PCRMINUS'
filter_out = merge(fail_SVID_plus_vs_minus_filter, fail_SVID_pairwise_minus_filter, by='VID', all=T)
filter_out[,4] = apply(filter_out[,c(2,3)],1,function(x){paste(x[!is.na(x)],collapse = ',')})
filter_out[,5] = apply(filter_out,1,function(x){strsplit(as.character(x[4]),',')[[1]][1]})

write.table(filter_out[,c(1,5)], opt$out_fail_filter, quote=F, sep='\t',col.names=F, row.names=F)

#integrate one vs. all comparison results to re-label format columns
print("write SVs with significant bias from one vs. all comparisons... ")
stat4_out = merge(stat4a, stat4b, by='VID', all=T)
stat4_out[is.na(stat4_out)] = 'pass'
stat4_out[,ncol(stat4_out)+1] = apply(stat4_out[,c(2:ncol(stat4_out))],1, function(x){paste(x[x!="pass"], collapse = ',')})
stat4_out_2 = stat4_out[,c(1,ncol(stat4_out))]
stat4_out_2 = stat4_out_2[stat4_out_2[,2]!="",]
write.table(stat4_out_2, opt$out_fail_format, quote=F, sep='\t',col.names=F, row.names=F)


