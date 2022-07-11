
library("optparse")

option_list = list(
        make_option(c( "--sample"), type="character", default=NULL, help="list of all samples", metavar="character"),
        make_option(c( "--loose_union"), type="character", default=NULL, help="loose union of all filters", metavar="character"),
        make_option(c( "--GQ_reacali"), type="character", default=NULL, help="GQrecalibrate filtering results", metavar="character"),
        make_option(c( "--boost_dyna"), type="character", default=NULL, help="boost with dynamic cutoffs filtering results", metavar="character"),
        make_option(c( "--boost_fix"), type="character", default=NULL, help="boost with fixed cutoffs filtering results", metavar="character"),
        make_option(c( "--minGQ10perc"), type="character", default=NULL, help="minGQ at 10%FDR filtering results", metavar="character"),
        make_option(c( "--boost_score"), type="character", default=NULL, help="per sample boost scores", metavar="character"),
        make_option(c( "--output"), type="character", default=NULL, help="name of bed file ", metavar="character"),
        make_option(c( "--SVID_anno"), type="character", default=NULL, help="annotation with all site-level training features", metavar="character"),
        make_option(c( "--SVID_coor"), type="character", default=NULL, help="annotation with SV coordinates", metavar="character"),
        make_option(c( "--SVID_NCR"), type="character", default=NULL, help="annotation with SVID and NCR", metavar="character"),
        make_option(c( "--ReCluster_anno"), type="character", default=NULL, help="annotation file between reclusterd and original SVID", metavar="character"),
        make_option(c( "--ReCluster_closest"), type="character", default=NULL, help="annotation file between reclusterd and original SVID", metavar="character"),
        make_option(c( "--ReCluster_sites"), type="character", default=NULL, help="coordiates of reclustered SVs", metavar="character")
 );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


data_polish<-function(dat){
  #colnames(dat)[c(18:23)] = c('RD_CN','RD_GQ','PE_GT','PE_GQ','SR_GT','SR_GQ')
  colnames_pos = match(c('GQ_reacali', 'boost_fix', 'boost_dyna', 'minGQ10.'), colnames(dat))
  for(i in colnames_pos){
    print(colnames(dat)[i])
    dat[,i] = as.character(dat[,i])
    dat[is.na(dat[,i]),][,i] ='no_info'
    dat[dat[,i]%in%c('0/1','1/1'),][,i]='info'
  }

  dat$EVIDENCE=as.character(dat$EVIDENCE)
  dat[is.na(dat$EVIDENCE),]$EVIDENCE = 'BAF'
  dat=dat[dat$SVTYPE%in%c('CPX','DEL','DUP','INS','INV'),]
  
  svtype_column_name = match(c('svtype','SVTYPE'), colnames(dat))
  colnames(dat)[svtype_column_name[1]] = 'SVTYPE'
  colnames(dat)[svtype_column_name[2]] = 'svtype'

  #dat = merge(dat, bs, by='SVID')
  if(nrow(dat[dat$svtype=="INS" & !dat$preds_builtin>-.6,])>0){
    dat[dat$svtype=="INS" & !dat$preds_builtin>-.6,]$boost_fix="no_info"  }
  if(nrow(dat[dat$svtype=="DUP" & !dat$preds_builtin>-.5,])>0){
    dat[dat$svtype=="DUP" & !dat$preds_builtin>-.5,]$boost_fix="no_info"     }
  dat$count_filters = apply(dat[,c('GQ_reacali', 'boost_fix', 'boost_dyna', 'minGQ10.')],1, function(x){length(x[x=="info"])})
  print(table(dat$count_filters))
  
  data_messy = dat[dat$svtype=="DEL" & dat$ALGORITHMS=="wham" & dat$EV=="SR" & dat$SVLEN<1000,]
  data_clean = dat[!dat$SVID%in%data_messy$SVID, ]
  data_clean = data_clean[data_clean$count_filters>0,]
  if(nrow(data_clean[is.na(data_clean$dnv_rate),])>0){
      data_clean[is.na(data_clean$dnv_rate),]$dnv_rate = 0
  }
  return(data_clean)
}

readin_anno<-function(annotation_file){
  anno=read.table(annotation_file, header=T)
  anno[,ncol(anno)+1] = anno$count_dnv / anno$count_children
  colnames(anno)[ncol(anno)] = 'dnv_rate'
  anno[,ncol(anno)+1] = 'over5Kb'
  colnames(anno)[ncol(anno)] = 'size_cate'
  anno[anno$SVLEN<5000,][,ncol(anno)] = '1to5Kb'
  anno[anno$SVLEN<1000,][,ncol(anno)] = '100to1Kb'
  anno[anno$SVLEN<100,][,ncol(anno)] = 'under100bp'
  return(anno)
}

readin_ncr<-function(SVID_ncr, SVID_anno, SVID_vs_closest_MEMBER){
  ncr=read.table(SVID_ncr, header=T)
  colnames(ncr)[1] = 'SVID'
  recluster1 = read.table(SVID_anno)
  recluster2 = read.table(SVID_vs_closest_MEMBER)
  ncr_uni = ncr[!ncr$SVID%in%recluster1[,1],]
  ncr_clu = ncr[ncr$SVID%in%recluster1[,1],]
  colnames(recluster2) = c('SVID2','SVID')
  ncr_clu=merge(ncr_clu, recluster2, by='SVID')
  ncr_clu=ncr_clu[,c(3,2)]
  colnames(ncr_clu)[1]='SVID'
  ncr2 = rbind(ncr_uni, ncr_clu)
  return(ncr2)
}

readin_cleanvcf<-function(cleanvcf){
  cleanvcf = read.table(cleanvcf, header=T)

  cleanvcf$PE_GQ = as.character(cleanvcf$PE_GQ)
  cleanvcf[cleanvcf$PE_GQ=="None",]$PE_GQ = -1
  cleanvcf$PE_GQ = as.integer(cleanvcf$PE_GQ)

  cleanvcf$SR_GQ = as.character(cleanvcf$SR_GQ)
  cleanvcf[cleanvcf$SR_GQ=="None",]$SR_GQ = -1
  cleanvcf$SR_GQ = as.integer(cleanvcf$SR_GQ)

  colnames(cleanvcf)[2]='cleanvcf'
  return(cleanvcf)
}

readin_sample_data<-function(loose_union, GQ_reacali, Boost_Dyna, Boost_Fix, minGQ10perc, boost_score){
  gq_recali = read.table(GQ_reacali, header=T)
  print(GQ_reacali)

  boost_dyna = read.table(Boost_Dyna, header=T)
  print(Boost_Dyna)

  boost_fix = read.table(Boost_Fix, header=T)
  print(Boost_Fix)

  minGQ_10 = read.table(minGQ10perc, header=T)
  print(minGQ10perc)


  overall_gt = readin_cleanvcf(loose_union)
  print(nrow(overall_gt))
  dat=merge(overall_gt, gq_recali[,c(1,2)], by='SVID', all=T)
  colnames(dat)[ncol(dat)] = 'GQ_reacali'
  print(nrow(dat))

  dat=merge(dat, boost_dyna[,c(1,2)], by='SVID', all=T)
  colnames(dat)[ncol(dat)] = 'boost_dyna'
  print(nrow(dat))

  dat=merge(dat, boost_fix[,c(1,2)], by='SVID', all=T)
  colnames(dat)[ncol(dat)] = 'boost_fix'
  print(nrow(dat))

  dat=merge(dat, minGQ_10[,c(1,2)], by='SVID', all=T)
  colnames(dat)[ncol(dat)] = 'minGQ10.'
  print(nrow(dat))

  bs=read.table(boost_score, header=T)
  dat=merge(dat, bs, by='SVID')
  print(nrow(dat))

  dat[,ncol(dat)+1] = apply(dat[,c("GQ_reacali","boost_dyna","boost_fix","minGQ10.")],1,function(x){length(x[!is.na(x)])})
  colnames(dat)[ncol(dat)] = 'count_filters'
  dat=dat[!is.na(dat$cleanvcf),]
  print(nrow(dat))

  return(dat)
}

#readin site level features: 
 #anno = readin_anno(opt$SVID_anno)
 #origin_coordinates = read.table(opt$SVID_coor, header=T)
 #ncr2 = readin_ncr(opt$SVID_NCR,opt$ReCluster_anno,opt$ReCluster_closest)
anno = readRDS(opt$SVID_anno)
origin_coordinates = readRDS(opt$SVID_coor)
ncr2 = readRDS(opt$SVID_NCR)

#readin sample level features: 
sample = read.table(opt$sample)
loose_union = read.table(opt$loose_union)
GQ_reacali = read.table(opt$GQ_reacali)
boost_dyna = read.table(opt$boost_dyna)
boost_fix = read.table(opt$boost_fix)
minGQ10perc = read.table(opt$minGQ10perc)
boost_score = read.table(opt$boost_score)
output = read.table(opt$output)
#filename_list is a table with sample name, cleanvcf, gq_recali, boost_dyna, boost_fix, minGQ_1, minGQ_5, minGQ_10, boost_score, output_filename in each column

for(i in 1:nrow(sample)){
  print(sample[i,1])
  SR_redun = read.table(opt$ReCluster_anno)
  colnames(SR_redun)= c('SVID','new_SVID')
  SR_recluster = read.table(opt$ReCluster_sites)
  colnames(SR_recluster)[4]='new_SVID'
  SR_redun=merge(SR_redun, SR_recluster[,c(1:4)], by='new_SVID')

  dat = readin_sample_data(as.character(loose_union[i,1]), as.character(GQ_reacali[i,1]), as.character(boost_dyna[i,1]), as.character(boost_fix[i,1]), as.character(minGQ10perc[i,1]), as.character(boost_score[i,1]))
  dat=merge(origin_coordinates[,c('SVID','X.chr','start','end')], dat, by='SVID')
  print(nrow(dat))

  redun_SV = dat[dat$SVID%in%SR_redun$SVID,]
  redun_SV = merge(redun_SV, SR_redun, by='SVID')
  kept_original_SVID=c()
  for(j in unique(redun_SV$new_SVID)){
    tmp = redun_SV[redun_SV$new_SVID==j,]
    tmp[,ncol(tmp)+1] = abs(tmp$V2 - tmp$start)
    kept = tmp[tmp[,ncol(tmp)]==min(tmp[,ncol(tmp)]),][1,]
    kept_original_SVID=c(kept_original_SVID, as.character(kept[1,1]))
  }

  redun_SV_V2 = redun_SV[redun_SV[,1]%in%kept_original_SVID,]
  dat_a = dat[!dat$SVID%in%redun_SV$SVID,c('SVID',"cleanvcf","EV","GQ","RD_CN","RD_GQ","PE_GT","PE_GQ","SR_GT","SR_GQ","GQ_reacali","boost_dyna","boost_fix","minGQ10.","preds_builtin","count_filters")]
  dat_b = redun_SV_V2[,c('new_SVID',"cleanvcf","EV","GQ","RD_CN","RD_GQ","PE_GT","PE_GQ","SR_GT","SR_GQ","GQ_reacali","boost_dyna","boost_fix","minGQ10.","preds_builtin","count_filters")]
  colnames(dat_b)[1]='SVID'

  #ryan_inte_SVID = read.table(as.character(filename_list[i,]$loose_union))
  #dat_a_2 = dat_a[dat_a$SVID%in%ryan_inte_SVID[,1],]
  #dat_b_2 = redun_SV_V2[redun_SV_V2$SVID%in%ryan_inte_SVID[,1],][,c('new_SVID',"cleanvcf","EV","GQ","RD_CN","RD_GQ","PE_GT","PE_GQ","SR_GT","SR_GQ","GQ_reacali","boost_dyna","boost_fix","minGQ1.","minGQ5.","minGQ10.","preds_builtin","count_filters")]
  #colnames(dat_b_2)[1]='SVID'

  dat_both_2=rbind(dat_a, dat_b)
  print(nrow(dat_both_2))
  out = merge(anno, dat_both_2, by='SVID')
  print(nrow(out))
  out = merge(out, ncr2, by='SVID')
  print(nrow(out))
  out=out[,c("SVID","SVTYPE","ALGORITHMS","EVIDENCE","NCR","SVLEN","svtype","cleanvcf", 
          "EV","GQ","PE_GQ", "PE_GT","RD_CN","RD_GQ", "SR_GQ","SR_GT","GQ_reacali", "boost_fix","boost_dyna","minGQ10.",
          "preds_builtin", "count_filters","size_cate","count_dnv", "count_children","count_non_ref","GC", "dnv_rate")]      
  out_polish = data_polish(out)

  write.table(out_polish, as.character(output[i,1]), quote=F, sep='\t', col.names=T, row.names=F)
}



