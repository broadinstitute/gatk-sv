#!/usr/bin/env Rscript
#script to train boost model and apply it on each sample

library(optparse)
library(data.table)
library(lightgbm)
library(pacman)

readin_site_features<-function(site_folder = './per_chr_anno/'){
  # read in the per-site traits across chromosomes
  list_file = list.files(site_folder)
  site_feature = read.table(paste(site_folder, list_file[1] ,sep = ''), header = T, comment.char = "")
  print(list_file[1])
  for(site_file in list_file[2:length(list_file)]){
    chr_feature = read.table(paste(site_folder, site_file, sep = ''), header = T, comment.char = "")
    print(site_file)
    site_feature = rbind(site_feature, chr_feature)
  }
  return(site_feature)
}

readin_sample_features<-function(sample_folder = './per_sample_anno/'){
  # read in training samples
  # sample_list = c('__hg00513__58c7ec','__hg00514_1__264528','__hg00512__766970','__hg00514__2b7cb6', '__hg00731__2f5d02','__hg00732__a5867c','__hg00733_1__388f2a','__hg00733__835b40','__na19238__c76ef5','__na19239__901548','__na19240__a01b77')
  file_list = list.files(sample_folder)
  sample_list = lapply(file_list,function(x){strsplit(as.character(x),'[.]')[[1]][1]})
  samp_feature = read.table(paste(sample_folder, file_list[1], sep=''), header = T, comment.char = "")
  #samp_feature[,ncol(samp_feature)+1] = sample_list[1]
  print(sample_list[1])
  #colnames(samp_feature)[ncol(samp_feature)] = 'sample'
  for(rec in c(2:length(file_list))){
    samp_feature_unit = read.table(paste(sample_folder, file_list[rec], sep = ''), header = T, comment.char = "")
    #samp_feature_unit[,ncol(samp_feature_unit)+1] = sample_list[rec]
    print(sample_list[rec])
    #colnames(samp_feature_unit)[ncol(samp_feature_unit)] = 'sample'
    samp_feature = rbind(samp_feature, samp_feature_unit)
  }
  return(samp_feature)
}

readin_training<-function(samp_feature, site_feature){
  training_set = merge(samp_feature, site_feature, by='SVID')
  #modify size category:
  training_set$size_cate = 's1_under250bp'
  training_set[training_set$SVLEN.x>250,]$size_cate = 's2_250bpto1Kb'
  training_set[training_set$SVLEN.x>1000,]$size_cate = 's3_1to5Kb'
  training_set[training_set$SVLEN.x>5000,]$size_cate = 's4_5to50Kb'
  training_set[training_set$SVLEN.x>50000,]$size_cate = 's5_over50Kb'
  
  for(i in c(1:ncol(training_set))){
    if(grepl('vs_',colnames(training_set)[i])){
      print(colnames(training_set)[i])
      if(nrow(training_set[training_set[,i]!='NO_OVR',])>0){
        training_set[,i]=as.character(training_set[,i])
        training_set[!is.na(training_set[,i]) & training_set[,i]!='NO_OVR',][,i]='OVR'
      } } }
  
  return(training_set)
}

readin_ssc_training<-function(){
  tp1 = read.table('per_sample_anno/gnomad_calls.ssc_samples.vs_array_calls.bed.gz')
  tp2 = read.table('per_sample_anno/gnomad_calls.ssc_samples.vs_exome_calls.bed.gz')
  fp = read.table('per_sample_anno/gnomad_calls.ssc_samples.original.vs_exome.wo_exome_calls.wo_array_calls.bed.gz')
}

readin_lg_cnv_training<-function(train_path = 'trainin_sets/array_training/'){
  ssc_hq_rare =   read.table(paste(train_path,'gnomad_calls.SSC_samples.lg_cnv.hq_rare.bed.gz', sep = ''), header = T, comment.char = "")
  ssc_hq_common = read.table(paste(train_path,'gnomad_calls.SSC_samples.lg_cnv.hq_common.bed.gz', sep = ''), header = T, comment.char = "")
  ssc_lq_rare =   read.table(paste(train_path,'gnomad_calls.SSC_samples.lg_cnv.lq_rare.bed.gz', sep = ''), header = T, comment.char = "")
  ssc_lq_common = read.table(paste(train_path,'gnomad_calls.SSC_samples.lg_cnv.lq_common.bed.gz', sep = ''), header = T, comment.char = "")
  ssc_hq_rare = ssc_hq_rare[,-ncol(ssc_hq_rare)]
  ssc_hq = rbind(ssc_hq_common, ssc_hq_rare)
  ssc_lq = rbind(ssc_lq_common, ssc_lq_rare)
  ssc_hq = merge(ssc_hq, site_feature, by='SVID')
  ssc_lq = merge(ssc_lq, site_feature, by='SVID')
  return(list(ssc_hq, ssc_lq))
}

train_filter_model_small<-function(data, feature_list = feature_list_oth, LD_SVID = LD_SVID){
  #set up a new column in the traning table to indicate number of supportive evidences
  #hq - high quality;  lq -  low quality
  data[,ncol(data)+1] = 0
  colnames(data)[ncol(data)] = 'piece_of_hq_evidences'
  data[,ncol(data)+1] = 0
  colnames(data)[ncol(data)] = 'piece_of_lq_evidences'
  
  #add external support to decide on truth and false training set
  # process by chromosome to avoid memory limits
  for(chr_name in unique(data$X.CHR)){
    print(chr_name)
    # add vapor support  
    data[data$X.CHR==chr_name & !is.na(data$VaPoR_GT) & data$VaPoR_GT%in%c('0/1','1/1'),]$piece_of_hq_evidences = 1+data[data$X.CHR==chr_name & !is.na(data$VaPoR_GT) & data$VaPoR_GT%in%c('0/1','1/1'),]$piece_of_hq_evidences 
    data[data$X.CHR==chr_name & !is.na(data$VaPoR_GT) & data$VaPoR_GT%in%c('0/0'),      ]$piece_of_lq_evidences = 1+data[data$X.CHR==chr_name & !is.na(data$VaPoR_GT) & data$VaPoR_GT%in%c('0/0'),      ]$piece_of_lq_evidences 
    # add PacBio support
    data[data$X.CHR==chr_name & data$vs_pacbio_ovr1a=='OVR',   ]$piece_of_hq_evidences = 1 + data[data$X.CHR==chr_name & data$vs_pacbio_ovr1a=='OVR',   ]$piece_of_hq_evidences
    data[data$X.CHR==chr_name & data$vs_pacbio_ovr1a=='NO_OVR',]$piece_of_lq_evidences = 1 + data[data$X.CHR==chr_name & data$vs_pacbio_ovr1a=='NO_OVR',]$piece_of_lq_evidences
    # add concordance support
    data[data$X.CHR==chr_name & data$concor_duplicates>1, ]$piece_of_hq_evidences = 1 + data[data$X.CHR==chr_name & data$concor_duplicates>1, ]$piece_of_hq_evidences
    data[data$X.CHR==chr_name & data$concor_duplicates==1,]$piece_of_lq_evidences = 1 + data[data$X.CHR==chr_name & data$concor_duplicates==1,]$piece_of_lq_evidences
    # add LD support:
    data[data$X.CHR==chr_name & data$SVID%in%LD_SVID[,1],]$piece_of_hq_evidences = 1+data[data$X.CHR==chr_name & data$SVID%in%LD_SVID[,1],]$piece_of_hq_evidences
  }
  
  #select high quality and low quality set: require two or more evidences
  hq = data[data$piece_of_hq_evidences>1,]
  lq = data[data$piece_of_lq_evidences>1,]
  print(c(nrow(hq), nrow(lq)))
  
  train_dat = rbind(hq,lq)[,feature_list]
  train_label=c(rep(1,nrow(hq)),rep(-1,nrow(lq)))
  train_set = cbind(train_dat,train_label)
  #train a multi-class light gbm model for deletions:
  set.seed(2)
  train_sample_all = sample(nrow(train_dat), nrow(train_dat))
  
  dtrain <- lgb.Dataset(data.matrix(train_dat), label = train_label)
  
  params <- list(   objective = "binary"   , metric = "auc"
                    , num_class = 1L   , learning_rate = 0.1    , min_data_in_leaf = 1L
                    , min_sum_hessian_in_leaf = 1.0, seed = 0 )
  
  model_builtin <- lgb.train(   params = params   , data = dtrain  , nrounds = 10L)
  
  tree_imp <- lgb.importance(model_builtin, percentage = TRUE)
  #lgb.plot.importance(tree_imp, top_n = 10L, measure = "Gain")
  
  return(model_builtin)
}

train_filter_model_medium<-function(data, feature_list = feature_list_cnv, LD_SVID = LD_SVID){
  #set up a new column in the traning table to indicate number of supportive evidences
  #hq - high quality;  lq -  low quality
  data[,ncol(data)+1] = 0
  colnames(data)[ncol(data)] = 'piece_of_hq_evidences'
  data[,ncol(data)+1] = 0
  colnames(data)[ncol(data)] = 'piece_of_lq_evidences'
  
  #add external support to decide on truth and false training set
  # process by chromosome to avoid memory limits
  for(chr_name in unique(data$X.CHR)){
    print(chr_name)
    # add vapor support
    data[data$X.CHR == chr_name & !is.na(data$VaPoR_GT) & data$VaPoR_GT%in%c('0/1','1/1'),]$piece_of_hq_evidences = 1+data[data$X.CHR == chr_name & !is.na(data$VaPoR_GT) & data$VaPoR_GT%in%c('0/1','1/1'),]$piece_of_hq_evidences 
    data[data$X.CHR == chr_name & !is.na(data$VaPoR_GT) & data$VaPoR_GT%in%c('0/0'),      ]$piece_of_lq_evidences = 1+data[data$X.CHR == chr_name & !is.na(data$VaPoR_GT) & data$VaPoR_GT%in%c('0/0'),      ]$piece_of_lq_evidences 
    # add PacBio support
    data[data$X.CHR == chr_name & data$vs_pacbio_ovr1a=='OVR',   ]$piece_of_hq_evidences = 1 + data[data$X.CHR == chr_name & data$vs_pacbio_ovr1a=='OVR',   ]$piece_of_hq_evidences
    data[data$X.CHR == chr_name & data$vs_pacbio_ovr1a=='NO_OVR',]$piece_of_lq_evidences = 1 + data[data$X.CHR == chr_name & data$vs_pacbio_ovr1a=='NO_OVR',]$piece_of_lq_evidences
    # add Bionano support
    data[data$X.CHR == chr_name & data$vs_bionano_ovr1a=='OVR',   ]$piece_of_hq_evidences = 1 + data[data$X.CHR == chr_name & data$vs_bionano_ovr1a=='OVR',   ]$piece_of_hq_evidences
    data[data$X.CHR == chr_name & data$vs_bionano_ovr1a=='NO_OVR',]$piece_of_lq_evidences = 1 + data[data$X.CHR == chr_name & data$vs_bionano_ovr1a=='NO_OVR',]$piece_of_lq_evidences
    # add concordance support
    data[data$X.CHR == chr_name & data$concor_duplicates>1, ]$piece_of_hq_evidences = 1 + data[data$X.CHR == chr_name & data$concor_duplicates>1, ]$piece_of_hq_evidences
    data[data$X.CHR == chr_name & data$concor_duplicates==1,]$piece_of_lq_evidences = 1 + data[data$X.CHR == chr_name & data$concor_duplicates==1,]$piece_of_lq_evidences
    # add LD support:
    data[data$X.CHR == chr_name & data$SVID%in%LD_SVID[,1],]$piece_of_hq_evidences = 1+data[data$X.CHR == chr_name & data$SVID%in%LD_SVID[,1],]$piece_of_hq_evidences
  }
  
  hq = data[data$piece_of_hq_evidences>1,]
  lq = data[data$piece_of_lq_evidences>2,]
  print(c(nrow(hq), nrow(lq)))
  
  train_dat = rbind(hq,lq)[,feature_list]
  train_label=c(rep(1,nrow(hq)),rep(-1,nrow(lq)))
  train_set = cbind(train_dat,train_label)
  #train a multi-class light gbm model for deletions:
  set.seed(2)
  train_sample_all = sample(nrow(train_dat), nrow(train_dat))
  
  dtrain <- lgb.Dataset(data.matrix(train_dat), label = train_label)
  
  params <- list(   objective = "binary"   , metric = "auc"
                    , num_class = 1L   , learning_rate = 0.1    , min_data_in_leaf = 1L
                    , min_sum_hessian_in_leaf = 1.0, seed=0 )
  
  model_builtin <- lgb.train(   params = params   , data = dtrain  , nrounds = 10L)
  
  tree_imp <- lgb.importance(model_builtin, percentage = TRUE)
  #lgb.plot.importance(tree_imp, top_n = 10L, measure = "Gain")
  
  return(model_builtin)
}

train_filter_model_large<-function(data, feature_list = feature_list_cnv, LD_SVID = LD_SVID, ssc_hq = ssc_hq, ssc_lq=ssc_lq){
  #set up a new column in the traning table to indicate number of supportive evidences
  #hq - high quality;  lq -  low quality
  data[,ncol(data)+1] = 0
  colnames(data)[ncol(data)] = 'piece_of_hq_evidences'
  data[,ncol(data)+1] = 0
  colnames(data)[ncol(data)] = 'piece_of_lq_evidences'
  # add vapor support  
  data[!is.na(data$VaPoR_GT) & data$VaPoR_GT%in%c('0/1','1/1'),]$piece_of_hq_evidences = 1+data[!is.na(data$VaPoR_GT) & data$VaPoR_GT%in%c('0/1','1/1'),]$piece_of_hq_evidences 
  data[!is.na(data$VaPoR_GT) & data$VaPoR_GT%in%c('0/0'),      ]$piece_of_lq_evidences = 1+data[!is.na(data$VaPoR_GT) & data$VaPoR_GT%in%c('0/0'),      ]$piece_of_lq_evidences 
  # add PacBio support
  data[data$vs_pacbio_ovr1a=='OVR',]$piece_of_hq_evidences = 1 + data[data$vs_pacbio_ovr1a=='OVR',]$piece_of_hq_evidences
  data[data$vs_pacbio_ovr1a=='NO_OVR',]$piece_of_lq_evidences = 1 + data[data$vs_pacbio_ovr1a=='NO_OVR',]$piece_of_lq_evidences
  # add Bionano support
  data[data$vs_bionano_ovr1a=='OVR',]$piece_of_hq_evidences = 1 + data[data$vs_bionano_ovr1a=='OVR',]$piece_of_hq_evidences
  data[data$vs_bionano_ovr1a=='NO_OVR',]$piece_of_lq_evidences = 1 + data[data$vs_bionano_ovr1a=='NO_OVR',]$piece_of_lq_evidences
  # add concordance support
  data[data$concor_duplicates>1,]$piece_of_hq_evidences = 1 + data[data$concor_duplicates>1,]$piece_of_hq_evidences
  data[data$concor_duplicates==1,]$piece_of_lq_evidences = 1 + data[data$concor_duplicates==1,]$piece_of_lq_evidences
  # add LD support:
  data[data$SVID%in%LD_SVID[,1],]$piece_of_hq_evidences = 1+data[data$SVID%in%LD_SVID[,1],]$piece_of_hq_evidences
  
  hq = data[data$piece_of_hq_evidences>1,]
  lq = data[data$piece_of_lq_evidences>2,]
  
  print(c(nrow(hq), nrow(lq)))
  train_dat = rbind(hq,lq)[,feature_list]
  train_label=c(rep(1,nrow(hq)),rep(-1,nrow(lq)))
  train_set = cbind(train_dat,train_label)
  
  train_dat_ssc = rbind(ssc_hq, ssc_lq)[,c(feature_list)]
  train_label_ssc = c(rep(1,nrow(ssc_hq)),rep(-1,nrow(ssc_lq)))
  train_set_ssc = cbind(train_dat_ssc,train_label_ssc)
  
  train_dat = rbind(train_dat, train_dat_ssc)
  train_label = c(train_label, train_label_ssc)
  #train a multi-class light gbm model for deletions:
  set.seed(2)
  train_sample_all = sample(nrow(train_dat), nrow(train_dat))
  
  dtrain <- lgb.Dataset(data.matrix(train_dat), label = train_label)
  
  params <- list(   objective = "binary"   , metric = "auc"
                    , num_class = 1L   , learning_rate = 0.1    , min_data_in_leaf = 1L
                    , min_sum_hessian_in_leaf = 1.0, seed=0 )
  
  model_builtin <- lgb.train(   params = params   , data = dtrain  , nrounds = 10L)
  
  tree_imp <- lgb.importance(model_builtin, percentage = TRUE)
  #lgb.plot.importance(tree_imp, top_n = 10L, measure = "Gain")
  
  return(model_builtin)
}

train_filter_model_INS<-function(data, feature_list = feature_list_oth, LD_SVID = LD_SVID){
  #set up a new column in the traning table to indicate number of supportive evidences
  #hq - high quality;  lq -  low quality
  data[,ncol(data)+1] = 0
  colnames(data)[ncol(data)] = 'piece_of_hq_evidences'
  data[,ncol(data)+1] = 0
  colnames(data)[ncol(data)] = 'piece_of_lq_evidences'
  
  #add external support to decide on truth and false training set
  # process by chromosome to avoid memory limits
  for(chr_name in unique(data$X.CHR)){
    print(chr_name)
    # add vapor support  
    data[data$X.CHR == chr_name & !is.na(data$VaPoR_GT) & data$VaPoR_GT%in%c('0/1','1/1'),]$piece_of_hq_evidences = 1+data[data$X.CHR == chr_name & !is.na(data$VaPoR_GT) & data$VaPoR_GT%in%c('0/1','1/1'),]$piece_of_hq_evidences 
    data[data$X.CHR == chr_name & !is.na(data$VaPoR_GT) & data$VaPoR_GT%in%c('0/0'),      ]$piece_of_lq_evidences = 1+data[data$X.CHR == chr_name & !is.na(data$VaPoR_GT) & data$VaPoR_GT%in%c('0/0'),      ]$piece_of_lq_evidences 
    # add PacBio support
    data[data$X.CHR == chr_name & data$vs_pacbio_ovr1a=='OVR',   ]$piece_of_hq_evidences = 1 + data[data$X.CHR == chr_name & data$vs_pacbio_ovr1a=='OVR',   ]$piece_of_hq_evidences
    data[data$X.CHR == chr_name & data$vs_pacbio_ovr1a=='NO_OVR',]$piece_of_lq_evidences = 1 + data[data$X.CHR == chr_name & data$vs_pacbio_ovr1a=='NO_OVR',]$piece_of_lq_evidences
    # add concordance support
    data[data$X.CHR == chr_name & data$concor_duplicates>1, ]$piece_of_hq_evidences = 1 + data[data$X.CHR == chr_name & data$concor_duplicates>1, ]$piece_of_hq_evidences
    data[data$X.CHR == chr_name & data$concor_duplicates==1,]$piece_of_lq_evidences = 1 + data[data$X.CHR == chr_name & data$concor_duplicates==1,]$piece_of_lq_evidences
    # add LD support:
    data[data$X.CHR == chr_name & data$SVID%in%LD_SVID[,1],]$piece_of_hq_evidences = 1+data[data$X.CHR == chr_name & data$SVID%in%LD_SVID[,1],]$piece_of_hq_evidences
  }
  
  hq = data[data$piece_of_hq_evidences>1,]
  lq = data[data$piece_of_lq_evidences>1,]
  print(c(nrow(hq), nrow(lq)))
  
  train_dat = rbind(hq,lq)[,feature_list]
  train_label=c(rep(1,nrow(hq)),rep(-1,nrow(lq)))
  train_set = cbind(train_dat,train_label)
  #train a multi-class light gbm model for deletions:
  set.seed(2)
  train_sample_all = sample(nrow(train_dat), nrow(train_dat))
  
  dtrain <- lgb.Dataset(data.matrix(train_dat), label = train_label)
  
  params <- list(   objective = "binary"   , metric = "auc"
                    , num_class = 1L   , learning_rate = 0.1    , min_data_in_leaf = 1L
                    , min_sum_hessian_in_leaf = 1.0, seed=0 )
  
  model_builtin <- lgb.train(   params = params   , data = dtrain  , nrounds = 10L)
  
  tree_imp <- lgb.importance(model_builtin, percentage = TRUE)
  #lgb.plot.importance(tree_imp, top_n = 10L, measure = "Gain")
  
  return(model_builtin)
}

train_del_dup_ins<-function(train, feature_list_cnv, feature_list_oth, LD_SVID, ssc_hq, ssc_lq){
  #feature_list_cnv=c("vs_raw_manta_ovr1a", "vs_raw_manta_ovr1b", "vs_raw_wham_ovr1a",  "vs_raw_wham_ovr1b",  "vs_raw_melt_ovr1a",  "vs_raw_melt_ovr1b","rd_median","rd_mean","rd_std","rd_median_le","rd_mean_le","rd_std_le","rd_median_ri","rd_mean_ri","rd_std_ri","PE_min","PE_max","SR_min","SR_max","size_cate","af_cate","ALGORITHMS","EVIDENCE","GC" ,"BothSideSupp","GT","GQ","RD_CN","RD_GQ","PE_GT","PE_GQ","SR_GT","SR_GQ","denovo_rate")
  #feature_list_oth=c("vs_raw_manta_ovr1a", "vs_raw_manta_ovr1b", "vs_raw_wham_ovr1a",  "vs_raw_wham_ovr1b",  "vs_raw_melt_ovr1a",  "vs_raw_melt_ovr1b","PE_min","PE_max","SR_min","SR_max","size_cate","af_cate","ALGORITHMS","EVIDENCE","GC" ,"BothSideSupp","GT","GQ","RD_CN","RD_GQ","PE_GT","PE_GQ","SR_GT","SR_GQ","denovo_rate")
  
  size = 's1_under250bp'
  train_recali_1a = train_filter_model_small(train[train$svtype=='DEL' & train$size_cate==size,],feature_list_oth,LD_SVID)
  train_recali_1b = train_filter_model_small(train[train$svtype=='DUP' & train$size_cate==size,],feature_list_oth,LD_SVID)
  
  size = 's2_250bpto1Kb'
  train_recali_2a = train_filter_model_small(train[train$svtype=='DEL' & train$size_cate==size,],feature_list_oth,LD_SVID)
  train_recali_2b = train_filter_model_small(train[train$svtype=='DUP' & train$size_cate==size,],feature_list_oth,LD_SVID)
  
  size = 's3_1to5Kb'
  train_recali_3a = train_filter_model_small(train[train$svtype=='DEL' & train$size_cate==size,],feature_list_oth,LD_SVID)
  train_recali_3b = train_filter_model_small(train[train$svtype=='DUP' & train$size_cate==size,],feature_list_oth,LD_SVID)
  
  size = 's4_5to50Kb'
  train_recali_4a = train_filter_model_medium(train[train$svtype=='DEL'& train$size_cate==size,], feature_list_cnv,LD_SVID)
  train_recali_4b = train_filter_model_medium(train[train$svtype=='DUP'& train$size_cate==size,], feature_list_cnv,LD_SVID)
  
  size = 's5_over50Kb'
  train_recali_5a = train_filter_model_large(train[train$svtype=='DEL'& train$size_cate==size,], feature_list_cnv,LD_SVID, ssc_hq, ssc_lq)
  train_recali_5b = train_filter_model_large(train[train$svtype=='DUP'& train$size_cate==size,], feature_list_cnv,LD_SVID, ssc_hq, ssc_lq)
  
  train_recali_1e = train_filter_model_INS(train[train$svtype=='INS',], feature_list_oth,LD_SVID)
  train_recali_1f = train_filter_model_INS(train[train$svtype%in%c('INS:ME:ALU')  ,], feature_list_oth,LD_SVID)
  train_recali_1g = train_filter_model_INS(train[train$svtype%in%c('INS:ME:LINE1','INS:ME:SVA','INS:ME') ,], feature_list_oth,LD_SVID)
  
  out=list(train_recali_1a,train_recali_1b,
           train_recali_2a,train_recali_2b,
           train_recali_3a,train_recali_3b,
           train_recali_4a,train_recali_4b,
           train_recali_5a,train_recali_5b,
           train_recali_1e,train_recali_1f,
           train_recali_1g)
  
  return(out)
}

test_del_dup_ins<-function(bgm_models, test, feature_list_cnv, feature_list_oth, cff=-0.8){
  models = bgm_models
  size='s1_under250bp'
  preds_builtin <- predict(models[[1]], data.matrix(test[test$svtype=='DEL' & test$size_cate==size ,][,feature_list_oth]), rawscore = TRUE, reshape = TRUE)
  train_recali_1a=cbind(test[test$svtype=='DEL' & test$size_cate==size ,], preds_builtin)
  
  preds_builtin <- predict(models[[2]], data.matrix(test[test$svtype=='DUP' & test$size_cate==size ,][,feature_list_oth]), rawscore = TRUE, reshape = TRUE)
  train_recali_1b=cbind(test[test$svtype=='DUP' & test$size_cate==size ,], preds_builtin)
  
  size='s2_250bpto1Kb'
  preds_builtin <- predict(models[[3]], data.matrix(test[test$svtype=='DEL' & test$size_cate==size ,][,feature_list_oth]), rawscore = TRUE, reshape = TRUE)
  train_recali_2a=cbind(test[test$svtype=='DEL' & test$size_cate==size ,], preds_builtin)
  
  preds_builtin <- predict(models[[4]], data.matrix(test[test$svtype=='DUP' & test$size_cate==size ,][,feature_list_oth]), rawscore = TRUE, reshape = TRUE)
  train_recali_2b=cbind(test[test$svtype=='DUP' & test$size_cate==size ,], preds_builtin)
  
  size='s3_1to5Kb'
  preds_builtin <- predict(models[[5]], data.matrix(test[test$svtype=='DEL' & test$size_cate==size ,][,feature_list_oth]), rawscore = TRUE, reshape = TRUE)
  train_recali_3a=cbind(test[test$svtype=='DEL' & test$size_cate==size ,], preds_builtin)
  
  preds_builtin <- predict(models[[6]], data.matrix(test[test$svtype=='DUP' & test$size_cate==size ,][,feature_list_oth]), rawscore = TRUE, reshape = TRUE)
  train_recali_3b=cbind(test[test$svtype=='DUP' & test$size_cate==size ,], preds_builtin)
  
  size = 's4_5to50Kb'
  preds_builtin <- predict(models[[7]], data.matrix(test[test$svtype=='DEL' & test$size_cate==size ,][,feature_list_cnv]), rawscore = TRUE, reshape = TRUE)
  train_recali_4a=cbind(test[test$svtype=='DEL' & test$size_cate==size ,], preds_builtin)
  
  preds_builtin <- predict(models[[8]], data.matrix(test[test$svtype=='DUP' & test$size_cate==size ,][,feature_list_cnv]), rawscore = TRUE, reshape = TRUE)
  train_recali_4b=cbind(test[test$svtype=='DUP' & test$size_cate==size ,], preds_builtin)
  
  size = 's5_over50Kb'
  preds_builtin <- predict(models[[9]], data.matrix(test[test$svtype=='DEL' & test$size_cate==size ,][,feature_list_cnv]), rawscore = TRUE, reshape = TRUE)
  train_recali_5a=cbind(test[test$svtype=='DEL' & test$size_cate==size ,], preds_builtin)
  
  preds_builtin <- predict(models[[10]], data.matrix(test[test$svtype=='DUP' & test$size_cate==size ,][,feature_list_cnv]), rawscore = TRUE, reshape = TRUE)
  train_recali_5b=cbind(test[test$svtype=='DUP' & test$size_cate==size ,], preds_builtin)
  
  
  preds_builtin <- predict(models[[11]], data.matrix(test[test$svtype=='INS'   ,][,feature_list_oth]), rawscore = TRUE, reshape = TRUE)
  train_recali_1e=cbind(test[test$svtype=='INS'   ,], preds_builtin)
  
  preds_builtin <- predict(models[[12]], data.matrix(test[test$svtype%in%c('INS:ME:ALU') ,][,feature_list_oth]), rawscore = TRUE, reshape = TRUE)
  train_recali_1f=cbind(test[test$svtype%in%c('INS:ME:ALU') ,], preds_builtin)
  
  preds_builtin <- predict(models[[13]], data.matrix(test[test$svtype%in%c('INS:ME:LINE1','INS:ME:SVA','INS:ME') ,][,feature_list_oth]), rawscore = TRUE, reshape = TRUE)
  train_recali_1g=cbind(test[test$svtype%in%c('INS:ME:LINE1','INS:ME:SVA','INS:ME') ,], preds_builtin)
  
  train_recali_1h = test[test$svtype%in%c("INV",'CPX','CTX','CNV','BND') ,]
  train_recali_1h[,ncol(train_recali_1h)+1]=1
  colnames(train_recali_1h)[ncol(train_recali_1h)]='preds_builtin'
  
  out=rbind(train_recali_1a,train_recali_1b,
            train_recali_2a,train_recali_2b,
            train_recali_3a,train_recali_3b,
            train_recali_4a,train_recali_4b,
            train_recali_5a,train_recali_5b,
            train_recali_1e,train_recali_1f,
            train_recali_1g,train_recali_1h)
  
  out[,ncol(out)+1]=0
  out[out$preds_builtin>cff,][,ncol(out)]=1
  colnames(out)[ncol(out)]='preds_out'
  
  return(out)
}

plot_roc_by_svtype<-function(roc_table_del, roc_table_dup, roc_table_ins, roc_table_mei, title_name){
  plot(c(0,1), c(0,1), frame.plot = F, type = 'n', xlab = 'False Positive Rate', ylab = 'True Positive Rate', main = title_name)
  for(i in c(0:5)){abline(h=i/5, col='grey', lty=1)}
  for(i in c(0:5)){abline(h=i/5-.1, col='grey', lty=2)}
  for(i in c(0:5)){abline(v=i/5, col='grey', lty=1)}
  for(i in c(0:5)){abline(v=i/5-.1, col='grey', lty=2)}
  lines(roc_table_del$fp, roc_table_del$tp, pch=20, type = 'b', lty=1, lwd=2, col=as.character(col_table[col_table[,1]=='DEL',2]))
  lines(roc_table_dup$fp, roc_table_dup$tp, pch=20, type = 'b', lty=1, lwd=2, col=as.character(col_table[col_table[,1]=='DUP',2]))
  lines(roc_table_ins$fp, roc_table_ins$tp, pch=20, type = 'b', lty=1, lwd=2, col=as.character(col_table[col_table[,1]=='INS',2]))
  lines(roc_table_mei$fp, roc_table_mei$tp, pch=20, type = 'b', lty=1, lwd=2, col=as.character(col_table[col_table[,1]=='MEI',2]))
  legend('bottomright', c('DEL','DUP','INS','MEI'), col=c(as.character(col_table[col_table[,1]=='DEL',2]),as.character(col_table[col_table[,1]=='DUP',2]),as.character(col_table[col_table[,1]=='INS',2]),as.character(col_table[col_table[,1]=='MEI',2])), bty = 'n', pch = 20, lwd=2, lty = 1, bg = 'white')
}

calcu_ROC_by_BS_vs_hgsv<-function(test_results){
  cff_list = seq(-4,4,by=.1)  
  out=data.frame('cff'=0,'tp'=0,'fp'=0)
  row_num = 0
  for(cff in cff_list){
    if(cff > min(test_results$preds_builtin) & cff<max(test_results$preds_builtin)){
      test_results$preds_out = 0
      test_results[test_results$preds_builtin>cff,]$preds_out=1  
      
      tp = nrow(test_results[test_results$vs_hgsv_ovr1a=='OVR' & test_results$preds_out ==1,])
      fp = nrow(test_results[test_results$vs_hgsv_ovr1a=='NO_OVR' & test_results$preds_out ==1,])
      tn = nrow(test_results[test_results$vs_hgsv_ovr1a=='NO_OVR' & test_results$preds_out ==0,])
      fn = nrow(test_results[test_results$vs_hgsv_ovr1a=='OVR' & test_results$preds_out ==0,])
      
      tpr = tp/(tp+fn)
      fpr = fp/(fp+tn)
      
      row_num=row_num+1
      out[row_num,1] = cff
      out[row_num,2] = tpr
      out[row_num,3] = fpr
    }
  }
  return(out)
}

plot_roc_by_svtype_all_sizes_vs_hgsv<-function(test_results,title_name='all size ranges'){
  roc_table_del = calcu_ROC_by_BS_vs_hgsv(test_results[test_results$svtype.x=='DEL',])
  roc_table_dup = calcu_ROC_by_BS_vs_hgsv(test_results[test_results$svtype.x=='DUP',])
  roc_table_ins = calcu_ROC_by_BS_vs_hgsv(test_results[test_results$svtype=='INS',])
  roc_table_mei = calcu_ROC_by_BS_vs_hgsv(test_results[test_results$svtype%in%c('INS:ME','INS:ME:ALU','INS:ME:LINE1','INS:ME:SVA'),])
  plot_roc_by_svtype(roc_table_del, roc_table_dup, roc_table_ins, roc_table_mei, title_name)
}

plot_roc_by_svtype_and_size_vs_hgsv<-function(test_results){
  roc_table_del_size1 = calcu_ROC_by_BS_vs_hgsv(test_results[test_results$svtype.x=='DEL' & test_results$size_cate%in%c('s1_under250bp','s2_250bpto1Kb'),])
  roc_table_dup_size1 = calcu_ROC_by_BS_vs_hgsv(test_results[test_results$svtype.x=='DUP' & test_results$size_cate%in%c('s1_under250bp','s2_250bpto1Kb'),])
  roc_table_ins_size1 = calcu_ROC_by_BS_vs_hgsv(test_results[test_results$svtype.x=='INS' & test_results$size_cate%in%c('s1_under250bp','s2_250bpto1Kb'),])
  
  roc_table_del_size2 = calcu_ROC_by_BS_vs_hgsv(test_results[test_results$svtype.x=='DEL' & test_results$size_cate%in%c('s3_1to5Kb'),])
  roc_table_dup_size2 = calcu_ROC_by_BS_vs_hgsv(test_results[test_results$svtype.x=='DUP' & test_results$size_cate%in%c('s3_1to5Kb'),])
  roc_table_ins_size2 = calcu_ROC_by_BS_vs_hgsv(test_results[test_results$svtype.x=='INS' & test_results$size_cate%in%c('s3_1to5Kb'),])
  
  roc_table_del_size3 = calcu_ROC_by_BS_vs_hgsv(test_results[test_results$svtype.x=='DEL' & test_results$size_cate%in%c('s4_5to50Kb'),])
  roc_table_dup_size3 = calcu_ROC_by_BS_vs_hgsv(test_results[test_results$svtype.x=='DUP' & test_results$size_cate%in%c('s4_5to50Kb'),])
  roc_table_ins_size3 = calcu_ROC_by_BS_vs_hgsv(test_results[test_results$svtype.x=='INS' & test_results$size_cate%in%c('s4_5to50Kb'),])
  
  roc_table_del_size4 = calcu_ROC_by_BS_vs_hgsv(test_results[test_results$svtype.x=='DEL' & test_results$size_cate%in%c('s5_over50Kb'),])
  roc_table_dup_size4 = calcu_ROC_by_BS_vs_hgsv(test_results[test_results$svtype.x=='DUP' & test_results$size_cate%in%c('s5_over50Kb'),])
  roc_table_ins_size4 = calcu_ROC_by_BS_vs_hgsv(test_results[test_results$svtype.x=='INS' & test_results$size_cate%in%c('s5_over50Kb'),])
  
  par(mfrow=c(2,2))
  plot_roc_by_svtype(roc_table_del_size1, roc_table_dup_size1, roc_table_ins_size1, 'under 1Kb')
  plot_roc_by_svtype(roc_table_del_size2, roc_table_dup_size2, roc_table_ins_size2, '1-5Kb')
  plot_roc_by_svtype(roc_table_del_size3, roc_table_dup_size3, roc_table_ins_size3, '5-50Kb')
  plot_roc_by_svtype(roc_table_del_size4, roc_table_dup_size4, roc_table_ins_size4, '>50Kb')
}

calcu_ROC_by_BS_vs_PacBio<-function(test_results){
  cff_list = seq(-4,4,by=.1)  
  out=data.frame('cff'=0,'tp'=0,'fp'=0)
  row_num = 0
  for(cff in cff_list){
    if(cff > min(test_results$preds_builtin) & cff<max(test_results$preds_builtin)){
      test_results$preds_out = 0
      test_results[test_results$preds_builtin>cff,]$preds_out=1  
      
      tp = nrow(test_results[test_results$vs_pacbio_ovr1a=='OVR' & test_results$preds_out ==1,])
      fp = nrow(test_results[test_results$vs_pacbio_ovr1a=='NO_OVR' & test_results$preds_out ==1,])
      tn = nrow(test_results[test_results$vs_pacbio_ovr1a=='NO_OVR' & test_results$preds_out ==0,])
      fn = nrow(test_results[test_results$vs_pacbio_ovr1a=='OVR' & test_results$preds_out ==0,])
      
      tpr = tp/(tp+fn)
      fpr = fp/(fp+tn)
      
      row_num=row_num+1
      out[row_num,1] = cff
      out[row_num,2] = tpr
      out[row_num,3] = fpr
    }
  }
  return(out)
}

plot_roc_by_svtype_all_sizes_vs_PacBio<-function(test_results,title_name='all size ranges'){
  roc_table_del = calcu_ROC_by_BS_vs_PacBio(test_results[test_results$svtype.x=='DEL',])
  roc_table_dup = calcu_ROC_by_BS_vs_PacBio(test_results[test_results$svtype.x=='DUP',])
  roc_table_ins = calcu_ROC_by_BS_vs_PacBio(test_results[test_results$svtype=='INS',])
  roc_table_mei = calcu_ROC_by_BS_vs_PacBio(test_results[test_results$svtype%in%c('INS:ME','INS:ME:ALU','INS:ME:LINE1','INS:ME:SVA'),])
  plot_roc_by_svtype(roc_table_del, roc_table_dup, roc_table_ins, roc_table_mei, title_name)
}

plot_roc_by_svtype_and_size_vs_PacBio<-function(test_results){
  roc_table_del_size1 = calcu_ROC_by_BS_vs_PacBio(test_results[test_results$svtype.x=='DEL' & test_results$size_cate%in%c('s1_under250bp','s2_250bpto1Kb'),])
  roc_table_dup_size1 = calcu_ROC_by_BS_vs_PacBio(test_results[test_results$svtype.x=='DUP' & test_results$size_cate%in%c('s1_under250bp','s2_250bpto1Kb'),])
  roc_table_ins_size1 = calcu_ROC_by_BS_vs_PacBio(test_results[test_results$svtype.x=='INS' & test_results$size_cate%in%c('s1_under250bp','s2_250bpto1Kb'),])
  
  roc_table_del_size2 = calcu_ROC_by_BS_vs_PacBio(test_results[test_results$svtype.x=='DEL' & test_results$size_cate%in%c('s3_1to5Kb'),])
  roc_table_dup_size2 = calcu_ROC_by_BS_vs_PacBio(test_results[test_results$svtype.x=='DUP' & test_results$size_cate%in%c('s3_1to5Kb'),])
  roc_table_ins_size2 = calcu_ROC_by_BS_vs_PacBio(test_results[test_results$svtype.x=='INS' & test_results$size_cate%in%c('s3_1to5Kb'),])
  
  roc_table_del_size3 = calcu_ROC_by_BS_vs_PacBio(test_results[test_results$svtype.x=='DEL' & test_results$size_cate%in%c('s4_5to50Kb'),])
  roc_table_dup_size3 = calcu_ROC_by_BS_vs_PacBio(test_results[test_results$svtype.x=='DUP' & test_results$size_cate%in%c('s4_5to50Kb'),])
  roc_table_ins_size3 = calcu_ROC_by_BS_vs_PacBio(test_results[test_results$svtype.x=='INS' & test_results$size_cate%in%c('s4_5to50Kb'),])
  
  roc_table_del_size4 = calcu_ROC_by_BS_vs_PacBio(test_results[test_results$svtype.x=='DEL' & test_results$size_cate%in%c('s5_over50Kb'),])
  roc_table_dup_size4 = calcu_ROC_by_BS_vs_PacBio(test_results[test_results$svtype.x=='DUP' & test_results$size_cate%in%c('s5_over50Kb'),])
  roc_table_ins_size4 = calcu_ROC_by_BS_vs_PacBio(test_results[test_results$svtype.x=='INS' & test_results$size_cate%in%c('s5_over50Kb'),])
  
  par(mfrow=c(2,2))
  plot_roc_by_svtype(roc_table_del_size1, roc_table_dup_size1, roc_table_ins_size1, 'under 1Kb')
  plot_roc_by_svtype(roc_table_del_size2, roc_table_dup_size2, roc_table_ins_size2, '1-5Kb')
  plot_roc_by_svtype(roc_table_del_size3, roc_table_dup_size3, roc_table_ins_size3, '5-50Kb')
  plot_roc_by_svtype(roc_table_del_size4, roc_table_dup_size4, roc_table_ins_size4, '>50Kb')
}

plot_upset_by_svtype_all_sizes<-function(test_results, boost_cff){
  bar0 = nrow(test_results[test_results$vs_hgsv_ovr1a=='NO_OVR' & test_results$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/0") & test_results$concor_duplicates==1,])
  bar1 = nrow(test_results[test_results$vs_hgsv_ovr1a=='OVR' & test_results$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/0") & test_results$concor_duplicates==1,])
  bar2 = nrow(test_results[test_results$vs_hgsv_ovr1a=='NO_OVR' & test_results$vs_pacbio_ovr1a=="OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/0") & test_results$concor_duplicates==1,])
  bar3 = nrow(test_results[test_results$vs_hgsv_ovr1a=='NO_OVR' & test_results$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/1", "1/1") & test_results$concor_duplicates==1,])
  bar4 = nrow(test_results[test_results$vs_hgsv_ovr1a=='NO_OVR' & test_results$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/0") & test_results$concor_duplicates>1,])
  
  bar12 = nrow(test_results[test_results$vs_hgsv_ovr1a=='OVR' & test_results$vs_pacbio_ovr1a=="OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/0") & test_results$concor_duplicates==1,])
  bar13 = nrow(test_results[test_results$vs_hgsv_ovr1a=='OVR' & test_results$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/1", "1/1") & test_results$concor_duplicates==1,])
  bar14 = nrow(test_results[test_results$vs_hgsv_ovr1a=='OVR' & test_results$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/0") & test_results$concor_duplicates>1,])
  bar23 = nrow(test_results[test_results$vs_hgsv_ovr1a=='NO_OVR' & test_results$vs_pacbio_ovr1a=="OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/1", "1/1") & test_results$concor_duplicates==1,])
  bar24 = nrow(test_results[test_results$vs_hgsv_ovr1a=='NO_OVR' & test_results$vs_pacbio_ovr1a=="OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/0") & test_results$concor_duplicates>1,])
  bar34 = nrow(test_results[test_results$vs_hgsv_ovr1a=='NO_OVR' & test_results$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/1", "1/1") & test_results$concor_duplicates>1,])
  
  bar123 = nrow(test_results[test_results$vs_hgsv_ovr1a=='OVR' & test_results$vs_pacbio_ovr1a=="OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/1", "1/1") & test_results$concor_duplicates==1,])
  bar124 = nrow(test_results[test_results$vs_hgsv_ovr1a=='OVR' & test_results$vs_pacbio_ovr1a=="OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/0") & test_results$concor_duplicates>1,])
  bar134 = nrow(test_results[test_results$vs_hgsv_ovr1a=='OVR' & test_results$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/1", "1/1") & test_results$concor_duplicates>1,])
  bar234 = nrow(test_results[test_results$vs_hgsv_ovr1a=='NO_OVR' & test_results$vs_pacbio_ovr1a=="OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/1", "1/1") & test_results$concor_duplicates>1,])
  bar1234 = nrow(test_results[test_results$vs_hgsv_ovr1a=='OVR' & test_results$vs_pacbio_ovr1a=="OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/1", "1/1") & test_results$concor_duplicates>1,])
  
  cff = boost_cff
  test_results$preds_out=0
  test_results[test_results$preds_builtin>cff,]$preds_out=1
  bar0b = nrow(test_results[test_results$preds_out==1 & test_results$vs_hgsv_ovr1a=='NO_OVR' & test_results$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/0") & test_results$concor_duplicates==1,])
  bar1b = nrow(test_results[test_results$preds_out==1 & test_results$vs_hgsv_ovr1a=='OVR' & test_results$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/0") & test_results$concor_duplicates==1,])
  bar2b = nrow(test_results[test_results$preds_out==1 & test_results$vs_hgsv_ovr1a=='NO_OVR' & test_results$vs_pacbio_ovr1a=="OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/0") & test_results$concor_duplicates==1,])
  bar3b = nrow(test_results[test_results$preds_out==1 & test_results$vs_hgsv_ovr1a=='NO_OVR' & test_results$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/1", "1/1") & test_results$concor_duplicates==1,])
  bar4b = nrow(test_results[test_results$preds_out==1 & test_results$vs_hgsv_ovr1a=='NO_OVR' & test_results$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/0") & test_results$concor_duplicates>1,])
  
  bar12b = nrow(test_results[test_results$preds_out==1 & test_results$vs_hgsv_ovr1a=='OVR' & test_results$vs_pacbio_ovr1a=="OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/0") & test_results$concor_duplicates==1,])
  bar13b = nrow(test_results[test_results$preds_out==1 & test_results$vs_hgsv_ovr1a=='OVR' & test_results$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/1", "1/1") & test_results$concor_duplicates==1,])
  bar14b = nrow(test_results[test_results$preds_out==1 & test_results$vs_hgsv_ovr1a=='OVR' & test_results$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/0") & test_results$concor_duplicates>1,])
  bar23b = nrow(test_results[test_results$preds_out==1 & test_results$vs_hgsv_ovr1a=='NO_OVR' & test_results$vs_pacbio_ovr1a=="OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/1", "1/1") & test_results$concor_duplicates==1,])
  bar24b = nrow(test_results[test_results$preds_out==1 & test_results$vs_hgsv_ovr1a=='NO_OVR' & test_results$vs_pacbio_ovr1a=="OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/0") & test_results$concor_duplicates>1,])
  bar34b = nrow(test_results[test_results$preds_out==1 & test_results$vs_hgsv_ovr1a=='NO_OVR' & test_results$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/1", "1/1") & test_results$concor_duplicates>1,])
  
  bar123b = nrow(test_results[test_results$preds_out==1 & test_results$vs_hgsv_ovr1a=='OVR' & test_results$vs_pacbio_ovr1a=="OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/1", "1/1") & test_results$concor_duplicates==1,])
  bar124b = nrow(test_results[test_results$preds_out==1 & test_results$vs_hgsv_ovr1a=='OVR' & test_results$vs_pacbio_ovr1a=="OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/0") & test_results$concor_duplicates>1,])
  bar134b = nrow(test_results[test_results$preds_out==1 & test_results$vs_hgsv_ovr1a=='OVR' & test_results$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/1", "1/1") & test_results$concor_duplicates>1,])
  bar234b = nrow(test_results[test_results$preds_out==1 & test_results$vs_hgsv_ovr1a=='NO_OVR' & test_results$vs_pacbio_ovr1a=="OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/1", "1/1") & test_results$concor_duplicates>1,])
  bar1234b = nrow(test_results[test_results$preds_out==1 & test_results$vs_hgsv_ovr1a=='OVR' & test_results$vs_pacbio_ovr1a=="OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/1", "1/1") & test_results$concor_duplicates>1,])
  
  bar_list = c(bar0, bar1,bar2,bar3,bar4, bar12, bar13, bar14, bar23, bar24, bar34, bar123, bar124, bar134, bar234, bar1234)
  bar_list_b = c(bar0b, bar1b,bar2b,bar3b,bar4b, bar12b, bar13b, bar14b, bar23b, bar24b, bar34b, bar123b, bar124b, bar134b, bar234b, bar1234b)
  y_range = range(bar1,bar2,bar3,bar4, bar12, bar13, bar14, bar23, bar24, bar34, bar123, bar124, bar134, bar234, bar1234)
  
  par(mar=c(0,6,1,2))
  par(fig=c(0,1,.65,1))
  plot(c(0,length(bar_list_b)), range(bar_list), frame.plot = F, type = 'n', xlab = '', ylab = 'Count of SVs / Genome', las=2, xaxt='n', mgp=c(4,1,0), main = 'SVs per sample Before / After Boost Model')
  bar_width = .4
  x_pos = .5
  bar_colors = c('grey', 'black')
  for(i in c(1:length(bar_list))){
    rect(x_pos - bar_width, 0,  x_pos + bar_width, bar_list[i], col = bar_colors[1] )
    rect(x_pos - bar_width, 0,  x_pos + bar_width, bar_list_b[i], col = bar_colors[2] )
    x_pos = x_pos + 1
  }
  legend('topright', c('From Module05', 'After Boost Filter'), title = 'Count of SVs', bty = 'n', col = bar_colors[c(1,2)], pch = 15, pt.cex = 2)
  
  par(fig=c(0,1,.3,.65), new=T)
  plot(c(0,length(bar_list_b)), range(0,1.1), frame.plot = F, type = 'n', xlab = '', ylab = 'Proportion of SVs \n / Genome', las=2, xaxt='n', mgp=c(3,1,0))
  bar_width = .4
  x_pos = .5
  bar_colors = c('grey', 'black')
  for(i in c(1:length(bar_list))){
    rect(x_pos - bar_width, 0,  x_pos + bar_width, 1, col = bar_colors[1] )
    rect(x_pos - bar_width, 0,  x_pos + bar_width, bar_list_b[i]/bar_list[i], col = bar_colors[2] )
    #text(x_pos, 1.05, bar_list[i], srt=45, cex = .6)
    x_pos = x_pos + 1
  }
  
  par(fig=c(0,1,0,.3), new=T)
  plot(c(0,length(bar_list_b)), c(0,4), frame.plot = F, type = 'n', xlab = '', ylab = '', xaxt='n', yaxt='n')
  axis(2,c(1:4)-.5, labels = c('Replicate \n Concordance','VaPoR \n Validated','PacBio \n Overlap','Illumina \n Overlap'), las=2, cex.axis=.8)
  for(i in c(1:length(bar_list_b))-.5){
    for(j in c(1:4)-.5){
      points(i,j, col='grey', pch=20, cex = 2)
    }
  }
  for(i in c(1:4)+.5){  points(i, 5-i, col='black', pch=20, cex = 2)}
  start_pos = 5.5
  for(i in c(3,4)-.5){  points(start_pos,i, col='black', pch=20, cex = 2)}
  for(i in c(2,4)-.5){  points(start_pos+1,i, col='black', pch=20, cex = 2)}
  for(i in c(1,4)-.5){  points(start_pos+2,i, col='black', pch=20, cex = 2)}
  for(i in c(2,3)-.5){  points(start_pos+3,i, col='black', pch=20, cex = 2)}
  for(i in c(1,3)-.5){  points(start_pos+4,i, col='black', pch=20, cex = 2)}
  for(i in c(1,2)-.5){  points(start_pos+5,i, col='black', pch=20, cex = 2)}
  for(i in c(2,3,4)-.5){  points(start_pos+6,i, col='black', pch=20, cex = 2)}
  for(i in c(1,3,4)-.5){  points(start_pos+7,i, col='black', pch=20, cex = 2)}
  for(i in c(1,2,4)-.5){  points(start_pos+8,i, col='black', pch=20, cex = 2)}
  for(i in c(1,2,3)-.5){  points(start_pos+9,i, col='black', pch=20, cex = 2)}
  for(i in c(1:4)-.5)  {  points(start_pos+10,i, col='black', pch=20, cex = 2)}
  lines(c(start_pos,start_pos),c(3:4)-.5)
  lines(c(start_pos+1,start_pos+1),c(2,4)-.5)
  lines(c(start_pos+2,start_pos+2),c(1,4)-.5)
  lines(c(start_pos+3,start_pos+3),c(2,3)-.5)
  lines(c(start_pos+4,start_pos+4),c(1,3)-.5)
  lines(c(start_pos+5,start_pos+5),c(1,2)-.5)
  lines(c(start_pos+6,start_pos+6),c(2,4)-.5)
  lines(c(start_pos+7,start_pos+7),c(1,4)-.5)
  lines(c(start_pos+8,start_pos+8),c(1,4)-.5)
  lines(c(start_pos+9,start_pos+9),c(1,3)-.5)
  lines(c(start_pos+10,start_pos+10),c(1,4)-.5)
  
}

plot_SVCounts_by_bs_score<-function(test_results, color, plot_type, title_name, digits=1000){
  plot(c(-4,2), c(0,as.integer(nrow(test_results)/digits)*digits+digits), frame.plot = F, type = 'n', xlab = 'Boost Score Cutoff', ylab = 'Count of SVs', yaxt='n')
  for(i in c(0:(as.integer(nrow(test_results)/digits)+1))){
    abline(h=i*digits, col='grey', lty=1)
    abline(h=i*digits-digits/2, col='grey', lty=2)
  }
  axis(2, c(0:(as.integer(nrow(test_results)/digits)+1))*digits, las=2)
  
  bs_cff_list = seq(-4,2, by = .2)
  bs_cff_table = data.frame('bs_cff' = 0, 'count_SVs' = 0)
  row_count = 0
  for(i in bs_cff_list){
    row_count=row_count+1
    bs_cff_table[row_count,1] = i
    bs_cff_table[row_count,2] = nrow(test_results[test_results$preds_builtin>i,])
  }
  lines(bs_cff_table[,1],bs_cff_table[,2], type = 'b', pch = plot_type, lty = 1, lwd=2, col = color)
  legend('topright', title_name, col = color, pch = plot_type, lty = 1, lwd = 2, bg = 'white', bty = 'n')
}

plot_SVCounts_by_bs_score_by_svtype<-function(test_results){
  par(mfrow=c(3,2))
  par(mar=c(4,4,0,0))
  plot_SVCounts_by_bs_score(test_results, 'black', 20, 'all SVs', 10000)
  plot_SVCounts_by_bs_score(test_results[test_results$svtype.x=='DEL',], as.character(col_table[col_table[,1]=='DEL',2]), 20, 'Deletions',10000)
  plot_SVCounts_by_bs_score(test_results[test_results$svtype.x=='DUP',], as.character(col_table[col_table[,1]=='DUP',2]), 20, 'Duplications')
  plot_SVCounts_by_bs_score(test_results[test_results$svtype.x=='INV',], as.character(col_table[col_table[,1]=='INV',2]), 20, 'Inversions')
  plot_SVCounts_by_bs_score(test_results[test_results$svtype.x=='INS',], as.character(col_table[col_table[,1]=='INS',2]), 20, 'Insertions')
  plot_SVCounts_by_bs_score(test_results[test_results$svtype%in%c('INS:ME','INS:ME:ALU','INS:ME:LINE1','INS:ME:SVA'),], as.character(col_table[col_table[,1]=='MEI',2]), 20, 'Insertions')
  
}

plot_SVCounts_by_bs_score_all_CV<-function(training_data, r1,r2,r3,r4,r5){
  training_data_V2= training_data[training_data$svtype!="BND",]
  bnd_SVID = training_data[training_data$svtype=="BND",]$SVID
  r1=r1[!r1$SVID%in%bnd_SVID,]
  r2=r2[!r2$SVID%in%bnd_SVID,]
  r3=r3[!r3$SVID%in%bnd_SVID,]
  r4=r4[!r4$SVID%in%bnd_SVID,]
  r5=r5[!r5$SVID%in%bnd_SVID,]
  
  sample_name_coord = read.table('../module10_benchmark/gnomad_HGSV.samples_PB_alignments')
  colnames(sample_name_coord)[c(1,2)]=c('gnomad_id','sample')
  sample_pop=read.table('../../HGSV_3/1KGP_2504/manifest/samples_continental_pop.tsv', sep = '\t', header = T, comment.char = "")
  sample_pop_mani = merge(sample_name_coord[,c(1,2)], sample_pop[,c('sample','Continental_group','pop', 'color')], by='sample')
  cff_list=seq(-4,1.6, by=.3)
  i=cff_list[1]
  tmp1 = r1[r1$preds_builtin>i,]
  tmp2 = r2[r2$preds_builtin>i,]
  tmp3 = r3[r3$preds_builtin>i,]
  tmp4 = r4[r4$preds_builtin>i,]
  tmp5 = r5[r5$preds_builtin>i,]
  stat=rbind(data.frame(table(tmp1$sample)),data.frame(table(tmp2$sample)),data.frame(table(tmp3$sample)),data.frame(table(tmp4$sample)),data.frame(table(tmp5$sample)))
  for(i in cff_list[c(2:length(cff_list))]){
    tmp1 = r1[r1$preds_builtin>i,]
    tmp2 = r2[r2$preds_builtin>i,]
    tmp3 = r3[r3$preds_builtin>i,]
    tmp4 = r4[r4$preds_builtin>i,]
    tmp5 = r5[r5$preds_builtin>i,]
    stat_tmp=rbind(data.frame(table(tmp1$sample)),
                   data.frame(table(tmp2$sample)),
                   data.frame(table(tmp3$sample)),
                   data.frame(table(tmp4$sample)),
                   data.frame(table(tmp5$sample)))
    stat = merge(stat, stat_tmp, by='Var1', all=T)
  }
  colnames(stat)[1]='gnomad_id'
  stat2 = merge(stat, sample_pop_mani, by='gnomad_id')
  
  par(fig=c(0,.8,0,1))
  par(mar=c(3,4,4,0))
  plot(c(-4,1.6),c(0,200000), frame.plot = F, type = 'n', xlab = "Boost Score Cutoff", ylab = 'Count of SVs / sample', yaxt='n', title = 'Count of SVs / sample after Boost', cex.axis=1.4)
  axis(2,c(0:10)*20000, labels = c(0:10)*20, las=2, cex.axis=1.4)
  for(i in c(0:10)){abline(h=i*20000, col='grey')}
  for(i in c(0:9)){abline(h=i*20000+10000, col='grey', lty=2)}
  for(i in c(1:nrow(stat2))){
    lines(cff_list, stat[i,c(1:length(cff_list))+1], col = stat2[i,]$color, type='b', pch=20, cex=.9, lty=1)
  }
  legend_table=unique(sample_pop_mani[,c('Continental_group','pop','color')])
  legend_table=legend_table[order(legend_table[,1]),]
  par(fig=c(.8,1,0,1), new=T)
  par(mar=c(3,0,4,0), new=T)
  plot(c(0,1), c(0,1), frame.plot = F, type = 'n', xlab = '', ylab = '', xaxt='n', yaxt='n')
  legend('center',paste(legend_table[,1],legend_table[,2],sep = '-'), col = legend_table$color, pch = 20, bty = 'n')
}

plot_upset_all_CV<-function(training_data, r1,r2,r3,r4,r5, boost_cff=-1){
  training_data_V2= training_data[training_data$svtype!="BND",]
  bnd_SVID = training_data[training_data$svtype=="BND",]$SVID
  r1_tmp=r1[!r1$SVID%in%bnd_SVID,]
  r2_tmp=r2[!r2$SVID%in%bnd_SVID,]
  r3_tmp=r3[!r3$SVID%in%bnd_SVID,]
  r4_tmp=r4[!r4$SVID%in%bnd_SVID,]
  r5_tmp=r5[!r5$SVID%in%bnd_SVID,]
  r1_tmp[,ncol(r1_tmp)+1] = 'r1'
  r2_tmp[,ncol(r2_tmp)+1] = 'r2'
  r3_tmp[,ncol(r3_tmp)+1] = 'r3'
  r4_tmp[,ncol(r4_tmp)+1] = 'r4'
  r5_tmp[,ncol(r5_tmp)+1] = 'r5'
  test_results_r1 =  merge(training_data_V2, r1_tmp, by=c('SVID','sample'))
  test_results_r2 =  merge(training_data_V2, r2_tmp, by=c('SVID','sample'))
  test_results_r3 =  merge(training_data_V2, r3_tmp, by=c('SVID','sample'))
  test_results_r4 =  merge(training_data_V2, r4_tmp, by=c('SVID','sample'))
  test_results_r5 =  merge(training_data_V2, r5_tmp, by=c('SVID','sample'))
  test_results=rbind(test_results_r1, test_results_r2, test_results_r3, test_results_r4, test_results_r5)
  
  bar0 = data.frame(table(test_results[test_results$vs_hgsv_ovr1a=='NO_OVR' & test_results$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/0") & test_results$concor_duplicates==1,][,c('sample','V5')]))
  bar1 = data.frame(table(test_results[test_results$vs_hgsv_ovr1a=='OVR' & test_results$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/0") & test_results$concor_duplicates==1,][,c('sample','V5')]))
  bar2 = data.frame(table(test_results[test_results$vs_hgsv_ovr1a=='NO_OVR' & test_results$vs_pacbio_ovr1a=="OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/0") & test_results$concor_duplicates==1,][,c('sample','V5')]))
  bar3 = data.frame(table(test_results[test_results$vs_hgsv_ovr1a=='NO_OVR' & test_results$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/1", "1/1") & test_results$concor_duplicates==1,][,c('sample','V5')]))
  bar4 = data.frame(table(test_results[test_results$vs_hgsv_ovr1a=='NO_OVR' & test_results$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/0") & test_results$concor_duplicates>1,][,c('sample','V5')]))
  
  bar12 = data.frame(table(test_results[test_results$vs_hgsv_ovr1a=='OVR' & test_results$vs_pacbio_ovr1a=="OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/0") & test_results$concor_duplicates==1,][,c('sample','V5')]))
  bar13 = data.frame(table(test_results[test_results$vs_hgsv_ovr1a=='OVR' & test_results$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/1", "1/1") & test_results$concor_duplicates==1,][,c('sample','V5')]))
  bar14 = data.frame(table(test_results[test_results$vs_hgsv_ovr1a=='OVR' & test_results$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/0") & test_results$concor_duplicates>1,][,c('sample','V5')]))
  bar23 = data.frame(table(test_results[test_results$vs_hgsv_ovr1a=='NO_OVR' & test_results$vs_pacbio_ovr1a=="OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/1", "1/1") & test_results$concor_duplicates==1,][,c('sample','V5')]))
  bar24 = data.frame(table(test_results[test_results$vs_hgsv_ovr1a=='NO_OVR' & test_results$vs_pacbio_ovr1a=="OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/0") & test_results$concor_duplicates>1,][,c('sample','V5')]))
  bar34 = data.frame(table(test_results[test_results$vs_hgsv_ovr1a=='NO_OVR' & test_results$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/1", "1/1") & test_results$concor_duplicates>1,][,c('sample','V5')]))
  
  bar123 = data.frame(table(test_results[test_results$vs_hgsv_ovr1a=='OVR' & test_results$vs_pacbio_ovr1a=="OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/1", "1/1") & test_results$concor_duplicates==1,][,c('sample','V5')]))
  bar124 = data.frame(table(test_results[test_results$vs_hgsv_ovr1a=='OVR' & test_results$vs_pacbio_ovr1a=="OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/0") & test_results$concor_duplicates>1,][,c('sample','V5')]))
  bar134 = data.frame(table(test_results[test_results$vs_hgsv_ovr1a=='OVR' & test_results$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/1", "1/1") & test_results$concor_duplicates>1,][,c('sample','V5')]))
  bar234 = data.frame(table(test_results[test_results$vs_hgsv_ovr1a=='NO_OVR' & test_results$vs_pacbio_ovr1a=="OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/1", "1/1") & test_results$concor_duplicates>1,][,c('sample','V5')]))
  bar1234 = data.frame(table(test_results[test_results$vs_hgsv_ovr1a=='OVR' & test_results$vs_pacbio_ovr1a=="OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/1", "1/1") & test_results$concor_duplicates>1,][,c('sample','V5')]))
  
  cff = boost_cff
  #test_results$preds_out=0
  #test_results[test_results$preds_builtin>cff,]$preds_out=1
  test_results_pass = test_results[test_results$preds_builtin>boost_cff,]
  bar0b = data.frame(table(test_results_pass[test_results_pass$vs_hgsv_ovr1a=='NO_OVR' & test_results_pass$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results_pass$VaPoR_GT) & test_results_pass$VaPoR_GT%in%c("0/0") & test_results_pass$concor_duplicates==1,][,c('sample','V5')]))
  bar1b = data.frame(table(test_results_pass[test_results_pass$vs_hgsv_ovr1a=='OVR' & test_results_pass$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results_pass$VaPoR_GT) & test_results_pass$VaPoR_GT%in%c("0/0") & test_results_pass$concor_duplicates==1,][,c('sample','V5')]))
  bar2b = data.frame(table(test_results_pass[test_results_pass$vs_hgsv_ovr1a=='NO_OVR' & test_results_pass$vs_pacbio_ovr1a=="OVR" & !is.na(test_results_pass$VaPoR_GT) & test_results_pass$VaPoR_GT%in%c("0/0") & test_results_pass$concor_duplicates==1,][,c('sample','V5')]))
  bar3b = data.frame(table(test_results_pass[test_results_pass$vs_hgsv_ovr1a=='NO_OVR' & test_results_pass$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results_pass$VaPoR_GT) & test_results_pass$VaPoR_GT%in%c("0/1", "1/1") & test_results_pass$concor_duplicates==1,][,c('sample','V5')]))
  bar4b = data.frame(table(test_results_pass[test_results_pass$vs_hgsv_ovr1a=='NO_OVR' & test_results_pass$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results_pass$VaPoR_GT) & test_results_pass$VaPoR_GT%in%c("0/0") & test_results_pass$concor_duplicates>1,][,c('sample','V5')]))
  
  bar12b = data.frame(table(test_results_pass[test_results_pass$vs_hgsv_ovr1a=='OVR' & test_results_pass$vs_pacbio_ovr1a=="OVR" & !is.na(test_results_pass$VaPoR_GT) & test_results_pass$VaPoR_GT%in%c("0/0") & test_results_pass$concor_duplicates==1,][,c('sample','V5')]))
  bar13b = data.frame(table(test_results_pass[test_results_pass$vs_hgsv_ovr1a=='OVR' & test_results_pass$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results_pass$VaPoR_GT) & test_results_pass$VaPoR_GT%in%c("0/1", "1/1") & test_results_pass$concor_duplicates==1,][,c('sample','V5')]))
  bar14b = data.frame(table(test_results_pass[test_results_pass$vs_hgsv_ovr1a=='OVR' & test_results_pass$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results_pass$VaPoR_GT) & test_results_pass$VaPoR_GT%in%c("0/0") & test_results_pass$concor_duplicates>1,][,c('sample','V5')]))
  bar23b = data.frame(table(test_results_pass[test_results_pass$vs_hgsv_ovr1a=='NO_OVR' & test_results_pass$vs_pacbio_ovr1a=="OVR" & !is.na(test_results_pass$VaPoR_GT) & test_results_pass$VaPoR_GT%in%c("0/1", "1/1") & test_results_pass$concor_duplicates==1,][,c('sample','V5')]))
  bar24b = data.frame(table(test_results_pass[test_results_pass$vs_hgsv_ovr1a=='NO_OVR' & test_results_pass$vs_pacbio_ovr1a=="OVR" & !is.na(test_results_pass$VaPoR_GT) & test_results_pass$VaPoR_GT%in%c("0/0") & test_results_pass$concor_duplicates>1,][,c('sample','V5')]))
  bar34b = data.frame(table(test_results_pass[test_results_pass$vs_hgsv_ovr1a=='NO_OVR' & test_results_pass$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results_pass$VaPoR_GT) & test_results_pass$VaPoR_GT%in%c("0/1", "1/1") & test_results_pass$concor_duplicates>1,][,c('sample','V5')]))
  
  bar123b = data.frame(table(test_results_pass[test_results_pass$vs_hgsv_ovr1a=='OVR' & test_results_pass$vs_pacbio_ovr1a=="OVR" & !is.na(test_results_pass$VaPoR_GT) & test_results_pass$VaPoR_GT%in%c("0/1", "1/1") & test_results_pass$concor_duplicates==1,][,c('sample','V5')]))
  bar124b = data.frame(table(test_results_pass[test_results_pass$vs_hgsv_ovr1a=='OVR' & test_results_pass$vs_pacbio_ovr1a=="OVR" & !is.na(test_results_pass$VaPoR_GT) & test_results_pass$VaPoR_GT%in%c("0/0") & test_results_pass$concor_duplicates>1,][,c('sample','V5')]))
  bar134b = data.frame(table(test_results_pass[test_results_pass$vs_hgsv_ovr1a=='OVR' & test_results_pass$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results_pass$VaPoR_GT) & test_results_pass$VaPoR_GT%in%c("0/1", "1/1") & test_results_pass$concor_duplicates>1,][,c('sample','V5')]))
  bar234b = data.frame(table(test_results_pass[test_results_pass$vs_hgsv_ovr1a=='NO_OVR' & test_results_pass$vs_pacbio_ovr1a=="OVR" & !is.na(test_results_pass$VaPoR_GT) & test_results_pass$VaPoR_GT%in%c("0/1", "1/1") & test_results_pass$concor_duplicates>1,][,c('sample','V5')]))
  bar1234b = data.frame(table(test_results_pass[test_results_pass$vs_hgsv_ovr1a=='OVR' & test_results_pass$vs_pacbio_ovr1a=="OVR" & !is.na(test_results_pass$VaPoR_GT) & test_results_pass$VaPoR_GT%in%c("0/1", "1/1") & test_results_pass$concor_duplicates>1,][,c('sample','V5')]))
  
  bar_list = c(median(bar0[bar0[,3]>0,][,3]), median(bar1[bar1[,3]>0,][,3]),median(bar2[bar2[,3]>0,][,3]),median(bar3[bar3[,3]>0,][,3]),median(bar4[bar4[,3]>0,][,3]), 
               median(bar12[bar12[,3]>0,][,3]), median(bar13[bar13[,3]>0,][,3]), median(bar14[bar14[,3]>0,][,3]), median(bar23[bar23[,3]>0,][,3]), median(bar24[bar24[,3]>0,][,3]), median(bar34[bar34[,3]>0,][,3]), 
               median(bar123[bar123[,3]>0,][,3]), median(bar124[bar124[,3]>0,][,3]), median(bar134[bar134[,3]>0,][,3]), median(bar234[bar234[,3]>0,][,3]), median(bar1234[bar1234[,3]>0,][,3]))
  bar_list_b = c(median(bar0b[bar0b[,3]>0,][,3]), median(bar1b[bar1b[,3]>0,][,3]),median(bar2b[bar2b[,3]>0,][,3]),median(bar3b[bar3b[,3]>0,][,3]),median(bar4b[bar4b[,3]>0,][,3]), 
                 median(bar12b[bar12b[,3]>0,][,3]), median(bar13b[bar13b[,3]>0,][,3]), median(bar14b[bar14b[,3]>0,][,3]), median(bar23b[bar23b[,3]>0,][,3]), median(bar24b[bar24b[,3]>0,][,3]), median(bar34b[bar34b[,3]>0,][,3]), 
                 median(bar123b[bar123b[,3]>0,][,3]), median(bar124b[bar124b[,3]>0,][,3]), median(bar134b[bar134b[,3]>0,][,3]), median(bar234b[bar234b[,3]>0,][,3]), median(bar1234b[bar1234b[,3]>0,][,3]))
  y_range = range(bar_list[c(2:length(bar_list))])
  
  par(mar=c(0,6,1,2))
  par(fig=c(0,1,.65,1))
  plot(c(0,length(bar_list_b)), range(bar_list), frame.plot = F, type = 'n', xlab = '', ylab = 'Count of SVs / Genome', las=2, xaxt='n', mgp=c(4,1,0), main = 'SVs per sample Before / After Boost Model')
  bar_width = .4
  x_pos = .5
  bar_colors = c('grey', 'black')
  for(i in c(1:length(bar_list))){
    rect(x_pos - bar_width, 0,  x_pos + bar_width, bar_list[i], col = bar_colors[1] )
    rect(x_pos - bar_width, 0,  x_pos + bar_width, bar_list_b[i], col = bar_colors[2] )
    x_pos = x_pos + 1
  }
  legend('topright', c('From Module05', 'After Boost Filter'), title = 'Count of SVs', bty = 'n', col = bar_colors[c(1,2)], pch = 15, pt.cex = 2)
  
  par(fig=c(0,1,.3,.65), new=T)
  plot(c(0,length(bar_list_b)), range(0,1.1), frame.plot = F, type = 'n', xlab = '', ylab = 'Proportion of SVs \n / Genome', las=2, xaxt='n', mgp=c(3,1,0))
  bar_width = .4
  x_pos = .5
  bar_colors = c('grey', 'black')
  for(i in c(1:length(bar_list))){
    rect(x_pos - bar_width, 0,  x_pos + bar_width, 1, col = bar_colors[1] )
    rect(x_pos - bar_width, 0,  x_pos + bar_width, bar_list_b[i]/bar_list[i], col = bar_colors[2] )
    #text(x_pos, 1.05, bar_list[i], srt=45, cex = .6)
    x_pos = x_pos + 1
  }
  
  par(fig=c(0,1,0,.3), new=T)
  plot(c(0,length(bar_list_b)), c(0,4), frame.plot = F, type = 'n', xlab = '', ylab = '', xaxt='n', yaxt='n')
  axis(2,c(1:4)-.5, labels = c('Replicate \n Concordance','VaPoR \n Validated','PacBio \n Overlap','Illumina \n Overlap'), las=2, cex.axis=.8)
  for(i in c(1:length(bar_list_b))-.5){
    for(j in c(1:4)-.5){
      points(i,j, col='grey', pch=20, cex = 2)
    }
  }
  for(i in c(1:4)+.5){  points(i, 5-i, col='black', pch=20, cex = 2)}
  start_pos = 5.5
  for(i in c(3,4)-.5){  points(start_pos,i, col='black', pch=20, cex = 2)}
  for(i in c(2,4)-.5){  points(start_pos+1,i, col='black', pch=20, cex = 2)}
  for(i in c(1,4)-.5){  points(start_pos+2,i, col='black', pch=20, cex = 2)}
  for(i in c(2,3)-.5){  points(start_pos+3,i, col='black', pch=20, cex = 2)}
  for(i in c(1,3)-.5){  points(start_pos+4,i, col='black', pch=20, cex = 2)}
  for(i in c(1,2)-.5){  points(start_pos+5,i, col='black', pch=20, cex = 2)}
  for(i in c(2,3,4)-.5){  points(start_pos+6,i, col='black', pch=20, cex = 2)}
  for(i in c(1,3,4)-.5){  points(start_pos+7,i, col='black', pch=20, cex = 2)}
  for(i in c(1,2,4)-.5){  points(start_pos+8,i, col='black', pch=20, cex = 2)}
  for(i in c(1,2,3)-.5){  points(start_pos+9,i, col='black', pch=20, cex = 2)}
  for(i in c(1:4)-.5)  {  points(start_pos+10,i, col='black', pch=20, cex = 2)}
  lines(c(start_pos,start_pos),c(3:4)-.5)
  lines(c(start_pos+1,start_pos+1),c(2,4)-.5)
  lines(c(start_pos+2,start_pos+2),c(1,4)-.5)
  lines(c(start_pos+3,start_pos+3),c(2,3)-.5)
  lines(c(start_pos+4,start_pos+4),c(1,3)-.5)
  lines(c(start_pos+5,start_pos+5),c(1,2)-.5)
  lines(c(start_pos+6,start_pos+6),c(2,4)-.5)
  lines(c(start_pos+7,start_pos+7),c(1,4)-.5)
  lines(c(start_pos+8,start_pos+8),c(1,4)-.5)
  lines(c(start_pos+9,start_pos+9),c(1,3)-.5)
  lines(c(start_pos+10,start_pos+10),c(1,4)-.5)
  
}

plot_upset_each_CV<-function(training_data, r1, boost_cff=-1){
  training_data_V2= training_data[training_data$svtype!="BND",]
  bnd_SVID = training_data[training_data$svtype=="BND",]$SVID
  r1_tmp=r1[!r1$SVID%in%bnd_SVID,]
  r1_tmp[,ncol(r1_tmp)+1] = 'r1'
  test_results_r1 =  merge(training_data_V2, r1_tmp, by=c('SVID','sample'))
  test_results=test_results_r1
  
  bar0 = data.frame(table(test_results[test_results$vs_hgsv_ovr1a=='NO_OVR' & test_results$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/0") & test_results$concor_duplicates==1,][,c('sample','V5')]))
  bar1 = data.frame(table(test_results[test_results$vs_hgsv_ovr1a=='OVR' & test_results$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/0") & test_results$concor_duplicates==1,][,c('sample','V5')]))
  bar2 = data.frame(table(test_results[test_results$vs_hgsv_ovr1a=='NO_OVR' & test_results$vs_pacbio_ovr1a=="OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/0") & test_results$concor_duplicates==1,][,c('sample','V5')]))
  bar3 = data.frame(table(test_results[test_results$vs_hgsv_ovr1a=='NO_OVR' & test_results$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/1", "1/1") & test_results$concor_duplicates==1,][,c('sample','V5')]))
  bar4 = data.frame(table(test_results[test_results$vs_hgsv_ovr1a=='NO_OVR' & test_results$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/0") & test_results$concor_duplicates>1,][,c('sample','V5')]))
  
  bar12 = data.frame(table(test_results[test_results$vs_hgsv_ovr1a=='OVR' & test_results$vs_pacbio_ovr1a=="OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/0") & test_results$concor_duplicates==1,][,c('sample','V5')]))
  bar13 = data.frame(table(test_results[test_results$vs_hgsv_ovr1a=='OVR' & test_results$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/1", "1/1") & test_results$concor_duplicates==1,][,c('sample','V5')]))
  bar14 = data.frame(table(test_results[test_results$vs_hgsv_ovr1a=='OVR' & test_results$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/0") & test_results$concor_duplicates>1,][,c('sample','V5')]))
  bar23 = data.frame(table(test_results[test_results$vs_hgsv_ovr1a=='NO_OVR' & test_results$vs_pacbio_ovr1a=="OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/1", "1/1") & test_results$concor_duplicates==1,][,c('sample','V5')]))
  bar24 = data.frame(table(test_results[test_results$vs_hgsv_ovr1a=='NO_OVR' & test_results$vs_pacbio_ovr1a=="OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/0") & test_results$concor_duplicates>1,][,c('sample','V5')]))
  bar34 = data.frame(table(test_results[test_results$vs_hgsv_ovr1a=='NO_OVR' & test_results$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/1", "1/1") & test_results$concor_duplicates>1,][,c('sample','V5')]))
  
  bar123 = data.frame(table(test_results[test_results$vs_hgsv_ovr1a=='OVR' & test_results$vs_pacbio_ovr1a=="OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/1", "1/1") & test_results$concor_duplicates==1,][,c('sample','V5')]))
  bar124 = data.frame(table(test_results[test_results$vs_hgsv_ovr1a=='OVR' & test_results$vs_pacbio_ovr1a=="OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/0") & test_results$concor_duplicates>1,][,c('sample','V5')]))
  bar134 = data.frame(table(test_results[test_results$vs_hgsv_ovr1a=='OVR' & test_results$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/1", "1/1") & test_results$concor_duplicates>1,][,c('sample','V5')]))
  bar234 = data.frame(table(test_results[test_results$vs_hgsv_ovr1a=='NO_OVR' & test_results$vs_pacbio_ovr1a=="OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/1", "1/1") & test_results$concor_duplicates>1,][,c('sample','V5')]))
  bar1234 = data.frame(table(test_results[test_results$vs_hgsv_ovr1a=='OVR' & test_results$vs_pacbio_ovr1a=="OVR" & !is.na(test_results$VaPoR_GT) & test_results$VaPoR_GT%in%c("0/1", "1/1") & test_results$concor_duplicates>1,][,c('sample','V5')]))
  
  cff = boost_cff
  #test_results$preds_out=0
  #test_results[test_results$preds_builtin>cff,]$preds_out=1
  test_results_pass = test_results[test_results$preds_builtin>boost_cff,]
  bar0b = data.frame(table(test_results_pass[test_results_pass$vs_hgsv_ovr1a=='NO_OVR' & test_results_pass$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results_pass$VaPoR_GT) & test_results_pass$VaPoR_GT%in%c("0/0") & test_results_pass$concor_duplicates==1,][,c('sample','V5')]))
  bar1b = data.frame(table(test_results_pass[test_results_pass$vs_hgsv_ovr1a=='OVR' & test_results_pass$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results_pass$VaPoR_GT) & test_results_pass$VaPoR_GT%in%c("0/0") & test_results_pass$concor_duplicates==1,][,c('sample','V5')]))
  bar2b = data.frame(table(test_results_pass[test_results_pass$vs_hgsv_ovr1a=='NO_OVR' & test_results_pass$vs_pacbio_ovr1a=="OVR" & !is.na(test_results_pass$VaPoR_GT) & test_results_pass$VaPoR_GT%in%c("0/0") & test_results_pass$concor_duplicates==1,][,c('sample','V5')]))
  bar3b = data.frame(table(test_results_pass[test_results_pass$vs_hgsv_ovr1a=='NO_OVR' & test_results_pass$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results_pass$VaPoR_GT) & test_results_pass$VaPoR_GT%in%c("0/1", "1/1") & test_results_pass$concor_duplicates==1,][,c('sample','V5')]))
  bar4b = data.frame(table(test_results_pass[test_results_pass$vs_hgsv_ovr1a=='NO_OVR' & test_results_pass$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results_pass$VaPoR_GT) & test_results_pass$VaPoR_GT%in%c("0/0") & test_results_pass$concor_duplicates>1,][,c('sample','V5')]))
  
  bar12b = data.frame(table(test_results_pass[test_results_pass$vs_hgsv_ovr1a=='OVR' & test_results_pass$vs_pacbio_ovr1a=="OVR" & !is.na(test_results_pass$VaPoR_GT) & test_results_pass$VaPoR_GT%in%c("0/0") & test_results_pass$concor_duplicates==1,][,c('sample','V5')]))
  bar13b = data.frame(table(test_results_pass[test_results_pass$vs_hgsv_ovr1a=='OVR' & test_results_pass$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results_pass$VaPoR_GT) & test_results_pass$VaPoR_GT%in%c("0/1", "1/1") & test_results_pass$concor_duplicates==1,][,c('sample','V5')]))
  bar14b = data.frame(table(test_results_pass[test_results_pass$vs_hgsv_ovr1a=='OVR' & test_results_pass$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results_pass$VaPoR_GT) & test_results_pass$VaPoR_GT%in%c("0/0") & test_results_pass$concor_duplicates>1,][,c('sample','V5')]))
  bar23b = data.frame(table(test_results_pass[test_results_pass$vs_hgsv_ovr1a=='NO_OVR' & test_results_pass$vs_pacbio_ovr1a=="OVR" & !is.na(test_results_pass$VaPoR_GT) & test_results_pass$VaPoR_GT%in%c("0/1", "1/1") & test_results_pass$concor_duplicates==1,][,c('sample','V5')]))
  bar24b = data.frame(table(test_results_pass[test_results_pass$vs_hgsv_ovr1a=='NO_OVR' & test_results_pass$vs_pacbio_ovr1a=="OVR" & !is.na(test_results_pass$VaPoR_GT) & test_results_pass$VaPoR_GT%in%c("0/0") & test_results_pass$concor_duplicates>1,][,c('sample','V5')]))
  bar34b = data.frame(table(test_results_pass[test_results_pass$vs_hgsv_ovr1a=='NO_OVR' & test_results_pass$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results_pass$VaPoR_GT) & test_results_pass$VaPoR_GT%in%c("0/1", "1/1") & test_results_pass$concor_duplicates>1,][,c('sample','V5')]))
  
  bar123b = data.frame(table(test_results_pass[test_results_pass$vs_hgsv_ovr1a=='OVR' & test_results_pass$vs_pacbio_ovr1a=="OVR" & !is.na(test_results_pass$VaPoR_GT) & test_results_pass$VaPoR_GT%in%c("0/1", "1/1") & test_results_pass$concor_duplicates==1,][,c('sample','V5')]))
  bar124b = data.frame(table(test_results_pass[test_results_pass$vs_hgsv_ovr1a=='OVR' & test_results_pass$vs_pacbio_ovr1a=="OVR" & !is.na(test_results_pass$VaPoR_GT) & test_results_pass$VaPoR_GT%in%c("0/0") & test_results_pass$concor_duplicates>1,][,c('sample','V5')]))
  bar134b = data.frame(table(test_results_pass[test_results_pass$vs_hgsv_ovr1a=='OVR' & test_results_pass$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results_pass$VaPoR_GT) & test_results_pass$VaPoR_GT%in%c("0/1", "1/1") & test_results_pass$concor_duplicates>1,][,c('sample','V5')]))
  bar234b = data.frame(table(test_results_pass[test_results_pass$vs_hgsv_ovr1a=='NO_OVR' & test_results_pass$vs_pacbio_ovr1a=="OVR" & !is.na(test_results_pass$VaPoR_GT) & test_results_pass$VaPoR_GT%in%c("0/1", "1/1") & test_results_pass$concor_duplicates>1,][,c('sample','V5')]))
  bar1234b = data.frame(table(test_results_pass[test_results_pass$vs_hgsv_ovr1a=='OVR' & test_results_pass$vs_pacbio_ovr1a=="OVR" & !is.na(test_results_pass$VaPoR_GT) & test_results_pass$VaPoR_GT%in%c("0/1", "1/1") & test_results_pass$concor_duplicates>1,][,c('sample','V5')]))
  
  bar_list = c(median(bar0[bar0[,3]>0,][,3]), median(bar1[bar1[,3]>0,][,3]),median(bar2[bar2[,3]>0,][,3]),median(bar3[bar3[,3]>0,][,3]),median(bar4[bar4[,3]>0,][,3]), 
               median(bar12[bar12[,3]>0,][,3]), median(bar13[bar13[,3]>0,][,3]), median(bar14[bar14[,3]>0,][,3]), median(bar23[bar23[,3]>0,][,3]), median(bar24[bar24[,3]>0,][,3]), median(bar34[bar34[,3]>0,][,3]), 
               median(bar123[bar123[,3]>0,][,3]), median(bar124[bar124[,3]>0,][,3]), median(bar134[bar134[,3]>0,][,3]), median(bar234[bar234[,3]>0,][,3]), median(bar1234[bar1234[,3]>0,][,3]))
  bar_min = c(min(bar0[bar0[,3]>0,][,3]), min(bar1[bar1[,3]>0,][,3]),min(bar2[bar2[,3]>0,][,3]),min(bar3[bar3[,3]>0,][,3]),min(bar4[bar4[,3]>0,][,3]), 
              min(bar12[bar12[,3]>0,][,3]), min(bar13[bar13[,3]>0,][,3]), min(bar14[bar14[,3]>0,][,3]), min(bar23[bar23[,3]>0,][,3]), min(bar24[bar24[,3]>0,][,3]), min(bar34[bar34[,3]>0,][,3]), 
              min(bar123[bar123[,3]>0,][,3]), min(bar124[bar124[,3]>0,][,3]), min(bar134[bar134[,3]>0,][,3]), min(bar234[bar234[,3]>0,][,3]), min(bar1234[bar1234[,3]>0,][,3]))
  bar_max = c(max(bar0[bar0[,3]>0,][,3]), max(bar1[bar1[,3]>0,][,3]),max(bar2[bar2[,3]>0,][,3]),max(bar3[bar3[,3]>0,][,3]),max(bar4[bar4[,3]>0,][,3]), 
              max(bar12[bar12[,3]>0,][,3]), max(bar13[bar13[,3]>0,][,3]), max(bar14[bar14[,3]>0,][,3]), max(bar23[bar23[,3]>0,][,3]), max(bar24[bar24[,3]>0,][,3]), max(bar34[bar34[,3]>0,][,3]), 
              max(bar123[bar123[,3]>0,][,3]), max(bar124[bar124[,3]>0,][,3]), max(bar134[bar134[,3]>0,][,3]), max(bar234[bar234[,3]>0,][,3]), max(bar1234[bar1234[,3]>0,][,3]))
  bar_list_b = c(median(bar0b[bar0b[,3]>0,][,3]), median(bar1b[bar1b[,3]>0,][,3]),median(bar2b[bar2b[,3]>0,][,3]),median(bar3b[bar3b[,3]>0,][,3]),median(bar4b[bar4b[,3]>0,][,3]), 
                 median(bar12b[bar12b[,3]>0,][,3]), median(bar13b[bar13b[,3]>0,][,3]), median(bar14b[bar14b[,3]>0,][,3]), median(bar23b[bar23b[,3]>0,][,3]), median(bar24b[bar24b[,3]>0,][,3]), median(bar34b[bar34b[,3]>0,][,3]), 
                 median(bar123b[bar123b[,3]>0,][,3]), median(bar124b[bar124b[,3]>0,][,3]), median(bar134b[bar134b[,3]>0,][,3]), median(bar234b[bar234b[,3]>0,][,3]), median(bar1234b[bar1234b[,3]>0,][,3]))
  y_range = range(bar_list[c(2:length(bar_list))])
  
  par(mar=c(0,6,1,2))
  par(fig=c(0,1,.65,1))
  plot(c(0,length(bar_list_b)),c(0,max(bar_max)), frame.plot = F, type = 'n', xlab = '', ylab = 'Count of SVs / Genome', las=2, xaxt='n', mgp=c(4,1,0), main = 'SVs per sample Before / After Boost Model')
  bar_width = .4
  x_pos = .5
  bar_colors = c('grey', 'black')
  for(i in c(1:length(bar_list))){
    rect(x_pos - bar_width, 0,  x_pos + bar_width, bar_list[i], col = bar_colors[1] )
    lines(c(x_pos - bar_width*.6,x_pos + bar_width*.6),c(bar_min[i], bar_min[i]), col='black', lwd=2)
    lines(c(x_pos - bar_width*.6,x_pos + bar_width*.6),c(bar_max[i], bar_max[i]), col='black', lwd=2)
    lines(c(x_pos,x_pos),c(bar_min[i], bar_max[i]), col='black', lwd=2)
    rect(x_pos - bar_width, 0,  x_pos + bar_width, bar_list_b[i], col = bar_colors[2] )
    x_pos = x_pos + 1
  }
  legend('topright', c('From Module05', 'After Boost Filter'), title = 'Count of SVs', bty = 'n', col = bar_colors[c(1,2)], pch = 15, pt.cex = 2)
  
  par(fig=c(0,1,.3,.65), new=T)
  plot(c(0,length(bar_list_b)), range(0,1.1), frame.plot = F, type = 'n', xlab = '', ylab = 'Proportion of SVs \n / Genome', las=2, xaxt='n', mgp=c(3,1,0))
  bar_width = .4
  x_pos = .5
  bar_colors = c('grey', 'black')
  for(i in c(1:length(bar_list))){
    rect(x_pos - bar_width, 0,  x_pos + bar_width, 1, col = bar_colors[1] )
    rect(x_pos - bar_width, 0,  x_pos + bar_width, bar_list_b[i]/bar_list[i], col = bar_colors[2] )
    #text(x_pos, 1.05, bar_list[i], srt=45, cex = .6)
    x_pos = x_pos + 1
  }
  
  par(fig=c(0,1,0,.3), new=T)
  plot(c(0,length(bar_list_b)), c(0,4), frame.plot = F, type = 'n', xlab = '', ylab = '', xaxt='n', yaxt='n')
  axis(2,c(1:4)-.5, labels = c('Replicate \n Concordance','VaPoR \n Validated','PacBio \n Overlap','Illumina \n Overlap'), las=2, cex.axis=.8)
  for(i in c(1:length(bar_list_b))-.5){
    for(j in c(1:4)-.5){
      points(i,j, col='grey', pch=20, cex = 2)
    }
  }
  for(i in c(1:4)+.5){  points(i, 5-i, col='black', pch=20, cex = 2)}
  start_pos = 5.5
  for(i in c(3,4)-.5){  points(start_pos,i, col='black', pch=20, cex = 2)}
  for(i in c(2,4)-.5){  points(start_pos+1,i, col='black', pch=20, cex = 2)}
  for(i in c(1,4)-.5){  points(start_pos+2,i, col='black', pch=20, cex = 2)}
  for(i in c(2,3)-.5){  points(start_pos+3,i, col='black', pch=20, cex = 2)}
  for(i in c(1,3)-.5){  points(start_pos+4,i, col='black', pch=20, cex = 2)}
  for(i in c(1,2)-.5){  points(start_pos+5,i, col='black', pch=20, cex = 2)}
  for(i in c(2,3,4)-.5){  points(start_pos+6,i, col='black', pch=20, cex = 2)}
  for(i in c(1,3,4)-.5){  points(start_pos+7,i, col='black', pch=20, cex = 2)}
  for(i in c(1,2,4)-.5){  points(start_pos+8,i, col='black', pch=20, cex = 2)}
  for(i in c(1,2,3)-.5){  points(start_pos+9,i, col='black', pch=20, cex = 2)}
  for(i in c(1:4)-.5)  {  points(start_pos+10,i, col='black', pch=20, cex = 2)}
  lines(c(start_pos,start_pos),c(3:4)-.5)
  lines(c(start_pos+1,start_pos+1),c(2,4)-.5)
  lines(c(start_pos+2,start_pos+2),c(1,4)-.5)
  lines(c(start_pos+3,start_pos+3),c(2,3)-.5)
  lines(c(start_pos+4,start_pos+4),c(1,3)-.5)
  lines(c(start_pos+5,start_pos+5),c(1,2)-.5)
  lines(c(start_pos+6,start_pos+6),c(2,4)-.5)
  lines(c(start_pos+7,start_pos+7),c(1,4)-.5)
  lines(c(start_pos+8,start_pos+8),c(1,4)-.5)
  lines(c(start_pos+9,start_pos+9),c(1,3)-.5)
  lines(c(start_pos+10,start_pos+10),c(1,4)-.5)
  
}

plot_fdr_vs_BScff<-function(training_data, r1,r2,r3,r4,r5){
  training_data_V2= training_data[training_data$svtype!="BND",]
  bnd_SVID = training_data[training_data$svtype=="BND",]$SVID
  r1_tmp=r1[!r1$SVID%in%bnd_SVID,]
  r2_tmp=r2[!r2$SVID%in%bnd_SVID,]
  r3_tmp=r3[!r3$SVID%in%bnd_SVID,]
  r4_tmp=r4[!r4$SVID%in%bnd_SVID,]
  r5_tmp=r5[!r5$SVID%in%bnd_SVID,]
  r1_tmp[,ncol(r1_tmp)+1] = 'r1'
  r2_tmp[,ncol(r2_tmp)+1] = 'r2'
  r3_tmp[,ncol(r3_tmp)+1] = 'r3'
  r4_tmp[,ncol(r4_tmp)+1] = 'r4'
  r5_tmp[,ncol(r5_tmp)+1] = 'r5'
  test_results_r1 =  merge(training_data_V2, r1_tmp, by=c('SVID','sample'))
  test_results_r2 =  merge(training_data_V2, r2_tmp, by=c('SVID','sample'))
  test_results_r3 =  merge(training_data_V2, r3_tmp, by=c('SVID','sample'))
  test_results_r4 =  merge(training_data_V2, r4_tmp, by=c('SVID','sample'))
  test_results_r5 =  merge(training_data_V2, r5_tmp, by=c('SVID','sample'))
  test_results=rbind(test_results_r1, test_results_r2, test_results_r3, test_results_r4, test_results_r5)
  
  out_table=data.frame('bs_cff'=0,'fdr'=0)
  num_table=0
  for(cff in seq(-3,0,by=.2)){
    print(cff)
    test_results_pass = test_results[test_results$preds_builtin>cff,]
    testx = data.frame(table(test_results_pass[,c('sample','V5')]))
    bar0b = data.frame(table(test_results_pass[test_results_pass$vs_hgsv_ovr1a=='NO_OVR' & test_results_pass$vs_pacbio_ovr1a=="NO_OVR" & !is.na(test_results_pass$VaPoR_GT) & test_results_pass$VaPoR_GT%in%c("0/0") & test_results_pass$concor_duplicates==1,][,c('sample','V5')]))
    stat = merge(testx, bar0b, by=c('sample','V5'))
    stat=stat[stat[,3]>0,]
    num_table=num_table+1
    out_table[num_table,1] = cff
    out_table[num_table,2] = sum(stat[,4])/sum(stat[,3])
  }
  
  par(fig=c(0,.7,0,1))
  plot(range(out_table[,1]), c(0,.35), frame.plot = F, type = 'n', xlab = '', ylab = '', yaxt='n', cex.axis=1.4)
  axis(2,c(0:7)/20, labels = paste(c(0:7)/20*100,'%', sep = ''), las=2, cex.axis=1.4)
  for(i in c(0:7)/20){abline(h=i, col='grey')}
  lines(out_table[,1], out_table[,2], pch=20, type = 'b', lty=1, lwd=2)
  return(out_table)
}

readin_site_features<-function(site_folder = './per_chr_anno/'){
  # read in the per-site traits across chromosomes
  list_file = list.files(site_folder)
  site_feature = read.table(paste(site_folder, list_file[1] ,sep = ''), header = T, comment.char = "")
  print(list_file[1])
  for(site_file in list_file[2:length(list_file)]){
    chr_feature = read.table(paste(site_folder, site_file, sep = ''), header = T, comment.char = "")
    print(site_file)
    site_feature = rbind(site_feature, chr_feature)
  }
  return(site_feature)
}

readin_training<-function(samp_feature, site_feature){
  training_set = merge(samp_feature, site_feature, by='SVID')
  #modify size category:
  training_set$size_cate = 's1_under250bp'
  training_set[training_set$SVLEN.x>250,]$size_cate = 's2_250bpto1Kb'
  training_set[training_set$SVLEN.x>1000,]$size_cate = 's3_1to5Kb'
  training_set[training_set$SVLEN.x>5000,]$size_cate = 's4_5to50Kb'
  training_set[training_set$SVLEN.x>50000,]$size_cate = 's5_over50Kb'
  
  for(i in c(1:ncol(training_set))){
    if(grepl('vs_',colnames(training_set)[i])){
      print(colnames(training_set)[i])
      if(nrow(training_set[training_set[,i]!='NO_OVR',])>0){
        training_set[,i]=as.character(training_set[,i])
        training_set[!is.na(training_set[,i]) & training_set[,i]!='NO_OVR',][,i]='OVR'
      } } }
  
  return(training_set)
}

test_del_dup_ins<-function(bgm_models, test, feature_list_cnv, feature_list_oth, cff=-0.8){
  models = bgm_models
  size='s1_under250bp'
  preds_builtin <- predict(models[[1]], data.matrix(test[test$svtype=='DEL' & test$size_cate==size ,][,feature_list_oth]), rawscore = TRUE, reshape = TRUE)
  train_recali_1a=cbind(test[test$svtype=='DEL' & test$size_cate==size ,], preds_builtin)
  
  preds_builtin <- predict(models[[2]], data.matrix(test[test$svtype=='DUP' & test$size_cate==size ,][,feature_list_oth]), rawscore = TRUE, reshape = TRUE)
  train_recali_1b=cbind(test[test$svtype=='DUP' & test$size_cate==size ,], preds_builtin)
  
  size='s2_250bpto1Kb'
  preds_builtin <- predict(models[[3]], data.matrix(test[test$svtype=='DEL' & test$size_cate==size ,][,feature_list_oth]), rawscore = TRUE, reshape = TRUE)
  train_recali_2a=cbind(test[test$svtype=='DEL' & test$size_cate==size ,], preds_builtin)
  
  preds_builtin <- predict(models[[4]], data.matrix(test[test$svtype=='DUP' & test$size_cate==size ,][,feature_list_oth]), rawscore = TRUE, reshape = TRUE)
  train_recali_2b=cbind(test[test$svtype=='DUP' & test$size_cate==size ,], preds_builtin)
  
  size='s3_1to5Kb'
  preds_builtin <- predict(models[[5]], data.matrix(test[test$svtype=='DEL' & test$size_cate==size ,][,feature_list_oth]), rawscore = TRUE, reshape = TRUE)
  train_recali_3a=cbind(test[test$svtype=='DEL' & test$size_cate==size ,], preds_builtin)
  
  preds_builtin <- predict(models[[6]], data.matrix(test[test$svtype=='DUP' & test$size_cate==size ,][,feature_list_oth]), rawscore = TRUE, reshape = TRUE)
  train_recali_3b=cbind(test[test$svtype=='DUP' & test$size_cate==size ,], preds_builtin)
  
  size = 's4_5to50Kb'
  preds_builtin <- predict(models[[7]], data.matrix(test[test$svtype=='DEL' & test$size_cate==size ,][,feature_list_cnv]), rawscore = TRUE, reshape = TRUE)
  train_recali_4a=cbind(test[test$svtype=='DEL' & test$size_cate==size ,], preds_builtin)
  
  preds_builtin <- predict(models[[8]], data.matrix(test[test$svtype=='DUP' & test$size_cate==size ,][,feature_list_cnv]), rawscore = TRUE, reshape = TRUE)
  train_recali_4b=cbind(test[test$svtype=='DUP' & test$size_cate==size ,], preds_builtin)
  
  size = 's5_over50Kb'
  preds_builtin <- predict(models[[9]], data.matrix(test[test$svtype=='DEL' & test$size_cate==size ,][,feature_list_cnv]), rawscore = TRUE, reshape = TRUE)
  train_recali_5a=cbind(test[test$svtype=='DEL' & test$size_cate==size ,], preds_builtin)
  
  preds_builtin <- predict(models[[10]], data.matrix(test[test$svtype=='DUP' & test$size_cate==size ,][,feature_list_cnv]), rawscore = TRUE, reshape = TRUE)
  train_recali_5b=cbind(test[test$svtype=='DUP' & test$size_cate==size ,], preds_builtin)
  
  
  preds_builtin <- predict(models[[11]], data.matrix(test[test$svtype=='INS'   ,][,feature_list_oth]), rawscore = TRUE, reshape = TRUE)
  train_recali_1e=cbind(test[test$svtype=='INS'   ,], preds_builtin)
  
  preds_builtin <- predict(models[[12]], data.matrix(test[test$svtype%in%c('INS:ME:ALU') ,][,feature_list_oth]), rawscore = TRUE, reshape = TRUE)
  train_recali_1f=cbind(test[test$svtype%in%c('INS:ME:ALU') ,], preds_builtin)
  
  preds_builtin <- predict(models[[13]], data.matrix(test[test$svtype%in%c('INS:ME:LINE1','INS:ME:SVA','INS:ME') ,][,feature_list_oth]), rawscore = TRUE, reshape = TRUE)
  train_recali_1g=cbind(test[test$svtype%in%c('INS:ME:LINE1','INS:ME:SVA','INS:ME') ,], preds_builtin)
  
  train_recali_1h = test[test$svtype%in%c("INV",'CPX','CTX','CNV','BND') ,]
  train_recali_1h[,ncol(train_recali_1h)+1]=1
  colnames(train_recali_1h)[ncol(train_recali_1h)]='preds_builtin'
  
  out=rbind(train_recali_1a,train_recali_1b,
            train_recali_2a,train_recali_2b,
            train_recali_3a,train_recali_3b,
            train_recali_4a,train_recali_4b,
            train_recali_5a,train_recali_5b,
            train_recali_1e,train_recali_1f,
            train_recali_1g,train_recali_1h)
  
  out[,ncol(out)+1]=0
  out[out$preds_builtin>cff,][,ncol(out)]=1
  colnames(out)[ncol(out)]='preds_out'
  
  return(out)
}

option_list = list(
  make_option(c( "--site_feature_path"), type="character", default=NULL, help="folder including the site-annotations across all chromosomes", metavar="character"),
  make_option(c( "--sample_feature_path"), type="character", default=NULL, help="folder including the site-annotations across all chromosomes", metavar="character"),
  make_option(c( "--LD_feature_path"), type="character", default=NULL, help="folder including the site-annotations across all chromosomes", metavar="character"),
  make_option(c( "--lg_cnv_path"), type="character", default=NULL, help="folder including the site-annotations across all chromosomes", metavar="character"),
  make_option(c( "--model_path"), type="character", default=NULL, help="folder including all the pretrained lgb models", metavar="character"),
  make_option(c( "--input_path"), type="character", default=NULL, help="folder including the annotated per-sample SVs", metavar="character"),
  make_option(c( "--output_path"), type="character", default=NULL, help="output path", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

train_site_path = opt$site_feature_path
train_sample_path = opt$sample_feature_path
train_LD_path = opt$LD_feature_path
train_lg_cnv = opt$lg_cnv_path
train_model_path = opt$model_path
test_sample_path = opt$input_path
output_path = opt$output_path

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# read in per sv site features
site_feature = readin_site_features(train_site_path)
site_feature[,ncol(site_feature)+1] = site_feature$count_dnv/site_feature$count_children
colnames(site_feature)[ncol(site_feature)] = 'dnv_rate'

# read in per-sample features 
samp_feature = readin_sample_features(train_sample_path)

# combine site-level and sample-level features for a complete training dataset
training_data = readin_training(samp_feature, site_feature)

#readin SVID with high LD with SNVs
LD_SVID_afr = read.table(paste(train_LD_path, 'gnomadV2_vs_V3.AFR.SVID', sep = ''))
LD_SVID_eur = read.table(paste(train_LD_path, 'gnomadV2_vs_V3.EUR.SVID', sep = ''))
LD_SVID = data.frame(unique(rbind(LD_SVID_afr, LD_SVID_eur)[,1]))

#readin training set for large CNVs on SSC samples
ssc_training = readin_lg_cnv_training(train_lg_cnv)
ssc_hq = ssc_training[[1]]
ssc_lq = ssc_training[[2]]

train_site_path = paste(paste(strsplit(train_site_path,'/')[[1]],collapse = '/'),'/',sep = '')
test_sample_path = paste(paste(strsplit(test_sample_path,'/')[[1]],collapse = '/'),'/',sep = '')
train_model_path = paste(paste(strsplit(train_model_path,'/')[[1]],collapse = '/'),'/',sep = '')
output_path = paste(paste(strsplit(output_path,'/')[[1]],collapse = '/'),'/',sep = '')
#read in per-site features
site_feature = readin_site_features(train_site_path)

#train bgm model on testing samples:
feature_list_cnv = c('svtype','size_cate','dnv_rate','GT','CN','CNQ','EV','GQ','PE_GQ','PE_GT','RD_CN','RD_GQ','SR_GQ','SR_GT','vs_raw_manta_ovr1a','vs_raw_manta_ovr1b','vs_raw_wham_ovr1a','vs_raw_wham_ovr1b','vs_raw_melt_ovr1a','vs_raw_melt_ovr1b','vs_raw_depth_ovr1a','vs_raw_depth_ovr1b','ALGORITHMS','EVIDENCE','sr_background_fail','sr_bothside_support','PASS','PESR_GT_OVERDISPERSION','UNRESOLVED','MULTIALLELIC','GC')
feature_list_oth = c('svtype','size_cate','dnv_rate','GT','EV','GQ','PE_GQ','PE_GT','SR_GQ','SR_GT','vs_raw_manta_ovr1a','vs_raw_manta_ovr1b','vs_raw_wham_ovr1a','vs_raw_wham_ovr1b','vs_raw_melt_ovr1a','vs_raw_melt_ovr1b','ALGORITHMS','EVIDENCE','sr_background_fail','sr_bothside_support','PASS','PESR_GT_OVERDISPERSION','UNRESOLVED','MULTIALLELIC','GC')
bgm_models = train_del_dup_ins(training_data,feature_list_cnv, feature_list_oth, LD_SVID, ssc_hq, ssc_lq)
saveRDS.lgb.Booster(bgm_models[[1]], file = paste(output_path, 'bgm_model.DEL_s1_under250bp.rds', sep = ''), raw=TRUE)
saveRDS.lgb.Booster(bgm_models[[2]], file = paste(output_path, 'bgm_model.DUP_s1_under250bp.rds', sep = ''), raw=TRUE)
saveRDS.lgb.Booster(bgm_models[[3]], file = paste(output_path, 'bgm_model.DEL_s2_250bpto1Kb.rds', sep = ''), raw=TRUE)
saveRDS.lgb.Booster(bgm_models[[4]], file = paste(output_path, 'bgm_model.DUP_s2_250bpto1Kb.rds', sep = ''), raw=TRUE)
saveRDS.lgb.Booster(bgm_models[[5]], file = paste(output_path, 'bgm_model.DEL_s3_1to5Kb.rds', sep = ''), raw=TRUE)
saveRDS.lgb.Booster(bgm_models[[6]], file = paste(output_path, 'bgm_model.DUP_s3_1to5Kb.rds', sep = ''), raw=TRUE)
saveRDS.lgb.Booster(bgm_models[[7]], file = paste(output_path, 'bgm_model.DEL_s4_5to50Kb.rds', sep = ''), raw=TRUE)
saveRDS.lgb.Booster(bgm_models[[8]], file = paste(output_path, 'bgm_model.DUP_s4_5to50Kb.rds', sep = ''), raw=TRUE)
saveRDS.lgb.Booster(bgm_models[[9]], file = paste(output_path, 'bgm_model.DEL_s5_over50Kb.rds', sep = ''), raw=TRUE)
saveRDS.lgb.Booster(bgm_models[[10]], file = paste(output_path, 'bgm_model.DUP_s5_over50Kb.rds', sep = ''), raw=TRUE)
saveRDS.lgb.Booster(bgm_models[[11]], file = paste(output_path, 'bgm_model.INS.rds', sep = ''), raw=TRUE)
saveRDS.lgb.Booster(bgm_models[[12]], file = paste(output_path, 'bgm_model.ALU.rds', sep = ''), raw=TRUE)
saveRDS.lgb.Booster(bgm_models[[13]], file = paste(output_path, 'bgm_model.MEI.rds', sep = ''), raw=TRUE)



#apply model:
list_files = paste(test_sample_path,list.files(test_sample_path, recursive=TRUE), sep='')
for(file_name in list_files){
  print(file_name)
  test_data=read.table(file_name, header = T, comment.char = "")
  test_sample = readin_training(test_data, site_feature)
  test_results = test_del_dup_ins(bgm_models,test_sample,feature_list_cnv, feature_list_oth,-1)
  write.table(test_results[,c('SVID','preds_builtin')], paste(output_path, strsplit(strsplit(file_name,'/')[[1]][length(strsplit(file_name,'/')[[1]])],'[.]')[[1]][1],'.boost_filtered', sep=''), quote = F, sep = '\t', col.names = T, row.names = F)
}



