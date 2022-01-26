library(data.table)
library(lightgbm)
library(pacman)

readin_site_features<-function(chr_list){
  # read in the per-site traits across chromosomes
  # eg of chr_list = c("chr18", "chr19", "chr20", "chr21", "chr22")
  site_feature = read.table(paste('per_chr_anno/gnomad-sv-v3.',chr_list[1],'.final_cleanup.info4.gz', sep = ''), header = T, comment.char = "")
  for(chr in chr_list[2:length(chr_list)]){
    chr_feature = read.table(paste('per_chr_anno/gnomad-sv-v3.',chr,'.final_cleanup.info4.gz', sep = ''), header = T, comment.char = "")
    site_feature = rbind(site_feature, chr_feature)
  }
  return(site_feature)
}

readin_sample_features<-function(sample_list){
  # read in training samples
  # sample_list = c('__hg00513__58c7ec','__hg00514_1__264528','__hg00512__766970','__hg00514__2b7cb6', '__hg00731__2f5d02','__hg00732__a5867c','__hg00733_1__388f2a','__hg00733__835b40','__na19238__c76ef5','__na19239__901548','__na19240__a01b77')
  samp_feature = read.table(paste('per_sample_anno/',sample_list[1],'.bed.gz', sep=''), header = T, comment.char = "")
  samp_feature[,ncol(samp_feature)+1] = sample_list[1]
  colnames(samp_feature)[ncol(samp_feature)] = 'sample'
  for(sample in sample_list[2:length(sample_list)]){
    samp_feature_unit = read.table(paste('per_sample_anno/',sample,'.bed.gz', sep = ''), header = T, comment.char = "")
    samp_feature_unit[,ncol(samp_feature_unit)+1] = sample
    colnames(samp_feature_unit)[ncol(samp_feature_unit)] = 'sample'
    samp_feature = rbind(samp_feature, samp_feature_unit)
  }
  return(samp_feature)
}

readin_training<-function(samp_feature, site_feature){
  training_set = merge(samp_feature, site_feature, by='SVID')
  #modify size category:
  training_set$size_cate = 's1_under250bp'
  training_set[training_set$length>249,]$size_cate = 's2_250bpto1Kb'
  training_set[training_set$length>999,]$size_cate = 's3_1to5Kb'
  training_set[training_set$length>4999,]$size_cate = 's4_5to50Kb'
  training_set[training_set$length>49999,]$size_cate = 's5_over50Kb'
  
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

train_filter_model_small<-function(data, feature_list = feature_list_oth){
  #set up a new column in the traning table to indicate number of supportive evidences
  #hq - high quality;  lq -  low quality
  data[,ncol(data)+1] = 0
  colnames(data)[ncol(data)] = 'piece_of_hq_evidences'
  data[,ncol(data)+1] = 0
  colnames(data)[ncol(data)] = 'piece_of_lq_evidences'
  # login vapor support  
  data[!is.na(data$VaPoR_GT) & data$VaPoR_GT%in%c('0/1','1/1'),]$piece_of_hq_evidences = 1+data[!is.na(data$VaPoR_GT) & data$VaPoR_GT%in%c('0/1','1/1'),]$piece_of_hq_evidences 
  data[!is.na(data$VaPoR_GT) & data$VaPoR_GT%in%c('0/0'),      ]$piece_of_lq_evidences = 1+data[!is.na(data$VaPoR_GT) & data$VaPoR_GT%in%c('0/0'),      ]$piece_of_lq_evidences 
  # login PacBio support
  data[data$vs_pacbio_ovr1a=='OVR',]$piece_of_hq_evidences = 1 + data[data$vs_pacbio_ovr1a=='OVR',]$piece_of_hq_evidences
  data[data$vs_pacbio_ovr1a=='NO_OVR',]$piece_of_lq_evidences = 1 + data[data$vs_pacbio_ovr1a=='NO_OVR',]$piece_of_lq_evidences
  # login concordance support
  data[data$concor_duplicates>1,]$piece_of_hq_evidences = 1 + data[data$concor_duplicates>1,]$piece_of_hq_evidences
  data[data$concor_duplicates==1,]$piece_of_lq_evidences = 1 + data[data$concor_duplicates==1,]$piece_of_lq_evidences
  
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
                    , min_sum_hessian_in_leaf = 1.0 )
  
  model_builtin <- lgb.train(   params = params   , data = dtrain  , nrounds = 10L)
  
  tree_imp <- lgb.importance(model_builtin, percentage = TRUE)
  #lgb.plot.importance(tree_imp, top_n = 10L, measure = "Gain")
  
  return(model_builtin)
}

train_filter_model_medium<-function(data, feature_list = feature_list_cnv){
  #set up a new column in the traning table to indicate number of supportive evidences
  #hq - high quality;  lq -  low quality
  data[,ncol(data)+1] = 0
  colnames(data)[ncol(data)] = 'piece_of_hq_evidences'
  data[,ncol(data)+1] = 0
  colnames(data)[ncol(data)] = 'piece_of_lq_evidences'
  # login vapor support  
  data[!is.na(data$VaPoR_GT) & data$VaPoR_GT%in%c('0/1','1/1'),]$piece_of_hq_evidences = 1+data[!is.na(data$VaPoR_GT) & data$VaPoR_GT%in%c('0/1','1/1'),]$piece_of_hq_evidences 
  data[!is.na(data$VaPoR_GT) & data$VaPoR_GT%in%c('0/0'),      ]$piece_of_lq_evidences = 1+data[!is.na(data$VaPoR_GT) & data$VaPoR_GT%in%c('0/0'),      ]$piece_of_lq_evidences 
  # login PacBio support
  data[data$vs_pacbio_ovr1a=='OVR',]$piece_of_hq_evidences = 1 + data[data$vs_pacbio_ovr1a=='OVR',]$piece_of_hq_evidences
  data[data$vs_pacbio_ovr1a=='NO_OVR',]$piece_of_lq_evidences = 1 + data[data$vs_pacbio_ovr1a=='NO_OVR',]$piece_of_lq_evidences
  # login Bionano support
  data[data$vs_bionano_ovr1a=='OVR',]$piece_of_hq_evidences = 1 + data[data$vs_bionano_ovr1a=='OVR',]$piece_of_hq_evidences
  data[data$vs_bionano_ovr1a=='NO_OVR',]$piece_of_lq_evidences = 1 + data[data$vs_bionano_ovr1a=='NO_OVR',]$piece_of_lq_evidences
  # login concordance support
  data[data$concor_duplicates>1,]$piece_of_hq_evidences = 1 + data[data$concor_duplicates>1,]$piece_of_hq_evidences
  data[data$concor_duplicates==1,]$piece_of_lq_evidences = 1 + data[data$concor_duplicates==1,]$piece_of_lq_evidences
  
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
                    , min_sum_hessian_in_leaf = 1.0 )
  
  model_builtin <- lgb.train(   params = params   , data = dtrain  , nrounds = 10L)
  
  tree_imp <- lgb.importance(model_builtin, percentage = TRUE)
  #lgb.plot.importance(tree_imp, top_n = 10L, measure = "Gain")
  
  return(model_builtin)
}

train_filter_model_large<-function(data, feature_list = feature_list_cnv){
  #set up a new column in the traning table to indicate number of supportive evidences
  #hq - high quality;  lq -  low quality
  data[,ncol(data)+1] = 0
  colnames(data)[ncol(data)] = 'piece_of_hq_evidences'
  data[,ncol(data)+1] = 0
  colnames(data)[ncol(data)] = 'piece_of_lq_evidences'
  # login vapor support  
  data[!is.na(data$VaPoR_GT) & data$VaPoR_GT%in%c('0/1','1/1'),]$piece_of_hq_evidences = 1+data[!is.na(data$VaPoR_GT) & data$VaPoR_GT%in%c('0/1','1/1'),]$piece_of_hq_evidences 
  data[!is.na(data$VaPoR_GT) & data$VaPoR_GT%in%c('0/0'),      ]$piece_of_lq_evidences = 1+data[!is.na(data$VaPoR_GT) & data$VaPoR_GT%in%c('0/0'),      ]$piece_of_lq_evidences 
  # login PacBio support
  data[data$vs_pacbio_ovr1a=='OVR',]$piece_of_hq_evidences = 1 + data[data$vs_pacbio_ovr1a=='OVR',]$piece_of_hq_evidences
  data[data$vs_pacbio_ovr1a=='NO_OVR',]$piece_of_lq_evidences = 1 + data[data$vs_pacbio_ovr1a=='NO_OVR',]$piece_of_lq_evidences
  # login Bionano support
  data[data$vs_bionano_ovr1a=='OVR',]$piece_of_hq_evidences = 1 + data[data$vs_bionano_ovr1a=='OVR',]$piece_of_hq_evidences
  data[data$vs_bionano_ovr1a=='NO_OVR',]$piece_of_lq_evidences = 1 + data[data$vs_bionano_ovr1a=='NO_OVR',]$piece_of_lq_evidences
  # login concordance support
  data[data$concor_duplicates>1,]$piece_of_hq_evidences = 1 + data[data$concor_duplicates>1,]$piece_of_hq_evidences
  data[data$concor_duplicates==1,]$piece_of_lq_evidences = 1 + data[data$concor_duplicates==1,]$piece_of_lq_evidences
  
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
                    , min_sum_hessian_in_leaf = 1.0 )
  
  model_builtin <- lgb.train(   params = params   , data = dtrain  , nrounds = 10L)
  
  tree_imp <- lgb.importance(model_builtin, percentage = TRUE)
  #lgb.plot.importance(tree_imp, top_n = 10L, measure = "Gain")
  
  return(model_builtin)
}

train_filter_model_INS<-function(data, feature_list = feature_list_oth){
  #set up a new column in the traning table to indicate number of supportive evidences
  #hq - high quality;  lq -  low quality
  data[,ncol(data)+1] = 0
  colnames(data)[ncol(data)] = 'piece_of_hq_evidences'
  data[,ncol(data)+1] = 0
  colnames(data)[ncol(data)] = 'piece_of_lq_evidences'
  # login vapor support  
  data[!is.na(data$VaPoR_GT) & data$VaPoR_GT%in%c('0/1','1/1'),]$piece_of_hq_evidences = 1+data[!is.na(data$VaPoR_GT) & data$VaPoR_GT%in%c('0/1','1/1'),]$piece_of_hq_evidences 
  data[!is.na(data$VaPoR_GT) & data$VaPoR_GT%in%c('0/0'),      ]$piece_of_lq_evidences = 1+data[!is.na(data$VaPoR_GT) & data$VaPoR_GT%in%c('0/0'),      ]$piece_of_lq_evidences 
  # login PacBio support
  data[data$vs_pacbio_ovr1a=='OVR',]$piece_of_hq_evidences = 1 + data[data$vs_pacbio_ovr1a=='OVR',]$piece_of_hq_evidences
  data[data$vs_pacbio_ovr1a=='NO_OVR',]$piece_of_lq_evidences = 1 + data[data$vs_pacbio_ovr1a=='NO_OVR',]$piece_of_lq_evidences
  # login concordance support
  data[data$concor_duplicates>1,]$piece_of_hq_evidences = 1 + data[data$concor_duplicates>1,]$piece_of_hq_evidences
  data[data$concor_duplicates==1,]$piece_of_lq_evidences = 1 + data[data$concor_duplicates==1,]$piece_of_lq_evidences
  
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
                    , min_sum_hessian_in_leaf = 1.0 )
  
  model_builtin <- lgb.train(   params = params   , data = dtrain  , nrounds = 10L)
  
  tree_imp <- lgb.importance(model_builtin, percentage = TRUE)
  #lgb.plot.importance(tree_imp, top_n = 10L, measure = "Gain")
  
  return(model_builtin)
}

train_del_dup_ins<-function(train){
  
  feature_list_cnv = c('svtype.y','size_cate','GT','CN','CNQ','EV','GQ','PE_GQ','PE_GT','RD_CN','RD_GQ','SR_GQ','SR_GT','vs_raw_manta_ovr1a','vs_raw_manta_ovr1b','vs_raw_wham_ovr1a','vs_raw_wham_ovr1b','vs_raw_melt_ovr1a','vs_raw_melt_ovr1b','vs_raw_depth_ovr1a','vs_raw_depth_ovr1b','ALGORITHMS','EVIDENCE','sr_background_fail','sr_bothside_support','PASS','PESR_GT_OVERDISPERSION','UNRESOLVED','MULTIALLELIC','GC')
  feature_list_oth = c('svtype.y','size_cate','GT','EV','GQ','PE_GQ','PE_GT','RD_CN','RD_GQ','SR_GQ','SR_GT','vs_raw_manta_ovr1a','vs_raw_manta_ovr1b','vs_raw_wham_ovr1a','vs_raw_wham_ovr1b','vs_raw_melt_ovr1a','vs_raw_melt_ovr1b','ALGORITHMS','EVIDENCE','sr_background_fail','sr_bothside_support','PASS','PESR_GT_OVERDISPERSION','UNRESOLVED','MULTIALLELIC','GC')
  
  #feature_list_cnv=c("vs_raw_manta_ovr1a", "vs_raw_manta_ovr1b", "vs_raw_wham_ovr1a",  "vs_raw_wham_ovr1b",  "vs_raw_melt_ovr1a",  "vs_raw_melt_ovr1b","rd_median","rd_mean","rd_std","rd_median_le","rd_mean_le","rd_std_le","rd_median_ri","rd_mean_ri","rd_std_ri","PE_min","PE_max","SR_min","SR_max","size_cate","af_cate","ALGORITHMS","EVIDENCE","GC" ,"BothSideSupp","GT","GQ","RD_CN","RD_GQ","PE_GT","PE_GQ","SR_GT","SR_GQ","denovo_rate")
  #feature_list_oth=c("vs_raw_manta_ovr1a", "vs_raw_manta_ovr1b", "vs_raw_wham_ovr1a",  "vs_raw_wham_ovr1b",  "vs_raw_melt_ovr1a",  "vs_raw_melt_ovr1b","PE_min","PE_max","SR_min","SR_max","size_cate","af_cate","ALGORITHMS","EVIDENCE","GC" ,"BothSideSupp","GT","GQ","RD_CN","RD_GQ","PE_GT","PE_GQ","SR_GT","SR_GQ","denovo_rate")

  size = 's1_under250bp'
  train_recali_1a = train_filter_model_small(train[train$svtype.y=='DEL' & train$size_cate==size,],feature_list_oth)
  train_recali_1b = train_filter_model_small(train[train$svtype.y=='DUP' & train$size_cate==size,],feature_list_oth)
  
  size = 's2_250bpto1Kb'
  train_recali_2a = train_filter_model_small(train[train$svtype.y=='DEL' & train$size_cate==size,],feature_list_oth)
  train_recali_2b = train_filter_model_small(train[train$svtype.y=='DUP' & train$size_cate==size,],feature_list_oth)
  
  size = 's3_1to5Kb'
  train_recali_3a = train_filter_model_small(train[train$svtype.y=='DEL' & train$size_cate==size,],feature_list_oth)
  train_recali_3b = train_filter_model_small(train[train$svtype.y=='DUP' & train$size_cate==size,],feature_list_oth)
  
  size = 's4_5to50Kb'
  train_recali_4a = train_filter_model_medium(train[train$svtype.y=='DEL'& train$size_cate==size,], feature_list_cnv)
  train_recali_4b = train_filter_model_medium(train[train$svtype.y=='DUP'& train$size_cate==size ,], feature_list_cnv)
  
  size = 's5_over50Kb'
  train_recali_5a = train_filter_model_large(train[train$svtype.y=='DEL'& train$size_cate==size,], feature_list_cnv)
  train_recali_5b = train_filter_model_large(train[train$svtype.y=='DUP'& train$size_cate==size ,], feature_list_cnv)
  
  train_recali_1e = train_filter_model_INS(train[train$svtype.y=='INS',], feature_list_oth)
  train_recali_1f = train_filter_model_INS(train[train$svtype.y%in%c('INS:ME:ALU')  ,], feature_list_oth)
  train_recali_1g = train_filter_model_INS(train[train$svtype.y%in%c('INS:ME:LINE1','INS:ME:SVA','INS:ME') ,], feature_list_oth)
  
  out=list(train_recali_1a,train_recali_1b,
           train_recali_2a,train_recali_2b,
           train_recali_3a,train_recali_3b,
           train_recali_4a,train_recali_4b,
           train_recali_5a,train_recali_5b,
           train_recali_1e,train_recali_1f,
           train_recali_1g)
  
  return(out)
}

test_del_dup_ins<-function(bgm_models, test, cff=-0.8){
  models = bgm_models

  feature_list_cnv = c('svtype.y','size_cate','GT','CN','CNQ','EV','GQ','PE_GQ','PE_GT','RD_CN','RD_GQ','SR_GQ','SR_GT','vs_raw_manta_ovr1a','vs_raw_manta_ovr1b','vs_raw_wham_ovr1a','vs_raw_wham_ovr1b','vs_raw_melt_ovr1a','vs_raw_melt_ovr1b','vs_raw_depth_ovr1a','vs_raw_depth_ovr1b','ALGORITHMS','EVIDENCE','sr_background_fail','sr_bothside_support','PASS','PESR_GT_OVERDISPERSION','UNRESOLVED','MULTIALLELIC','GC')
  feature_list_oth = c('svtype.y','size_cate','GT','EV','GQ','PE_GQ','PE_GT','RD_CN','RD_GQ','SR_GQ','SR_GT','vs_raw_manta_ovr1a','vs_raw_manta_ovr1b','vs_raw_wham_ovr1a','vs_raw_wham_ovr1b','vs_raw_melt_ovr1a','vs_raw_melt_ovr1b','ALGORITHMS','EVIDENCE','sr_background_fail','sr_bothside_support','PASS','PESR_GT_OVERDISPERSION','UNRESOLVED','MULTIALLELIC','GC')
  
  size='s1_under250bp'
  preds_builtin <- predict(models[[1]], data.matrix(test[test$svtype.y=='DEL' & test$size_cate==size ,][,feature_list_oth]), rawscore = TRUE, reshape = TRUE)
  train_recali_1a=cbind(test[test$svtype.y=='DEL' & test$size_cate==size ,], preds_builtin)
  
  preds_builtin <- predict(models[[2]], data.matrix(test[test$svtype.y=='DUP' & test$size_cate==size ,][,feature_list_oth]), rawscore = TRUE, reshape = TRUE)
  train_recali_1b=cbind(test[test$svtype.y=='DUP' & test$size_cate==size ,], preds_builtin)
  
  size='s2_250bpto1Kb'
  preds_builtin <- predict(models[[3]], data.matrix(test[test$svtype.y=='DEL' & test$size_cate==size ,][,feature_list_oth]), rawscore = TRUE, reshape = TRUE)
  train_recali_2a=cbind(test[test$svtype.y=='DEL' & test$size_cate==size ,], preds_builtin)
  
  preds_builtin <- predict(models[[4]], data.matrix(test[test$svtype.y=='DUP' & test$size_cate==size ,][,feature_list_oth]), rawscore = TRUE, reshape = TRUE)
  train_recali_2b=cbind(test[test$svtype.y=='DUP' & test$size_cate==size ,], preds_builtin)
  
  size='s3_1to5Kb'
  preds_builtin <- predict(models[[5]], data.matrix(test[test$svtype.y=='DEL' & test$size_cate==size ,][,feature_list_oth]), rawscore = TRUE, reshape = TRUE)
  train_recali_3a=cbind(test[test$svtype.y=='DEL' & test$size_cate==size ,], preds_builtin)
  
  preds_builtin <- predict(models[[6]], data.matrix(test[test$svtype.y=='DUP' & test$size_cate==size ,][,feature_list_oth]), rawscore = TRUE, reshape = TRUE)
  train_recali_3b=cbind(test[test$svtype.y=='DUP' & test$size_cate==size ,], preds_builtin)
  
  size = 's4_5to50Kb'
  preds_builtin <- predict(models[[7]], data.matrix(test[test$svtype.y=='DEL' & test$size_cate==size ,][,feature_list_cnv]), rawscore = TRUE, reshape = TRUE)
  train_recali_4a=cbind(test[test$svtype.y=='DEL' & test$size_cate==size ,], preds_builtin)
  
  preds_builtin <- predict(models[[8]], data.matrix(test[test$svtype.y=='DUP' & test$size_cate==size ,][,feature_list_cnv]), rawscore = TRUE, reshape = TRUE)
  train_recali_4b=cbind(test[test$svtype.y=='DUP' & test$size_cate==size ,], preds_builtin)

  size = 's5_over50Kb'
  preds_builtin <- predict(models[[9]], data.matrix(test[test$svtype.y=='DEL' & test$size_cate==size ,][,feature_list_cnv]), rawscore = TRUE, reshape = TRUE)
  train_recali_5a=cbind(test[test$svtype.y=='DEL' & test$size_cate==size ,], preds_builtin)
  
  preds_builtin <- predict(models[[10]], data.matrix(test[test$svtype.y=='DUP' & test$size_cate==size ,][,feature_list_cnv]), rawscore = TRUE, reshape = TRUE)
  train_recali_5b=cbind(test[test$svtype.y=='DUP' & test$size_cate==size ,], preds_builtin)
  
  
  preds_builtin <- predict(models[[11]], data.matrix(test[test$svtype.y=='INS'   ,][,feature_list_oth]), rawscore = TRUE, reshape = TRUE)
  train_recali_1e=cbind(test[test$svtype.y=='INS'   ,], preds_builtin)
  
  preds_builtin <- predict(models[[12]], data.matrix(test[test$svtype.y%in%c('INS:ME:ALU') ,][,feature_list_oth]), rawscore = TRUE, reshape = TRUE)
  train_recali_1f=cbind(test[test$svtype.y%in%c('INS:ME:ALU') ,], preds_builtin)
  
  preds_builtin <- predict(models[[13]], data.matrix(test[test$svtype.y%in%c('INS:ME:LINE1','INS:ME:SVA','INS:ME') ,][,feature_list_oth]), rawscore = TRUE, reshape = TRUE)
  train_recali_1g=cbind(test[test$svtype.y%in%c('INS:ME:LINE1','INS:ME:SVA','INS:ME') ,], preds_builtin)
  
  train_recali_1h = test[test$svtype.y%in%c("INV",'CPX','CTX','CNV','BND') ,]
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
  roc_table_ins = calcu_ROC_by_BS_vs_hgsv(test_results[test_results$svtype.y=='INS',])
  roc_table_mei = calcu_ROC_by_BS_vs_hgsv(test_results[test_results$svtype.y%in%c('INS:ME','INS:ME:ALU','INS:ME:LINE1','INS:ME:SVA'),])
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
  roc_table_ins = calcu_ROC_by_BS_vs_PacBio(test_results[test_results$svtype.y=='INS',])
  roc_table_mei = calcu_ROC_by_BS_vs_PacBio(test_results[test_results$svtype.y%in%c('INS:ME','INS:ME:ALU','INS:ME:LINE1','INS:ME:SVA'),])
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
  
  par(mar=c(0,4,2,2))
  par(fig=c(0,1,.65,1))
  plot(c(0,length(bar_list_b)), range(bar_list), frame.plot = F, type = 'n', xlab = '', ylab = 'Count of SVs', las=2, xaxt='n')
  bar_width = .4
  x_pos = .5
  bar_colors = c('grey', 'black')
  for(i in c(1:length(bar_list))){
    rect(x_pos - bar_width, 0,  x_pos + bar_width, bar_list[i], col = bar_colors[1] )
    rect(x_pos - bar_width, 0,  x_pos + bar_width, bar_list_b[i], col = bar_colors[2] )
    x_pos = x_pos + 1
  }
  
  par(fig=c(0,1,.3,.65), new=T)
  plot(c(0,length(bar_list_b)), range(0,1.1), frame.plot = F, type = 'n', xlab = '', ylab = 'Proportion of SVs', las=2, xaxt='n')
  bar_width = .4
  x_pos = .5
  bar_colors = c('grey', 'black')
  for(i in c(1:length(bar_list))){
    rect(x_pos - bar_width, 0,  x_pos + bar_width, 1, col = bar_colors[1] )
    rect(x_pos - bar_width, 0,  x_pos + bar_width, bar_list_b[i]/bar_list[i], col = bar_colors[2] )
    text(x_pos, 1.05, bar_list[i], srt=45)
    x_pos = x_pos + 1
  }
  
  par(fig=c(0,1,0,.3), new=T)
  plot(c(0,length(bar_list_b)), c(0,4), frame.plot = F, type = 'n', xlab = '', ylab = '', xaxt='n', yaxt='n')
  axis(2,c(1:4)-.5, labels = c('concor','VaPoR','PB','1KGP'), las=2)
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
  plot_SVCounts_by_bs_score(test_results[test_results$svtype.y%in%c('INS:ME','INS:ME:ALU','INS:ME:LINE1','INS:ME:SVA'),], as.character(col_table[col_table[,1]=='MEI',2]), 20, 'Insertions')
  
}
  
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
setwd('~/Google Drive/My Drive/Talkowski_Lab/gnomAD_V3/module07/')
col_table=read.table('SV_colors.txt', comment.char = "")
col_table[nrow(col_table)+1,]=c('MEI','purple')
# read in per sv site features
training_chr_list = c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22', 'chrY')
site_feature = readin_site_features(training_chr_list)

# read in per-sample features 
training_sample_list = c('__hg00513__58c7ec','__hg00512__766970','__hg00514__2b7cb6', '__hg00731__2f5d02','__hg00732__a5867c','__hg00733__835b40','__na19238__c76ef5','__na19239__901548','__na19240__a01b77')
samp_feature = readin_sample_features(training_sample_list)

# combine site-level and sample-level features for a complete training dataset
training_data = readin_training(samp_feature, site_feature)

# train bgm model
bgm_models = train_del_dup_ins(training_data)

#apply bgm model on testing samples:
for(sample_name in c('__h_ij_hg00514_hg00514_1__d4092a', '__hg00514_1__264528', '__h_ij_hg00733_hg00733_2__c7cf50', '__hg00733_1__388f2a')){
  test_sample = readin_training(read.table(paste('per_sample_anno/',sample_name,'.bed.gz', sep=''), header = T, comment.char = ""), site_feature)
  test_results = test_del_dup_ins(bgm_models,test_sample,-1)
  write.table(test_results, paste('per_sample_anno/',sample_name,'.boost_filtered', sep=''), quote = F, sep = '\t', col.names = T, row.names = F)
  }

png('boost_qc_plots/count_SVs_by_bs.png')
plot_SVCounts_by_bs_score(test_results)
dev.off()

png('boost_qc_plots/ROC_all_sizes.png')
plot_roc_by_svtype_all_sizes_vs_hgsv(test_results)
dev.off()

plot_roc_by_svtype_all_sizes_vs_PacBio(test_results)

png('boost_qc_plots/ROC_by_size.png')
plot_roc_by_svtype_and_size_vs_hgsv(test_results)
dev.off()

for(bs_cff in c(-2,-1,0)){
  png(paste('boost_qc_plots/upsetplot','.bs',bs_cff,'.png', sep = ''))
  plot_upset_by_svtype_all_sizes(test_results, bs_cff)
  dev.off()
  for(svtype in c('DEL','DUP','INS')){
    print(c(svtype, bs_cff))
    png(paste('boost_qc_plots/upsetplot.',svtype,'.bs',bs_cff,'.png', sep = ''))
    plot_upset_by_svtype_all_sizes(test_results[test_results$svtype.x==svtype,], bs_cff)
    dev.off()
  }
}


hist(test_results$preds_builtin, xlab='boost score', main = 'distribution of boost score')
hist(test_results[test_results$vs_hgsv_ovr1a=="OVR",]$preds_builtin, col='black', add=T)
legend('topright',c('all SVs','overlap with 1KGP calls'), pch = 15, pt.cex = 2, bty = 'n', col = c('grey','black'))
table(test_results[,c('preds_out','vs_hgsv_ovr1a','svtype.x')])
table(test_results[,c('preds_out','vs_pacbio_ovr1a','svtype.x')])
table(test_results[,c('preds_out','concor_duplicates','svtype.x')])

lg_cnv_hq1=read.table('per_sample_anno/gnomad_calls.ssc_samples.vs_array_calls.bed.gz')
lg_cnv_hq2=read.table('per_sample_anno/gnomad_calls.ssc_samples.vs_exome_calls.bed.gz')
lg_cnv_lq = read.table('per_sample_anno/gnomad_calls.ssc_samples.original.vs_exome.wo_exome_calls.wo_array_calls.bed.gz')
stat1 = data.frame(table(lg_cnv_hq1[,4]))
stat2 = data.frame(table(lg_cnv_hq2[,4]))
stat3 = data.frame(table(lg_cnv_lq[,4]))
par(mfrow=c(2,1))
hist(stat1[,2], col='red', xlim = c(0,7000), main = 'AC of high qual hg CNVs')
hist(stat2[,2], col='blue', add=T, breaks = 20)
legend('topright', c('vs_array','vs_exome'), col = c('red', 'blue'), bty = 'n', pch = 15)
hist(stat3[,2], col='grey',  breaks = 50, xlim = c(0,7000), main = 'AC of high qual lg CNVs')





