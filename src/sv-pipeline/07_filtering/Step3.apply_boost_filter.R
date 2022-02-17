#!/usr/bin/env Rscript
#script to apply the pre-trained boost model on each sample

library(data.table)
library(lightgbm)
library(pacman)
library(optparse)

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
  make_option(c( "--model_path"), type="character", default=NULL, help="folder including all the pretrained lgb models", metavar="character"),
  make_option(c( "--input_path"), type="character", default=NULL, help="folder including the annotated per-sample SVs", metavar="character")
  
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

train_site_path = opt$site_feature_path
test_sample_path = opt$input_path
train_model_path = opt$model_path

train_site_path = paste(paste(strsplit(train_site_path,'/')[[1]],collapse = '/'),'/',sep = '')
test_sample_path = paste(paste(strsplit(test_sample_path,'/')[[1]],collapse = '/'),'/',sep = '')
train_model_path = paste(paste(strsplit(train_model_path,'/')[[1]],collapse = '/'),'/',sep = '')
#read in per-site features
site_feature = readin_site_features(train_site_path)

#apply bgm model on testing samples:
feature_list_cnv = c('svtype','size_cate','GT','CN','CNQ','EV','GQ','PE_GQ','PE_GT','RD_CN','RD_GQ','SR_GQ','SR_GT','vs_raw_manta_ovr1a','vs_raw_manta_ovr1b','vs_raw_wham_ovr1a','vs_raw_wham_ovr1b','vs_raw_melt_ovr1a','vs_raw_melt_ovr1b','vs_raw_depth_ovr1a','vs_raw_depth_ovr1b','ALGORITHMS','EVIDENCE','sr_background_fail','sr_bothside_support','PASS','PESR_GT_OVERDISPERSION','UNRESOLVED','MULTIALLELIC','GC')
feature_list_oth = c('svtype','size_cate','GT','EV','GQ','PE_GQ','PE_GT','SR_GQ','SR_GT','vs_raw_manta_ovr1a','vs_raw_manta_ovr1b','vs_raw_wham_ovr1a','vs_raw_wham_ovr1b','vs_raw_melt_ovr1a','vs_raw_melt_ovr1b','ALGORITHMS','EVIDENCE','sr_background_fail','sr_bothside_support','PASS','PESR_GT_OVERDISPERSION','UNRESOLVED','MULTIALLELIC','GC')

#readin lgb models
bgm_models = list(readRDS.lgb.Booster(paste(train_model_path,'bgm_model.DEL_s1_under250bp.rds',sep = '')),
                  readRDS.lgb.Booster(paste(train_model_path,'bgm_model.DUP_s1_under250bp.rds',sep = '')),
                  readRDS.lgb.Booster(paste(train_model_path,'bgm_model.DEL_s2_250bpto1Kb.rds',sep = '')),
                  readRDS.lgb.Booster(paste(train_model_path,'bgm_model.DUP_s2_250bpto1Kb.rds',sep = '')),
                  readRDS.lgb.Booster(paste(train_model_path,'bgm_model.DEL_s3_1to5Kb.rds',sep = '')),
                  readRDS.lgb.Booster(paste(train_model_path,'bgm_model.DUP_s3_1to5Kb.rds',sep = '')),
                  readRDS.lgb.Booster(paste(train_model_path,'bgm_model.DEL_s4_5to50Kb.rds',sep = '')),
                  readRDS.lgb.Booster(paste(train_model_path,'bgm_model.DUP_s4_5to50Kb.rds',sep = '')),
                  readRDS.lgb.Booster(paste(train_model_path,'bgm_model.DEL_s5_over50Kb.rds',sep = '')),
                  readRDS.lgb.Booster(paste(train_model_path,'bgm_model.DUP_s5_over50Kb.rds',sep = '')),
                  readRDS.lgb.Booster(paste(train_model_path,'bgm_model.INS.rds',sep = '')),
                  readRDS.lgb.Booster(paste(train_model_path,'bgm_model.ALU.rds',sep = '')),
                  readRDS.lgb.Booster(paste(train_model_path,'bgm_model.MEI.rds',sep = '')))

#apply model:
list_files = paste(test_sample_path,list.files(test_sample_path), sep='')
for(file_name in list_files){
  print(file_name)
  test_data=read.table(file_name, header = T, comment.char = "")
  test_sample = readin_training(test_data, site_feature)
  test_results = test_del_dup_ins(bgm_models,test_sample,feature_list_cnv, feature_list_oth,-1)
  write.table(test_results[,c('SVID','preds_builtin')], paste(strsplit(file_name,'[.]')[[1]][1],'.boost_filtered', sep=''), quote = F, sep = '\t', col.names = T, row.names = F)
}



