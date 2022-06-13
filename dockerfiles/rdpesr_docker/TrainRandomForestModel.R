#！R
#Rscript to train random forest model

data_readin<-function(filename, inheri){
  dat=read.table(filename, header = T)
  #colnames(dat)[c(18:23)] = c('RD_CN','RD_GQ','PE_GT','PE_GQ','SR_GT','SR_GQ')
  for(i in c(24:29)){
    dat[,i] = as.character(dat[,i])
    dat[is.na(dat[,i]),][,i] ='no_info'
    dat[dat[,i]%in%c('0/1','1/1'),][,i]='info'
  }
  
  dat_dedun = cluster_sites(dat[,c("X.chrom", "start", "end", "SVID", "svtype")])
  dat=dat[dat$SVID%in%dat_dedun$SVID,]
  
  dat$EVIDENCE=as.character(dat$EVIDENCE)
  dat[is.na(dat$EVIDENCE),]$EVIDENCE = 'BAF'
  dat=dat[dat$svtype%in%c('CPX','DEL','DUP','INS','INV'),]
  dat[dat$PE_GQ=="None",]$PE_GQ = -1
  dat[dat$SR_GQ=="None",]$SR_GQ = -1
  dat$PE_GQ = as.integer(as.character(dat$PE_GQ))
  dat$SR_GQ = as.integer(as.character(dat$SR_GQ))
  dat$cleanvcf_vs_ref = as.character(dat$cleanvcf_vs_ref)
  dat[dat$cleanvcf_vs_ref!="NO_OVR",]$cleanvcf_vs_ref="OVR"
  dat=dat[,c("SVID","SVTYPE","ALGORITHMS","EVIDENCE","NCR","SVLEN","svtype" ,"cleanvcf_vs_ref","cleanvcf","EV","GQ","PE_GQ","PE_GT","RD_CN","RD_GQ","SR_GQ","SR_GT","GQ_reacali","boost_fix","boost_dyna","minGQ10.","preds_builtin")]
  
  #dat = merge(dat, bs, by='SVID')
  dat[dat$svtype=="INS" & !dat$preds_builtin>-.6,]$boost_fix="no_info"
  dat[dat$svtype=="DUP" & !dat$preds_builtin>-.5,]$boost_fix="no_info"
  dat[,ncol(dat)+1] = apply(dat,1, function(x){length(x[x=="info"])})
  colnames(dat)[ncol(dat)] = 'count_filters'
  
  data_messy = dat[dat$svtype=="DEL" & dat$ALGORITHMS=="wham" & dat$EV=="SR" & dat$SVLEN<1000,]
  data_clean = dat[!dat$SVID%in%data_messy$SVID, ]
  data_clean = data_clean[data_clean$count_filters>0,]
  data_clean[,ncol(data_clean)+1] = 'over5Kb'
  data_clean[data_clean$SVLEN<5000,][,ncol(data_clean)] = '1to5Kb'
  data_clean[data_clean$SVLEN<1000,][,ncol(data_clean)] = '100to1Kb'
  data_clean[data_clean$SVLEN<100,][,ncol(data_clean)] = 'under100bp'
  colnames(data_clean)[ncol(data_clean)] = 'size_cate'
  
  data_clean = merge(data_clean, inheri, by='SVID')
  data_clean[,ncol(data_clean)+1] = data_clean$count_dnv / data_clean$count_children
  colnames(data_clean)[ncol(data_clean)] = 'dnv_rate'
  data_clean[is.na(data_clean$dnv_rate),]$dnv_rate = 0
  return(data_clean)
}

vapor_readin<-function(filename){
  vapor=read.table(filename, header = T)
  vapor = vapor[vapor[,1]!="CHR",]
  vapor = vapor[!is.na(vapor$VaPoR_GT),]
  vapor$VaPoR_QS = as.double(vapor$VaPoR_QS)
  vapor$VaPoR_GS = as.double(vapor$VaPoR_GS)
  vapor$VaPoR_GQ = as.double(vapor$VaPoR_GQ)
  vapor$VaPoR_GT = as.character(vapor$VaPoR_GT)
  vapor[vapor$VaPoR_GS>0,]$VaPoR_GT = '0/1'
  return(vapor)
}

sample_readin<-function(sample, inheri){
  dat_sample = data_readin(paste('./per_sample_inte_vs_pb/per_sample_bed_integrated/',sample,'.vs_pb.inte.gz', sep = ''), inheri)
  vapor_sample = vapor_readin(paste('./vapor_results/cleanvcf_',sample,'.bed.gz', sep=''))
  dat_sample=merge(dat_sample, vapor_sample[,c("SVID", "VaPoR_GT")], by='SVID', all=T)
  dat_sample=dat_sample[!is.na(dat_sample$SVTYPE),]
  return(dat_sample)
}

train_data_readin<-function(train_samples, inheri){
  train_dat = sample_readin(train_samples[1], inheri)
  if(length(train_samples)>1){
    for(sample in train_samples[2:length(train_samples)]){
      tmp_dat = sample_readin(sample, inheri)
      train_dat = rbind(train_dat, tmp_dat)
    }
  }
  return(train_dat) 
}

train_model_on_PBSV<-function(traindata, svtype,size_cate){
  #train RF model on PacBio SV results
  colnames_to_remove = c('SVID','cleanvcf_vs_ref','VaPoR_GT')
  colnums = match(colnames_to_remove, colnames(traindata))
  if(size_cate!="all"){
    traindata6_dup<-traindata[traindata$SVTYPE==svtype & traindata$size_cate==size_cate  ,-c(colnums)]
    trainclass<-traindata[traindata$SVTYPE==svtype & traindata$size_cate==size_cate,c('cleanvcf_vs_ref')]
  }
  if(size_cate=="all"){
    traindata6_dup<-traindata[traindata$SVTYPE==svtype ,-c(colnums)]
    trainclass<-traindata[traindata$SVTYPE==svtype ,c('cleanvcf_vs_ref')]
  }
  trained_model = train(traindata6_dup, trainclass, method='rf')
  return(trained_model)
}

train_model_on_VaPoR<-function(traindata, svtype,size_cate){
  #train RF model on PacBio SV results
  traindata = traindata[!is.na(traindata$VaPoR_GT),]
  colnames_to_remove = c('SVID','cleanvcf_vs_ref','VaPoR_GT')
  colnums = match(colnames_to_remove, colnames(traindata))
  if(size_cate!="all"){
    traindata6_dup<-traindata[traindata$SVTYPE==svtype & traindata$size_cate==size_cate  ,-c(colnums)]
    trainclass<-traindata[traindata$SVTYPE==svtype & traindata$size_cate==size_cate,c('VaPoR_GT')]
  }
  if(size_cate=="all"){
    traindata6_dup<-traindata[traindata$SVTYPE==svtype ,-c(colnums)]
    trainclass<-traindata[traindata$SVTYPE==svtype ,c("VaPoR_GT")]
  }
  trainclass=as.character(trainclass)
  trainclass[trainclass=="1/1"]='0/1'
  trained_model = train(traindata6_dup, trainclass, method='rf')
  return(trained_model)
}

integrate_RF_results_PBSV<-function(train_test_del_sz1, cff){
  train_test_del_sz1[,ncol(train_test_del_sz1)+1] = 'NO_OVR'
  colnames(train_test_del_sz1)[ncol(train_test_del_sz1)] = "fitted"
  train_test_del_sz1[train_test_del_sz1$OVR>cff,]$fitted = 'OVR'
  return(train_test_del_sz1)
}

integrate_RF_results_VaPoR<-function(train_test_del_sz1, cff){
  train_test_del_sz1[,ncol(train_test_del_sz1)+1] = 'NO_OVR'
  colnames(train_test_del_sz1)[ncol(train_test_del_sz1)] = "fitted"
  train_test_del_sz1[train_test_del_sz1$`0/1`>cff,]$fitted = 'OVR'
  return(train_test_del_sz1)
}

apply_model_in_PBSV<-function(sample_name){
  test_1=sample_readin(sample_name,inheri)
  test_del_sz1 <- predict(PBSV_model_del_sz1, test_1[test_1$SVTYPE=="DEL" & test_1$size_cate=='under100bp',], type='prob')
  test_dup_sz1 <- predict(PBSV_model_dup_sz1, test_1[test_1$SVTYPE=="DUP" & test_1$size_cate=='under100bp',], type='prob')
  test_ins_sz1 <- predict(PBSV_model_ins_sz1, test_1[test_1$SVTYPE=="INS" & test_1$size_cate=='under100bp',], type='prob')
  test_del_sz2 <- predict(PBSV_model_del_sz2, test_1[test_1$SVTYPE=="DEL" & test_1$size_cate=='100to1Kb',], type='prob')
  test_dup_sz2 <- predict(PBSV_model_dup_sz2, test_1[test_1$SVTYPE=="DUP" & test_1$size_cate=='100to1Kb',], type='prob')
  test_ins_sz2 <- predict(PBSV_model_ins_sz2, test_1[test_1$SVTYPE=="INS" & test_1$size_cate=='100to1Kb',], type='prob')
  test_del_sz3 <- predict(PBSV_model_del_sz3, test_1[test_1$SVTYPE=="DEL" & test_1$size_cate=='1to5Kb',], type='prob')
  test_dup_sz3 <- predict(PBSV_model_dup_sz3, test_1[test_1$SVTYPE=="DUP" & test_1$size_cate=='1to5Kb',], type='prob')
  test_ins_sz3 <- predict(PBSV_model_ins_sz3, test_1[test_1$SVTYPE=="INS" & test_1$size_cate=='1to5Kb',], type='prob')
  test_del_sz4 <- predict(PBSV_model_del_sz4, test_1[test_1$SVTYPE=="DEL" & test_1$size_cate=='over5Kb',], type='prob')
  test_dup_sz4 <- predict(PBSV_model_dup_sz4, test_1[test_1$SVTYPE=="DUP" & test_1$size_cate=='over5Kb',], type='prob')
  test_ins_sz4 <- predict(PBSV_model_ins_sz4, test_1[test_1$SVTYPE=="INS" & test_1$size_cate=='over5Kb',], type='prob')
  test_alu <- predict(PBSV_model_ins_alu, test_1[test_1$SVTYPE=="INS:ME:ALU" ,], type='prob')
  test_l1 <- predict(PBSV_model_ins_l1, test_1[test_1$SVTYPE=="INS:ME:LINE1" ,], type='prob')
  test_sva <- predict(PBSV_model_ins_sva, test_1[test_1$SVTYPE=="INS:ME:SVA" ,], type='prob')
  train_test_del_sz1 = cbind(test_1[test_1$SVTYPE=="DEL" & test_1$size_cate=='under100bp',], test_del_sz1)
  train_test_dup_sz1 = cbind(test_1[test_1$SVTYPE=="DUP" & test_1$size_cate=='under100bp',], test_dup_sz1)
  train_test_ins_sz1 = cbind(test_1[test_1$SVTYPE=="INS" & test_1$size_cate=='under100bp',], test_ins_sz1)
  train_test_del_sz2 = cbind(test_1[test_1$SVTYPE=="DEL" & test_1$size_cate=='100to1Kb',], test_del_sz2)
  train_test_dup_sz2 = cbind(test_1[test_1$SVTYPE=="DUP" & test_1$size_cate=='100to1Kb',], test_dup_sz2)
  train_test_ins_sz2 = cbind(test_1[test_1$SVTYPE=="INS" & test_1$size_cate=='100to1Kb',], test_ins_sz2)
  train_test_del_sz3 = cbind(test_1[test_1$SVTYPE=="DEL" & test_1$size_cate=='1to5Kb',], test_del_sz3)
  train_test_dup_sz3 = cbind(test_1[test_1$SVTYPE=="DUP" & test_1$size_cate=='1to5Kb',], test_dup_sz3)
  train_test_ins_sz3 = cbind(test_1[test_1$SVTYPE=="INS" & test_1$size_cate=='1to5Kb',], test_ins_sz3)
  train_test_del_sz4 = cbind(test_1[test_1$SVTYPE=="DEL" & test_1$size_cate=='over5Kb',], test_del_sz4)
  train_test_dup_sz4 = cbind(test_1[test_1$SVTYPE=="DUP" & test_1$size_cate=='over5Kb',], test_dup_sz4)
  train_test_ins_sz4 = cbind(test_1[test_1$SVTYPE=="INS" & test_1$size_cate=='over5Kb',], test_ins_sz4)
  train_test_alu = cbind(test_1[test_1$SVTYPE=="INS:ME:ALU" ,], test_alu)
  train_test_l1 = cbind(test_1[test_1$SVTYPE=="INS:ME:LINE1" ,], test_l1)
  train_test_sva = cbind(test_1[test_1$SVTYPE=="INS:ME:SVA" ,], test_sva)
  
  pdf(paste('./per_sample_inte_vs_pb/plots/RF_predict_vs_PBSV_Supp', sample_name, 'pdf', sep='.'))
  par(mfrow=c(3,4))
  par(mar=rep(2,4))
  cff_list=c(.5,.4,.5,.5,.5,.5,.5,.5,.6,.5,.5,.5)
  boxplot(train_test_del_sz1$OVR~train_test_del_sz1$cleanvcf_vs_ref, frame.plot=F, main = 'DEL <100bp', ylim=c(0,1), col=c('blue','red'), pch=20, xaxt='n')
  abline(h=cff_list[1])
  train_test_del_sz1 = integrate_RF_results(train_test_del_sz1, cff_list[1])
  boxplot(train_test_del_sz2$OVR~train_test_del_sz2$cleanvcf_vs_ref, frame.plot=F, main = 'DEL 100bp-1Kb', ylim=c(0,1), col=c('blue','red'), pch=20, xaxt='n')
  abline(h=cff_list[2])
  train_test_del_sz2 = integrate_RF_results(train_test_del_sz2, cff_list[2])
  boxplot(train_test_del_sz3$OVR~train_test_del_sz3$cleanvcf_vs_ref, frame.plot=F, main = 'DEL 1-5Kb', ylim=c(0,1), col=c('blue','red'), pch=20, xaxt='n')
  abline(h=cff_list[3])
  train_test_del_sz3 = integrate_RF_results(train_test_del_sz3, cff_list[3])
  boxplot(train_test_del_sz4$OVR~train_test_del_sz4$cleanvcf_vs_ref, frame.plot=F, main = 'DEL >5Kb', ylim=c(0,1), col=c('blue','red'), pch=20, xaxt='n')
  abline(h=cff_list[4])
  train_test_del_sz4 = integrate_RF_results(train_test_del_sz4, cff_list[4])
  boxplot(train_test_dup_sz1$OVR~train_test_dup_sz1$cleanvcf_vs_ref, frame.plot=F, main = 'DUP <100bp', ylim=c(0,1), col=c('blue','red'), pch=20, xaxt='n')
  abline(h=cff_list[5])
  train_test_dup_sz1 = integrate_RF_results(train_test_dup_sz1, cff_list[5])
  boxplot(train_test_dup_sz2$OVR~train_test_dup_sz2$cleanvcf_vs_ref, frame.plot=F, main = 'DUP 100bp-1Kb', ylim=c(0,1), col=c('blue','red'), pch=20, xaxt='n')
  abline(h=cff_list[6])
  train_test_dup_sz2 = integrate_RF_results(train_test_dup_sz2, cff_list[6])
  boxplot(train_test_dup_sz3$OVR~train_test_dup_sz3$cleanvcf_vs_ref, frame.plot=F, main = 'DUP 1-5Kb', ylim=c(0,1), col=c('blue','red'), pch=20, xaxt='n')
  abline(h=cff_list[7])
  train_test_dup_sz3 = integrate_RF_results(train_test_dup_sz3, cff_list[7])
  boxplot(train_test_dup_sz4$OVR~train_test_dup_sz4$cleanvcf_vs_ref, frame.plot=F, main = 'DUP >5Kb', ylim=c(0,1), col=c('blue','red'), pch=20, xaxt='n')
  abline(h=cff_list[8])
  train_test_dup_sz4 = integrate_RF_results(train_test_dup_sz4, cff_list[8])
  boxplot(train_test_ins_sz1$OVR~train_test_ins_sz1$cleanvcf_vs_ref, frame.plot=F, main = 'INS <100bp', ylim=c(0,1), col=c('blue','red'), pch=20, xaxt='n')
  abline(h=cff_list[9])
  train_test_ins_sz1 = integrate_RF_results(train_test_ins_sz1, cff_list[9])
  boxplot(train_test_ins_sz2$OVR~train_test_ins_sz2$cleanvcf_vs_ref, frame.plot=F, main = 'INS 100bp-1Kb', ylim=c(0,1), col=c('blue','red'), pch=20, xaxt='n')
  abline(h=cff_list[10])
  train_test_ins_sz2 = integrate_RF_results(train_test_ins_sz2, cff_list[10])
  boxplot(train_test_ins_sz3$OVR~train_test_ins_sz3$cleanvcf_vs_ref, frame.plot=F, main = 'INS 1-5Kb', ylim=c(0,1), col=c('blue','red'), pch=20, xaxt='n')
  abline(h=cff_list[11])
  train_test_ins_sz3 = integrate_RF_results(train_test_ins_sz3, cff_list[11])
  boxplot(train_test_ins_sz4$OVR~train_test_ins_sz4$cleanvcf_vs_ref, frame.plot=F, main = 'INS >5Kb', ylim=c(0,1), col=c('blue','red'), pch=20, xaxt='n')
  abline(h=cff_list[12])
  train_test_ins_sz4 = integrate_RF_results(train_test_ins_sz4, cff_list[12])
  dev.off()
  
  train_test_alu = integrate_RF_results(train_test_alu, .5)
  train_test_l1 = integrate_RF_results(train_test_l1, .5)
  train_test_sva = integrate_RF_results(train_test_sva, .5)
  out_dat = rbind(train_test_del_sz1, train_test_del_sz2, train_test_del_sz3, train_test_del_sz4,
                  train_test_dup_sz1, train_test_dup_sz2, train_test_dup_sz3, train_test_dup_sz4,
                  train_test_ins_sz1, train_test_ins_sz2, train_test_ins_sz3, train_test_ins_sz4,
                  train_test_alu, train_test_l1, train_test_sva)
  return(out_dat)
}

apply_model_in_VaPoR<-function(sample_name){
  test_1=sample_readin(sample_name,inheri)
  test_del_sz1 <- predict(VaPoR_model_del_sz1, test_1[test_1$SVTYPE=="DEL" & test_1$size_cate=='under100bp',], type='prob')
  test_dup_sz1 <- predict(VaPoR_model_dup_sz1, test_1[test_1$SVTYPE=="DUP" & test_1$size_cate=='under100bp',], type='prob')
  test_ins_sz1 <- predict(VaPoR_model_ins_sz1, test_1[test_1$SVTYPE=="INS" & test_1$size_cate=='under100bp',], type='prob')
  test_del_sz2 <- predict(VaPoR_model_del_sz2, test_1[test_1$SVTYPE=="DEL" & test_1$size_cate=='100to1Kb',], type='prob')
  test_dup_sz2 <- predict(VaPoR_model_dup_sz2, test_1[test_1$SVTYPE=="DUP" & test_1$size_cate=='100to1Kb',], type='prob')
  test_ins_sz2 <- predict(VaPoR_model_ins_sz2, test_1[test_1$SVTYPE=="INS" & test_1$size_cate=='100to1Kb',], type='prob')
  test_del_sz3 <- predict(VaPoR_model_del_sz3, test_1[test_1$SVTYPE=="DEL" & test_1$size_cate=='1to5Kb',], type='prob')
  test_dup_sz3 <- predict(VaPoR_model_dup_sz3, test_1[test_1$SVTYPE=="DUP" & test_1$size_cate=='1to5Kb',], type='prob')
  test_ins_sz3 <- predict(VaPoR_model_ins_sz3, test_1[test_1$SVTYPE=="INS" & test_1$size_cate=='1to5Kb',], type='prob')
  test_del_sz4 <- predict(VaPoR_model_del_sz4, test_1[test_1$SVTYPE=="DEL" & test_1$size_cate=='over5Kb',], type='prob')
  test_dup_sz4 <- predict(VaPoR_model_dup_sz4, test_1[test_1$SVTYPE=="DUP" & test_1$size_cate=='over5Kb',], type='prob')
  test_ins_sz4 <- predict(VaPoR_model_ins_sz4, test_1[test_1$SVTYPE=="INS" & test_1$size_cate=='over5Kb',], type='prob')
  test_alu <- predict(VaPoR_model_ins_alu, test_1[test_1$SVTYPE=="INS:ME:ALU" ,], type='prob')
  test_l1 <- predict(VaPoR_model_ins_l1, test_1[test_1$SVTYPE=="INS:ME:LINE1" ,], type='prob')
  test_sva <- predict(VaPoR_model_ins_sva, test_1[test_1$SVTYPE=="INS:ME:SVA" ,], type='prob')

  train_test_del_sz1 = cbind(test_1[test_1$SVTYPE=="DEL" & test_1$size_cate=='under100bp',], test_del_sz1)
  train_test_dup_sz1 = cbind(test_1[test_1$SVTYPE=="DUP" & test_1$size_cate=='under100bp',], test_dup_sz1)
  train_test_ins_sz1 = cbind(test_1[test_1$SVTYPE=="INS" & test_1$size_cate=='under100bp',], test_ins_sz1)
  train_test_del_sz2 = cbind(test_1[test_1$SVTYPE=="DEL" & test_1$size_cate=='100to1Kb',], test_del_sz2)
  train_test_dup_sz2 = cbind(test_1[test_1$SVTYPE=="DUP" & test_1$size_cate=='100to1Kb',], test_dup_sz2)
  train_test_ins_sz2 = cbind(test_1[test_1$SVTYPE=="INS" & test_1$size_cate=='100to1Kb',], test_ins_sz2)
  train_test_del_sz3 = cbind(test_1[test_1$SVTYPE=="DEL" & test_1$size_cate=='1to5Kb',], test_del_sz3)
  train_test_dup_sz3 = cbind(test_1[test_1$SVTYPE=="DUP" & test_1$size_cate=='1to5Kb',], test_dup_sz3)
  train_test_ins_sz3 = cbind(test_1[test_1$SVTYPE=="INS" & test_1$size_cate=='1to5Kb',], test_ins_sz3)
  train_test_del_sz4 = cbind(test_1[test_1$SVTYPE=="DEL" & test_1$size_cate=='over5Kb',], test_del_sz4)
  train_test_dup_sz4 = cbind(test_1[test_1$SVTYPE=="DUP" & test_1$size_cate=='over5Kb',], test_dup_sz4)
  train_test_ins_sz4 = cbind(test_1[test_1$SVTYPE=="INS" & test_1$size_cate=='over5Kb',], test_ins_sz4)
  train_test_alu = cbind(test_1[test_1$SVTYPE=="INS:ME:ALU" ,], test_alu)
  train_test_l1 = cbind(test_1[test_1$SVTYPE=="INS:ME:LINE1" ,], test_l1)
  train_test_sva = cbind(test_1[test_1$SVTYPE=="INS:ME:SVA" ,], test_sva)
  
  cff_list=c(.5,.4,.5,.5,.6,.6,.4,.2,.5,.5,.2,.5)
 
 pdf(paste('./per_sample_inte_vs_pb/plots/RF_predict_vs_VaPoR_Supp', sample_name, 'pdf', sep='.'))
  par(mfrow=c(3,4))
  par(mar=rep(2,4))
  boxplot(train_test_del_sz1$`0/1`~train_test_del_sz1$VaPoR_GT, frame.plot=F, main = 'DEL <100bp', ylim=c(0,1), col=c('blue','red'), pch=20, xaxt='n')
  abline(h=cff_list[1])
  boxplot(train_test_del_sz2$`0/1`~train_test_del_sz2$VaPoR_GT, frame.plot=F, main = 'DEL 100bp-1Kb', ylim=c(0,1), col=c('blue','red'), pch=20, xaxt='n')
  abline(h=cff_list[2])
  boxplot(train_test_del_sz3$`0/1`~train_test_del_sz3$VaPoR_GT, frame.plot=F, main = 'DEL 1-5Kb', ylim=c(0,1), col=c('blue','red'), pch=20, xaxt='n')
  abline(h=cff_list[3])
  boxplot(train_test_del_sz4$`0/1`~train_test_del_sz4$VaPoR_GT, frame.plot=F, main = 'DEL >5Kb', ylim=c(0,1), col=c('blue','red'), pch=20, xaxt='n')
  abline(h=cff_list[4])
  boxplot(train_test_dup_sz1$`0/1`~train_test_dup_sz1$VaPoR_GT, frame.plot=F, main = 'DUP <100bp', ylim=c(0,1), col=c('blue','red'), pch=20, xaxt='n')
  abline(h=cff_list[5])
  boxplot(train_test_dup_sz2$`0/1`~train_test_dup_sz2$VaPoR_GT, frame.plot=F, main = 'DUP 100bp-1Kb', ylim=c(0,1), col=c('blue','red'), pch=20, xaxt='n')
  abline(h=cff_list[6])
  boxplot(train_test_dup_sz3$`0/1`~train_test_dup_sz3$VaPoR_GT, frame.plot=F, main = 'DUP 1-5Kb', ylim=c(0,1), col=c('blue','red'), pch=20, xaxt='n')
  abline(h=cff_list[7])
  boxplot(train_test_dup_sz4$`0/1`~train_test_dup_sz4$VaPoR_GT, frame.plot=F, main = 'DUP >5Kb', ylim=c(0,1), col=c('blue','red'), pch=20, xaxt='n')
  abline(h=cff_list[8])
  boxplot(train_test_ins_sz1$`0/1`~train_test_ins_sz1$VaPoR_GT, frame.plot=F, main = 'INS <100bp', ylim=c(0,1), col=c('blue','red'), pch=20, xaxt='n')
  abline(h=cff_list[9])
  boxplot(train_test_ins_sz2$`0/1`~train_test_ins_sz2$VaPoR_GT, frame.plot=F, main = 'INS 100bp-1Kb', ylim=c(0,1), col=c('blue','red'), pch=20, xaxt='n')
  abline(h=cff_list[10])
  boxplot(train_test_ins_sz3$`0/1`~train_test_ins_sz3$VaPoR_GT, frame.plot=F, main = 'INS 1-5Kb', ylim=c(0,1), col=c('blue','red'), pch=20, xaxt='n')
  abline(h=cff_list[11])
  boxplot(train_test_ins_sz4$`0/1`~train_test_ins_sz4$VaPoR_GT, frame.plot=F, main = 'INS >5Kb', ylim=c(0,1), col=c('blue','red'), pch=20, xaxt='n')
  abline(h=cff_list[12])
  dev.off()
  
  train_test_del_sz1 = integrate_RF_results_vapor(train_test_del_sz1, cff_list[1])
  train_test_del_sz2 = integrate_RF_results_vapor(train_test_del_sz2, cff_list[2])
  train_test_del_sz3 = integrate_RF_results_vapor(train_test_del_sz3, cff_list[3])
  train_test_del_sz4 = integrate_RF_results_vapor(train_test_del_sz4, cff_list[4])
  train_test_dup_sz1 = integrate_RF_results_vapor(train_test_dup_sz1, cff_list[5])
  train_test_dup_sz2 = integrate_RF_results_vapor(train_test_dup_sz2, cff_list[6])
  train_test_dup_sz3 = integrate_RF_results_vapor(train_test_dup_sz3, cff_list[7])
  train_test_dup_sz4 = integrate_RF_results_vapor(train_test_dup_sz4, cff_list[8])
  train_test_ins_sz1 = integrate_RF_results_vapor(train_test_ins_sz1, cff_list[9])
  train_test_ins_sz2 = integrate_RF_results_vapor(train_test_ins_sz2, cff_list[10])
  train_test_ins_sz3 = integrate_RF_results_vapor(train_test_ins_sz3, cff_list[11])
  train_test_ins_sz4 = integrate_RF_results_vapor(train_test_ins_sz4, cff_list[12])
  
  train_test_alu = integrate_RF_results_vapor(train_test_alu, .5)
  train_test_l1  = integrate_RF_results_vapor(train_test_l1, .5)
  train_test_sva = integrate_RF_results_vapor(train_test_sva, .5)
  out_dat = rbind(train_test_del_sz1, train_test_del_sz2, train_test_del_sz3, train_test_del_sz4,
                  train_test_dup_sz1, train_test_dup_sz2, train_test_dup_sz3, train_test_dup_sz4,
                  train_test_ins_sz1, train_test_ins_sz2, train_test_ins_sz3, train_test_ins_sz4,
                  train_test_alu, train_test_l1, train_test_sva)
  out_dat=out_dat[!is.na(out_dat$VaPoR_GT),]
  return(out_dat)
}

calcu_stat_PBSV<-function(test_1){
  tmp_stat=data.frame('svtype'=0,'size_cate'=0,'TP' = 0, 'FP' = 0,'FN'=0,'TN'=0, 'FPR' = 0)
  tmp_nrow = 1
  for(i in c('DEL', 'DUP', 'INS')){
    for(j in unique(test_1$size_cate)){
      tmp = test_1[test_1$SVTYPE==i & test_1$size_cate==j,]
      tmp_stat[tmp_nrow, c(1,2)] = c(i,j)
      tmp_stat[tmp_nrow, c(3:ncol(tmp_stat))] = c(
        nrow(tmp[tmp$fitted=="OVR" & tmp$cleanvcf_vs_ref=="OVR",]),
        nrow(tmp[tmp$fitted=="OVR" & tmp$cleanvcf_vs_ref=="NO_OVR",]),
        nrow(tmp[tmp$fitted=="NO_OVR" & tmp$cleanvcf_vs_ref=="OVR",]),
        nrow(tmp[tmp$fitted=="NO_OVR" & tmp$cleanvcf_vs_ref=="NO_OVR",]),
        nrow(tmp[tmp$fitted=="OVR" & tmp$cleanvcf_vs_ref=="NO_OVR",])/nrow(tmp[tmp$fitted=="OVR",]))
      tmp_nrow = tmp_nrow+1
    }
  }
  for(i in c('INS:ME:ALU','INS:ME:LINE1','INS:ME:SVA')){
    tmp = test_1[test_1$SVTYPE==i,]
    tmp_stat[tmp_nrow, c(1,2)] = c(i,'all')
    tmp_stat[tmp_nrow, c(3:ncol(tmp_stat))] = c(
      nrow(tmp[tmp$fitted=="OVR" & tmp$cleanvcf_vs_ref=="OVR",]),
      nrow(tmp[tmp$fitted=="OVR" & tmp$cleanvcf_vs_ref=="NO_OVR",]),
      nrow(tmp[tmp$fitted=="NO_OVR" & tmp$cleanvcf_vs_ref=="OVR",]),
      nrow(tmp[tmp$fitted=="NO_OVR" & tmp$cleanvcf_vs_ref=="NO_OVR",]),
      nrow(tmp[tmp$fitted=="OVR" & tmp$cleanvcf_vs_ref=="NO_OVR",])/nrow(tmp[tmp$fitted=="OVR",]))
    tmp_nrow = tmp_nrow+1
  }
  return(tmp_stat)
}

calcu_stat_VaPoR<-function(test_1){
  tmp_stat=data.frame('svtype'=0,'size_cate'=0,'TP' = 0, 'FP' = 0,'FN'=0,'TN'=0, 'FPR' = 0)
  tmp_nrow = 1
  for(i in c('DEL', 'DUP', 'INS')){
    for(j in unique(test_1$size_cate)){
      tmp = test_1[test_1$SVTYPE==i & test_1$size_cate==j,]
      tmp_stat[tmp_nrow, c(1,2)] = c(i,j)
      tmp_stat[tmp_nrow, c(3:ncol(tmp_stat))] = c(
        nrow(tmp[tmp$fitted=="OVR" & tmp$VaPoR_GT=="0/1",]),
        nrow(tmp[tmp$fitted=="OVR" & tmp$VaPoR_GT=="0/0",]),
        nrow(tmp[tmp$fitted=="NO_OVR" & tmp$VaPoR_GT=="0/1",]),
        nrow(tmp[tmp$fitted=="NO_OVR" & tmp$VaPoR_GT=="0/0",]),
        nrow(tmp[tmp$fitted=="OVR" & tmp$VaPoR_GT=="0/0",])/nrow(tmp[tmp$fitted=="OVR",]))
      tmp_nrow = tmp_nrow+1
    }
  }
  for(i in c('INS:ME:ALU','INS:ME:LINE1','INS:ME:SVA')){
    tmp = test_1[test_1$SVTYPE==i,]
    tmp_stat[tmp_nrow, c(1,2)] = c(i,'all')
    tmp_stat[tmp_nrow, c(3:ncol(tmp_stat))] = c(
      nrow(tmp[tmp$fitted=="OVR" & tmp$VaPoR_GT=="0/1",]),
      nrow(tmp[tmp$fitted=="OVR" & tmp$VaPoR_GT=="0/0",]),
      nrow(tmp[tmp$fitted=="NO_OVR" & tmp$VaPoR_GT=="0/1",]),
      nrow(tmp[tmp$fitted=="NO_OVR" & tmp$VaPoR_GT=="0/0",]),
      nrow(tmp[tmp$fitted=="OVR" & tmp$VaPoR_GT=="0/0",])/nrow(tmp[tmp$fitted=="OVR",]))
    tmp_nrow = tmp_nrow+1
  }
  return(tmp_stat)
}


#!/usr/bin/env Rscript
library("optparse")

option_list = list(
        make_option(c( "--inheri"), type="character", default=NULL, help="inheritance information", metavar="character"),
        make_option(c( "--genomic_context"), type="character", default=NULL, help="genomic context information", metavar="character"),
        make_option(c( "--train"), type="character", default=NULL, help="R project including training data", metavar="character"),
        make_option(c( "--output_PBSV"), type="character", default=NULL, help="name of PBSV trained model", metavar="character"),
        make_option(c( "--output_VaPoR"), type="character", default=NULL, help="name of VaPoR trained model", metavar="character"),
        make_option(c( "--svtype"), type="character", default=NULL, help="svtype", metavar="character"),
        make_option(c( "--size_cate"), type="character", default=NULL, help="size category", metavar="character")
 );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


library(caret)


inheri = read.table(opt$inheri, header = T)
gc=read.table(opt$genomic_context, header=T)
inheri = merge(inheri, gc, by='SVID')
load(opt$train)

PBSV_model = train_model_on_PBSV(train_dat, opt$svtype, opt$size_cate)

VaPoR_model = train_model_on_VaPoR(train_dat, opt$svtype, opt$size_cate)

saveRDS(PBSV_model, opt$output_PBSV)
saveRDS(VaPoR_model, opt$output_VaPoR)


