library("optparse")

option_list = list(
        make_option(c( "--input"), type="character", default=NULL, help="input file including all filtering features", metavar="character"),

        make_option(c( "--PBSV_model_del_sz1"), type="character", default=NULL, help="pre trained random forest models on deletions under 100bp", metavar="character"),
        make_option(c( "--PBSV_model_del_sz2"), type="character", default=NULL, help="pre trained random forest models on deletions between 100bp to 1Kb", metavar="character"),
        make_option(c( "--PBSV_model_del_sz3"), type="character", default=NULL, help="pre trained random forest models on deletions between 1Kb and 5Kb", metavar="character"),
        make_option(c( "--PBSV_model_del_sz4"), type="character", default=NULL, help="pre trained random forest models on deletions over5Kb", metavar="character"),
        make_option(c( "--PBSV_model_dup_sz1"), type="character", default=NULL, help="pre trained random forest models on duplications under 100bp", metavar="character"),
        make_option(c( "--PBSV_model_dup_sz2"), type="character", default=NULL, help="pre trained random forest models on duplications between 100bp to 1Kb", metavar="character"),
        make_option(c( "--PBSV_model_dup_sz3"), type="character", default=NULL, help="pre trained random forest models on duplications between 1Kb and 5Kb", metavar="character"),
        make_option(c( "--PBSV_model_dup_sz4"), type="character", default=NULL, help="pre trained random forest models on duplications over5Kb", metavar="character"),
        make_option(c( "--PBSV_model_ins_sz1"), type="character", default=NULL, help="pre trained random forest models on insertions under 100bp", metavar="character"),
        make_option(c( "--PBSV_model_ins_sz2"), type="character", default=NULL, help="pre trained random forest models on insertions between 100bp to 1Kb", metavar="character"),
        make_option(c( "--PBSV_model_ins_sz3"), type="character", default=NULL, help="pre trained random forest models on insertions between 1Kb and 5Kb", metavar="character"),
        make_option(c( "--PBSV_model_ins_sz4"), type="character", default=NULL, help="pre trained random forest models on insertions over5Kb", metavar="character"),
        make_option(c( "--PBSV_model_ins_alu"), type="character", default=NULL, help="pre trained random forest models on ALU", metavar="character"),
        make_option(c( "--PBSV_model_ins_l1"), type="character", default=NULL, help="pre trained random forest models on LINE1", metavar="character"),
        make_option(c( "--PBSV_model_ins_sva"), type="character", default=NULL, help="pre trained random forest models on SVA", metavar="character"),

        make_option(c( "--VaPoR_model_del_sz1"), type="character", default=NULL, help="pre trained random forest models on deletions under 100bp", metavar="character"),
        make_option(c( "--VaPoR_model_del_sz2"), type="character", default=NULL, help="pre trained random forest models on deletions between 100bp to 1Kb", metavar="character"),
        make_option(c( "--VaPoR_model_del_sz3"), type="character", default=NULL, help="pre trained random forest models on deletions between 1Kb and 5Kb", metavar="character"),
        make_option(c( "--VaPoR_model_del_sz4"), type="character", default=NULL, help="pre trained random forest models on deletions over5Kb", metavar="character"),
        make_option(c( "--VaPoR_model_dup_sz1"), type="character", default=NULL, help="pre trained random forest models on duplications under 100bp", metavar="character"),
        make_option(c( "--VaPoR_model_dup_sz2"), type="character", default=NULL, help="pre trained random forest models on duplications between 100bp to 1Kb", metavar="character"),
        make_option(c( "--VaPoR_model_dup_sz3"), type="character", default=NULL, help="pre trained random forest models on duplications between 1Kb and 5Kb", metavar="character"),
        make_option(c( "--VaPoR_model_dup_sz4"), type="character", default=NULL, help="pre trained random forest models on duplications over5Kb", metavar="character"),
        make_option(c( "--VaPoR_model_ins_sz1"), type="character", default=NULL, help="pre trained random forest models on insertions under 100bp", metavar="character"),
        make_option(c( "--VaPoR_model_ins_sz2"), type="character", default=NULL, help="pre trained random forest models on insertions between 100bp to 1Kb", metavar="character"),
        make_option(c( "--VaPoR_model_ins_sz3"), type="character", default=NULL, help="pre trained random forest models on insertions between 1Kb and 5Kb", metavar="character"),
        make_option(c( "--VaPoR_model_ins_sz4"), type="character", default=NULL, help="pre trained random forest models on insertions over5Kb", metavar="character"),
        make_option(c( "--VaPoR_model_ins_alu"), type="character", default=NULL, help="pre trained random forest models on ALU", metavar="character"),
        make_option(c( "--VaPoR_model_ins_l1"), type="character", default=NULL, help="pre trained random forest models on LINE1", metavar="character"),
        make_option(c( "--VaPoR_model_ins_sva"), type="character", default=NULL, help="pre trained random forest models on SVA", metavar="character"));

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

apply_test_data<-function(test_dat, PBSV_model_del_sz1,PBSV_model_del_sz2,PBSV_model_del_sz3,PBSV_model_del_sz4,PBSV_model_dup_sz1,PBSV_model_dup_sz2,PBSV_model_dup_sz3,PBSV_model_dup_sz4,PBSV_model_ins_sz1,PBSV_model_ins_sz2,PBSV_model_ins_sz3,PBSV_model_ins_sz4,PBSV_model_ins_alu,PBSV_model_ins_l1,PBSV_model_ins_sva)
{
  test_del_sz1 <- predict(PBSV_model_del_sz1, test_dat[test_dat$SVTYPE=="DEL" & test_dat$size_cate=='under100bp',], type='prob')
  test_del_sz2 <- predict(PBSV_model_del_sz2, test_dat[test_dat$SVTYPE=="DEL" & test_dat$size_cate=='100to1Kb',], type='prob')
  test_del_sz3 <- predict(PBSV_model_del_sz3, test_dat[test_dat$SVTYPE=="DEL" & test_dat$size_cate=='1to5Kb',], type='prob')
  test_del_sz4 <- predict(PBSV_model_del_sz4, test_dat[test_dat$SVTYPE=="DEL" & test_dat$size_cate=='over5Kb',], type='prob')
  
  test_dup_sz1 <- predict(PBSV_model_dup_sz1, test_dat[test_dat$SVTYPE=="DUP" & test_dat$size_cate=='under100bp',], type='prob')
  test_dup_sz2 <- predict(PBSV_model_dup_sz2, test_dat[test_dat$SVTYPE=="DUP" & test_dat$size_cate=='100to1Kb',], type='prob')
  test_dup_sz3 <- predict(PBSV_model_dup_sz3, test_dat[test_dat$SVTYPE=="DUP" & test_dat$size_cate=='1to5Kb',], type='prob')
  test_dup_sz4 <- predict(PBSV_model_dup_sz4, test_dat[test_dat$SVTYPE=="DUP" & test_dat$size_cate=='over5Kb',], type='prob')
  
  test_ins_sz1 <- predict(PBSV_model_ins_sz1, test_dat[test_dat$SVTYPE=="INS" & test_dat$size_cate=='under100bp',], type='prob')
  test_ins_sz2 <- predict(PBSV_model_ins_sz2, test_dat[test_dat$SVTYPE=="INS" & test_dat$size_cate=='100to1Kb',], type='prob')
  test_ins_sz3 <- predict(PBSV_model_ins_sz3, test_dat[test_dat$SVTYPE=="INS" & test_dat$size_cate=='1to5Kb',], type='prob')
  test_ins_sz4 <- predict(PBSV_model_ins_sz4, test_dat[test_dat$SVTYPE=="INS" & test_dat$size_cate=='over5Kb',], type='prob')
  
  test_alu <- predict(PBSV_model_ins_alu, test_dat[test_dat$SVTYPE=="INS:ME:ALU" ,], type='prob')
  test_l1 <- predict(PBSV_model_ins_l1, test_dat[test_dat$SVTYPE=="INS:ME:LINE1" ,], type='prob')
  test_sva <- predict(PBSV_model_ins_sva, test_dat[test_dat$SVTYPE=="INS:ME:SVA" ,], type='prob')
  
  
  train_test_del_sz1 = cbind(test_dat[test_dat$SVTYPE=="DEL" & test_dat$size_cate=='under100bp',], test_del_sz1)
  train_test_dup_sz1 = cbind(test_dat[test_dat$SVTYPE=="DUP" & test_dat$size_cate=='under100bp',], test_dup_sz1)
  train_test_ins_sz1 = cbind(test_dat[test_dat$SVTYPE=="INS" & test_dat$size_cate=='under100bp',], test_ins_sz1)
  train_test_del_sz2 = cbind(test_dat[test_dat$SVTYPE=="DEL" & test_dat$size_cate=='100to1Kb',], test_del_sz2)
  train_test_dup_sz2 = cbind(test_dat[test_dat$SVTYPE=="DUP" & test_dat$size_cate=='100to1Kb',], test_dup_sz2)
  train_test_ins_sz2 = cbind(test_dat[test_dat$SVTYPE=="INS" & test_dat$size_cate=='100to1Kb',], test_ins_sz2)
  train_test_del_sz3 = cbind(test_dat[test_dat$SVTYPE=="DEL" & test_dat$size_cate=='1to5Kb',], test_del_sz3)
  train_test_dup_sz3 = cbind(test_dat[test_dat$SVTYPE=="DUP" & test_dat$size_cate=='1to5Kb',], test_dup_sz3)
  train_test_ins_sz3 = cbind(test_dat[test_dat$SVTYPE=="INS" & test_dat$size_cate=='1to5Kb',], test_ins_sz3)
  train_test_del_sz4 = cbind(test_dat[test_dat$SVTYPE=="DEL" & test_dat$size_cate=='over5Kb',], test_del_sz4)
  train_test_dup_sz4 = cbind(test_dat[test_dat$SVTYPE=="DUP" & test_dat$size_cate=='over5Kb',], test_dup_sz4)
  train_test_ins_sz4 = cbind(test_dat[test_dat$SVTYPE=="INS" & test_dat$size_cate=='over5Kb',], test_ins_sz4)
  train_test_alu = cbind(test_dat[test_dat$SVTYPE=="INS:ME:ALU" ,], test_alu)
  train_test_l1 = cbind(test_dat[test_dat$SVTYPE=="INS:ME:LINE1" ,], test_l1)
  train_test_sva = cbind(test_dat[test_dat$SVTYPE=="INS:ME:SVA" ,], test_sva)
  
  out_dat = rbind(train_test_del_sz1, train_test_del_sz2, train_test_del_sz3, train_test_del_sz4,
                  train_test_dup_sz1, train_test_dup_sz2, train_test_dup_sz3, train_test_dup_sz4,
                  train_test_ins_sz1, train_test_ins_sz2, train_test_ins_sz3, train_test_ins_sz4,
                  train_test_alu, train_test_l1, train_test_sva)
 
  return(out_dat) 
}


PBSV_model_del_sz1 = readRDS(opt$PBSV_model_del_sz1)
PBSV_model_del_sz2 = readRDS(opt$PBSV_model_del_sz2)
PBSV_model_del_sz3 = readRDS(opt$PBSV_model_del_sz3)
PBSV_model_del_sz4 = readRDS(opt$PBSV_model_del_sz4)

PBSV_model_dup_sz1 = readRDS(opt$PBSV_model_dup_sz1)
PBSV_model_dup_sz2 = readRDS(opt$PBSV_model_dup_sz2)
PBSV_model_dup_sz3 = readRDS(opt$PBSV_model_dup_sz3)
PBSV_model_dup_sz4 = readRDS(opt$PBSV_model_dup_sz4)

PBSV_model_ins_sz1 = readRDS(opt$PBSV_model_ins_sz1)
PBSV_model_ins_sz2 = readRDS(opt$PBSV_model_ins_sz2)
PBSV_model_ins_sz3 = readRDS(opt$PBSV_model_ins_sz3)
PBSV_model_ins_sz4 = readRDS(opt$PBSV_model_ins_sz4)

PBSV_model_ins_alu = readRDS(opt$PBSV_model_ins_alu)
PBSV_model_ins_l1 = readRDS(opt$PBSV_model_ins_l1)
PBSV_model_ins_sva = readRDS(opt$PBSV_model_ins_sva)


VaPoR_model_del_sz1 = readRDS(opt$VaPoR_model_del_sz1)
VaPoR_model_del_sz2 = readRDS(opt$VaPoR_model_del_sz2)
VaPoR_model_del_sz3 = readRDS(opt$VaPoR_model_del_sz3)
VaPoR_model_del_sz4 = readRDS(opt$VaPoR_model_del_sz4)

VaPoR_model_dup_sz1 = readRDS(opt$VaPoR_model_dup_sz1)
VaPoR_model_dup_sz2 = readRDS(opt$VaPoR_model_dup_sz2)
VaPoR_model_dup_sz3 = readRDS(opt$VaPoR_model_dup_sz3)
VaPoR_model_dup_sz4 = readRDS(opt$VaPoR_model_dup_sz4)

VaPoR_model_ins_sz1 = readRDS(opt$VaPoR_model_ins_sz1)
VaPoR_model_ins_sz2 = readRDS(opt$VaPoR_model_ins_sz2)
VaPoR_model_ins_sz3 = readRDS(opt$VaPoR_model_ins_sz3)
VaPoR_model_ins_sz4 = readRDS(opt$VaPoR_model_ins_sz4)

VaPoR_model_ins_alu = readRDS(opt$VaPoR_model_ins_alu)
VaPoR_model_ins_l1 = readRDS(opt$VaPoR_model_ins_l1)
VaPoR_model_ins_sva = readRDS(opt$VaPoR_model_ins_sva)

input_list = read.table(opt$input)
input_list[,2] = apply(input_list,1,function(x){tail(strsplit(as.character(x[1]),'[.]')[[1]],3)[1]})
for(i in 1:nrow(input_list)){
  test_dat=read.table(as.character(input_list[i,1]), header = T)
  test_dat[test_dat$SVTYPE=="INS:ME",]$SVTYPE="INS"
  pbsv_test = apply_test_data(test_dat, PBSV_model_del_sz1,PBSV_model_del_sz2,PBSV_model_del_sz3,PBSV_model_del_sz4,PBSV_model_dup_sz1,PBSV_model_dup_sz2,PBSV_model_dup_sz3,PBSV_model_dup_sz4,PBSV_model_ins_sz1,PBSV_model_ins_sz2,PBSV_model_ins_sz3,PBSV_model_ins_sz4,PBSV_model_ins_alu,PBSV_model_ins_l1,PBSV_model_ins_sva)
  vapor_test = apply_test_data(test_dat, VaPoR_model_del_sz1,VaPoR_model_del_sz2,VaPoR_model_del_sz3,VaPoR_model_del_sz4,VaPoR_model_dup_sz1,VaPoR_model_dup_sz2,VaPoR_model_dup_sz3,VaPoR_model_dup_sz4,VaPoR_model_ins_sz1,VaPoR_model_ins_sz2,VaPoR_model_ins_sz3,VaPoR_model_ins_sz4,VaPoR_model_ins_alu,VaPoR_model_ins_l1,VaPoR_model_ins_sva)
  out = merge(pbsv_test[,c('SVID','NO_OVR','OVR')], vapor_test[,c('SVID','0/0','0/1')], by='SVID')
  colnames(out)=c('SVID','PBSV_neg','PBSV_pos','VaPoR_neg','VaPoR_pos')
  write.table(out, paste('RF_results.', as.character(input_list[i,2]), '.tsv', sep=''), quote=F, sep = '\t', col.names = T, row.names = F)

}



