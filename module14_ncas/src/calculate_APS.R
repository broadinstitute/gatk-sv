add_size_cate<-function(dat){
  sv_size_seg = c(c(5:10)*10, c(1:50)*20,c(11:50)*100,c(6:20)*1000,c(3:10)*10000, c(2:10)*100000)
  dat[,ncol(dat)+1] = 'over5Kb'
  colnames(dat)[ncol(dat)] = 'size_cate'
  dat[dat$SVLEN<5000,]$size_cate = '400bpto5Kb'
  dat[dat$SVLEN<400,]$size_cate = 'under400bp'
  return(dat)
}
add_size_decile_and_ps_by_all_sizes<-function(data, cate_name, decile = 100){
  size = data[,c("name","SVLEN")]
  size = size[order(size[,2]),]
  size[,3] = c(1:nrow(size))
  size[,4] = as.integer(c(1:nrow(size))/(nrow(size)+1)*decile)/decile
  colnames(size)[c(3,4)]=c('size_rank','size_decile')
  data = merge(data, size[,-2], by = c('name'))
  data[,ncol(data)+1] = 0
  colnames(data)[ncol(data)] = 'singleton_prop'

  for(i in unique(data$size_decile)){
    singleton_prop = nrow(data[data$size_decile==i & data$AC==1,]) / nrow(data[data$size_decile==i,])
    data[data$size_decile==i,]$singleton_prop = singleton_prop
  }
  out = data[,c("name",'size_rank', "size_decile", 'singleton_prop')]
  out[,ncol(out)+1] = cate_name
  colnames(out)[ncol(out)] = 'aps_cate'
  return(out)
}
add_size_decile_and_ps_by_unique_sizes<-function(data, cate_name, decile = 100){
  size = data[,c("name","SVLEN")]
  size = size[order(size[,2]),]
  unique_size = unique(size[,c("SVLEN","SVLEN")])
  unique_size[,2] = c(1:nrow(unique_size))/(nrow(unique_size)+1)
  unique_size[,3] = as.integer(unique_size[,2]*decile)/decile
  colnames(unique_size) = c("SVLEN",'size_rank', 'size_decile')
  size = merge(size, unique_size, by='SVLEN')
  data = merge(data, size, by = c('name'))
  data[,ncol(data)+1] = 0
  colnames(data)[ncol(data)] = 'singleton_prop'
  
  for(i in unique(data$size_decile)){
    singleton_prop = nrow(data[data$size_decile==i & data$AC==1,]) / nrow(data[data$size_decile==i,])
    data[data$size_decile==i,]$singleton_prop = singleton_prop
  }
  out = data[,c("name",'size_rank', "size_decile", 'singleton_prop')]
  out[,ncol(out)+1] = cate_name
  colnames(out)[ncol(out)] = 'aps_cate'
  return(out)
}
add_ps_annotations<-function(train_data,test_data){
  pass_ps_size_ranges = calculate_size_ranges(train_data)
  test_data = test_data[!is.na(test_data$aps_cate),]
  
  out_data = test_data[1,]
  for(i in unique(test_data$aps_cate)){
    if(!is.na(i)){
      print(c(i))
      tmp_data = test_data[test_data$aps_cate==i,]
      tmp = pass_ps_size_ranges[pass_ps_size_ranges$aps_cate==i,]
      tmp = tmp[order(tmp$size_decile),]
      for(j in 1:nrow(tmp)){
        tmp_data[!tmp_data$SVLEN<tmp[j,]$min_size,]$size_decile = tmp[j,]$size_decile
      }
      out_data = rbind(out_data, tmp_data)
    }
  }
  out_data = out_data[-1,]
  return(out_data)
}
aps_table.by_gene_list <-function(pass_aps_intragenic_Testing, gene_list){
  
  out_table = data.frame('func'=0,'aps'=0)
  row_num = 1
  out_table[row_num,1] = 'lof'
  out_table[row_num,2] = mean(pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$LOF) & pass_aps_intragenic_Testing$LOF%in%gene_list,]$APS)
  row_num = row_num+1
  out_table[row_num,1] = 'lof_del'
  out_table[row_num,2] = mean(pass_aps_intragenic_Testing[pass_aps_intragenic_Testing$SVTYPE=="DEL" & !is.na(pass_aps_intragenic_Testing$LOF) & pass_aps_intragenic_Testing$LOF%in%gene_list,]$APS)
  row_num = row_num+1
  out_table[row_num,1] = 'lof_ins'
  out_table[row_num,2] = mean(pass_aps_intragenic_Testing[pass_aps_intragenic_Testing$SVTYPE=="INS" & !is.na(pass_aps_intragenic_Testing$LOF) & pass_aps_intragenic_Testing$LOF%in%gene_list,]$APS)
  row_num = row_num+1
  out_table[row_num,1] = 'lof_dup'
  out_table[row_num,2] = mean(pass_aps_intragenic_Testing[pass_aps_intragenic_Testing$SVTYPE=="DUP" & !is.na(pass_aps_intragenic_Testing$LOF) & pass_aps_intragenic_Testing$LOF%in%gene_list,]$APS)
  row_num = row_num+1
  out_table[row_num,1] = 'lof_oth'
  out_table[row_num,2] = mean(pass_aps_intragenic_Testing[!pass_aps_intragenic_Testing$SVTYPE%in%c("DEL",'DUP','INS') & !is.na(pass_aps_intragenic_Testing$LOF) & pass_aps_intragenic_Testing$LOF%in%gene_list,]$APS)
  row_num = row_num+1
  out_table[row_num,1] = 'cg'
  out_table[row_num,2] = mean(pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$COPY_GAIN) & pass_aps_intragenic_Testing$COPY_GAIN%in%gene_list,]$APS)
  row_num = row_num+1
  out_table[row_num,1] = 'ied'
  out_table[row_num,2] = mean(pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$INTRAGENIC_EXON_DUP) & pass_aps_intragenic_Testing$INTRAGENIC_EXON_DUP%in%gene_list,]$APS)
  row_num = row_num+1
  out_table[row_num,1] = 'ped'
  out_table[row_num,2] = mean(pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PARTIAL_EXON_DUP) & pass_aps_intragenic_Testing$PARTIAL_EXON_DUP%in%gene_list,]$APS)
  row_num = row_num+1
  out_table[row_num,1] = 'tss'
  out_table[row_num,2] = mean(pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$TSS_DUP) & pass_aps_intragenic_Testing$TSS_DUP%in%gene_list,]$APS)
  row_num = row_num+1
  out_table[row_num,1] = 'dp'
  out_table[row_num,2] = mean(pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$DUP_PARTIAL) & pass_aps_intragenic_Testing$DUP_PARTIAL%in%gene_list,]$APS)
  row_num = row_num+1
  out_table[row_num,1] = 'inv'
  out_table[row_num,2] = mean(pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$INV_SPAN) & pass_aps_intragenic_Testing$INV_SPAN%in%gene_list,]$APS)
  row_num = row_num+1
  out_table[row_num,1] = 'intronic'
  out_table[row_num,2] = mean(pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$INTRONIC) & pass_aps_intragenic_Testing$INTRONIC%in%gene_list,]$APS)
  row_num = row_num+1
  out_table[row_num,1] = 'nearest_tss'
  out_table[row_num,2] = mean(pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$NEAREST_TSS) & pass_aps_intragenic_Testing$NEAREST_TSS%in%gene_list,]$APS)
  row_num = row_num+1
  out_table[row_num,1] = 'promoter'
  out_table[row_num,2] = mean(pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PROMOTER) & pass_aps_intragenic_Testing$PROMOTER%in%gene_list,]$APS)
  row_num = row_num+1
  out_table[row_num,1] = 'promoter_del'
  out_table[row_num,2] = mean(pass_aps_intragenic_Testing[pass_aps_intragenic_Testing$SVTYPE=="DEL" & !is.na(pass_aps_intragenic_Testing$PROMOTER) & pass_aps_intragenic_Testing$PROMOTER%in%gene_list,]$APS)
  row_num = row_num+1
  out_table[row_num,1] = 'promoter_ins'
  out_table[row_num,2] = mean(pass_aps_intragenic_Testing[pass_aps_intragenic_Testing$SVTYPE=="INS" & !is.na(pass_aps_intragenic_Testing$PROMOTER) & pass_aps_intragenic_Testing$PROMOTER%in%gene_list,]$APS)
  row_num = row_num+1
  out_table[row_num,1] = 'promoter_dup'
  out_table[row_num,2] = mean(pass_aps_intragenic_Testing[pass_aps_intragenic_Testing$SVTYPE=="DUP" & !is.na(pass_aps_intragenic_Testing$PROMOTER) & pass_aps_intragenic_Testing$PROMOTER%in%gene_list,]$APS)
  row_num = row_num+1
  out_table[row_num,1] = 'promoter_oth'
  out_table[row_num,2] = mean(pass_aps_intragenic_Testing[!pass_aps_intragenic_Testing$SVTYPE%in%c("DEL",'DUP','INS') & !is.na(pass_aps_intragenic_Testing$PROMOTER) & pass_aps_intragenic_Testing$PROMOTER%in%gene_list,]$APS)
  row_num = row_num+1
  out_table[row_num,1] = 'utr'
  out_table[row_num,2] = mean(pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$UTR) & pass_aps_intragenic_Testing$UTR%in%gene_list,]$APS)
  row_num = row_num+1
  out_table[row_num,1] = 'utr_del'
  out_table[row_num,2] = mean(pass_aps_intragenic_Testing[pass_aps_intragenic_Testing$SVTYPE=="DEL" & !is.na(pass_aps_intragenic_Testing$UTR) & pass_aps_intragenic_Testing$UTR%in%gene_list,]$APS)
  row_num = row_num+1
  out_table[row_num,1] = 'utr_ins'
  out_table[row_num,2] = mean(pass_aps_intragenic_Testing[pass_aps_intragenic_Testing$SVTYPE=="INS" & !is.na(pass_aps_intragenic_Testing$UTR) & pass_aps_intragenic_Testing$UTR%in%gene_list,]$APS)
  row_num = row_num+1
  out_table[row_num,1] = 'utr_dup'
  out_table[row_num,2] = mean(pass_aps_intragenic_Testing[pass_aps_intragenic_Testing$SVTYPE=="DUP" & !is.na(pass_aps_intragenic_Testing$UTR) & pass_aps_intragenic_Testing$UTR%in%gene_list,]$APS)
  row_num = row_num+1
  out_table[row_num,1] = 'utr_oth'
  out_table[row_num,2] = mean(pass_aps_intragenic_Testing[!pass_aps_intragenic_Testing$SVTYPE%in%c("DEL",'DUP','INS') & !is.na(pass_aps_intragenic_Testing$UTR) & pass_aps_intragenic_Testing$UTR%in%gene_list,]$APS)
  return(out_table)
} 
build_marginal_APS_model.by_all_SVs<-function(dat_2, func_anno, aps_model_prefix, output_file_name){
  #extract SVID of intergenic SVs
  integenic = func_anno[func_anno$PREDICTED_INTERGENIC=="True",]
  intronic = func_anno[!is.na(func_anno$PREDICTED_INTRONIC),]
  #split intergenic SVs as training data
  pass_intergenic_Training = add_size_cate(dat_2[dat_2$name%in%integenic$name,])
  pass_inv_training = dat_2[!is.na(dat_2$PREDICTED_INTRONIC),]
  pass_inv_training = add_size_cate(pass_inv_training[pass_inv_training$SVTYPE=="INV",])
  pass_intergenic_Training = rbind(pass_intergenic_Training, pass_inv_training)
  #split intragenic SVs as testing data
  pass_intragenic_Testing = add_size_cate(dat_2)
  #annotate intergenic training data with the categories and size grids used for marginal APS calculation
  pass_ps_intergenic_Training = reorganize_supp_evi.all_sizes.V2(pass_intergenic_Training)
  #annotate intragenic training data with the categories and size grids used for marginal APS calculation
  pass_ps_intragenic_Testing = reorganize_supp_evi.all_sizes.V2(pass_intragenic_Testing)
  #correct the size grids for testing data according to training data:
  pass_ps_intragenic_Testing = add_ps_annotations(pass_ps_intergenic_Training, pass_ps_intragenic_Testing)
  
  #build the marginal APS model
  aps_models = calculate_marginal_APS_model(pass_ps_intergenic_Training, aps_model_prefix)
  #annotate the testing data with APS
  pass_aps_intragenic_Testing= merge(pass_ps_intragenic_Testing, aps_models, by=c('size_decile','aps_cate'))
  pass_aps_intragenic_Testing[,ncol(pass_aps_intragenic_Testing)+1] = 0-pass_aps_intragenic_Testing$fitted_PS
  colnames(pass_aps_intragenic_Testing)[ncol(pass_aps_intragenic_Testing)] = 'APS'
  pass_aps_intragenic_Testing[pass_aps_intragenic_Testing$AC==1,]$APS = 1-pass_aps_intragenic_Testing[pass_aps_intragenic_Testing$AC==1,]$fitted_PS
  #annotate the testing data with function predictions
  write.table(pass_aps_intragenic_Testing[,c("name","singleton_prop", "PS", "fitted_PS", "residual_PS", "APS" )], output_file_name, quote = F, sep = '\t', col.names = T, row.names = F)
  return(pass_aps_intragenic_Testing[,c("name","singleton_prop", "PS", "fitted_PS", "residual_PS", "APS" )])
}
build_marginal_APS_model.by_uniq_SVs<-function(dat_2, func_anno, aps_model_prefix, output_file_name){
  #extract SVID of intergenic SVs
  integenic = func_anno[func_anno$PREDICTED_INTERGENIC=="True",]
  intronic = func_anno[!is.na(func_anno$PREDICTED_INTRONIC),]
  #split intergenic SVs as training data
  pass_intergenic_Training = add_size_cate(dat_2[dat_2$name%in%integenic$name,])
  pass_inv_training = dat_2[!is.na(dat_2$PREDICTED_INTRONIC),]
  pass_inv_training = add_size_cate(pass_inv_training[pass_inv_training$SVTYPE=="INV",])
  pass_intergenic_Training = rbind(pass_intergenic_Training, pass_inv_training)
  #split intragenic SVs as testing data
  pass_intragenic_Testing = add_size_cate(dat_2)
  #annotate intergenic training data with the categories and size grids used for marginal APS calculation
  pass_ps_intergenic_Training = reorganize_supp_evi.unique_sizes.V2(pass_intergenic_Training)
  #annotate intragenic training data with the categories and size grids used for marginal APS calculation
  pass_ps_intragenic_Testing = reorganize_supp_evi.unique_sizes.V2(pass_intragenic_Testing)
  #correct the size grids for testing data according to training data:
  pass_ps_intragenic_Testing = add_ps_annotations(pass_ps_intergenic_Training, pass_ps_intragenic_Testing)
  
  #build the marginal APS model
  aps_models = calculate_marginal_APS_model(pass_ps_intergenic_Training, aps_model_prefix)
  #annotate the testing data with APS
  pass_aps_intragenic_Testing= merge(pass_ps_intragenic_Testing, aps_models, by=c('size_decile','aps_cate'))
  pass_aps_intragenic_Testing[,ncol(pass_aps_intragenic_Testing)+1] = 0-pass_aps_intragenic_Testing$fitted_PS
  colnames(pass_aps_intragenic_Testing)[ncol(pass_aps_intragenic_Testing)] = 'APS'
  pass_aps_intragenic_Testing[pass_aps_intragenic_Testing$AC==1,]$APS = 1-pass_aps_intragenic_Testing[pass_aps_intragenic_Testing$AC==1,]$fitted_PS
  #annotate the testing data with function predictions
  pass_aps_intragenic_Testing = merge(pass_aps_intragenic_Testing,func_anno[,c("name", "PREDICTED_BREAKEND_EXONIC","PREDICTED_COPY_GAIN","PREDICTED_DUP_PARTIAL","PREDICTED_INTERGENIC","PREDICTED_INTRAGENIC_EXON_DUP","PREDICTED_INTRONIC","PREDICTED_INV_SPAN","PREDICTED_LOF","PREDICTED_MSV_EXON_OVERLAP", "PREDICTED_NEAREST_TSS","PREDICTED_PARTIAL_EXON_DUP","PREDICTED_PROMOTER", "PREDICTED_TSS_DUP","PREDICTED_UTR"   )], by='name')
  write.table(pass_aps_intragenic_Testing[,c("name","singleton_prop", "PS", "fitted_PS", "residual_PS", "APS" )], output_file_name, quote = F, sep = '\t', col.names = T, row.names = F)
  return(pass_aps_intragenic_Testing[,c("name","singleton_prop", "PS", "fitted_PS", "residual_PS", "APS" )])
}
#calculate APS of SVs of a specific group
calculate_APS<-function(pass_aps_intragenic_Testing, target_SVID, permutations = 100, sample_prop = .9){
  target_SVs = pass_aps_intragenic_Testing[pass_aps_intragenic_Testing$name%in%target_SVID,]
  APS = mean(target_SVs$APS)
  out_stat=data.frame('permutate'=0, 'APS' = APS)
  for(i in c(1:permutations)){
    print(paste('permutation ', i, ' ...', sep=''))
    sample_SVs = target_SVs[sample(c(1:nrow(target_SVs)), as.integer(nrow(target_SVs)*sample_prop)),]
    sample_APS = mean(sample_SVs$APS)
    out_stat[nrow(out_stat)+1,1] = i
    out_stat[nrow(out_stat),] = sample_APS
  }
  return(out_stat)
}
#calculate APS for SVs annotated for certain functional interuptions
calculate_APS_across_functional_annotations<-function(pass_aps_intragenic_Testing){
  
  plof_all = calculate_APS(pass_aps_intragenic_Testing, pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_LOF),]$name)
  plof_DEL = calculate_APS(pass_aps_intragenic_Testing, pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_LOF) & pass_aps_intragenic_Testing$SVTYPE=="DEL",]$name)
  plof_INS = calculate_APS(pass_aps_intragenic_Testing, pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_LOF) & pass_aps_intragenic_Testing$SVTYPE=="INS",]$name)
  plof_DUP = calculate_APS(pass_aps_intragenic_Testing, pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_LOF) & pass_aps_intragenic_Testing$SVTYPE=="DUP",]$name)
  plof_oth = calculate_APS(pass_aps_intragenic_Testing, pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_LOF) & pass_aps_intragenic_Testing$SVTYPE%in%c("INV","CPX"),]$name)
  cg = calculate_APS(pass_aps_intragenic_Testing, pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_COPY_GAIN),]$name)
  ied = calculate_APS(pass_aps_intragenic_Testing, pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_INTRAGENIC_EXON_DUP),]$name)
  ped = calculate_APS(pass_aps_intragenic_Testing, pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_PARTIAL_EXON_DUP),]$name)
  tss = calculate_APS(pass_aps_intragenic_Testing, pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_TSS_DUP),]$name)
  dp = calculate_APS(pass_aps_intragenic_Testing, pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_DUP_PARTIAL),]$name)
  dp = calculate_APS(pass_aps_intragenic_Testing, pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_DUP_PARTIAL),]$name)
  inv = calculate_APS(pass_aps_intragenic_Testing, pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_INV_SPAN),]$name)
  promoter = calculate_APS(pass_aps_intragenic_Testing, pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_PROMOTER),]$name)
  intron = calculate_APS(pass_aps_intragenic_Testing, pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_INTRONIC),]$name)
  utr = calculate_APS(pass_aps_intragenic_Testing, pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_UTR),]$name)
  intergenic = calculate_APS(pass_aps_intragenic_Testing, pass_aps_intragenic_Testing[pass_aps_intragenic_Testing$PREDICTED_INTERGENIC=="True",]$name)

  all = calculate_APS(pass_aps_intragenic_Testing, pass_aps_intragenic_Testing$name)
  
  
  aps_table = data.frame('func'=0,'aps'=0, 'lower_CI'=0,'higher_CI'=0)
  aps_table[1,] = c('ALL',all[1,2], quantile(all[2:101,2], .05), quantile(all[2:101,2], .95))
  aps_table[nrow(aps_table)+1,] = c('plof_all',plof_all[1,2], quantile(plof_all[2:101,2], .05), quantile(plof_all[2:101,2], .95))
  aps_table[nrow(aps_table)+1,] = c('plof_ins',plof_INS[1,2], quantile(plof_INS[2:101,2], .05), quantile(plof_INS[2:101,2], .95))
  aps_table[nrow(aps_table)+1,] = c('ied',ied[1,2], quantile(ied[2:101,2], .05), quantile(ied[2:101,2], .95))
  aps_table[nrow(aps_table)+1,] = c('plof_del',plof_DEL[1,2], quantile(plof_DEL[2:101,2], .05), quantile(plof_DEL[2:101,2], .95))
  aps_table[nrow(aps_table)+1,] = c('ped',ped[1,2], quantile(ped[2:101,2], .05), quantile(ped[2:101,2], .95))
  aps_table[nrow(aps_table)+1,] = c('plof_oth',plof_oth[1,2], quantile(plof_oth[2:101,2], .05), quantile(plof_oth[2:101,2], .95))
  aps_table[nrow(aps_table)+1,] = c('plof_dup',plof_DUP[1,2], quantile(plof_DUP[2:101,2], .05), quantile(plof_DUP[2:101,2], .95))
  aps_table[nrow(aps_table)+1,] = c('inv' ,inv[1,2], quantile(inv[2:101,2], .05), quantile(inv[2:101,2], .95))
  aps_table[nrow(aps_table)+1,] = c('tss',tss[1,2], quantile(tss[2:101,2], .05), quantile(tss[2:101,2], .95))
  aps_table[nrow(aps_table)+1,] = c('dp' ,dp[1,2], quantile(dp[2:101,2], .05), quantile(dp[2:101,2], .95))
  aps_table[nrow(aps_table)+1,] = c('promotor' ,promoter[1,2], quantile(promoter[2:101,2], .05), quantile(promoter[2:101,2], .95))
  aps_table[nrow(aps_table)+1,] = c('intron' ,intron[1,2], quantile(intron[2:101,2], .05), quantile(intron[2:101,2], .95))
  aps_table[nrow(aps_table)+1,] = c('cg' ,cg[1,2],  quantile(cg[2:101,2], .05), quantile(cg[2:101,2], .95))
  aps_table[nrow(aps_table)+1,] = c('intergenic' ,intergenic[1,2],  quantile(intergenic[2:101,2], .05), quantile(intergenic[2:101,2], .95))
  aps_table[nrow(aps_table)+1,] = c('utr' ,utr[1,2],  quantile(utr[2:101,2], .05), quantile(utr[2:101,2], .95))
  
  return(aps_table)
}
#calculate marginal APS for each SV
calculate_marginal_APS_model<-function(pass_ps_intergenic_Training, pdf_prefix = 'APS_models.'){
  aps_models = data.frame('xvalues'=0,'yvalues'=0,'fitted_y'=0,'V4'=0)
  for(i in sort(unique(pass_ps_intergenic_Training$aps_cate))){
    print(paste(pdf_prefix,i, '.pdf', sep = ''))
    pdf(paste(pdf_prefix,i, '.pdf', sep = ''), width = 6, height = 3)
    par(mfrow=c(1,2))
    par(mar=c(3,3,2,1))
    
    band_width = 6
    if(i=='DEL,US_RM,others,400bpto5Kb' | i=='DEL,SD_SR,others,400bpto5Kb' | i=='INS,US_RM,others,over250bp'){band_width = 6}
    #if(i=="DEL,SD_SR,others,300bpto5Kb"  ){band_width=20}
    #if(i=="DEL,SD_SR,others,under300bp"| "DEL,SD_SR,sr_only,under300bp"){band_width=40}
    if(!is.na(i)){
      tmp = pass_ps_intergenic_Training[!is.na(pass_ps_intergenic_Training$aps_cate) & pass_ps_intergenic_Training$aps_cate==i,]

      tmp_stat = unique(tmp[,c("size_decile", "singleton_prop")])
      tmp_stat[,3] = 0
      colnames(tmp_stat)[3] = 'median_size'
      tmp_stat = tmp_stat[order(tmp_stat$size_decile),]
      for(j in unique(tmp_stat[,1])){
        size_list = tmp[tmp$size_decile==j,]$SVLEN
        tmp_stat[tmp_stat[,1]==j,3] = median(size_list) 
      }
      
      xvalues = tmp_stat$size_decile

      #xvalues = log10(tmp_stat$median_size)
      yvalues = tmp_stat$singleton_prop
      y_fitted_model = roll_apply_nlm(xvalues, yvalues, band_width,i)
      
      svtype = strsplit(as.character(i),',')[[1]][1]
      # Plot the chart with new data by fitting it to a prediction from 100 data points.
      plot(xvalues, yvalues,  col = colors[colors[1]==svtype,2], pch=19, xaxt='n', main=i, cex.main = 1, cex = 1, las=2)
      axis(1, xvalues[c(1,round(length(xvalues)/2),length(xvalues))], labels = tmp_stat$median_size[c(1,round(length(xvalues)/2),length(xvalues))],mgp=c(.5,.3,0))
      new.data <- data.frame(xvalues = seq(min(xvalues),max(xvalues),len = 100))
      lines(y_fitted_model$xvalues, y_fitted_model$fitted_y, lwd=2, col='black')
      axis(1, log10(c(1:5)*100), labels = paste(c(1:5)*100,'bp', sep=''), mgp = c(.5,.5,0))
      axis(1, log10(c(500,c(1:5)*1000)), labels = c('500bp',paste(c(1:5),'Kb', sep='')),  mgp = c(.5,.5,0))
      axis(1, log10(c(0:10)*5000), labels = paste(c(0:10)*5,'Kb', sep=''), mgp = c(.5,.5,0))
      
      residual = y_fitted_model$yvalues-y_fitted_model$fitted_y
      plot(xvalues, residual,  col = colors[colors[1]==svtype,2], pch=17, xaxt='n', main=i, ylim=c(-.3,.3), cex.main = 1, cex = 1, las=2)
      axis(1, xvalues[c(1,round(length(xvalues)/2),length(xvalues))], labels = tmp_stat$median_size[c(1,round(length(xvalues)/2),length(xvalues))],mgp=c(.5,.3,0))
      abline(h=0, col='black', lwd=2)
      axis(1, log10(c(1:5)*100), labels = paste(c(1:5)*100,'bp', sep=''),  mgp = c(.5,.5,0))
      axis(1, log10(c(500,c(1:5)*1000)), labels = c('500bp',paste(c(1:5),'Kb', sep='')),  mgp = c(.5,.5,0))
      axis(1, log10(c(0:10)*5000), labels = paste(c(0:10)*5,'Kb', sep=''),  mgp = c(.5,.5,0))
    }
    aps_models = rbind(aps_models, y_fitted_model)
    dev.off()
  }
  aps_models[,5] = aps_models[,2]-aps_models[,3]
  colnames(aps_models) = c('size_decile','PS','fitted_PS','aps_cate','residual_PS')
  return(aps_models)
}
#calculate the size ranges of each SV category used in APS calculation
calculate_size_ranges<-function(data){
  out = data.frame('aps_cate' = 0, 'size_decile' = 0, 'min_size' = 0, 'max_size' = 0)
  row_out = 0
  for(i in unique(data$aps_cate)){
    tmp1 = data[!is.na(data$aps_cate) & data$aps_cate==i,]
    if(!is.na(i)){
      print(c(i))
      for(j in unique(tmp1$size_decile)){
        tmp2 = tmp1[tmp1$size_decile==j, ]
        size_range = range(tmp2$SVLEN)
        row_out = row_out+1
        out[row_out,1] = i
        out[row_out,2] = j
        out[row_out,3] = size_range[1]
        out[row_out,4] = size_range[2]
      }
    }
  }
  out = out[-1,]
  return(out)
}
#extract the actual singleton proportion for each type of SVs
calculate_overall_singleton_proportion<-function(data){
  
  singleton_prop=data.frame('function'=0,'singleton_prop'=0,'mean_PS'=0,'mean_fitted_PS'=0,'mean_APS'=0)
  row_num = 1
  singleton_prop[row_num,1] = 'ALL'
  singleton_prop[row_num,2] = nrow(data[data$AC==1,])/nrow(data)
  singleton_prop[row_num,3] = mean(data$PS)
  singleton_prop[row_num,4] = mean(data$fitted_PS)  
  singleton_prop[row_num,5] = mean(data$APS)
  row_num = row_num+1
  singleton_prop[row_num,1] = 'pLoF_del'
  singleton_prop[row_num,2] = nrow(data[data$SVTYPE=="DEL" & !is.na(data$PREDICTED_LOF) & data$AC==1,])/nrow(data[data$SVTYPE=="DEL" & !is.na(data$PREDICTED_LOF),])
  singleton_prop[row_num,3] = mean(data[data$SVTYPE=="DEL" & !is.na(data$PREDICTED_LOF),]$PS)
  singleton_prop[row_num,4] = mean(data[data$SVTYPE=="DEL" & !is.na(data$PREDICTED_LOF),]$fitted_PS)  
  singleton_prop[row_num,5] = mean(data[data$SVTYPE=="DEL" & !is.na(data$PREDICTED_LOF),]$APS)
  row_num = row_num+1
  singleton_prop[row_num,1] = 'pLoF_ins'
  singleton_prop[row_num,2] = nrow(data[data$SVTYPE=="INS" & !is.na(data$PREDICTED_LOF) & data$AC==1,])/nrow(data[data$SVTYPE=="INS" & !is.na(data$PREDICTED_LOF),])
  singleton_prop[row_num,3] = mean(data[data$SVTYPE=="INS" & !is.na(data$PREDICTED_LOF),]$PS)
  singleton_prop[row_num,4] = mean(data[data$SVTYPE=="INS" & !is.na(data$PREDICTED_LOF),]$fitted_PS)  
  singleton_prop[row_num,5] = mean(data[data$SVTYPE=="INS" & !is.na(data$PREDICTED_LOF),]$APS)
  row_num = row_num+1
  singleton_prop[row_num,1] = 'pLoF_oth'
  singleton_prop[row_num,2] = nrow(data[!data$SVTYPE%in%c('DEL',"INS") & !is.na(data$PREDICTED_LOF) & data$AC==1,])/nrow(data[!data$SVTYPE%in%c('DEL',"INS") & !is.na(data$PREDICTED_LOF),])
  singleton_prop[row_num,3] = mean(data[!data$SVTYPE%in%c('DEL',"INS") & !is.na(data$PREDICTED_LOF),]$PS)
  singleton_prop[row_num,4] = mean(data[!data$SVTYPE%in%c('DEL',"INS") & !is.na(data$PREDICTED_LOF),]$fitted_PS)  
  singleton_prop[row_num,5] = mean(data[!data$SVTYPE%in%c('DEL',"INS") & !is.na(data$PREDICTED_LOF),]$APS)
  row_num = row_num+1
  singleton_prop[row_num,1] = 'CG'
  singleton_prop[row_num,2] = nrow(data[!is.na(data$PREDICTED_COPY_GAIN) & data$AC==1,])/nrow(data[!is.na(data$PREDICTED_COPY_GAIN),])
  singleton_prop[row_num,3] = mean(data[!is.na(data$PREDICTED_COPY_GAIN),]$PS)
  singleton_prop[row_num,4] = mean(data[!is.na(data$PREDICTED_COPY_GAIN),]$fitted_PS)  
  singleton_prop[row_num,5] = mean(data[!is.na(data$PREDICTED_COPY_GAIN),]$APS)
  row_num = row_num+1
  singleton_prop[row_num,1] = 'IED'
  singleton_prop[row_num,2] = nrow(data[!is.na(data$PREDICTED_INTRAGENIC_EXON_DUP) & data$AC==1,])/nrow(data[!is.na(data$PREDICTED_INTRAGENIC_EXON_DUP),])
  singleton_prop[row_num,3] = mean(data[!is.na(data$PREDICTED_INTRAGENIC_EXON_DUP),]$PS)
  singleton_prop[row_num,4] = mean(data[!is.na(data$PREDICTED_INTRAGENIC_EXON_DUP),]$fitted_PS)  
  singleton_prop[row_num,5] = mean(data[!is.na(data$PREDICTED_INTRAGENIC_EXON_DUP),]$APS)
  row_num = row_num+1
  singleton_prop[row_num,1] = 'PED'
  singleton_prop[row_num,2] = nrow(data[!is.na(data$PREDICTED_PARTIAL_EXON_DUP) & data$AC==1,])/nrow(data[!is.na(data$PREDICTED_PARTIAL_EXON_DUP),])
  singleton_prop[row_num,3] = mean(data[!is.na(data$PREDICTED_PARTIAL_EXON_DUP),]$PS)
  singleton_prop[row_num,4] = mean(data[!is.na(data$PREDICTED_PARTIAL_EXON_DUP),]$fitted_PS)  
  singleton_prop[row_num,5] = mean(data[!is.na(data$PREDICTED_PARTIAL_EXON_DUP),]$APS)
  row_num = row_num+1
  singleton_prop[row_num,1] = 'TSS'
  singleton_prop[row_num,2] = nrow(data[!is.na(data$PREDICTED_TSS_DUP) & data$AC==1,])/nrow(data[!is.na(data$PREDICTED_TSS_DUP),])
  singleton_prop[row_num,3] = mean(data[!is.na(data$PREDICTED_TSS_DUP),]$PS)
  singleton_prop[row_num,4] = mean(data[!is.na(data$PREDICTED_TSS_DUP),]$fitted_PS)  
  singleton_prop[row_num,5] = mean(data[!is.na(data$PREDICTED_TSS_DUP),]$APS)
  row_num = row_num+1
  singleton_prop[row_num,1] = 'DP'
  singleton_prop[row_num,2] = nrow(data[!is.na(data$PREDICTED_DUP_PARTIAL) & data$AC==1,])/nrow(data[!is.na(data$PREDICTED_DUP_PARTIAL),])
  singleton_prop[row_num,3] = mean(data[!is.na(data$PREDICTED_DUP_PARTIAL),]$PS)
  singleton_prop[row_num,4] = mean(data[!is.na(data$PREDICTED_DUP_PARTIAL),]$fitted_PS)  
  singleton_prop[row_num,5] = mean(data[!is.na(data$PREDICTED_DUP_PARTIAL),]$APS)
  #row_num = row_num+1
  #singleton_prop[row_num,1] = 'MSV'
  #singleton_prop[row_num,2] = nrow(data[!is.na(data$PREDICTED_MSV_EXON_OVERLAP) & data$AC==1,])/nrow(data[!is.na(data$PREDICTED_MSV_EXON_OVERLAP),])
  #singleton_prop[row_num,3] = mean(data[!is.na(data$PREDICTED_MSV_EXON_OVERLAP),]$PS)
  #singleton_prop[row_num,4] = mean(data[!is.na(data$PREDICTED_MSV_EXON_OVERLAP),]$fitted_PS)  
  #singleton_prop[row_num,5] = mean(data[!is.na(data$PREDICTED_MSV_EXON_OVERLAP),]$APS)
  row_num = row_num+1
  singleton_prop[row_num,1] = 'INV'
  singleton_prop[row_num,2] = nrow(data[!is.na(data$PREDICTED_INV_SPAN) & data$AC==1,])/nrow(data[!is.na(data$PREDICTED_INV_SPAN),])
  singleton_prop[row_num,3] = mean(data[!is.na(data$PREDICTED_INV_SPAN),]$PS)
  singleton_prop[row_num,4] = mean(data[!is.na(data$PREDICTED_INV_SPAN),]$fitted_PS)  
  singleton_prop[row_num,5] = mean(data[!is.na(data$PREDICTED_INV_SPAN),]$APS)
  row_num = row_num+1
  singleton_prop[row_num,1] = 'utr_del'
  singleton_prop[row_num,2] = nrow(data[data$SVTYPE=="DEL" & !is.na(data$PREDICTED_UTR) & data$AC==1,])/nrow(data[data$SVTYPE=="DEL" & !is.na(data$PREDICTED_UTR),])
  singleton_prop[row_num,3] = mean(data[data$SVTYPE=="DEL" & !is.na(data$PREDICTED_UTR),]$PS)
  singleton_prop[row_num,4] = mean(data[data$SVTYPE=="DEL" & !is.na(data$PREDICTED_UTR),]$fitted_PS)  
  singleton_prop[row_num,5] = mean(data[data$SVTYPE=="DEL" & !is.na(data$PREDICTED_UTR),]$APS)
  row_num = row_num+1
  singleton_prop[row_num,1] = 'utr_ins'
  singleton_prop[row_num,2] = nrow(data[data$SVTYPE=="INS" & !is.na(data$PREDICTED_UTR) & data$AC==1,])/nrow(data[data$SVTYPE=="INS" & !is.na(data$PREDICTED_UTR),])
  singleton_prop[row_num,3] = mean(data[data$SVTYPE=="INS" & !is.na(data$PREDICTED_UTR),]$PS)
  singleton_prop[row_num,4] = mean(data[data$SVTYPE=="INS" & !is.na(data$PREDICTED_UTR),]$fitted_PS)  
  singleton_prop[row_num,5] = mean(data[data$SVTYPE=="INS" & !is.na(data$PREDICTED_UTR),]$APS)
  row_num = row_num+1
  singleton_prop[row_num,1] = 'utr_oth'
  singleton_prop[row_num,2] = nrow(data[!data$SVTYPE%in%c('DEL',"INS") & !is.na(data$PREDICTED_UTR) & data$AC==1,])/nrow(data[!data$SVTYPE%in%c('DEL',"INS") & !is.na(data$PREDICTED_UTR),])
  singleton_prop[row_num,3] = mean(data[!data$SVTYPE%in%c('DEL',"INS") & !is.na(data$PREDICTED_UTR),]$PS)
  singleton_prop[row_num,4] = mean(data[!data$SVTYPE%in%c('DEL',"INS") & !is.na(data$PREDICTED_UTR),]$fitted_PS)  
  singleton_prop[row_num,5] = mean(data[!data$SVTYPE%in%c('DEL',"INS") & !is.na(data$PREDICTED_UTR),]$APS)
  row_num = row_num+1
  singleton_prop[row_num,1] = 'promoter_del'
  singleton_prop[row_num,2] = nrow(data[data$SVTYPE=="DEL" & !is.na(data$PREDICTED_PROMOTER) & data$AC==1,])/nrow(data[data$SVTYPE=="DEL" & !is.na(data$PREDICTED_PROMOTER),])
  singleton_prop[row_num,3] = mean(data[data$SVTYPE=="DEL" & !is.na(data$PREDICTED_PROMOTER),]$PS)
  singleton_prop[row_num,4] = mean(data[data$SVTYPE=="DEL" & !is.na(data$PREDICTED_PROMOTER),]$fitted_PS)  
  singleton_prop[row_num,5] = mean(data[data$SVTYPE=="DEL" & !is.na(data$PREDICTED_PROMOTER),]$APS)
  row_num = row_num+1
  singleton_prop[row_num,1] = 'promoter_ins'
  singleton_prop[row_num,2] = nrow(data[data$SVTYPE=="INS" & !is.na(data$PREDICTED_PROMOTER) & data$AC==1,])/nrow(data[data$SVTYPE=="INS" & !is.na(data$PREDICTED_PROMOTER),])
  singleton_prop[row_num,3] = mean(data[data$SVTYPE=="INS" & !is.na(data$PREDICTED_PROMOTER),]$PS)
  singleton_prop[row_num,4] = mean(data[data$SVTYPE=="INS" & !is.na(data$PREDICTED_PROMOTER),]$fitted_PS)  
  singleton_prop[row_num,5] = mean(data[data$SVTYPE=="INS" & !is.na(data$PREDICTED_PROMOTER),]$APS)
  row_num = row_num+1
  singleton_prop[row_num,1] = 'promoter_oth'
  singleton_prop[row_num,2] = nrow(data[!data$SVTYPE%in%c('DEL',"INS") & !is.na(data$PREDICTED_PROMOTER) & data$AC==1,])/nrow(data[!data$SVTYPE%in%c('DEL',"INS") & !is.na(data$PREDICTED_PROMOTER),])
  singleton_prop[row_num,3] = mean(data[!data$SVTYPE%in%c('DEL',"INS") & !is.na(data$PREDICTED_PROMOTER),]$PS)
  singleton_prop[row_num,4] = mean(data[!data$SVTYPE%in%c('DEL',"INS") & !is.na(data$PREDICTED_PROMOTER),]$fitted_PS)  
  singleton_prop[row_num,5] = mean(data[!data$SVTYPE%in%c('DEL',"INS") & !is.na(data$PREDICTED_PROMOTER),]$APS)
  row_num = row_num+1
  singleton_prop[row_num,1] = 'intronic_del'
  singleton_prop[row_num,2] = nrow(data[data$SVTYPE=="DEL" & !is.na(data$PREDICTED_INTRONIC) & data$AC==1,])/nrow(data[data$SVTYPE=="DEL" & !is.na(data$PREDICTED_INTRONIC),])
  singleton_prop[row_num,3] = mean(data[data$SVTYPE=="DEL" & !is.na(data$PREDICTED_INTRONIC),]$PS)
  singleton_prop[row_num,4] = mean(data[data$SVTYPE=="DEL" & !is.na(data$PREDICTED_INTRONIC),]$fitted_PS)  
  singleton_prop[row_num,5] = mean(data[data$SVTYPE=="DEL" & !is.na(data$PREDICTED_INTRONIC),]$APS)
  row_num = row_num+1
  singleton_prop[row_num,1] = 'intronic_ins'
  singleton_prop[row_num,2] = nrow(data[data$SVTYPE=="INS" & !is.na(data$PREDICTED_INTRONIC) & data$AC==1,])/nrow(data[data$SVTYPE=="INS" & !is.na(data$PREDICTED_INTRONIC),])
  singleton_prop[row_num,3] = mean(data[data$SVTYPE=="INS" & !is.na(data$PREDICTED_INTRONIC),]$PS)
  singleton_prop[row_num,4] = mean(data[data$SVTYPE=="INS" & !is.na(data$PREDICTED_INTRONIC),]$fitted_PS)  
  singleton_prop[row_num,5] = mean(data[data$SVTYPE=="INS" & !is.na(data$PREDICTED_INTRONIC),]$APS)
  row_num = row_num+1
  singleton_prop[row_num,1] = 'intronic_oth'
  singleton_prop[row_num,2] = nrow(data[!data$SVTYPE%in%c('DEL',"INS") & !is.na(data$PREDICTED_INTRONIC) & data$AC==1,])/nrow(data[!data$SVTYPE%in%c('DEL',"INS") & !is.na(data$PREDICTED_INTRONIC),])
  singleton_prop[row_num,3] = mean(data[!data$SVTYPE%in%c('DEL',"INS") & !is.na(data$PREDICTED_INTRONIC),]$PS)
  singleton_prop[row_num,4] = mean(data[!data$SVTYPE%in%c('DEL',"INS") & !is.na(data$PREDICTED_INTRONIC),]$fitted_PS)  
  singleton_prop[row_num,5] = mean(data[!data$SVTYPE%in%c('DEL',"INS") & !is.na(data$PREDICTED_INTRONIC),]$APS)
  return(singleton_prop)
}
calculate_decile_table <- function(data){
  decile_table = data.frame('decile'=0,'aps_loeuf'=0, 'aps_pHaplo'=0, 'aps_pTriplo'=0)
  for(i in c(1:10)){
    decile_table[i,1] = i
    decile_table[i,2] = mean(data[data$ptv_oe_dec==i,]$APS)
    decile_table[i,3] = mean(data[!data$phaplo_oe_cent > i*10  & data$phaplo_oe_cent > (i-1)*10,]$APS)
    decile_table[i,4] = mean(data[!data$ptriplo_oe_cent > i*10 & data$ptriplo_oe_cent > (i-1)*10,]$APS)
  }
  return(decile_table)
}
calcu_sd_from_SVcounts<-function(sv_count){        
  slope = -0.5087
  intersect = -0.3093
  y = slope * log10(sv_count) + intersect
  return(10^y)
}
merge_SVs_with_loeuf_phaplo_ptriplo <- function(pass_aps_intragenic_Testing, snv_data, function_name){
  data = pass_aps_intragenic_Testing[,c("name","size_decile","aps_cate", "X.chrom","start","end","svtype","samples","AC", "ALGORITHMS","AN","CHR2", "CPX_INTERVALS","CPX_TYPE","END", "END2","EVIDENCE","SVLEN", "SVTYPE","AF","N_BI_GENOS","N_HOMREF","N_HET","N_HOMALT", "genomic_context","FILTER","size_cate", "EVIDENCE_reorganized","size_rank","singleton_prop","PS","fitted_PS","residual_PS","APS", 
                                        function_name)]
  data = data[!is.na(data[,ncol(data)]),]
  data[,ncol(data)+1]=data[,ncol(data)]
  colnames(data)[ncol(data)] = 'gene'
  data[,ncol(data)+1]=apply(data, 1, function(x){length(strsplit(as.character(x[match('gene', colnames(data))]),',')[[1]])})
  colnames(data)[ncol(data)] = 'count_genes'
  for(i in 1:nrow(data)){
    if(data[i,]$count_genes>1){
      print(i)
      gene_list = strsplit(as.character(data[i,]$gene), ',')[[1]]
      gene_snv_info = snv_data[snv_data$gene%in%gene_list, ]
      if(nrow(gene_snv_info)>0){
        picked_gene = gene_snv_info[gene_snv_info$LOEUF==min(gene_snv_info$LOEUF), ]
        data[i,]$gene = picked_gene$gene[1]
      }
      if(nrow(gene_snv_info)==0){
        data[i,]$gene = gene_list[1]
      }
    }
  }
  
  data=merge(data, snv_data, by='gene')
  return(data)
}
plot_aps_across_func_anno<-function(aps_table, pdf_file){
  pdf(pdf_file)
  par(mar=c(6,4,4,4))
  plot(c(0,nrow(aps_table)), c(-.07,.14), frame.plot = F,type = 'n', xlab = '', ylab = 'APS', xaxt='n', las=2, yaxt='n')
  axis(1,c(1:nrow(aps_table)), aps_table$func, las=2, cex.axis=1.5)
  axis(2,c(-3:6)/50, labels = c(-3:6)/50, las=2, cex.axis=1.5)
  abline(h=c(-3:6)/50, col='grey', lty=2)
  abline(h=0, col='grey')
  points(c(1:nrow(aps_table)), aps_table$aps, pch=20, cex=2.5)
  for(i in c(1:nrow(aps_table))){
    lines(c(i,i), aps_table[i,c(3,4)])
  }
  dev.off()
}
plot_aps_vs_constraint<-function(data, title){
  plot(c(1,10), range(data[,c(2:4)]), frame.plot = F, type = 'n', xlab = 'decile', ylab = 'APS', las=2, main = title)
  abline(h= c(-4:4)/20, lty=2, col='grey')
  lines(data[,1], data$aps_loeuf, col = 'red', pch=20, type = 'b', lwd=2)
  lines(data[,1], data$aps_pHaplo, col = 'darkgreen', pch=17, type = 'b', lwd=2)
  lines(data[,1], data$aps_pTriplo, col = 'blue', pch=15,  type = 'b', lwd=2)
  legend('bottomleft', c('LOEUF', 'pHaplo', 'Triplo'), col = c('red','darkgreen','blue'), pch=c(20,17,15), lwd=2, bty = 'n')
}
plot_aps_acros_gene_constraint<-function(pass_aps_intragenic_Testing, snv_data, prefix = 'aps_across_functions.by_all_SVLEN.'){
  
  pass_aps_intragenic_Testing = select_most_constraint_gene_all_functions(pass_aps_intragenic_Testing, snv_data)
  
  aps_table.constraint_0 = aps_table.by_gene_list(pass_aps_intragenic_Testing, snv_data$gene)
  aps_table.constraint_1 = aps_table.by_gene_list(pass_aps_intragenic_Testing, snv_data[snv_data$ptv_oe_dec%in%c(1,2),]$gene)
  aps_table.constraint_2 = aps_table.by_gene_list(pass_aps_intragenic_Testing, snv_data[snv_data$ptv_oe_dec%in%c(3,4),]$gene)
  aps_table.constraint_3 = aps_table.by_gene_list(pass_aps_intragenic_Testing, snv_data[snv_data$ptv_oe_dec%in%c(5,6),]$gene)
  aps_table.constraint_4 = aps_table.by_gene_list(pass_aps_intragenic_Testing, snv_data[snv_data$ptv_oe_dec%in%c(7,8),]$gene)
  aps_table.constraint_5 = aps_table.by_gene_list(pass_aps_intragenic_Testing, snv_data[snv_data$ptv_oe_dec%in%c(9,10),]$gene)
  
  aps_table.constraint = merge(aps_table.constraint_0, aps_table.constraint_1, by='func')
  aps_table.constraint = merge(aps_table.constraint,   aps_table.constraint_2, by='func')
  aps_table.constraint = merge(aps_table.constraint,   aps_table.constraint_3, by='func')
  aps_table.constraint = merge(aps_table.constraint,   aps_table.constraint_4, by='func')
  aps_table.constraint = merge(aps_table.constraint,   aps_table.constraint_5, by='func')
  colnames(aps_table.constraint) = c('func','all_genes','constraint_1_2','constraint_3_4','constraint_5_6','constraint_7_8','constraint_9_10')
  write.table(aps_table.constraint, paste(prefix, 'by_constraint_level.tsv', sep='.'), quote = F, sep = '\t', col.names = T, row.names = F)
  
  pdf(paste(prefix, 'by_constraint_level.V1.pdf', sep = '.'))
  par(fig=c(0,1,.3,1))
  par(mar=c(6,4,2,2))
  plot(c(1,23),c(-.2,.2), frame.plot = F, type = 'n', xlab = '', ylab = 'APS', xaxt='n', las=2)
  abline(h=0)
  axis(1,c(1:nrow(aps_table.constraint_0)), aps_table.constraint_0[,1], las=2)
  points(c(1:nrow(aps_table.constraint_1)), aps_table.constraint_1[,2], col='red', pch=20, type='b')
  points(c(1:nrow(aps_table.constraint_2)), aps_table.constraint_2[,2], col='orange', pch=20, type='b')
  points(c(1:nrow(aps_table.constraint_3)), aps_table.constraint_3[,2], col='brown', pch=20, type='b')
  points(c(1:nrow(aps_table.constraint_4)), aps_table.constraint_4[,2], col='darkgreen', pch=20, type='b')
  points(c(1:nrow(aps_table.constraint_5)), aps_table.constraint_5[,2], col='blue', pch=20, type='b')
  legend("bottomright",c('1-2', '5-6','9-10'), col=c('red','brown','blue'), bty = 'n', pch=20, title = 'LOEUF decile')
  dev.off()
  
  pdf(paste(prefix, 'by_constraint_level.V2.pdf', sep = '.'))
  par(mfrow=c(3,2))
  par(mar=rep(2,4))
  
  plot(c(1,23),c(-.2,.2), frame.plot = F, type = 'n', xlab = '', ylab = 'APS', xaxt='n', las=2)
  abline(h=0)
  axis(1,c(1:nrow(aps_table.constraint_0)), aps_table.constraint_0[,1], las=2)
  points(c(1:nrow(aps_table.constraint_0)), aps_table.constraint_0[,2], col='black', pch=20, type='b')
  
  plot(c(1,23),c(-.2,.2), frame.plot = F, type = 'n', xlab = '', ylab = 'APS', xaxt='n', las=2)
  abline(h=0)
  axis(1,c(1:nrow(aps_table.constraint_1)), aps_table.constraint_1[,1], las=2)
  points(c(1:nrow(aps_table.constraint_1)), aps_table.constraint_1[,2], col='red', pch=20, type='b')
  
  plot(c(1,23),c(-.2,.2), frame.plot = F, type = 'n', xlab = '', ylab = 'APS', xaxt='n', las=2)
  abline(h=0)
  axis(1,c(1:nrow(aps_table.constraint_1)), aps_table.constraint_1[,1], las=2)
  points(c(1:nrow(aps_table.constraint_2)), aps_table.constraint_2[,2], col='orange', pch=20, type='b')
  
  plot(c(1,23),c(-.2,.2), frame.plot = F, type = 'n', xlab = '', ylab = 'APS', xaxt='n', las=2)
  abline(h=0)
  axis(1,c(1:nrow(aps_table.constraint_1)), aps_table.constraint_1[,1], las=2)
  points(c(1:nrow(aps_table.constraint_3)), aps_table.constraint_3[,2], col='brown', pch=20, type='b')
  
  plot(c(1,23),c(-.2,.2), frame.plot = F, type = 'n', xlab = '', ylab = 'APS', xaxt='n', las=2)
  abline(h=0)
  axis(1,c(1:nrow(aps_table.constraint_1)), aps_table.constraint_1[,1], las=2)
  points(c(1:nrow(aps_table.constraint_4)), aps_table.constraint_4[,2], col='darkgreen', pch=20, type='b')
  
  plot(c(1,23),c(-.2,.2), frame.plot = F, type = 'n', xlab = '', ylab = 'APS', xaxt='n', las=2)
  abline(h=0)
  axis(1,c(1:nrow(aps_table.constraint_1)), aps_table.constraint_1[,1], las=2)
  points(c(1:nrow(aps_table.constraint_5)), aps_table.constraint_5[,2], col='blue', pch=20, type='b')
  legend("bottomright",c('1-2', '5-6','9-10'), col=c('red','brown','blue'), bty = 'n', pch=20, title = 'LOEUF decile')
  
  dev.off()
}
plot_SV_length<-function(pass_aps_intragenic_Testing, pdf_name = 'function_anno.distribution_of_SV_length.pdf'){
  ied = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_INTRAGENIC_EXON_DUP),]
  ped = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_PARTIAL_EXON_DUP),]
  tss = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_TSS_DUP),]
  dp = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_DUP_PARTIAL),]
  cg = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_COPY_GAIN),]
  lof = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_LOF),]
  promoter = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_PROMOTER),]
  utr = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_UTR),]
  inv = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_INV_SPAN),]
  intronic = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_INTRONIC),]
  intergenic = pass_aps_intragenic_Testing[pass_aps_intragenic_Testing$PREDICTED_INTERGENIC=="True",]
  all = pass_aps_intragenic_Testing
  
  pdf(pdf_name)
  par(mfrow=c(3,4))
  hist(log10(ied$SVLEN), main = 'ied')
  hist(log10(ped$SVLEN), main = 'ped')
  hist(log10(tss$SVLEN), main = 'tss')
  hist(log10(dp$SVLEN), main = 'dp')
  hist(log10(lof$SVLEN), main = 'lof')
  hist(log10(cg$SVLEN), main = 'cg')
  hist(log10(promoter$SVLEN), main = 'promoter')
  hist(log10(utr$SVLEN), main = 'utr')
  hist(log10(inv$SVLEN), main = 'inv')
  hist(log10(intronic$SVLEN), main = 'intronic')
  hist(log10(intergenic$SVLEN), main = 'intergenic')
  hist(log10(all$SVLEN), main = 'all')
  dev.off()
}
plot_ps_residual<-function(pass_aps_intragenic_Testing, pdf_name = 'function_anno.distribution_of_marginal_PS_residual.pdf'){
  ied = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_INTRAGENIC_EXON_DUP),]
  ped = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_PARTIAL_EXON_DUP),]
  tss = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_TSS_DUP),]
  dp = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_DUP_PARTIAL),]
  cg = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_COPY_GAIN),]
  lof = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_LOF),]
  promoter = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_PROMOTER),]
  utr = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_UTR),]
  inv = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_INV_SPAN),]
  intronic = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_INTRONIC),]
  intergenic = pass_aps_intragenic_Testing[pass_aps_intragenic_Testing$PREDICTED_INTERGENIC=="True",]
  all = pass_aps_intragenic_Testing
  
  pdf(pdf_name)
  par(mfrow=c(3,4))
  hist(ied$residual_PS, main = 'ied')
  hist(ped$residual_PS, main = 'ped')
  hist(tss$residual_PS, main = 'tss')
  hist(dp$residual_PS, main = 'dp')
  hist(lof$residual_PS, main = 'lof')
  hist(cg$residual_PS, main = 'cg')
  hist(promoter$residual_PS, main = 'promoter')
  hist(utr$residual_PS, main = 'utr')
  hist(inv$residual_PS, main = 'inv')
  hist(intronic$residual_PS, main = 'intronic')
  hist(intergenic$residual_PS, main = 'intergenic')
  hist(all$residual_PS, main = 'all')
  dev.off()
}
plot_ps<-function(pass_aps_intragenic_Testing, pdf_name = 'function_anno.distribution_of_marginal_PS.pdf'){
  ied = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_INTRAGENIC_EXON_DUP),]
  ped = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_PARTIAL_EXON_DUP),]
  tss = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_TSS_DUP),]
  dp = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_DUP_PARTIAL),]
  cg = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_COPY_GAIN),]
  lof = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_LOF),]
  promoter = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_PROMOTER),]
  utr = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_UTR),]
  inv = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_INV_SPAN),]
  intronic = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_INTRONIC),]
  intergenic = pass_aps_intragenic_Testing[pass_aps_intragenic_Testing$PREDICTED_INTERGENIC=="True",]
  all = pass_aps_intragenic_Testing
  
  pdf(pdf_name)
  par(mfrow=c(3,4))
  hist(ied$PS, main = 'ied')
  hist(ped$PS, main = 'ped')
  hist(tss$PS, main = 'tss')
  hist(dp$PS, main = 'dp')
  hist(lof$PS, main = 'lof')
  hist(cg$PS, main = 'cg')
  hist(promoter$PS, main = 'promoter')
  hist(utr$PS, main = 'utr')
  hist(inv$PS, main = 'inv')
  hist(intronic$PS, main = 'intronic')
  hist(intergenic$PS, main = 'intergenic')
  hist(all$PS, main = 'all')
  dev.off()
}
plot_ps_fitted<-function(pass_aps_intragenic_Testing, pdf_name = 'function_anno.distribution_of_marginal_PS_fitted.pdf'){
  ied = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_INTRAGENIC_EXON_DUP),]
  ped = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_PARTIAL_EXON_DUP),]
  tss = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_TSS_DUP),]
  dp = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_DUP_PARTIAL),]
  cg = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_COPY_GAIN),]
  lof = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_LOF),]
  promoter = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_PROMOTER),]
  utr = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_UTR),]
  inv = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_INV_SPAN),]
  intronic = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_INTRONIC),]
  intergenic = pass_aps_intragenic_Testing[pass_aps_intragenic_Testing$PREDICTED_INTERGENIC=="True",]
  all = pass_aps_intragenic_Testing
  
  pdf(pdf_name)
  par(mfrow=c(3,4))
  hist(ied$fitted_PS, main = 'ied')
  hist(ped$fitted_PS, main = 'ped')
  hist(tss$fitted_PS, main = 'tss')
  hist(dp$fitted_PS, main = 'dp')
  hist(lof$fitted_PS, main = 'lof')
  hist(cg$fitted_PS, main = 'cg')
  hist(promoter$fitted_PS, main = 'promoter')
  hist(utr$fitted_PS, main = 'utr')
  hist(inv$fitted_PS, main = 'inv')
  hist(intronic$fitted_PS, main = 'intronic')
  hist(intergenic$fitted_PS, main = 'intergenic')
  hist(all$fitted_PS, main = 'all')
  dev.off()
}
plot_aps<-function(pass_aps_intragenic_Testing, pdf_name = 'function_anno.distribution_of_marginal_APS.pdf'){
  ied = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_INTRAGENIC_EXON_DUP),]
  ped = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_PARTIAL_EXON_DUP),]
  tss = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_TSS_DUP),]
  dp = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_DUP_PARTIAL),]
  cg = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_COPY_GAIN),]
  lof = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_LOF),]
  promoter = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_PROMOTER),]
  utr = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_UTR),]
  inv = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_INV_SPAN),]
  intronic = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_INTRONIC),]
  intergenic = pass_aps_intragenic_Testing[pass_aps_intragenic_Testing$PREDICTED_INTERGENIC=="True",]
  all = pass_aps_intragenic_Testing
  
  pdf(pdf_name)
  par(mfrow=c(3,4))
  hist(ied$APS, main = 'ied')
  hist(ped$APS, main = 'ped')
  hist(tss$APS, main = 'tss')
  hist(dp$APS, main = 'dp')
  hist(lof$APS, main = 'lof')
  hist(cg$APS, main = 'cg')
  hist(promoter$APS, main = 'promoter')
  hist(utr$APS, main = 'utr')
  hist(inv$APS, main = 'inv')
  hist(intronic$APS, main = 'intronic')
  hist(intergenic$APS, main = 'intergenic')
  hist(all$APS, main = 'all')
  dev.off()
}
remove_overlap_func_annotations<-function(data){
  SVID = data[!is.na(data$PREDICTED_LOF),]$name
  if(nrow(data[!is.na(data$PREDICTED_INTRAGENIC_EXON_DUP) & data$name%in%SVID,])>0){
    data[!is.na(data$PREDICTED_INTRAGENIC_EXON_DUP) & data$name%in%SVID,]$PREDICTED_INTRAGENIC_EXON_DUP = NA
  }
  
  SVID = c(SVID, data[!is.na(data$PREDICTED_INTRAGENIC_EXON_DUP),]$name)
  if(nrow(data[!is.na(data$PREDICTED_PARTIAL_EXON_DUP) & data$name%in%SVID,])>0){
    data[!is.na(data$PREDICTED_PARTIAL_EXON_DUP) & data$name%in%SVID,]$PREDICTED_PARTIAL_EXON_DUP = NA
  }
  
  SVID = c(SVID, data[!is.na(data$PREDICTED_PARTIAL_EXON_DUP),]$name)
  if(nrow(data[!is.na(data$PREDICTED_COPY_GAIN) & data$name%in%SVID,])>0){
    data[!is.na(data$PREDICTED_COPY_GAIN) & data$name%in%SVID,]$PREDICTED_COPY_GAIN = NA
  }
  
  SVID = c(SVID, data[!is.na(data$PREDICTED_COPY_GAIN),]$name)
  if(nrow(data[!is.na(data$PREDICTED_TSS_DUP) & data$name%in%SVID,])>0){
    data[!is.na(data$PREDICTED_TSS_DUP) & data$name%in%SVID,]$PREDICTED_TSS_DUP = NA
  }
  
  SVID = c(SVID, data[!is.na(data$PREDICTED_TSS_DUP),]$name)
  if(nrow(data[!is.na(data$PREDICTED_MSV_EXON_OVERLAP) & data$name%in%SVID,])>0){
    data[!is.na(data$PREDICTED_MSV_EXON_OVERLAP) & data$name%in%SVID,]$PREDICTED_MSV_EXON_OVERLAP = NA
  }
  
  SVID = c(SVID, data[!is.na(data$PREDICTED_MSV_EXON_OVERLAP),]$name)
  if(nrow(data[!is.na(data$PREDICTED_DUP_PARTIAL) & data$name%in%SVID,])>0){
    data[!is.na(data$PREDICTED_DUP_PARTIAL) & data$name%in%SVID,]$PREDICTED_DUP_PARTIAL = NA
  }
  
  SVID = c(SVID, data[!is.na(data$PREDICTED_DUP_PARTIAL),]$name)
  if(nrow(data[!is.na(data$PREDICTED_BREAKEND_EXONIC) & data$name%in%SVID,])>0){
    data[!is.na(data$PREDICTED_BREAKEND_EXONIC) & data$name%in%SVID,]$PREDICTED_BREAKEND_EXONIC = NA
  }
  
  SVID = c(SVID, data[!is.na(data$PREDICTED_BREAKEND_EXONIC),]$name)
  if(nrow(data[!is.na(data$PREDICTED_UTR) & data$name%in%SVID,])>0){
    data[!is.na(data$PREDICTED_UTR) & data$name%in%SVID,]$PREDICTED_UTR = NA
  }
  
  SVID = c(SVID, data[!is.na(data$PREDICTED_UTR),]$name)
  if(nrow(data[!is.na(data$PREDICTED_PROMOTER) & data$name%in%SVID,])>0){
    data[!is.na(data$PREDICTED_PROMOTER) & data$name%in%SVID,]$PREDICTED_PROMOTER = NA
  }
  
  SVID = c(SVID, data[!is.na(data$PREDICTED_PROMOTER),]$name)
  if(nrow(data[!is.na(data$PREDICTED_INTRONIC) & data$name%in%SVID,])>0){
    data[!is.na(data$PREDICTED_INTRONIC) & data$name%in%SVID,]$PREDICTED_INTRONIC = NA
  }
  
  SVID = c(SVID, data[!is.na(data$PREDICTED_INTRONIC),]$name)
  if(nrow(data[!is.na(data$PREDICTED_INV_SPAN) & data$name%in%SVID,])>0){
    data[!is.na(data$PREDICTED_INV_SPAN) & data$name%in%SVID,]$PREDICTED_INV_SPAN = NA
  }
  
  SVID = c(SVID, data[!is.na(data$data$PREDICTED_INV_SPAN),]$name)
  return(data)
}
#add the aps_cate info
reorganize_supp_evi.all_sizes<-function(pass){
  #re-organize SV supportive evidences
  if(nrow(pass[is.na(pass$EVIDENCE),])>0){
    pass[is.na(pass$EVIDENCE),]$EVIDENCE = "BAF"
  }
  
  pass[,ncol(pass)+1] = 'depth_only'
  pass[grepl("PE",pass$EVIDENCE),][,ncol(pass)] = 'others'
  pass[grepl("SR",pass$EVIDENCE),][,ncol(pass)] = 'others'
  pass[pass$EVIDENCE=="SR",][,ncol(pass)] = 'sr_only'
  colnames(pass)[ncol(pass)] = 'EVIDENCE_reorganized'
  # not differentiate supports INS, INV and CPX, since there are too few depth or sr only calls
  pass[pass$SVTYPE%in%c('INV','CPX','INS'),]$EVIDENCE_reorganized = 'others'
  # for deletions and duplication, merge depth -only in simple repeated regions into others
  pass[pass$SVTYPE%in%c('DEL','DUP') & pass$genomic_context=="SR" & pass$EVIDENCE_reorganized=='depth_only',]$EVIDENCE_reorganized = 'others'
  # merge sr-only duplication in SD into others
  pass[pass$SVTYPE%in%c('DUP') & pass$genomic_context=="SD" & pass$EVIDENCE_reorganized=='sr_only',]$EVIDENCE_reorganized = 'others'
  
  #pass[pass$genomic_context%in%c('SR'),]$genomic_context = 'SR'
  #pass[pass$genomic_context%in%c('SD','US','RM'),]$genomic_context = 'US_RM_SD'
  
  pass[pass$EVIDENCE_reorganized=='depth_only',]$size_cate='all_sizes'
  pass[pass$SVTYPE=="DEL" & pass$EVIDENCE_reorganized=='sr_only' ,]$size_cate='over400bp'
  pass[pass$SVTYPE=="DEL" & pass$EVIDENCE_reorganized=='sr_only' & pass$SVLEN<400,]$size_cate='250to400bp'
  pass[pass$SVTYPE=="DEL" & pass$EVIDENCE_reorganized=='sr_only' & pass$SVLEN<250,]$size_cate='under250bp'
  pass[pass$SVTYPE=="DUP" & pass$EVIDENCE_reorganized=='sr_only',]$size_cate='over250bp'
  pass[pass$SVTYPE=="DUP" & pass$EVIDENCE_reorganized=='sr_only' & pass$SVLEN<250,]$size_cate='under250bp'
  pass[pass$SVTYPE=='INS',]$size_cate='over250bp'
  pass[pass$SVTYPE=='INS' & pass$SVLEN<250,]$size_cate='under250bp'
  
  pass[pass$svtype=="INS:ME:ALU" & pass$SVLEN>281,]$SVLEN = 281
  #alu_anno =  add_size_decile_and_ps_by_all_sizes(pass[pass$svtype=="INS:ME:ALU",],    'ALU',   min(as.integer(length(unique(pass[pass$svtype=="INS:ME:ALU",]$SVLEN))/10),  as.integer(nrow(pass[pass$svtype=="INS:ME:ALU",])/100), 100))
  #line1_anno = add_size_decile_and_ps_by_all_sizes(pass[pass$svtype=="INS:ME:LINE1",], 'LINE1', min(as.integer(length(unique(pass[pass$svtype=="INS:ME:LINE1",]$SVLEN))/10),as.integer(nrow(pass[pass$svtype=="INS:ME:LINE1",])/100), 100))
  #sva_anno =  add_size_decile_and_ps_by_all_sizes(pass[pass$svtype=="INS:ME:SVA",],    'SVA',   min(as.integer(length(unique(pass[pass$svtype=="INS:ME:SVA",]$SVLEN))/10),  as.integer(nrow(pass[pass$svtype=="INS:ME:SVA",])/100), 100))
  #inv_anno =  add_size_decile_and_ps_by_all_sizes(pass[pass$SVTYPE=="INV",],'INV', min(as.integer(length(unique(pass[pass$SVTYPE=="INV",]$SVLEN))/10), as.integer(nrow(pass[pass$SVTYPE=="INV",])/100),100))
  #cpx_anno =  add_size_decile_and_ps_by_all_sizes(pass[pass$SVTYPE=="CPX",],'CPX', min(as.integer(length(unique(pass[pass$SVTYPE=="CPX",]$SVLEN))/10), as.integer(nrow(pass[pass$SVTYPE=="CPX",])/100),100))
  
  alu_anno =  add_size_decile_and_ps_by_all_sizes(pass[pass$svtype=="INS:ME:ALU",],    'ALU',   min( as.integer(nrow(pass[pass$svtype=="INS:ME:ALU",])/100), 100))
  line1_anno = add_size_decile_and_ps_by_all_sizes(pass[pass$svtype=="INS:ME:LINE1",], 'LINE1', min(as.integer(nrow(pass[pass$svtype=="INS:ME:LINE1",])/100), 100))
  sva_anno =  add_size_decile_and_ps_by_all_sizes(pass[pass$svtype=="INS:ME:SVA",],    'SVA',   min(  as.integer(nrow(pass[pass$svtype=="INS:ME:SVA",])/100), 100))
  inv_anno =  add_size_decile_and_ps_by_all_sizes(pass[pass$SVTYPE=="INV",],'INV', min( as.integer(nrow(pass[pass$SVTYPE=="INV",])/100),100))
  cpx_anno =  add_size_decile_and_ps_by_all_sizes(pass[pass$SVTYPE=="CPX",],'CPX', min( as.integer(nrow(pass[pass$SVTYPE=="CPX",])/100),100))
  all_anno = rbind(inv_anno, cpx_anno, alu_anno, line1_anno, sva_anno)
  
  for(svtype in c('INS', 'DEL','DUP')){
    for(gc in unique(pass$genomic_context)){
      for(evi in unique(pass$EVIDENCE_reorganized)){
        for(size in unique(pass$size_cate)){
          tmp = pass[pass$SVTYPE==svtype & !pass$svtype%in%c('INS:ME:ALU','INS:ME:LINE1','INS:ME:SVA') & pass$genomic_context==gc & pass$EVIDENCE_reorganized==evi & pass$size_cate==size, ]
          if(nrow(tmp)>0){
            if(evi ==  "depth_only"){
              #tmp_anno = add_size_decile_and_ps_by_all_sizes(tmp, paste(svtype, gc, evi, size, sep=','), min(as.integer(length(unique(tmp$SVLEN))/10), as.integer(nrow(tmp)/200),50))
              tmp_anno = add_size_decile_and_ps_by_all_sizes(tmp, paste(svtype, gc, evi, size, sep=','), min(as.integer(nrow(tmp)/100),50))
            }
            if(evi !=  "depth_only"){
              #tmp_anno = add_size_decile_and_ps_by_all_sizes(tmp, paste(svtype, gc, evi, size, sep=','), min(as.integer(length(unique(tmp$SVLEN))/10), as.integer(nrow(tmp)/200),100))
              tmp_anno = add_size_decile_and_ps_by_all_sizes(tmp, paste(svtype, gc, evi, size, sep=','), min(as.integer(nrow(tmp)/100),100))
            }
            all_anno = rbind(all_anno, tmp_anno)
          }
        }
      }
    }
  }
  
  out = merge(pass, all_anno, by='name',all=T)
  return(out)
}
reorganize_supp_evi.all_sizes.V2<-function(pass){
  #re-organize SV supportive evidences
  if(nrow(pass[is.na(pass$EVIDENCE),])>0){
    pass[is.na(pass$EVIDENCE),]$EVIDENCE = "BAF"
  }
  
  pass[,ncol(pass)+1] = 'depth_only'
  pass[grepl("PE",pass$EVIDENCE),][,ncol(pass)] = 'others'
  pass[grepl("SR",pass$EVIDENCE),][,ncol(pass)] = 'others'
  pass[pass$EVIDENCE=="SR",][,ncol(pass)] = 'sr_only'
  colnames(pass)[ncol(pass)] = 'EVIDENCE_reorganized'
  # not differentiate supports INS, INV and CPX, since there are too few depth or sr only calls
  pass[pass$SVTYPE%in%c('INV','CPX','INS'),]$EVIDENCE_reorganized = 'others'
  # for deletions and duplication, merge depth-only in simple repeated regions into others
  pass[pass$SVTYPE%in%c('DEL','DUP') & pass$genomic_context=="SR" & pass$EVIDENCE_reorganized=='depth_only',]$EVIDENCE_reorganized = 'others'
  # merge sr-only duplication in SD into others
  pass[pass$SVTYPE%in%c('DUP') & pass$genomic_context=="SD" & pass$EVIDENCE_reorganized=='sr_only',]$EVIDENCE_reorganized = 'others'
  
  #pass[pass$genomic_context%in%c('SR'),]$genomic_context = 'SR'
  #pass[pass$genomic_context%in%c('SD','US','RM'),]$genomic_context = 'US_RM_SD'
  
  pass[pass$EVIDENCE_reorganized=='depth_only',]$size_cate='all_sizes'
  pass[pass$EVIDENCE_reorganized=='depth_only' & pass$genomic_context%in%c('US','RM'),]$genomic_context='US_RM'
  
  pass[pass$SVTYPE%in%c('DEL',"DUP") & pass$EVIDENCE_reorganized=='sr_only',]$size_cate='over250bp'
  pass[pass$SVTYPE%in%c('DEL',"DUP") & pass$EVIDENCE_reorganized=='sr_only' & pass$SVLEN<250,]$size_cate='under250bp'
  pass[pass$SVTYPE%in%c('DEL',"DUP") & pass$EVIDENCE_reorganized=='sr_only' & pass$genomic_context%in%c('US','RM'),]$genomic_context = 'US_RM'
  
  pass[pass$SVTYPE%in%c('DEL',"DUP") & pass$EVIDENCE_reorganized=='others' & pass$genomic_context%in%c('US','RM','SD'),]$size_cate='over5kb'
  pass[pass$SVTYPE%in%c('DEL',"DUP") & pass$EVIDENCE_reorganized=='others' & pass$genomic_context%in%c('US','RM','SD') & pass$SVLEN<5000,]$size_cate='1to5Kb'
  pass[pass$SVTYPE%in%c('DEL',"DUP") & pass$EVIDENCE_reorganized=='others' & pass$genomic_context%in%c('US','RM','SD') & pass$SVLEN<1000,]$size_cate='400to1000bp'
  pass[pass$SVTYPE%in%c('DEL',"DUP") & pass$EVIDENCE_reorganized=='others' & pass$genomic_context%in%c('US','RM','SD') & pass$SVLEN<400,]$size_cate='under400bp'
  pass[pass$SVTYPE%in%c('DEL') & pass$EVIDENCE_reorganized=='others' & pass$genomic_context%in%c('US','RM','SD'),]$genomic_context = 'US_RM_SD'
  pass[pass$SVTYPE%in%c('DUP') & pass$EVIDENCE_reorganized=='others' & pass$genomic_context%in%c('US','RM'),]$genomic_context = 'US_RM'
  
  pass[pass$SVTYPE%in%c('DEL',"DUP") & pass$EVIDENCE_reorganized=='others' & pass$genomic_context%in%c('SR'),]$size_cate='over800'
  pass[pass$SVTYPE%in%c('DEL',"DUP") & pass$EVIDENCE_reorganized=='others' & pass$genomic_context%in%c('SR') & pass$SVLEN<800,]$size_cate='250to800bp'
  pass[pass$SVTYPE%in%c('DEL',"DUP") & pass$EVIDENCE_reorganized=='others' & pass$genomic_context%in%c('SR') & pass$SVLEN<250,]$size_cate='under250bp'
  
  pass[pass$SVTYPE=="DUP" & pass$genomic_context=="SD",]$size_cate = 'all_sizes'
  pass[pass$SVTYPE=="DEL" & pass$EVIDENCE_reorganized=='sr_only' & pass$size_cate == 'over250bp' & pass$genomic_context%in%c('SD','SR'),]$genomic_context = 'SD_SR'
  
  pass[pass$SVTYPE=='INS',]$size_cate='over250bp'
  pass[pass$SVTYPE=='INS' & pass$SVLEN<250,]$size_cate='under250bp'
  pass[pass$SVTYPE=='INS' & pass$genomic_context%in%c('US','RM'),]$genomic_context = 'US_RM'
  pass[pass$SVTYPE=='INS' & pass$genomic_context%in%c('SD','SR'),]$size_cate = 'all_sizes'
  #pass[pass$SVTYPE=='INS' & pass$genomic_context%in%c('SD','SD'),]$genomic_context = 'SD_SR'
  
  pass[pass$svtype=="INS:ME:ALU" & pass$SVLEN>281,]$SVLEN = 281

  alu_anno =   add_size_decile_and_ps_by_all_sizes(pass[pass$genomic_context%in%c('US_RM') & pass$svtype=="INS:ME:ALU",],    'ALU,clean',  min(as.integer(nrow(pass[pass$genomic_context%in%c('US_RM') & pass$svtype=="INS:ME:ALU",])/100), 100))
  line1_anno = add_size_decile_and_ps_by_all_sizes(pass[pass$genomic_context%in%c('US_RM') & pass$svtype=="INS:ME:LINE1",], 'LINE1,clean', min(as.integer(nrow(pass[pass$genomic_context%in%c('US_RM') & pass$svtype=="INS:ME:LINE1",])/100), 100))
  sva_anno =   add_size_decile_and_ps_by_all_sizes(pass[pass$genomic_context%in%c('US_RM') & pass$svtype=="INS:ME:SVA",],    'SVA,clean',  min(as.integer(nrow(pass[pass$genomic_context%in%c('US_RM') & pass$svtype=="INS:ME:SVA",])/100), 100))
 
  alu_anno.SD =   add_size_decile_and_ps_by_all_sizes(pass[pass$genomic_context%in%c('SD') & pass$svtype=="INS:ME:ALU",],    'ALU,SD',  min(as.integer(nrow(pass[pass$genomic_context%in%c('SD') & pass$svtype=="INS:ME:ALU",])/100), 100))
  line1_anno.SD = add_size_decile_and_ps_by_all_sizes(pass[pass$genomic_context%in%c('SD') & pass$svtype=="INS:ME:LINE1",], 'LINE1,SD', min(as.integer(nrow(pass[pass$genomic_context%in%c('SD') & pass$svtype=="INS:ME:LINE1",])/100), 100))
  sva_anno.SD =   add_size_decile_and_ps_by_all_sizes(pass[pass$genomic_context%in%c('SD') & pass$svtype=="INS:ME:SVA",],    'SVA,SD',  min(as.integer(nrow(pass[pass$genomic_context%in%c('SD') & pass$svtype=="INS:ME:SVA",])/100), 100))
 
  alu_anno.SR =   add_size_decile_and_ps_by_all_sizes(pass[pass$genomic_context%in%c('SR') & pass$svtype=="INS:ME:ALU",],    'ALU,SR',  min(as.integer(nrow(pass[pass$genomic_context%in%c('SR') & pass$svtype=="INS:ME:ALU",])/100), 100))
  line1_anno.SR = add_size_decile_and_ps_by_all_sizes(pass[pass$genomic_context%in%c('SR') & pass$svtype=="INS:ME:LINE1",], 'LINE1,SR', min(as.integer(nrow(pass[pass$genomic_context%in%c('SR') & pass$svtype=="INS:ME:LINE1",])/100), 100))
  sva_anno.SR =   add_size_decile_and_ps_by_all_sizes(pass[pass$genomic_context%in%c('SR') & pass$svtype=="INS:ME:SVA",],    'SVA,SR',  min(as.integer(nrow(pass[pass$genomic_context%in%c('SR') & pass$svtype=="INS:ME:SVA",])/100), 100))
  
  inv_anno =   add_size_decile_and_ps_by_all_sizes(pass[pass$SVTYPE=="INV",],'INV', min( as.integer(nrow(pass[pass$SVTYPE=="INV",])/100),100))
  cpx_anno =   add_size_decile_and_ps_by_all_sizes(pass[pass$SVTYPE=="CPX",],'CPX', min( as.integer(nrow(pass[pass$SVTYPE=="CPX",])/100),100))
  all_anno = rbind(inv_anno, cpx_anno, alu_anno, line1_anno, sva_anno , alu_anno.SD, line1_anno.SD, sva_anno.SD, alu_anno.SR, line1_anno.SR, sva_anno.SR)
  
  for(svtype in c('INS', 'DEL','DUP')){
    for(gc in unique(pass$genomic_context)){
      for(evi in unique(pass$EVIDENCE_reorganized)){
        for(size in unique(pass$size_cate)){
          tmp = pass[pass$SVTYPE==svtype & !pass$svtype%in%c('INS:ME:ALU','INS:ME:LINE1','INS:ME:SVA') & pass$genomic_context==gc & pass$EVIDENCE_reorganized==evi & pass$size_cate==size, ]
          if(nrow(tmp)>0){
            if(evi ==  "depth_only"){
              #tmp_anno = add_size_decile_and_ps_by_all_sizes(tmp, paste(svtype, gc, evi, size, sep=','), min(as.integer(length(unique(tmp$SVLEN))/10), as.integer(nrow(tmp)/200),50))
              tmp_anno = add_size_decile_and_ps_by_all_sizes(tmp, paste(svtype, gc, evi, size, sep=','), min(as.integer(nrow(tmp)/100),50))
            }
            if(evi !=  "depth_only"){
              #tmp_anno = add_size_decile_and_ps_by_all_sizes(tmp, paste(svtype, gc, evi, size, sep=','), min(as.integer(length(unique(tmp$SVLEN))/10), as.integer(nrow(tmp)/200),100))
              tmp_anno = add_size_decile_and_ps_by_all_sizes(tmp, paste(svtype, gc, evi, size, sep=','), min(as.integer(nrow(tmp)/100),100))
            }
            all_anno = rbind(all_anno, tmp_anno)
          }
        }
      }
    }
  }
  
  out = merge(pass, all_anno, by='name',all=T)
  return(out)
}
reorganize_supp_evi.unique_sizes<-function(pass){
  #re-organize SV supportive evidences
  if(nrow(pass[is.na(pass$EVIDENCE),])>0){
    pass[is.na(pass$EVIDENCE),]$EVIDENCE = "BAF"
  }
  
  pass[,ncol(pass)+1] = 'depth_only'
  pass[grepl("PE",pass$EVIDENCE),][,ncol(pass)] = 'others'
  pass[grepl("SR",pass$EVIDENCE),][,ncol(pass)] = 'others'
  pass[pass$EVIDENCE=="SR",][,ncol(pass)] = 'sr_only'
  colnames(pass)[ncol(pass)] = 'EVIDENCE_reorganized'
  # not differentiate supports INS, INV and CPX, since there are too few depth or sr only calls
  pass[pass$SVTYPE%in%c('INV','CPX','INS'),]$EVIDENCE_reorganized = 'others'
  # for deletions and duplication, merge depth -only in simple repeated regions into others
  pass[pass$SVTYPE%in%c('DEL','DUP') & pass$genomic_context=="SR" & pass$EVIDENCE_reorganized=='depth_only',]$EVIDENCE_reorganized = 'others'
  # merge sr-only duplication in SD into others
  pass[pass$SVTYPE%in%c('DUP') & pass$genomic_context=="SD" & pass$EVIDENCE_reorganized=='sr_only',]$EVIDENCE_reorganized = 'others'
  
  pass[pass$genomic_context%in%c('SR'),]$genomic_context = 'SR'
  pass[pass$genomic_context%in%c('SD','US','RM'),]$genomic_context = 'US_RM_SD'
  
  pass[pass$EVIDENCE_reorganized=='depth_only',]$size_cate='all_sizes'
  pass[pass$SVTYPE=="DEL" & pass$EVIDENCE_reorganized=='sr_only' ,]$size_cate='over400bp'
  pass[pass$SVTYPE=="DEL" & pass$EVIDENCE_reorganized=='sr_only' & pass$SVLEN<400,]$size_cate='250to400bp'
  pass[pass$SVTYPE=="DEL" & pass$EVIDENCE_reorganized=='sr_only' & pass$SVLEN<250,]$size_cate='under250bp'
  pass[pass$SVTYPE=="DUP" & pass$EVIDENCE_reorganized=='sr_only',]$size_cate='over250bp'
  pass[pass$SVTYPE=="DUP" & pass$EVIDENCE_reorganized=='sr_only' & pass$SVLEN<250,]$size_cate='under250bp'
  pass[pass$SVTYPE=='INS',]$size_cate='over250bp'
  pass[pass$SVTYPE=='INS' & pass$SVLEN<250,]$size_cate='under250bp'
  
  pass[pass$svtype=="INS:ME:ALU" & pass$SVLEN>281,]$SVLEN = 281
  alu_anno =  add_size_decile_and_ps_by_unique_sizes(pass[pass$svtype=="INS:ME:ALU",],    'ALU',   min(as.integer(length(unique(pass[pass$svtype=="INS:ME:ALU",]$SVLEN))/10),  as.integer(nrow(pass[pass$svtype=="INS:ME:ALU",])/100), 100))
  line1_anno = add_size_decile_and_ps_by_unique_sizes(pass[pass$svtype=="INS:ME:LINE1",], 'LINE1', min(as.integer(length(unique(pass[pass$svtype=="INS:ME:LINE1",]$SVLEN))/10),as.integer(nrow(pass[pass$svtype=="INS:ME:LINE1",])/100), 100))
  sva_anno =  add_size_decile_and_ps_by_unique_sizes(pass[pass$svtype=="INS:ME:SVA",],    'SVA',   min(as.integer(length(unique(pass[pass$svtype=="INS:ME:SVA",]$SVLEN))/10),  as.integer(nrow(pass[pass$svtype=="INS:ME:SVA",])/100), 100))
  inv_anno =  add_size_decile_and_ps_by_unique_sizes(pass[pass$SVTYPE=="INV",],           'INV',   min(as.integer(length(unique(pass[pass$SVTYPE=="INV",]$SVLEN))/10), as.integer(nrow(pass[pass$SVTYPE=="INV",])/100),100))
  cpx_anno =  add_size_decile_and_ps_by_unique_sizes(pass[pass$SVTYPE=="CPX",],           'CPX',   min(as.integer(length(unique(pass[pass$SVTYPE=="CPX",]$SVLEN))/10), as.integer(nrow(pass[pass$SVTYPE=="CPX",])/100),100))
  all_anno = rbind(inv_anno, cpx_anno, alu_anno, line1_anno, sva_anno)
  
  for(svtype in c('INS', 'DEL','DUP')){
    for(gc in unique(pass$genomic_context)){
      for(evi in unique(pass$EVIDENCE_reorganized)){
        for(size in unique(pass$size_cate)){
          tmp = pass[pass$SVTYPE==svtype & !pass$svtype%in%c('INS:ME:ALU','INS:ME:LINE1','INS:ME:SVA') & pass$genomic_context==gc & pass$EVIDENCE_reorganized==evi & pass$size_cate==size, ]
          if(nrow(tmp)>0){
            if(evi ==  "depth_only"){
              tmp_anno = add_size_decile_and_ps_by_unique_sizes(tmp, paste(svtype, gc, evi, size, sep=','), min(as.integer(length(unique(tmp$SVLEN))/10), as.integer(nrow(tmp)/200),50))
            }
            if(evi !=  "depth_only"){
              tmp_anno = add_size_decile_and_ps_by_unique_sizes(tmp, paste(svtype, gc, evi, size, sep=','), min(as.integer(length(unique(tmp$SVLEN))/10), as.integer(nrow(tmp)/200),100))
            }
            all_anno = rbind(all_anno, tmp_anno)
          }
        }
      }
    }
  }
  
  out = merge(pass, all_anno, by='name',all=T)
  return(out)
}
reorganize_supp_evi.unique_sizes.V2<-function(pass){
  #re-organize SV supportive evidences
  if(nrow(pass[is.na(pass$EVIDENCE),])>0){
    pass[is.na(pass$EVIDENCE),]$EVIDENCE = "BAF"
  }
  
  pass[,ncol(pass)+1] = 'depth_only'
  pass[grepl("PE",pass$EVIDENCE),][,ncol(pass)] = 'others'
  pass[grepl("SR",pass$EVIDENCE),][,ncol(pass)] = 'others'
  pass[pass$EVIDENCE=="SR",][,ncol(pass)] = 'sr_only'
  colnames(pass)[ncol(pass)] = 'EVIDENCE_reorganized'
  # not differentiate supports INS, INV and CPX, since there are too few depth or sr only calls
  pass[pass$SVTYPE%in%c('INV','CPX','INS'),]$EVIDENCE_reorganized = 'others'
  # for deletions and duplication, merge depth-only in simple repeated regions into others
  pass[pass$SVTYPE%in%c('DEL','DUP') & pass$genomic_context=="SR" & pass$EVIDENCE_reorganized=='depth_only',]$EVIDENCE_reorganized = 'others'
  # merge sr-only duplication in SD into others
  pass[pass$SVTYPE%in%c('DUP') & pass$genomic_context=="SD" & pass$EVIDENCE_reorganized=='sr_only',]$EVIDENCE_reorganized = 'others'
  
  #pass[pass$genomic_context%in%c('SR'),]$genomic_context = 'SR'
  #pass[pass$genomic_context%in%c('SD','US','RM'),]$genomic_context = 'US_RM_SD'
  
  pass[pass$EVIDENCE_reorganized=='depth_only',]$size_cate='all_sizes'
  pass[pass$EVIDENCE_reorganized=='depth_only' & pass$genomic_context%in%c('US','RM'),]$genomic_context='US_RM'
  
  pass[pass$SVTYPE%in%c('DEL',"DUP") & pass$EVIDENCE_reorganized=='sr_only',]$size_cate='over250bp'
  pass[pass$SVTYPE%in%c('DEL',"DUP") & pass$EVIDENCE_reorganized=='sr_only' & pass$SVLEN<250,]$size_cate='under250bp'
  pass[pass$SVTYPE%in%c('DEL',"DUP") & pass$EVIDENCE_reorganized=='sr_only' & pass$genomic_context%in%c('US','RM'),]$genomic_context = 'US_RM'
  
  pass[pass$SVTYPE%in%c('DEL',"DUP") & pass$EVIDENCE_reorganized=='others' & pass$genomic_context%in%c('US','RM','SD'),]$size_cate='over5kb'
  pass[pass$SVTYPE%in%c('DEL',"DUP") & pass$EVIDENCE_reorganized=='others' & pass$genomic_context%in%c('US','RM','SD') & pass$SVLEN<5000,]$size_cate='1to5Kb'
  pass[pass$SVTYPE%in%c('DEL',"DUP") & pass$EVIDENCE_reorganized=='others' & pass$genomic_context%in%c('US','RM','SD') & pass$SVLEN<1000,]$size_cate='400to1000bp'
  pass[pass$SVTYPE%in%c('DEL',"DUP") & pass$EVIDENCE_reorganized=='others' & pass$genomic_context%in%c('US','RM','SD') & pass$SVLEN<400,]$size_cate='under400bp'
  pass[pass$SVTYPE%in%c('DEL') & pass$EVIDENCE_reorganized=='others' & pass$genomic_context%in%c('US','RM','SD'),]$genomic_context = 'US_RM_SD'
  pass[pass$SVTYPE%in%c('DUP') & pass$EVIDENCE_reorganized=='others' & pass$genomic_context%in%c('US','RM'),]$genomic_context = 'US_RM'
  
  pass[pass$SVTYPE%in%c('DEL',"DUP") & pass$EVIDENCE_reorganized=='others' & pass$genomic_context%in%c('SR'),]$size_cate='over800'
  pass[pass$SVTYPE%in%c('DEL',"DUP") & pass$EVIDENCE_reorganized=='others' & pass$genomic_context%in%c('SR') & pass$SVLEN<800,]$size_cate='250to800bp'
  pass[pass$SVTYPE%in%c('DEL',"DUP") & pass$EVIDENCE_reorganized=='others' & pass$genomic_context%in%c('SR') & pass$SVLEN<250,]$size_cate='under250bp'
  
  pass[pass$SVTYPE=="DUP" & pass$genomic_context=="SD",]$size_cate = 'all_sizes'
  pass[pass$SVTYPE=="DEL" & pass$EVIDENCE_reorganized=='sr_only' & pass$size_cate == 'over250bp' & pass$genomic_context%in%c('SD','SR'),]$genomic_context = 'SD_SR'
  
  pass[pass$SVTYPE=='INS',]$size_cate='over250bp'
  pass[pass$SVTYPE=='INS' & pass$SVLEN<250,]$size_cate='under250bp'
  pass[pass$SVTYPE=='INS' & pass$genomic_context%in%c('US','RM'),]$genomic_context = 'US_RM'
  pass[pass$SVTYPE=='INS' & pass$genomic_context%in%c('SD','SR'),]$size_cate = 'all_sizes'
  
  pass[pass$svtype=="INS:ME:ALU" & pass$SVLEN>281,]$SVLEN = 281
  
  alu_anno =   add_size_decile_and_ps_by_unique_sizes(pass[pass$genomic_context%in%c('US_RM') & pass$svtype=="INS:ME:ALU",],    'ALU,clean',  min(as.integer(nrow(pass[pass$genomic_context%in%c('US_RM') & pass$svtype=="INS:ME:ALU",])/100), 100))
  line1_anno = add_size_decile_and_ps_by_unique_sizes(pass[pass$genomic_context%in%c('US_RM') & pass$svtype=="INS:ME:LINE1",], 'LINE1,clean', min(as.integer(nrow(pass[pass$genomic_context%in%c('US_RM') & pass$svtype=="INS:ME:LINE1",])/100), 100))
  sva_anno =   add_size_decile_and_ps_by_unique_sizes(pass[pass$genomic_context%in%c('US_RM') & pass$svtype=="INS:ME:SVA",],    'SVA,clean',  min(as.integer(nrow(pass[pass$genomic_context%in%c('US_RM') & pass$svtype=="INS:ME:SVA",])/100), 100))
  
  alu_anno.SD =   add_size_decile_and_ps_by_unique_sizes(pass[pass$genomic_context%in%c('SD') & pass$svtype=="INS:ME:ALU",],    'ALU,SD',  min(as.integer(nrow(pass[pass$genomic_context%in%c('SD') & pass$svtype=="INS:ME:ALU",])/100), 100))
  line1_anno.SD = add_size_decile_and_ps_by_unique_sizes(pass[pass$genomic_context%in%c('SD') & pass$svtype=="INS:ME:LINE1",], 'LINE1,SD', min(as.integer(nrow(pass[pass$genomic_context%in%c('SD') & pass$svtype=="INS:ME:LINE1",])/100), 100))
  sva_anno.SD =   add_size_decile_and_ps_by_unique_sizes(pass[pass$genomic_context%in%c('SD') & pass$svtype=="INS:ME:SVA",],    'SVA,SD',  min(as.integer(nrow(pass[pass$genomic_context%in%c('SD') & pass$svtype=="INS:ME:SVA",])/100), 100))
  
  alu_anno.SR =   add_size_decile_and_ps_by_unique_sizes(pass[pass$genomic_context%in%c('SR') & pass$svtype=="INS:ME:ALU",],    'ALU,SR',  min(as.integer(nrow(pass[pass$genomic_context%in%c('SR') & pass$svtype=="INS:ME:ALU",])/100), 100))
  line1_anno.SR = add_size_decile_and_ps_by_unique_sizes(pass[pass$genomic_context%in%c('SR') & pass$svtype=="INS:ME:LINE1",], 'LINE1,SR', min(as.integer(nrow(pass[pass$genomic_context%in%c('SR') & pass$svtype=="INS:ME:LINE1",])/100), 100))
  sva_anno.SR =   add_size_decile_and_ps_by_unique_sizes(pass[pass$genomic_context%in%c('SR') & pass$svtype=="INS:ME:SVA",],    'SVA,SR',  min(as.integer(nrow(pass[pass$genomic_context%in%c('SR') & pass$svtype=="INS:ME:SVA",])/100), 100))
  
  inv_anno =   add_size_decile_and_ps_by_unique_sizes(pass[pass$SVTYPE=="INV",],'INV', min( as.integer(nrow(pass[pass$SVTYPE=="INV",])/100),100))
  cpx_anno =   add_size_decile_and_ps_by_unique_sizes(pass[pass$SVTYPE=="CPX",],'CPX', min( as.integer(nrow(pass[pass$SVTYPE=="CPX",])/100),100))
  all_anno = rbind(inv_anno, cpx_anno, alu_anno, line1_anno, sva_anno , alu_anno.SD, line1_anno.SD, sva_anno.SD, alu_anno.SR, line1_anno.SR, sva_anno.SR)
  
  for(svtype in c('INS', 'DEL','DUP')){
    for(gc in unique(pass$genomic_context)){
      for(evi in unique(pass$EVIDENCE_reorganized)){
        for(size in unique(pass$size_cate)){
          print(c(svtype, gc, evi, size))
          tmp = pass[pass$SVTYPE==svtype & !pass$svtype%in%c('INS:ME:ALU','INS:ME:LINE1','INS:ME:SVA') & pass$genomic_context==gc & pass$EVIDENCE_reorganized==evi & pass$size_cate==size, ]
          if(nrow(tmp)>0){
            if(evi ==  "depth_only"){
              #tmp_anno = add_size_decile_and_ps_by_unique_sizes(tmp, paste(svtype, gc, evi, size, sep=','), min(as.integer(length(unique(tmp$SVLEN))/10), as.integer(nrow(tmp)/200),50))
              tmp_anno = add_size_decile_and_ps_by_unique_sizes(tmp, paste(svtype, gc, evi, size, sep=','), min(as.integer(nrow(tmp)/100),50))
            }
            if(evi !=  "depth_only"){
              #tmp_anno = add_size_decile_and_ps_by_unique_sizes(tmp, paste(svtype, gc, evi, size, sep=','), min(as.integer(length(unique(tmp$SVLEN))/10), as.integer(nrow(tmp)/200),100))
              tmp_anno = add_size_decile_and_ps_by_unique_sizes(tmp, paste(svtype, gc, evi, size, sep=','), min(as.integer(nrow(tmp)/100),100))
            }
            all_anno = rbind(all_anno, tmp_anno)
          }
        }
      }
    }
  }
  
  out = merge(pass, all_anno, by='name',all=T)
  return(out)
}
#apply lm model
roll_apply_nlm<-function(xvalues, yvalues, bandwidth, aps_label){
  fitted_model = data.frame('xvalues' = 0, 'yvalues' = 0, 'fitted_y' = 0)
  for(i in 1:length(xvalues)){
    tmp_x = xvalues[c(max(0,i-bandwidth):min(length(xvalues), i+bandwidth))]
    tmp_y = yvalues[c(max(0,i-bandwidth):min(length(xvalues), i+bandwidth))]
    
    model <- lm(tmp_y ~ tmp_x)
    # model <- nls(tmp_y ~ b1*tmp_x^2 + b2*tmp_x + b3,start = list(b1=1, b2=1, b3=1))
    #model <- nls(tmp_y ~ b1 * (1 + b2 * tmp_x) ^ b3,start = list(b1 =1,b2 = 1, b3=1))
    fitted_y = predict(model)
    
    i_index = match(xvalues[i], tmp_x)
    fitted_model[i,1] = tmp_x[i_index]
    fitted_model[i,2] = tmp_y[i_index]
    fitted_model[i,3] = fitted_y[i_index]
    
  }
  fitted_model[,4] = aps_label
  return(fitted_model)
}
#Process SNV gene data
read.snvdata <- function(SNVdata.in, gene.metadata.in){
  #Read & clean
  snv.data <- read.table(SNVdata.in, header=T, sep=',')  # sep=','
  metadata <- read.table(gene.metadata.in, header=T, comment.char="")  # comment.char=""
  merged <- merge(x=snv.data, y=metadata, by.x="gene_gencodeV33", by.y="gene", sort=F)
  merged <- merged[which(!(merged$chromosome %in% c("chrX", "chrY"))), ]
  #Assign oe deciles
  #   merged$mis_oe_dec <- ceiling(10*rank(merged$oe_mis_upper)/(nrow(merged)+1))  # mis_oe --> oe_mis_upper
  merged$ptv_oe_dec <- ceiling(10*rank(merged$LOEUF)/(nrow(merged)+1))  # ptv_oe --> LOEUF
  #Assign oe to 40 bins
  #   merged$mis_oe_binrank <- ceiling(40*rank(merged$oe_mis_upper)/(nrow(merged)+1))
  merged$ptv_oe_binrank <- ceiling(40*rank(merged$LOEUF)/(nrow(merged)+1))
  #Assign oe percentiles
  #   merged$mis_oe_cent <- ceiling(100*rank(merged$oe_mis_upper)/(nrow(merged)+1))
  merged$ptv_oe_cent <- ceiling(100*rank(merged$LOEUF)/(nrow(merged)+1))
  #Return formatted data
  return(merged)
}
###Sets sv types & colors
SetSvTypesAndColors <- function(svtypes.file = NULL){
  if(!is.null(svtypes.file)){
    svtypes <- read.table(svtypes.file, sep="\t", header=F, comment.char="", check.names=F)
    svtypes <- as.data.frame(apply(svtypes, 2, as.character))
    colnames(svtypes) <- c("svtype", "color")
  }else{
    require(RColorBrewer, quietly=T)
    svtypes.v <- unique(dat$SVTYPE)
    svtypes.c <- brewer.pal(length(svtypes.v), "Dark2")
    svtypes <- data.frame("svtype"=svtypes.v,
                          "color"=svtypes.c)
  }
  return(svtypes)
}
#for SVs that interrupt mutliple genes, select the most constraint one (lowest LOEUF)
select_most_constraint_gene_single_function<-function(pass_aps_intragenic_Testing, snv_data, function_name){
  function_col = match(function_name, colnames(pass_aps_intragenic_Testing))
  tmp = pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing[,function_col]),c("name", function_name)]
  tmp[,3] = apply(tmp,1,function(x){length(strsplit(as.character(x[2]),',')[[1]])})
  tmp[,4] = tmp[,2]
  colnames(tmp)[4] = gsub('PREDICTED_','',function_name)
  for(i in 1:nrow(tmp)){
    if(tmp[i,3]>1){
      gene_list = strsplit(as.character(tmp[i,2]),',')[[1]]
      gene_constraint = snv_data[snv_data$gene%in%gene_list,]
      gene_most_constraint = gene_constraint[gene_constraint$LOEUF==min(gene_constraint$LOEUF),]$gene[1]
      tmp[i,4] = gene_most_constraint
    }
  }
  pass_aps_intragenic_Testing = merge(pass_aps_intragenic_Testing, tmp[,c(1,4)], by='name', all=T)
  return(pass_aps_intragenic_Testing)
}
select_most_constraint_gene_all_functions<-function(pass_aps_intragenic_Testing, snv_data){
  pass_aps_intragenic_Testing = select_most_constraint_gene_single_function(pass_aps_intragenic_Testing, snv_data, 'PREDICTED_BREAKEND_EXONIC')
  pass_aps_intragenic_Testing = select_most_constraint_gene_single_function(pass_aps_intragenic_Testing, snv_data, 'PREDICTED_LOF')
  pass_aps_intragenic_Testing = select_most_constraint_gene_single_function(pass_aps_intragenic_Testing, snv_data, 'PREDICTED_COPY_GAIN')
  pass_aps_intragenic_Testing = select_most_constraint_gene_single_function(pass_aps_intragenic_Testing, snv_data, 'PREDICTED_TSS_DUP')
  pass_aps_intragenic_Testing = select_most_constraint_gene_single_function(pass_aps_intragenic_Testing, snv_data, 'PREDICTED_DUP_PARTIAL')
  pass_aps_intragenic_Testing = select_most_constraint_gene_single_function(pass_aps_intragenic_Testing, snv_data, 'PREDICTED_INTRAGENIC_EXON_DUP')
  pass_aps_intragenic_Testing = select_most_constraint_gene_single_function(pass_aps_intragenic_Testing, snv_data, 'PREDICTED_PARTIAL_EXON_DUP')
  pass_aps_intragenic_Testing = select_most_constraint_gene_single_function(pass_aps_intragenic_Testing, snv_data, 'PREDICTED_INV_SPAN')
  pass_aps_intragenic_Testing = select_most_constraint_gene_single_function(pass_aps_intragenic_Testing, snv_data, 'PREDICTED_PROMOTER')
  pass_aps_intragenic_Testing = select_most_constraint_gene_single_function(pass_aps_intragenic_Testing, snv_data, 'PREDICTED_UTR')
  pass_aps_intragenic_Testing = select_most_constraint_gene_single_function(pass_aps_intragenic_Testing, snv_data, 'PREDICTED_NEAREST_TSS')
  pass_aps_intragenic_Testing = select_most_constraint_gene_single_function(pass_aps_intragenic_Testing, snv_data, 'PREDICTED_INTRONIC')
  return(pass_aps_intragenic_Testing)
}
aps_by_size <-function(pass_aps_intragenic_Testing){
  stat=data.frame(matrix(c('s1_under250bp','s2_250to400bp','s3_400bpto1Kb','s4_1to5Kb','s5_5to100Kb','s6_over100Kb')))
  colnames(stat)[1] = 'size_cate' 
  
  i=1
  tmp = pass_aps_intragenic_Testing[pass_aps_intragenic_Testing$SVLEN<250,]
  stat[i,2] = mean(tmp$APS)
  
  i=i+1
  tmp = pass_aps_intragenic_Testing[!pass_aps_intragenic_Testing$SVLEN<250 & pass_aps_intragenic_Testing$SVLEN<400,]
  stat[i,2] = mean(tmp$APS)
  
  i=i+1
  tmp = pass_aps_intragenic_Testing[!pass_aps_intragenic_Testing$SVLEN<400 & pass_aps_intragenic_Testing$SVLEN<1000,]
  stat[i,2] = mean(tmp$APS)
  
  i=i+1
  tmp = pass_aps_intragenic_Testing[!pass_aps_intragenic_Testing$SVLEN<1000 & pass_aps_intragenic_Testing$SVLEN<5000,]
  stat[i,2] = mean(tmp$APS)
  
  i=i+1
  tmp = pass_aps_intragenic_Testing[!pass_aps_intragenic_Testing$SVLEN<5000 & pass_aps_intragenic_Testing$SVLEN<100000,]
  stat[i,2] = mean(tmp$APS)
  
  i=i+1
  tmp = pass_aps_intragenic_Testing[!pass_aps_intragenic_Testing$SVLEN<100000,]
  stat[i,2] = mean(tmp$APS)
  
  return(stat)
}
aps_permutation<-function(pass_aps_intragenic_Testing,random_shuffle = 1000, sample_prop = .9){
  aps_shuffle = data.frame("all"=0, "PREDICTED_INTERGENIC" = 0,
                           "PREDICTED_LOF" = 0, "PREDICTED_LOF_DEL" = 0, "PREDICTED_LOF_INS" = 0, "PREDICTED_LOF_oth" = 0, 
                           "PREDICTED_COPY_GAIN"=0,"PREDICTED_INTRAGENIC_EXON_DUP"=0,"PREDICTED_PARTIAL_EXON_DUP"=0,
                           "PREDICTED_DUP_PARTIAL"=0, "PREDICTED_TSS_DUP"=0, "PREDICTED_INV_SPAN"=0,"PREDICTED_COMPLEX_CODING_DISRUPTION"=0,
                           "PREDICTED_NEAREST_TSS"=0,"PREDICTED_PROMOTER"=0,"PREDICTED_UTR"=0,"PREDICTED_INTRONIC"=0)
  pass_aps_intragenic_functions = pass_aps_intragenic_Testing[,c('name','SVTYPE','APS',colnames(aps_shuffle)[colnames(aps_shuffle)%in%colnames(pass_aps_intragenic_Testing)])]
  sample_size = as.integer(nrow(pass_aps_intragenic_functions) *sample_prop)
  for(i in c(1:random_shuffle)){
    print(i)
    sample_index=sample(c(1:nrow(pass_aps_intragenic_functions)), sample_size)
    tmp = pass_aps_intragenic_functions[sample_index,]
    
    j=1 
    aps_shuffle[i,j] = mean(tmp$APS)
    j=2 
    aps_shuffle[i,j] = mean(tmp[tmp$PREDICTED_INTERGENIC=='True',]$APS)
    j=3
    aps_shuffle[i,j] = mean(tmp[!is.na(tmp$PREDICTED_LOF),]$APS)
    j=4 
    aps_shuffle[i,j] = mean(tmp[tmp$SVTYPE=="DEL" & !is.na(tmp$PREDICTED_LOF),]$APS)
    j=5 
    aps_shuffle[i,j] = mean(tmp[tmp$SVTYPE=="INS" & !is.na(tmp$PREDICTED_LOF),]$APS)
    j=6 
    aps_shuffle[i,j] = mean(tmp[!tmp$SVTYPE%in%c("DEL","INS") & !is.na(tmp$PREDICTED_LOF),]$APS)
    
    for(j in c(7:ncol(aps_shuffle))){
      tmp2 = tmp[!is.na(tmp[,colnames(tmp)==colnames(aps_shuffle)[j]]),]
      if(nrow(tmp2) >0){
        aps_shuffle[i,j] = mean(tmp2$APS)
      }
      if(nrow(tmp2)==0){
        aps_shuffle[i,j] = NA
      }
    }
  }
  return(aps_shuffle)
}
plot_aps_permutation<-function(pass_aps_intragenic_Testing, aps_shuffle){
  
  mean_aps = data.frame(aps_shuffle[1,])
  mean_aps[1,1] = mean(pass_aps_intragenic_Testing$APS)
  mean_aps[1,2] = mean(pass_aps_intragenic_Testing[pass_aps_intragenic_Testing$PREDICTED_INTERGENIC=="True",]$APS)
  mean_aps[1,3] = mean(pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing$PREDICTED_LOF),]$APS)
  mean_aps[1,4] = mean(pass_aps_intragenic_Testing[pass_aps_intragenic_Testing$SVTYPE=="DEL" & !is.na(pass_aps_intragenic_Testing$PREDICTED_LOF),]$APS)
  mean_aps[1,5] = mean(pass_aps_intragenic_Testing[pass_aps_intragenic_Testing$SVTYPE=="INS" & !is.na(pass_aps_intragenic_Testing$PREDICTED_LOF),]$APS)
  mean_aps[1,6] = mean(pass_aps_intragenic_Testing[!pass_aps_intragenic_Testing$SVTYPE%in%c("DEL","INS") & !is.na(pass_aps_intragenic_Testing$PREDICTED_LOF),]$APS)
  for(i in c(7:ncol(aps_shuffle))){
    mean_aps[1,i] = mean(pass_aps_intragenic_Testing[!is.na(pass_aps_intragenic_Testing[,colnames(pass_aps_intragenic_Testing)==colnames(aps_shuffle)[i]]),]$APS)
  }
  
  par(mar=c(8,4,4,4))
  plot(c(0,ncol(aps_shuffle)), c(-.2,.2), frame.plot = F, type = 'n', xlab = '', ylab = 'APS', xaxt='n', las=2)
  for(i in c(-2:2)){abline(h=i/10, col='grey')}
  for(i in c(-2:2)){abline(h=(i-.5)/10, col='grey', lty=2)}
  axis(1,c(1:ncol(aps_shuffle)),gsub('PREDICTED_','',colnames(aps_shuffle)), las=2)
  for(i in c(1:ncol(aps_shuffle))){
    aps_range = range(aps_shuffle[,i])
    points(i, mean_aps[1,i], cex=1.5, pch = 20)
    lines(rep(i,2), aps_range, lwd = 1, col = 'black')
  }
}
#readin and organize SNV data including loeuf, triplo- and haplo-sensitivity data
snv_data_readin<-function(SNVdata.in, gene.metadata.in, mane_genes.in, rcnv.in){
  snv.data <- read.snvdata(SNVdata.in, gene.metadata.in)
  mane.genes <- read.table(mane_genes.in, header=F)
  genes <- intersect(sort(unique(as.character(snv.data$gene_gencodeV33))), sort(unique(as.character(mane.genes$V1))))
  
  rcnv_dat=read.table(rcnv.in, header = T, comment.char = "", sep = '\t')
  colnames(rcnv_dat)[1] = 'gene'
  rcnv_dat$ptriplo_oe_cent <- ceiling(100*rank(1-rcnv_dat$pTriplo)/(nrow(rcnv_dat)+1))
  rcnv_dat$phaplo_oe_cent <- ceiling(100*rank(1-rcnv_dat$pHaplo)/(nrow(rcnv_dat)+1))
  
  snv.data = merge(snv.data, rcnv_dat, by='gene')
  return(snv.data)
}
calcu_aps_vs_gene_table<-function(aps_sv,aps_sd_vs_svcount){
  aps_sv_gene_interrupt = data.frame('interrupt' = 0,'aps'=0,'sv_count' = 0)
  aps_sv_gene_interrupt[nrow(aps_sv_gene_interrupt),1] = 'all'
  aps_sv_gene_interrupt[nrow(aps_sv_gene_interrupt),2] = mean(aps_sv$marginal_aps.all_SVLEN)
  aps_sv_gene_interrupt[nrow(aps_sv_gene_interrupt),3] = nrow(aps_sv)
  lof = aps_sv[!is.na(aps_sv$PREDICTED_LOF),]
  aps_sv_gene_interrupt[nrow(aps_sv_gene_interrupt)+1,1] = 'lof'
  aps_sv_gene_interrupt[nrow(aps_sv_gene_interrupt),2] = mean(lof$marginal_aps.all_SVLEN)
  aps_sv_gene_interrupt[nrow(aps_sv_gene_interrupt),3] = nrow(lof)
  copy_gain = aps_sv[!is.na(aps_sv$PREDICTED_COPY_GAIN),]
  aps_sv_gene_interrupt[nrow(aps_sv_gene_interrupt)+1,1] = 'cg'
  aps_sv_gene_interrupt[nrow(aps_sv_gene_interrupt),2] = mean(copy_gain$marginal_aps.all_SVLEN)
  aps_sv_gene_interrupt[nrow(aps_sv_gene_interrupt),3] = nrow(copy_gain)
  ied = aps_sv[!is.na(aps_sv$PREDICTED_INTRAGENIC_EXON_DUP),]
  aps_sv_gene_interrupt[nrow(aps_sv_gene_interrupt)+1,1] = 'ied'
  aps_sv_gene_interrupt[nrow(aps_sv_gene_interrupt),2] = mean(ied$marginal_aps.all_SVLEN)
  aps_sv_gene_interrupt[nrow(aps_sv_gene_interrupt),3] = nrow(ied)
  ped = aps_sv[!is.na(aps_sv$PREDICTED_PARTIAL_EXON_DUP),]
  aps_sv_gene_interrupt[nrow(aps_sv_gene_interrupt)+1,1] = 'ped'
  aps_sv_gene_interrupt[nrow(aps_sv_gene_interrupt),2] = mean(ped$marginal_aps.all_SVLEN)
  aps_sv_gene_interrupt[nrow(aps_sv_gene_interrupt),3] = nrow(ped)
  tss_dup = aps_sv[!is.na(aps_sv$PREDICTED_TSS_DUP),]
  aps_sv_gene_interrupt[nrow(aps_sv_gene_interrupt)+1,1] = 'tss_dup'
  aps_sv_gene_interrupt[nrow(aps_sv_gene_interrupt),2] = mean(tss_dup$marginal_aps.all_SVLEN)
  aps_sv_gene_interrupt[nrow(aps_sv_gene_interrupt),3] = nrow(tss_dup)
  dup_partial = aps_sv[!is.na(aps_sv$PREDICTED_DUP_PARTIAL),]
  aps_sv_gene_interrupt[nrow(aps_sv_gene_interrupt)+1,1] = 'dup_partial'
  aps_sv_gene_interrupt[nrow(aps_sv_gene_interrupt),2] = mean(dup_partial$marginal_aps.all_SVLEN)
  aps_sv_gene_interrupt[nrow(aps_sv_gene_interrupt),3] = nrow(dup_partial)
  utr=  aps_sv[!is.na(aps_sv$PREDICTED_UTR),]
  aps_sv_gene_interrupt[nrow(aps_sv_gene_interrupt)+1,1] = 'utr'
  aps_sv_gene_interrupt[nrow(aps_sv_gene_interrupt),2] = mean(utr$marginal_aps.all_SVLEN)
  aps_sv_gene_interrupt[nrow(aps_sv_gene_interrupt),3] = nrow(utr)
  promoter = aps_sv[!is.na(aps_sv$PREDICTED_PROMOTER),]
  aps_sv_gene_interrupt[nrow(aps_sv_gene_interrupt)+1,1] = 'promoter'
  aps_sv_gene_interrupt[nrow(aps_sv_gene_interrupt),2] = mean(promoter$marginal_aps.all_SVLEN)
  aps_sv_gene_interrupt[nrow(aps_sv_gene_interrupt),3] = nrow(promoter)
  intronic=  aps_sv[!is.na(aps_sv$PREDICTED_INTRONIC),]
  aps_sv_gene_interrupt[nrow(aps_sv_gene_interrupt)+1,1] = 'intronic'
  aps_sv_gene_interrupt[nrow(aps_sv_gene_interrupt),2] = mean(intronic$marginal_aps.all_SVLEN)
  aps_sv_gene_interrupt[nrow(aps_sv_gene_interrupt),3] = nrow(intronic)
  intergenic = aps_sv[aps_sv$PREDICTED_INTERGENIC=="True",]
  aps_sv_gene_interrupt[nrow(aps_sv_gene_interrupt)+1,1] = 'intergenic'
  aps_sv_gene_interrupt[nrow(aps_sv_gene_interrupt),2] = mean(intergenic$marginal_aps.all_SVLEN)
  aps_sv_gene_interrupt[nrow(aps_sv_gene_interrupt),3] = nrow(intergenic)
  
  aps_sv_gene_interrupt[,ncol(aps_sv_gene_interrupt)+1] = sapply(aps_sv_gene_interrupt$sv_count, function(x){10^(log10(x)*aps_sd_vs_svcount$coefficients[2] + aps_sd_vs_svcount$coefficients[1])})
  colnames(aps_sv_gene_interrupt)[ncol(aps_sv_gene_interrupt)] = 'aps_sd'
  return(aps_sv_gene_interrupt)
}
plot_aps_vs_gene<-function(aps_sv_gene_interrupt){
  par(mar=c(8,4,4,4))
  plot(c(0,nrow(aps_sv_gene_interrupt)),c(-0.04,0.08), frame.plot = F, type = 'n', xlab = '', ylab = 'APS', yaxt='n', xaxt='n')
  axis(2,c(-2:4)/50, las=2)
  abline(h=c(-2:4)/50, col='grey')
  abline(h=0)
  aps_sv_gene_interrupt_noncoding = aps_sv_gene_interrupt[aps_sv_gene_interrupt$interrupt%in%c('intronic','intergenic','all'),]
  aps_sv_gene_interrupt_coding = aps_sv_gene_interrupt[!aps_sv_gene_interrupt$interrupt%in%c('all','intergenic','intronic'),]
  aps_sv_gene_interrupt_coding = aps_sv_gene_interrupt_coding[c(1,3,7,4,6,5,2,8),]
  aps_sv_gene_interrupt = rbind(aps_sv_gene_interrupt_coding, aps_sv_gene_interrupt_noncoding)
  for(i in 1:nrow(aps_sv_gene_interrupt)){
    points(i,aps_sv_gene_interrupt[i,]$aps, pch=20, cex=1.5)
    lines(rep(i,2), c(aps_sv_gene_interrupt[i,]$aps - 1.96*aps_sv_gene_interrupt[i,]$aps_sd,  aps_sv_gene_interrupt[i,]$aps + 1.96*aps_sv_gene_interrupt[i,]$aps_sd))
  }
  axis(1,c(1:nrow(aps_sv_gene_interrupt)), labels = aps_sv_gene_interrupt$interrupt, las=2)
}
plot_aps_vs_gene.vertical<-function(aps_sv_gene_interrupt){
  par(mar=c(4,4,4,4))
  plot(c(-0.1,0.1), c(0,-nrow(aps_sv_gene_interrupt)),frame.plot = F, type = 'n', ylab = '', xlab = 'APS', yaxt='n', xaxt='n')
  axis(1,c(-5:5)/50, las=2)
  abline(v=c(-5:5)/50, col='grey')
  abline(b=0)
  aps_sv_gene_interrupt_noncoding = aps_sv_gene_interrupt[aps_sv_gene_interrupt$interrupt%in%c('intronic','intergenic','all'),]
  aps_sv_gene_interrupt_coding = aps_sv_gene_interrupt[!aps_sv_gene_interrupt$interrupt%in%c('all','intergenic','intronic'),]
  aps_sv_gene_interrupt_coding = aps_sv_gene_interrupt_coding[c(1,3,7,4,6,5,2,8),]
  aps_sv_gene_interrupt = rbind(aps_sv_gene_interrupt_coding, aps_sv_gene_interrupt_noncoding)
  for(i in 1:nrow(aps_sv_gene_interrupt)){
    points(aps_sv_gene_interrupt[i,]$aps, -i,pch=20, cex=3)
    lines( c(aps_sv_gene_interrupt[i,]$aps - 1.96*aps_sv_gene_interrupt[i,]$aps_sd,  aps_sv_gene_interrupt[i,]$aps + 1.96*aps_sv_gene_interrupt[i,]$aps_sd),-rep(i,2))
  }
  axis(2,-c(1:nrow(aps_sv_gene_interrupt)), labels = aps_sv_gene_interrupt$interrupt, las=2)
}
#calculate singleton proportion across svtype, genomic context, SV size, and evidence
calcu_singleton_proportion.by_svtype<-function(dat_2){
  sp_by_svtype = data.frame(table(dat_2$SVTYPE))
  sp_by_svtype[,3] = sapply(sp_by_svtype[,1], function(x){nrow(dat_2[dat_2$AC==1 & dat_2$SVTYPE==x,])})
  sp_by_svtype[,4] = sp_by_svtype[,3] /sp_by_svtype[,2]
  colnames(sp_by_svtype)=c('SVTYPE','all_SVs','singletons','singleton_proportion')
  return(sp_by_svtype)
}
calcu_singleton_proportion.by_genomic_context<-function(dat_2){
  sp_by_gc = data.frame(table(dat_2$genomic_context))
  sp_by_gc[,3] = sapply(sp_by_gc[,1], function(x){nrow(dat_2[dat_2$AC==1 & dat_2$genomic_context==x,])})
  sp_by_gc[,4] = sp_by_gc[,3] /sp_by_gc[,2]
  colnames(sp_by_gc)=c('genomic_context','all_SVs','singletons','singleton_proportion')
  return(sp_by_gc)
}
calcu_singleton_proportion.by_size<-function(dat_2){
  dat_2[,ncol(dat_2)+1] = 's0_under_100bp'
  colnames(dat_2)[ncol(dat_2)] = 'size_range'
  dat_2[dat_2$SVLEN>100,]$size_range = 's1_100to250bp'
  dat_2[dat_2$SVLEN>250,]$size_range = 's2_250to500bp'
  dat_2[dat_2$SVLEN>500,]$size_range = 's3_500bpto5Kb'
  dat_2[dat_2$SVLEN>5000,]$size_range = 's4_over5Kb'
  
  sp_by_size = data.frame(table(dat_2$size_range))
  sp_by_size[,3] = sapply(sp_by_size[,1], function(x){nrow(dat_2[dat_2$AC==1 & dat_2$size_range==x,])})
  sp_by_size[,4] = sp_by_size[,3] /sp_by_size[,2]
  colnames(sp_by_size)=c('size_range','all_SVs','singletons','singleton_proportion')
  return(sp_by_size)
}
calcu_singleton_proportion.by_evidence<-function(dat_2){
  sp_by_evi = data.frame('evidence' = 0,'all_SVs' = 0,'singletons' = 0,'singleton_proportion' = 0)
  sp_by_evi[1,1] = 'SR_only'
  sp_by_evi[1,2] = nrow(dat_2[dat_2$EVIDENCE=='SR',])
  sp_by_evi[1,3] = nrow(dat_2[dat_2$AC==1 & dat_2$EVIDENCE=='SR',])
  
  sp_by_evi[2,1] = 'RD_only'
  sp_by_evi[2,2] = nrow(dat_2[dat_2$EVIDENCE=='RD',])
  sp_by_evi[2,3] = nrow(dat_2[dat_2$AC==1 & dat_2$EVIDENCE=='RD',])
  
  sp_by_evi[3,1] = 'other'
  sp_by_evi[3,2] = nrow(dat_2[!dat_2$EVIDENCE%in%c('SR','RD'),])
  sp_by_evi[3,3] = nrow(dat_2[dat_2$AC==1 & !dat_2$EVIDENCE%in%c('SR','RD'),])
  
  sp_by_evi[,4] = sp_by_evi[,3] / sp_by_evi[,2]
  return(sp_by_evi)
}
plot_singleton_proportion.by_svtype<-function(sp_by_svtype, sv_color, cex_num = 2){
  plot(c(.5,nrow(sp_by_svtype)),c(0.3,1), frame.plot = F, type = 'n', xlab = '', ylab = 'Singleton Proportion', xaxt='n', las=2)
  abline(h=c(0:10)/10, col='grey')
  svtype_list = c('DEL','DUP', 'INS','INV', 'CPX', 'CTX')
  axis(1,c(1:length(svtype_list)), svtype_list, cex.axis = 1.5)
  for(i in 1:length(svtype_list)){
    points(i, sp_by_svtype[sp_by_svtype$SVTYPE==svtype_list[i],]$singleton_proportion, pch=20, col = sv_color[sv_color[,1] == svtype_list[i],2], cex = cex_num)
  }
}
plot_singleton_proportion.by_genomic_context<-function(sp_by_gc, gc_color, cex_num = 2){
  plot(c(.5,nrow(sp_by_gc)),c(0.1,0.5), frame.plot = F, type = 'n', xlab = '', ylab = 'Singleton Proportion', xaxt='n', las=2)
  abline(h=c(0:10)/10, col='grey')
  gc_list = c('US','RM','SD','SR')
  axis(1,c(1:length(gc_list)), c('unique\nsequences','repeat\nmasks','segmental\nduplicates','simple\nrepeats'), cex.axis = 1.5)
  for(i in 1:length(gc_list)){
    points(i, sp_by_gc[sp_by_gc$genomic_context==gc_list[i],]$singleton_proportion, pch=20, cex = cex_num, col = gc_color[gc_color[,1] == gc_list[i],2])
  }
}
plot_singleton_proportion.by_size<-function(sp_by_size, cex_num = 2){
  plot(c(.5,nrow(sp_by_size)),c(.3,.45), frame.plot = F, type = 'n', xlab = '', ylab = 'Singleton Proportion', xaxt='n', las=2)
  abline(h=c(0:10)/20, col='grey')
  axis(1,c(1:nrow(sp_by_size)), c('100-250bp','250-500bp','500bp-5Kb','>5Kb'))
  for(i in 1:nrow(sp_by_size)){
    points(i, sp_by_size[i,]$singleton_proportion, pch=20, cex = cex_num)
  }
}
calcu_gc_color <-function(){
  colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  gc_color = data.frame('genomc_context' = 0,'color' = 0)
  gc_color[1,1] = 'US'
  gc_color[2,1] = 'RM'
  gc_color[3,1] = 'SD'
  gc_color[4,1] = 'SR'
  gc_color[,2] = colorBlindGrey8[c(3,6,7,8)]
  return(gc_color)
}
calcu_ps_by_size_table<-function(tmp3, bin_count, svtype, gc, evi){
  bin_size = nrow(tmp3)/bin_count
  tmp3[,ncol(tmp3)+1] = as.integer((c(1:nrow(tmp3))-1)/bin_size)
  colnames(tmp3)[ncol(tmp3)] = 'size_bin'
  ps_by_size = data.frame(table(tmp3$size_bin))
  ps_by_size[,3] = sapply(ps_by_size[,1], function(x){nrow(tmp3[tmp3$size_bin==x & tmp3$AC==1,]) / nrow(tmp3[tmp3$size_bin==x,])})
  ps_by_size[,4] = sapply(ps_by_size[,1], function(x){min(tmp3[tmp3$size_bin==x,]$SVLEN)})
  ps_by_size[,5] = sapply(ps_by_size[,1], function(x){max(tmp3[tmp3$size_bin==x,]$SVLEN)})
  ps_by_size[,ncol(ps_by_size)+1] = svtype
  ps_by_size[,ncol(ps_by_size)+1] = gc
  ps_by_size[,ncol(ps_by_size)+1] = evi
  return(ps_by_size)
}
plot_ps_across_cate<-function(dat_2){
  svtype_list = c('DEL','DUP')
  gc_list = unique(dat_2$genomic_context)
  evidence_list = c('sr_only','depth_only','other')
  pdf('../../calcu_aps/singleton_prop.by_svtype_gc_size.legend.pdf')
  plot(c(0,1),c(0,1),frame.plot = F, type = 'n', xlab = '', ylab = '', xaxt='n', yaxt='n')
  legend('center',c('US','RM','SD','SR'), col = c("#56B4E9", "#0072B2", "#F0E442", "#D55E00"), pch=20, lwd=2, bty='n')
  dev.off()
  pdf('../../calcu_aps/singleton_prop.by_svtype_gc_size.DEL_DUP.pdf', height = 4, width = 8)
  par(mfrow=c(length(svtype_list), length(evidence_list)))
  par(mar=c(1,4,2,4))
  for(svtype in svtype_list){
    tmp1 = dat_2[dat_2$SVTYPE==svtype,]
    for(evi in evidence_list){
      if(evi == 'sr_only'){
        tmp3 = tmp1[tmp1$EVIDENCE=='SR',]
      }
      if(evi == 'depth_only'){
        tmp3 = tmp1[tmp1$EVIDENCE=='RD',]
      }
      if(evi == 'other'){
        tmp3 = tmp1[!tmp1$EVIDENCE%in%c('SR','RD'),]
      }
      
      tmp3 = tmp3[!is.na(tmp3$SVLEN),]
      tmp3 = tmp3[order(tmp3$SVLEN),]
      
      bin_count = 30
      colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
      
      gc = 'US'
      tmp2 = tmp3[tmp3$genomic_context==gc,]
      ps_by_size = calcu_ps_by_size_table(tmp2, bin_count, svtype, gc, evi)
      for(gc in c('RM','SD','SR')){
        tmp2 = tmp3[tmp3$genomic_context==gc,]
        ps_by_size.tmp = calcu_ps_by_size_table(tmp2, bin_count, svtype, gc, evi)
        ps_by_size = rbind(ps_by_size,ps_by_size.tmp)}
      ps_by_size[,ncol(ps_by_size)+1] = "#56B4E9"
      ps_by_size[ps_by_size[,7] == 'RM',][,ncol(ps_by_size)] = "#0072B2"
      ps_by_size[ps_by_size[,7] == 'SD',][,ncol(ps_by_size)] = "#F0E442"
      ps_by_size[ps_by_size[,7] == 'SR',][,ncol(ps_by_size)] = "#D55E00"
      colnames(ps_by_size)[ncol(ps_by_size)] = 'color'
      
      plot(c(0,bin_count), range(ps_by_size[,3]),frame.plot = F, type = 'n', xlab = 'size', ylab = 'PS', xaxt='n', las=2, main = paste(svtype, evi, sep=','))
      axis(1, ps_by_size[,4], las=2)
      for(i in c('US','RM','SD','SR')){
        lines(ps_by_size[ps_by_size[,7] == i,1], ps_by_size[ps_by_size[,7] == i,3], col = ps_by_size[ps_by_size[,7] == i,]$color, pch=20, type = 'b')
      }
    }
  }
  dev.off()
  pdf('../../calcu_aps/singleton_prop.by_svtype_gc_size.non_DEL_DUP.pdf', height = 4, width = 8)
  par(mfrow=c(2,3))
  par(mar=c(1,4,2,4))
  for(svtype in c('INS')){
    tmp1 = dat_2[dat_2$svtype==svtype,]
    for(evi in c('sr_only','other')){
      if(evi == 'sr_only'){
        tmp3 = tmp1[tmp1$EVIDENCE=='SR',]
      }
      if(evi == 'other'){
        tmp3 = tmp1[!tmp1$EVIDENCE%in%c('SR'),]
      }
      
      tmp3 = tmp3[!is.na(tmp3$SVLEN),]
      tmp3 = tmp3[order(tmp3$SVLEN),]
      
      bin_count = 30
      colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
      
      gc = 'US'
      tmp2 = tmp3[tmp3$genomic_context==gc,]
      ps_by_size = calcu_ps_by_size_table(tmp2, bin_count, svtype, gc, evi)
      for(gc in c('RM','SD','SR')){
        tmp2 = tmp3[tmp3$genomic_context==gc,]
        ps_by_size.tmp = calcu_ps_by_size_table(tmp2, bin_count, svtype, gc, evi)
        ps_by_size = rbind(ps_by_size,ps_by_size.tmp)}
      ps_by_size[,ncol(ps_by_size)+1] = "#56B4E9"
      ps_by_size[ps_by_size[,7] == 'RM',][,ncol(ps_by_size)] = "#0072B2"
      ps_by_size[ps_by_size[,7] == 'SD',][,ncol(ps_by_size)] = "#F0E442"
      ps_by_size[ps_by_size[,7] == 'SR',][,ncol(ps_by_size)] = "#D55E00"
      colnames(ps_by_size)[ncol(ps_by_size)] = 'color'
      
      plot(c(0,bin_count), range(ps_by_size[,3]),frame.plot = F, type = 'n', xlab = 'size', ylab = 'PS', xaxt='n', las=2, main = paste(svtype, evi, sep=','))
      axis(1, ps_by_size[,4], las=2)
      for(i in c('US','RM','SD','SR')){
        lines(ps_by_size[ps_by_size[,7] == i,1], ps_by_size[ps_by_size[,7] == i,3], col = ps_by_size[ps_by_size[,7] == i,]$color, pch=20, type = 'b')
      }
    }
  }
  for(svtype in c('INS:ME:ALU','INS:ME:LINE1','INS:ME:SVA')){
    tmp1 = dat_2[dat_2$svtype==svtype,]
    tmp3 = tmp1
    tmp3 = tmp3[!is.na(tmp3$SVLEN),]
    tmp3 = tmp3[order(tmp3$SVLEN),]
    
    bin_count = 30
    colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    
    gc = 'US'
    tmp2 = tmp3[tmp3$genomic_context==gc,]
    ps_by_size = calcu_ps_by_size_table(tmp2, bin_count, svtype, gc, evi)
    for(gc in c('RM','SD','SR')){
      tmp2 = tmp3[tmp3$genomic_context==gc,]
      ps_by_size.tmp = calcu_ps_by_size_table(tmp2, bin_count, svtype, gc, evi)
      ps_by_size = rbind(ps_by_size,ps_by_size.tmp)}
    ps_by_size[,ncol(ps_by_size)+1] = "#56B4E9"
    ps_by_size[ps_by_size[,7] == 'RM',][,ncol(ps_by_size)] = "#0072B2"
    ps_by_size[ps_by_size[,7] == 'SD',][,ncol(ps_by_size)] = "#F0E442"
    ps_by_size[ps_by_size[,7] == 'SR',][,ncol(ps_by_size)] = "#D55E00"
    colnames(ps_by_size)[ncol(ps_by_size)] = 'color'
    
    plot(c(0,bin_count), range(ps_by_size[,3]),frame.plot = F, type = 'n', xlab = 'size', ylab = 'PS', xaxt='n', las=2, main = svtype)
    axis(1, ps_by_size[,4], las=2)
    for(i in c('US','RM','SD','SR')){
      lines(ps_by_size[ps_by_size[,7] == i,1], ps_by_size[ps_by_size[,7] == i,3], col = ps_by_size[ps_by_size[,7] == i,]$color, pch=20, type = 'b')
    }
  }
  dev.off()
}


#!R
library("optparse")

option_list = list(
        make_option(c("--SV_color"), type="character", default=NULL, help="color for each type of SVs", metavar="character"),
        make_option(c("--genomic_context"), type="character", default=NULL, help="SV sites information", metavar="character"),
        make_option(c("--SV_function_prediction"), type="character", default=NULL, help="SV sites with predicted functional consequences", metavar="character"),
        make_option(c('-s', "--SV_info"), type="character", default=NULL, help="SV sites information", metavar="character"),
        make_option(c('-p', "--prefix"), type="character", default=NULL, help="prefix for output file", metavar="character"),
        make_option(c('-f', "--output_folder"), type="character", default=NULL, help="folder for output files", metavar="character")
 );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

SV_color = opt$SV_color
SV_info = opt$SV_info
genomic_context = opt$genomic_context
output_folder = opt$output_folder
prefix = opt$prefix
SV_function_prediction = opt$SV_function_prediction
#readin the SV type and colors
print('read in SV colors ...')
colors = read.table(SV_color, header = F, comment.char = "")
svtypes <- SetSvTypesAndColors(SV_color)

#readin the SV sites with AF info
print('read in SV sites ...')
dat=read.table(SV_info, header = T, sep = '\t', comment.char = "")
print('read in SVID vs. genomic context ...')
gc=read.table(genomic_context, header = F)
colnames(gc) = c('name', 'genomic_context')
print('merge genomic context into SV info ...')
dat = merge(dat, gc, by='name')
## restrict to PASS autosomal SVs without small depth only DUPs
print('restrict to autosomal PASS events ...')
dat_2 = dat[dat$FILTER=="PASS" & !is.na(dat$AF) & !dat$X.chrom%in%c('chrX','chrY') & !is.na(dat$SVLEN),]
print('remove small depth only DUPs ...')
sm_depth_only_dup = dat[dat$SVTYPE=="DUP" & dat$ALGORITHMS=="depth" & dat$SVLEN<20000,]
dat_2 = dat_2[!dat_2$name%in%sm_depth_only_dup$name,] #remove depth only small DUPs

#explore confounding factors of SV freq
print('plot correlation of SV frequencies vs. sv type, size and genomic context ...')
sp_by_svtype = calcu_singleton_proportion.by_svtype(dat_2)
sp_by_gc = calcu_singleton_proportion.by_genomic_context(dat_2)
sp_by_size = calcu_singleton_proportion.by_size(dat_2)
sp_by_evi = calcu_singleton_proportion.by_evidence(dat_2)
gc_color = calcu_gc_color()
pdf(paste(prefix, 'singleton_prop.by_svtype_gc_size.pdf', sep='.'), height = 4,width = 12)
par(mfrow=c(1,3))
par(mar=c(3,6,3,3))
plot_singleton_proportion.by_svtype(sp_by_svtype, colors, 5)
plot_singleton_proportion.by_genomic_context(sp_by_gc, gc_color, 5)
plot_singleton_proportion.by_size(sp_by_size[2:nrow(sp_by_size),], 5)
dev.off()

#readin the SV sites with functional annotation info
func_anno=read.table(SV_function_prediction, header=T, comment.char="")
func_anno= func_anno[,c("name", "PREDICTED_BREAKEND_EXONIC", "PREDICTED_COPY_GAIN", "PREDICTED_DUP_PARTIAL", "PREDICTED_INTERGENIC", "PREDICTED_INTRAGENIC_EXON_DUP", "PREDICTED_INTRONIC", "PREDICTED_INV_SPAN", "PREDICTED_LOF", "PREDICTED_MSV_EXON_OVERLAP", "PREDICTED_NEAREST_TSS", "PREDICTED_PARTIAL_EXON_DUP", "PREDICTED_PROMOTER", "PREDICTED_TSS_DUP", "PREDICTED_UTR")]
func_anno = remove_overlap_func_annotations(func_anno)
if(!file.exists(output_folder)){
  dir.create(output_folder)
}
#build APS model, sort SVs by all SV length
print('build APS model based on the length of all SVs ...')
pass_aps_intragenic_Testing.all_SVs = build_marginal_APS_model.by_all_SVs(dat_2,  func_anno, aps_model_prefix = paste(output_folder,'APS_models.all_SVLEN.', sep='/'), output_file_name = paste(prefix, '.marginal_APS_annotated.by_all_SVLEN.bed', sep='') )
#build APS model, sort SVs by unique SV length
print('build APS model based on the length of unique SVs ...')
pass_aps_intragenic_Testing.uniq_SVs = build_marginal_APS_model.by_uniq_SVs(dat_2, func_anno, aps_model_prefix = paste(output_folder, 'APS_models.uniq_SVLEN.', sep='/'), output_file_name = paste(prefix, '.marginal_APS_annotated.by_unique_SVLEN.bed', sep='') )
#combine both models to generate final output:
print('combine APS from both models ...')
aps_both = merge(pass_aps_intragenic_Testing.all_SVs[,c("name", "APS")], pass_aps_intragenic_Testing.uniq_SVs[,c("name", "APS")], by='name')
colnames(aps_both) = c('name','marginal_aps.all_SVLEN','marginal_aps.uniq_SVLEN')
write.table(aps_both, paste(prefix, 'SVID_aps', sep='.'), quote = F, sep = '\t', col.names = T, row.names = F)


#calculate sd of aps by sv count
print('calculate deviation of SVs by SV count ...')
sv_count_list = c(c(1:10)*10, c(2:10)*100, c(2:10)*1000,c(2:10)*10000,c(2:10)*100000)
aps_sd_stat=data.frame(sv_count_list)
for(i in c(2:5)){aps_sd_stat[,i] = 0}
colnames(aps_sd_stat)[c(2:5)] = c('APS.all_SVLEN','APS_sd.all_SVLEN','APS.uniq_SVLEN','APS_sd.uniq_SVLEN')
permu_round = 1000
for(i in 1:nrow(aps_sd_stat)){
  print(i)
  tmp.all_SVLEN = c()
  tmp.uniq_SVLEN = c()
  for(j in c(1:permu_round)){
    tmp = aps_both[sample(1:nrow(aps_both), aps_sd_stat[i,1]),]
    tmp.all_SVLEN = c(tmp.all_SVLEN, mean(tmp$marginal_aps.all_SVLEN))
    tmp.uniq_SVLEN = c(tmp.uniq_SVLEN, mean(tmp$marginal_aps.uniq_SVLEN))
  }
  aps_sd_stat[i,2] = mean(tmp.all_SVLEN)
  aps_sd_stat[i,3] = sd(tmp.all_SVLEN)
  aps_sd_stat[i,4] = mean(tmp.uniq_SVLEN)
  aps_sd_stat[i,5] = sd(tmp.uniq_SVLEN)
  
}
write.table(aps_sd_stat, paste(prefix, 'aps_deviation.tsv', sep='.'), quote = F, sep = '\t', col.names = T, row.names = F)
#aps_sd_vs_svcount=lm(log10(aps_sd_stat$APS_sd.all_SVLEN) ~ log10(aps_sd_stat$sv_count_list))










