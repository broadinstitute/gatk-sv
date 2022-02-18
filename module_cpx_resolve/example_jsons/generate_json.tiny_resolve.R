setwd('~/Google Drive/Talkowski_Lab/gnomAD_V3/module10_cpx_unresolved_manual/')
mani=read.table('../manifests_final/gnomAD-SV_v3.master_sample_metadata.pre_00c_qc.gatksv_sample_id.callers', header = T)
mani_batch=read.table('../manifests_final/module00c_01_03_04_04rerun_04b_results.batch.tsv', header = T)
header=read.table('Step1_TinyResolve/tiny_resolve.header.json', quote = "", sep = '\t')
for(i in unique(mani$batch)){
  mani_tmp = mani[mani$batch==i,]
  mani_tmp[,ncol(mani_tmp)+1] = mani_batch[mani_batch$batch==i,]$pe
  tmp = header
  tmp[,1]=as.character(tmp[,1])
  tmp[nrow(tmp)+1,1]=paste('    "TinyResolve.samples" : [',paste(paste('"',mani_tmp[,1],'"',sep = ''), collapse = ',\n        '),'],',sep = '')
  tmp[nrow(tmp)+1,1]=paste('    "TinyResolve.manta_vcfs" : [',paste(paste('"',mani_tmp$std_manta,'"',sep = ''), collapse = ',\n        '),'],',sep = '')
  tmp[nrow(tmp)+1,1]=paste('    "TinyResolve.discfile" : [',paste(paste('"',mani_tmp[,ncol(mani_tmp)],'"',sep = ''), collapse = ',\n        '),']',sep = '')
  tmp[nrow(tmp)+1,1]='}'
  write.table(tmp, paste('tiny_resolve.',i,'.json',sep = ''), quote = F, sep = '\t', col.names = F, row.names = F)
}


