#!R
#plot count of SVs stat

#!/usr/bin/env Rscript
library("optparse")

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="input file, should include 5 columns representing sample_id, count_SVs, chr, svtype and size_range", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="name of output folder")
  );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

input=opt$input
output=opt$output

stat=read.table(input, header = T)

size_list=c('all_sizes','50to500bp','500bpto5Kb', '5Kbto1Mb', 'over1Mb')
chr_list=c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY')
stat$CHR=factor(stat$CHR, levels = chr_list)
svtype_list=unique(stat$SVTYPE)

for(svtype in svtype_list){
  pdf(paste(paste(output,'Count_SVs_Per_sample.boxplot',sep='/'),gsub('.stat','',tail(strsplit(input,'/')[[1]],1)),svtype,'pdf',sep = '.'))
  par(mfrow=c(5,1))
  par(mar=c(3,4,3,2))
  for(size in size_list){
    tmp = stat[stat$SIZE_Range==size & stat$SVTYPE==svtype,]
    if(nrow(tmp)>0){
      boxplot(tmp$Count_SV ~ tmp$CHR, pch=20, col='blue', main=size, frame.plot=F, ylab=paste('Count - ',svtype), xaxt='n')
      axis(1,c(1:length(chr_list)), labels = chr_list, las=2)
    }
    else{
      boxplot(tmp$Count_SV ~ tmp$CHR, pch=20, col='blue', main=size, frame.plot=F, ylim=c(0,1), ylab=paste('Count - ',svtype))
      axis(1,c(1:length(chr_list)), labels = chr_list, las=2)
    }
  }
  dev.off()

  pdf(paste(paste(output, 'Count_SVs_Per_sample.lineplot', sep='/'),gsub('.stat','',tail(strsplit(input,'/')[[1]],1)),svtype,'pdf',sep = '.'))
  par(mfrow=c(5,1))
  par(mar=c(3,4,3,2))
  for(size in size_list){
    tmp = stat[stat$SIZE_Range==size & stat$SVTYPE==svtype,]
    plot(c(1,length(chr_list)), c(0,max(tmp$Count_SV)), xlab = '', ylab = paste('Count -', svtype), xaxt='n', type = 'n', frame.plot = F)
    axis(1,c(1:length(chr_list)), labels = chr_list, las=2)

    chr_count = data.frame(matrix(ncol = 2, nrow = length(chr_list)))
    colnames(chr_count)=c('CHR','Count_SV')
    rownames(chr_count)=chr_list
    for(sample in unique(stat$Sample_ID)){
      tmp2=tmp[tmp$Sample_ID==sample,c('CHR','Count_SV')]
      if(nrow(tmp2)<nrow(chr_count)){
        tmp3=chr_count[!rownames(chr_count)%in%tmp2$CHR,]
        for(k in 1:nrow(tmp3)){
          tmp2[nrow(tmp2)+1,1]=rownames(tmp3)[k]
          tmp2[nrow(tmp2),2]=0
        }
      }
      tmp2[,1]=factor(tmp2[,1], levels = chr_list)
      tmp2=tmp2[order(tmp2[,1]),]
      lines(c(1:length(chr_list)), tmp2[,2], type = 'b', pch=20, cex=.5)
    }
  }
  dev.off()
}



