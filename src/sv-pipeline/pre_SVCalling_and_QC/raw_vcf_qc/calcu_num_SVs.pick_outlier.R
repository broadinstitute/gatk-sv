#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
#
# Copyright © 2018 Xuefang Zhao <XZHAO12@mgh.harvard.edu>
# Distributed under terms of the MIT license.

library("optparse")

generate_outlier_matrics<-function(dat){
    samples = unique(dat[,4])
    outlier = data.frame('#CHROM'=0,'SVTYPE'=0,'Mean'=0,'Median'=0,'STD'=0,'Outlier_Sample'=0,'Outlier_Number'=0,'Outlier_Cate'=0)
    rec=0
    for(i in unique(dat[,1])){
        for(j in unique(dat[,2])){
            tmp = dat[dat[,1]==i & dat[,2]==j, ]
            outlier_name = c()
            outlier_nums = c()
            if(nrow(tmp)< length(samples)){
                for(x in samples[!samples%in%tmp[,4]]){
                    outlier_name = c(outlier_name, x)
                    outlier_nums = c(outlier_nums, 0)
                }
            outlier_tmp = tmp[abs(tmp[,3] - median(tmp[,3])) > 5*sd(tmp[,3]), ]
            if( nrow(outlier_tmp[!is.na(outlier_tmp[,1]),])>0 ){
                outlier_name = c(outlier_name, outlier_tmp[,4])
                outlier_nums = c(outlier_nums, outlier_tmp[,3])
            }
                print(length(outlier_name))
                print(length(outlier_nums))
            for(k in 1:length(outlier_name)){
                print(c(k, outlier_name[k]))
                rec=rec+1
                outlier[rec,1]=i
                outlier[rec,2]=j
                outlier[rec,3]=mean(tmp[,3])
                outlier[rec,4]=median(tmp[,3])
                outlier[rec,5]=sd(tmp[,3])
                outlier[rec,6]=outlier_name[k]
                outlier[rec,7]=outlier_nums[k]
                outlier[rec,8]='low'
                if (outlier_nums[k] > median(tmp[,3])){
                    outlier[rec,8]='high'
                }
            }
            }
        }
    }
    return(outlier)
    }

write_output<-function(outlier, outname){
    colnames(outlier)[1] = '#CHROM'
    out1 = outlier[outlier[,8]=='low',]
    out2 = outlier[outlier[,8]=='high',]
    write.table(out1, paste(outname, 'low', sep='.'), quote=F, sep='\t', col.names=T, row.names=F )
    write.table(out1, paste(outname, 'high', sep='.'), quote=F, sep='\t', col.names=T, row.names=F )
}

option_list = list(
            make_option(c("-s", "--statfile"), type="character", default=NULL,
                  help="name of stats concatinated from all samples", metavar="character"),
            make_option(c("-o", "--outname"), type="character", default=NULL,
                  help="name of output file", metavar="character")
     );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

dat=read.table(opt$statfile)
outlier = generate_outlier_matrics(dat)
write_output(outlier, opt$outname)



