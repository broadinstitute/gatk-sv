#!R
#script to calculate cwas scores

library("optparse")

option_list = list(
        make_option(c('-i', "--input"), type="character", default=NULL, help="input, list of files to integrate", metavar="character"),
        make_option(c('-o', "--output"), type="character", default=NULL, help="name of output file", metavar="character")
 );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


integrate_sv_vs_nc<-function(sv_vs_nc_files){
        nc_element = gsub('.stat','',gsub('gnomAD_SV_v3.vs.','',sv_vs_nc_files[1]))
        sv_vs_nc=read.table(sv_vs_nc_files[1])
        colnames(sv_vs_nc)=c(nc_element,'name')
        for(i in c(2:length(sv_vs_nc_files))){
                print(i)
                nc_tmp = gsub('.stat','',gsub('gnomAD_SV_v3.vs.','',sv_vs_nc_files[i]))
                sv_tm = read.table(sv_vs_nc_files[i])
                colnames(sv_tm)=c(nc_tmp,'name')
                sv_vs_nc = merge(sv_vs_nc, sv_tm, by='name', all=T)
        }
        return(sv_vs_nc)
}


sv_vs_nc_files = read.table(opt$input)
sv_vs_nc_files_a = sv_vs_nc_files[,1]
print(sv_vs_nc_files_a[1])
sv_vs_nc_a = integrate_sv_vs_nc(sv_vs_nc_files_a)
write.table(sv_vs_nc_a, opt$output, quote=F, sep='\t', col.names=T, row.names=F)


