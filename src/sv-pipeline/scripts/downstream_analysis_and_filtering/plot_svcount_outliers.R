library(optparse)

getJitter <- function(x,dens){
  jit <- dens$y[head(which(dens$x>=x),1)]
  return(jit)
}
jitterX <- function(y,dens){
  return(jitter(1,amount=getJitter(y,dens)))
}

###List of command-line options
option_list <- list(
  make_option(c("-i","--input"), type="character", 
              help="stat table with count of SVs per sample"),
  make_option(c("-s","--svcount"), type="character", 
              help="stat table with count of overall SVs sites"),
  make_option(c("-o","--outputpath"), type="character", 
              help="output folder for all plots"),  
  make_option(c("-m","--max"), type="integer",
              help="max number of IQR to plot")
 )


###Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog INFILE OUTDIR",
                                option_list=option_list),
                   positional_arguments=TRUE)

opts <- args$options
dat=read.table(opts$input, header = T)
all_counts=read.table(opts$svcount)
all_counts=all_counts[all_counts[,1]>0,]
output_prefix = gsub('.txt','',opts$input)

svtype=unique(dat$svtype)
x_max= opts$max
for(i in svtype){
  pdf(paste(opts$outputpath,'/',output_prefix,'.',i,'.pdf',sep=''))

  svtype_counts = dat[dat$svtype==i,]
  quantile_list=quantile(svtype_counts$count)
  iqr = IQR(svtype_counts$count)

  par(fig=c(0,1,.6,1))
  par(mar=c(2,4,2,4))
  plot(c(0,x_max), c(max(min(svtype_counts$count), quantile_list[2]-x_max*iqr), min(max(svtype_counts$count),quantile_list[4]+x_max*iqr) ), frame.plot = F, type = 'n', xlab = 'n * IQR', ylab = 'Count of SVs', las=2, xaxt='n', main=i)
  axis(1)
  
  abline(h=quantile_list[3], col='grey', lty=1,lwd=2)
  for(n in c(0:x_max)){
    abline(h=quantile_list[4]+n*iqr, col='grey', lty=2)
    abline(h=quantile_list[2]-n*iqr, col='grey', lty=2)
  }
  
  sample_counts = data.frame('iqr'=0,'sample_counts'=0)
  rec=0
  for(j in 1:x_max){
    cutoffs = c(quantile_list[2]-j*iqr, quantile_list[4]+j*iqr)
    counts_outside_cutoffs=svtype_counts[svtype_counts$count<cutoffs[2] & svtype_counts$count>cutoffs[1],]
    
    vals=counts_outside_cutoffs$count
    dens <- density(vals, na.rm=T)
    df_x = sapply(vals,function(y){jitterX(y,dens)})
    df_y = vals
    dens$y <- 0.35*dens$y/max(dens$y,na.rm=T)
    polygon(x=c(-dens$y,rev(dens$y))+j,y=c(dens$x,rev(dens$x)),
            border=adjustcolor("darkgreen",alpha=0.3),col=NA, yaxt='n', lwd=2)
    points(df_x + j-1,df_y,pch=21,lwd=1.5,cex=0.1,col="darkgreen")
    
    rec=rec+1
    sample_counts[rec,1]=j
    sample_counts[rec,2]=nrow(counts_outside_cutoffs)
  }
  
  par(fig=c(0,1,.3,.6), new=T)
  counts_for_svtype = all_counts[all_counts[,4]==i,]
  plot(c(0,x_max), c(min(counts_for_svtype[,3]),max(counts_for_svtype[,3])), frame.plot = F, type = 'n', xlab = 'n * IQR', ylab = '# SV Sites', las=2, xaxt='n')
  axis(1)
  lines(counts_for_svtype[,1],counts_for_svtype[,3], pch=18, lwd=2, type='b', cex=2)
  
  par(fig=c(0,1,0,.3),new=T)
  par(mar=c(4,4,2,4), new=T)
  plot(c(0,x_max), c(min(counts_for_svtype[,2]),max(counts_for_svtype[,2])), frame.plot = F, type = 'n', xlab = 'n * IQR', ylab = '# Outlier Samples', las=2, xaxt='n')
  lines(counts_for_svtype[,1],counts_for_svtype[,2], pch=18, lwd=2, type='b', cex=2)
  axis(1)
  lines(counts_for_svtype[,1],counts_for_svtype[,2], pch=18, lwd=2, type='b', cex=2)
  dev.off()
}

