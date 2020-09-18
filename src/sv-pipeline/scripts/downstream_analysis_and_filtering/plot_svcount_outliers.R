getJitter <- function(x,dens){
  jit <- dens$y[head(which(dens$x>=x),1)]
  return(jit)
}
jitterX <- function(y,dens){
  return(jitter(1,amount=getJitter(y,dens)))
}


library(optparse)
###List of command-line options
option_list <- list(
  make_option(c("-i","--input"), type="character", 
              help="stat table with count of SVs per sample"),
  make_option(c("-m","--max"), type="integer",
              help="max number of IQR to plot")
 )


###Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog INFILE OUTDIR",
                                option_list=option_list),
                   positional_arguments=TRUE)

opts <- args$options
dat=read.table(opts$input, header = T)
output_prefix = gsub('.txt','',opts$input)

svtype=unique(dat$svtype)
x_max= opts$max
for(i in svtype){

  tmp = dat[dat$svtype==i,]
  quantil_list=quantile(tmp$count)
  iqr = IQR(tmp$count)
  
  pdf(paste(output_prefix,'.',i,'.pdf',sep=''))
  par(fig=c(0,1,.4,1))
  par(mar=c(2,4,2,4))
  plot(c(0,x_max), c(max(min(tmp$count), quantil_list[2]-x_max*iqr), min(max(tmp$count),quantil_list[4]+x_max*iqr) ), frame.plot = F, type = 'n', xlab = 'n * IQR', ylab = 'Count of SVs', las=2, xaxt='n', main=i)
  axis(1)
  
  abline(h=quantil_list[3], col='grey', lty=1,lwd=2)
  for(n in c(0:x_max)){
    abline(h=quantil_list[4]+n*iqr, col='grey', lty=2)
    abline(h=quantil_list[2]-n*iqr, col='grey', lty=2)
  }
  
  sample_counts = data.frame('iqr'=0,'sample_counts'=0)
  rec=0
  for(j in c(1:x_max)){
    cff = c(quantil_list[2]-j*iqr, quantil_list[4]+j*iqr)
    tmp2=tmp[tmp$count<cff[2] & tmp$count>cff[1],]
    
    vals=tmp2$count
    dens <- density(vals, na.rm=T)
    df_x = sapply(vals,function(y){jitterX(y,dens)})
    df_y = vals
    dens$y <- 0.35*dens$y/max(dens$y,na.rm=T)
    polygon(x=c(-dens$y,rev(dens$y))+j,y=c(dens$x,rev(dens$x)),
            border=adjustcolor("darkgreen",alpha=0.3),col=NA, yaxt='n', lwd=2)
    points(df_x + j-1,df_y,pch=21,lwd=1.5,cex=0.1,col="darkgreen")
    
    rec=rec+1
    sample_counts[rec,1]=j
    sample_counts[rec,2]=nrow(tmp2)
  }
  
  par(fig=c(0,1,0,.4),new=T)
  par(mar=c(4,4,2,4), new=T)
  plot(c(0,x_max), c(min(sample_counts[,2]),max(sample_counts[,2])), frame.plot = F, type = 'n', xlab = 'n * IQR', ylab = '# samples kept', las=2, xaxt='n')
  axis(1)
  lines(sample_counts[,1],sample_counts[,2], pch=18, lwd=2, type='b', cex=2)
  dev.off()
}

