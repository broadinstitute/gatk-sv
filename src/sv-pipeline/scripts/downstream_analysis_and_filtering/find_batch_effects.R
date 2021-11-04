#!/usr/bin/env Rscript

# Talkowski SV pipeline downstream analysis helper script

# Make list of all nonredundant pairs of batches from an input list of batches


###Set global parameters
options(stringsAsFactors=F,scipen=1000)
allpops <- c("AFR","ASN","EUR","HSP")
svtypes <- c("DEL","DUP","MCNV","INS","INV","CPX")


###################
###HELPER FUNCTIONS
###################
#Allele frequency correlation plot between datasets
plot.AFcorr <- function(dat,batch1=NULL,batch2=NULL,title=NULL,axlims){
  #Prepare plot
  AF.pairs <- data.frame(dat$b1.AF,dat$b2.AF)
  #Artificially assign non-zero AFs for sites that appear in 0 samples
  AF.pairs[which(AF.pairs[,1]==0),1] <- 1/900
  AF.pairs[which(AF.pairs[,2]==0),2] <- 1/900
  logscale.all <- log10(as.numeric(unlist(sapply(c(0:9),function(i){(1:9)*(10^i)}))))
  logscale.major <- 0:9
  major.labels <- sapply(logscale.major,function(i){expression(paste(i^"th"))})
  par(mar=c(3.2,1.5,1.5,3.2))
  plot(x=log10(c(min(AF.pairs[,1],na.rm=T),1)),
       y=log10(c(min(AF.pairs[,2],na.rm=T),1)),
       type="n",xaxt="n",yaxt="n",xlab="",ylab="",
       xlim=log10(axlims),ylim=log10(axlims))
  axis(1,at=-logscale.all,labels=NA,tck=-0.015,lwd=0.7)
  axis(1,at=-logscale.major,labels=NA,tck=-0.03,lwd=1.1)
  mtext(1,text=bquote("log"[10] ~ "("*.(batch1) ~ "AF)"),line=2)
  axis(4,at=-logscale.all,labels=NA,tck=-0.015,lwd=0.7)
  axis(4,at=-logscale.major,labels=NA,tck=-0.03,lwd=1.1)
  mtext(4,text=bquote("log"[10] ~ "("*.(batch2) ~ "AF)"),line=2)
  sapply(-logscale.major,function(i){
    # axis(1,at=i,labels=bquote('10'^.(i)),tick=F,line=-0.6,cex.axis=0.8)
    # axis(4,at=i,labels=bquote('10'^.(i)),tick=F,line=-0.4,cex.axis=0.8,las=2)
    axis(1,at=i,labels=i,tick=F,line=-0.6,cex.axis=0.8)
    axis(4,at=i,labels=i,tick=F,line=-0.4,cex.axis=0.8,las=2)
  })
  abline(h=log10(1/900),v=log10(1/900),lty=3)
  axis(1,at=log10(1/900),labels="AC=0",line=-0.8,cex.axis=0.7,tick=F)
  axis(4,at=log10(1/900),labels="AC=0",line=-0.6,cex.axis=0.7,tick=F,las=2)
  mtext(3,text=title,line=0.1,font=2)
  
  #Add points
  pt.cex <- 0.4
  alpha <- 0.25
  points(x=log10(AF.pairs[which(dat$bonf.p>=0.05),1]),
         y=log10(AF.pairs[which(dat$bonf.p>=0.05),2]),
         pch=19,cex=pt.cex,lwd=0,
         col=adjustcolor("gray50",alpha=alpha))
  points(x=log10(AF.pairs[which(dat$bonf.p<0.05),1]),
         y=log10(AF.pairs[which(dat$bonf.p<0.05),2]),
         pch=19,cex=pt.cex,lwd=0,
         col=adjustcolor("red",alpha=alpha))
  
  #Add stats
  abline(lm(AF.pairs[,1] ~ AF.pairs[,2]),col="gray10",lty=2)
  AB.cor <- format(round(cor(AF.pairs[,1],AF.pairs[,2])^2,3),nsmall=3)
  text(x=par("usr")[1],y=par("usr")[4]-(0.085*(par("usr")[4]-par("usr")[3])),
       labels=bquote(italic(R)^2 == .(AB.cor)),cex=1.4,pos=4)
}


###Read command-line arguments
args <- commandArgs(trailingOnly=T)
infile <- as.character(args[1])
batch1 <- as.character(args[2])
batch2 <- as.character(args[3])
OUTPREFIX <- as.character(args[4])

# #Dev parameters:
# infile <- "~/scratch/gnomAD_AF_table.merged.txt.gz"
# infile <- "~/scratch/gnomAD_v2_SV_MASTER.gnomAD_v2_SV_PCRPLUS_Q1_batch_1_vs_gnomAD_v2_SV_PCRMINUS_Q2_batch_1.AF_comparison_table.txt.gz"
# batch1 <- "gnomAD_v2_SV_PCRPLUS_Q1_batch_1"
# batch2 <- "gnomAD_v2_SV_PCRMINUS_Q2_batch_1"
# OUTPREFIX <- "~/scratch/gnomAD_v2_SV_MASTER"


###Process input data
dat <- read.table(infile,header=T,sep="\t",comment.char="")
dat$bonf.p <- p.adjust(dat$chisq.p,method="bonferroni")
write.table(dat,paste(OUTPREFIX,".",batch1,"_vs_",batch2,".freq_table_wBonferroni.txt",sep=""),
            col.names=T,row.names=F,sep="\t",quote=F)


###Write list of significant batch effect variants
bad.vars <- dat$VID[which(dat$bonf.p<0.05)]
write.table(bad.vars,paste(OUTPREFIX,".",batch1,"_vs_",batch2,".batch_effect_variants.txt",sep=""),
            col.names=F,row.names=F,sep="\t",quote=F)


###Plot AF correlations, one per SVTYPE
axlims <- c(1/900,1)
png(paste(OUTPREFIX,".",batch1,"_vs_",batch2,".AF_correlation_scatterplot.ALL.png",sep=""),
    height=6*300,width=6*300,res=400)
plot.AFcorr(dat=dat,batch1=batch1,batch2=batch2,title="All SV",axlims=axlims)
dev.off()
sapply(svtypes,function(svtype){
  subdat <- dat[grep(svtype,dat$VID,fixed=T),]
  if(nrow(subdat)>0){
    png(paste(OUTPREFIX,".",batch1,"_vs_",batch2,".AF_correlation_scatterplot.",svtype,".png",sep=""),
        height=6*300,width=6*300,res=400)
    plot.AFcorr(dat=subdat,batch1=batch1,batch2=batch2,title=svtype,axlims=axlims)
    dev.off()
  }
})

