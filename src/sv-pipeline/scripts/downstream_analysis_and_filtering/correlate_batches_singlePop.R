#!/usr/bin/env Rscript

# Talkowski SV pipeline downstream analysis helper script

# Compute & plot AF correlation coefficient


###Set global parameters
options(stringsAsFactors=F,scipen=1000)
svtypes <- c("DEL","DUP","INS","MCNV","INV","CPX","BND")
allpops <- c("AFR","ASN","EUR","HSP")


###################
###HELPER FUNCTIONS
###################
#Compute correlation coefficient for a single ancestry for a single pair of batches
correlate.batchpair.byPop <- function(dat,batch1,batch2,pop){
  #Subset data for each batch:
  #1: restrict to sites with >0 AC in at least one batch
  #2: exclude sites on chrX or chrY (if any exist)
  b1.dat <- dat[,c(1:3,grep(batch1,colnames(dat),fixed=T))]
  b2.dat <- dat[,c(1:3,grep(batch2,colnames(dat),fixed=T))]
  allosome.VIDs.toexclude <- unique(c(grep("_X_",b1.dat$VID,fixed=T),grep("_Y_",b1.dat$VID,fixed=T)))
  if(length(allosome.VIDs.toexclude)>0){
    b1.dat <- b1.dat[-allosome.VIDs.toexclude,]
    b2.dat <- b2.dat[-allosome.VIDs.toexclude,]
  }
  b1.AC <- b1.dat[,grep(paste(pop,"AC",sep="_"),colnames(b1.dat),fixed=T)]
  b1.AN <- b1.dat[,grep(paste(pop,"AN",sep="_"),colnames(b1.dat),fixed=T)]
  b2.AC <- b2.dat[,grep(paste(pop,"AC",sep="_"),colnames(b2.dat),fixed=T)]
  b2.AN <- b2.dat[,grep(paste(pop,"AN",sep="_"),colnames(b2.dat),fixed=T)]
  nonzero.sites <- which(b1.AC>0 | b2.AC>0)
  b1.AF <- b1.AC[nonzero.sites]/b1.AN[nonzero.sites]
  b2.AF <- b2.AC[nonzero.sites]/b2.AN[nonzero.sites]
  if(length(b1.AF)>0 & length(b2.AF)>0){
    return(cor(b1.AF,b2.AF,use="pairwise.complete.obs")^2)
  }else{
    return(NA)
  }
}
#Compute correlation coefficient for a single ancestry for all pairs of batches (returns 36x36 matrix)
correlate.allbatches.byPop <- function(dat,batches.list,pop){
  res <- sapply(batches.list$batch,function(batch1){
    sapply(batches.list$batch,function(batch2){
      correlate.batchpair.byPop(dat=dat,batch1=batch1,batch2=batch2,pop=pop)
    })
  })
}


#################
###PLOT FUNCTIONS
#################
#Plot heatmap of single AF correlation matrix
plot.heat <- function(mat,title=NULL,cex.labels=1){
  #Prep plot parameters
  n.batches <- nrow(mat)
  col.pal <- colorRampPalette(c("#440154","#365C8C","#25A584","#FDE725"))(101)
  mat <- floor(mat*100)+1
  #Prep plot area
  par(mar=c(0.1,4.5,4.5,0.1))
  plot(x=c(0,-n.batches),y=c(0,-n.batches),type="n",
       xaxt="n",xlab="",xaxs="i",yaxt="n",ylab="",yaxs="i")
  #Plot heatmap
  sapply(1:n.batches,function(i){
    rect(xleft=-i+1,xright=-i,
         ybottom=seq(-n.batches,-1),
         ytop=seq(-n.batches+1,0),
         col=col.pal[mat[,i]],border=col.pal[mat[,i]])
  })
  #Add batch labels
  revised.labels <- rownames(mat)
  revised.labels <- gsub("gnomAD_v2_SV_","",revised.labels)
  revised.labels <- gsub("PCRPLUS","PCR+",revised.labels)
  revised.labels <- gsub("PCRMINUS","PCR-",revised.labels)
  revised.labels <- gsub("_Q"," ",revised.labels)
  revised.labels <- gsub("_batch_",".",revised.labels)
  label.colors <- rep("#E01925",times=n.batches)
  label.colors[grep("PCR+",revised.labels,fixed=T)] <- "#1F6AAB"
  sapply(1:n.batches,function(i){
    axis(2,at=-n.batches+i-0.5,tick=F,line=-0.8,las=2,cex.axis=cex.labels*0.7,
         labels=revised.labels[i],col.axis=label.colors[i])
    axis(3,at=-i+0.5,tick=F,line=-0.8,las=2,cex.axis=cex.labels*0.7,
         labels=revised.labels[i],col.axis=label.colors[i])
  })
  axis(3,at=par("usr")[1],hadj=1,line=-1,tick=F,
       labels=bquote("R"^2 * ""["AF"]))
  mtext(3,line=3,text=title)
  box()
}
#Plot heatmap key
plot.heat.key <- function(){
  col.pal <- colorRampPalette(c("#440154","#365C8C","#25A584","#FDE725"))(101)
  par(mar=c(0.5,2.5,2.5,2.5))
  plot(x=c(0,1),y=c(0,100),type="n",
       xaxt="n",xlab="",xaxs="i",yaxt="n",ylab="",yaxs="i")
  rect(xleft=0,xright=1,
       ybottom=seq(0,99),ytop=seq(1,100),
       col=col.pal,border=col.pal)
  abline(h=seq(0,100,20),col="white",lwd=0.5,lty=2)
  axis(4,at=seq(0,100,20),tck=-0.02,labels=NA)
  axis(4,at=seq(0,100,20),tick=F,las=2,
       labels=seq(0,1,0.2),line=-0.4,cex.axis=0.8)
  mtext(3,text=bquote("R"^2 * ""["AF"]),line=0.1)
  box()
}
#Generate jitter residuals for a vector of values based on their density
sina.jitter <- function(vals){
  d <- density(vals)
  dv <- approx(d$x,d$y,xout=vals)
  dv <- dv$y/max(dv$y)
  dv.j <- sapply(1:length(vals),function(i){
    jitter(vals[i],amount=dv[i])
  })
  return(dv.j)
}
#Generate sina points to add to existing plot
sina.plot <- function(vals,y.at,color,width=0.25){
  j <- (width*sina.jitter(vals=rep(0,times=length(vals))))+y.at
  points(x=vals,y=j,pch=21,cex=0.25,col=color,bg=adjustcolor(color,alpha=0.3))
  segments(x0=median(vals),x1=median(vals),
           y0=y.at-width,y1=y.at+width)
  points(x=median(vals),y=y.at,pch=18)
}
#Plot sina's of AF correlations per batch
plot.vert.sina <- function(mat,colors,title=NULL){
  #Prep plot area
  n.batches <- nrow(mat)
  par(bty="n",mar=c(0.5,3.5,3,1))
  plot(x=c(0,1),y=c(0,-n.batches),type="n",
       xaxt="n",xlab="",yaxt="n",ylab="")
  abline(v=seq(0,1,0.2),col="gray80")
  axis(3,at=axTicks(1),labels=NA)
  axis(3,at=axTicks(1),tick=F,line=-0.5,cex.axis=0.8,
       labels=round(axTicks(1),1))
  mtext(3,line=0.75,text=bquote("R"^2 * ""["AF"]),cex=0.8)
  mtext(3,line=2,text=title)
  #Add sina plots
  sapply(1:ncol(mat),function(i){
    sina.plot(vals=mat[,i],y.at=-i+0.5,color=colors[i])
    axis(4,at=-i+0.5,cex.axis=0.7,labels=format(round(median(mat[,i]),2),nsmall=2),line=-1.5,tick=F,las=2)
  })
  #Add batch labels
  revised.labels <- rownames(mat)
  revised.labels <- gsub("gnomAD_v2_SV_","",revised.labels)
  revised.labels <- gsub("PCRPLUS","PCR+",revised.labels)
  revised.labels <- gsub("PCRMINUS","PCR-",revised.labels)
  revised.labels <- gsub("_Q"," ",revised.labels)
  revised.labels <- gsub("_batch_",".",revised.labels)
  sapply(1:ncol(mat),function(i){
    axis(2,at=-i+0.5,tick=F,line=-1.25,las=2,cex.axis=0.7,
         labels=revised.labels[i],col.axis=colors[i])
  })
}



###Read command-line arguments
args <- commandArgs(trailingOnly=T)
batches.list.in <- as.character(args[1])
freq.table.in <- as.character(args[2])
pop <- as.character(args[3])
OUTPREFIX <- as.character(args[4])

# #Dev parameters
# batches.list.in <- "~/scratch/gnomAD_v2_SV_MASTER.batches.list"
# freq.table.in <- "~/scratch/gnomAD_v2_SV_MASTER.merged_AF_table.txt.gz"
# pop <- "ASN"
# OUTPREFIX <- "~/scratch/gnomAD_v2_SV_MASTER.ASN"


###Process input data
batches.list <- data.frame("batch"=read.table(batches.list.in,header=F,sep="\t")[,1])
dat <- read.table(freq.table.in,header=T,sep="\t",comment.char="")
colnames(dat)[1:3] <- c("VID","SVLEN","SVTYPE")
#One correlation matrix per svtype
cor.mats <- lapply(svtypes,function(svtype){
  correlate.allbatches.byPop(dat=dat[which(dat$SVTYPE==svtype),],batches.list=batches.list,pop=pop)
})
all.mat <- correlate.allbatches.byPop(dat=dat,batches.list=batches.list,pop=pop)
cor.mats <- c(list(all.mat),cor.mats)
names(cor.mats) <- c("ALL",svtypes)


###Write correlation coefficient matrices to file
sapply(1:length(cor.mats),function(i){
  m <- cor.mats[[i]]
  m <- cbind(data.frame("batch"=rownames(m)),m)
  rownames(m) <- NULL
  write.table(m,paste(OUTPREFIX,pop,names(cor.mats)[i],"R2_matrix.txt",sep="."),
              col.names=T,row.names=F,sep="\t",quote=F)
})


###Plot correlation coefficient matrix heatmaps
#One heatmap per SVTYPE
sapply(1:length(cor.mats),function(i){
  pdf(paste(OUTPREFIX,pop,names(cor.mats)[i],"R2_matrix_heatmap.pdf",sep="."),
      height=5,width=(7/6)*5)
  layout(matrix(c(1,1,1,3,2,4),byrow=F,nrow=3),
         heights=c(2,2,2),widths=c(6,1))
  plot.heat(mat=cor.mats[[i]],title=names(cor.mats)[i])
  plot.heat.key()
  dev.off()
})
#One strip of heatmaps across all SVTYPES
pdf(paste(OUTPREFIX,pop,"R2_matrix_heatmap_strip.pdf",sep="."),
    height=3,width=3*length(cor.mats)*(((4*length(cor.mats))+1)/(4*length(cor.mats))))
layout(matrix(seq(1,length(cor.mats)+1),byrow=F,nrow=1),
       widths=c(rep(4,length(cor.mats)),1))
sapply(1:length(cor.mats),function(i){
  plot.heat(mat=cor.mats[[i]],title=names(cor.mats)[i],cex.labels=0.7)
})
plot.heat.key()
dev.off()


###Plot sina's of cross-batch correlation coefficients per batch
perbatch.colors <- rep("#E01925",times=nrow(batches.list))
perbatch.colors[grep("PCRPLUS",rownames(cor.mats[[1]]),fixed=T)] <- "#1F6AAB"
sapply(1:length(cor.mats),function(i){
  pdf(paste(OUTPREFIX,pop,names(cor.mats)[i],"perBatch_R2_sina_plot.pdf",sep="."),
      height=5,width=6)
  plot.vert.sina(mat=cor.mats[[i]],
                 colors=perbatch.colors,
                 title=names(cor.mats)[[i]])
  dev.off()
})

