#!/usr/bin/env Rscript

# Helper script to determine categories of SV that
# don't meet a minimum Mendelian violation rate (MVR)


###Set master parameters
options(stringsAsFactors=F,scipen=1000)


###################
###HELPER FUNCTIONS
###################
#Clean input data
cleanDat <- function(dat){
  #Drop variant IDs
  if(colnames(dat)[1] %in% c("#VID","VID","X.VID")){
    dat <- dat[,-1]
  }
  #Drop variants with no trios considered
  dat <- dat[which(dat$PCRMINUS_TRIOS_CONSIDERED>0),]
  # #Bin by size
  # dat$SVLEN.category <- NA
  # dat$SVLEN.category[which(dat$SVLEN<1000)] <- "SMALL"
  # dat$SVLEN.category[which(dat$SVLEN>=1000 & dat$SVLEN<5000)] <- "MEDIUM"
  # dat$SVLEN.category[which(dat$SVLEN>=5000)] <- "LARGE"
  dat$DNR <- dat$PCRMINUS_APPARENT_DE_NOVO/dat$PCRMINUS_TRIOS_CONSIDERED
  dat$MVR <- 1-(dat$PCRMINUS_MENDELIAN/(dat$PCRMINUS_TRIOS_CONSIDERED-dat$PCRMINUS_INCOMPLETE_GENOTYPE_TRIOS-dat$PCRMINUS_NO_VARIANT_TRIOS))
  #Return
  return(dat)
}
#Classify data as high or low MVR
classifyDatByMVR <- function(dat,cutoff=0.05,deNovo.only=F){
  if(deNovo.only==T){
    dat$MVR <- dat$APPARENT_DE_NOVO/dat$TRIOS_CONSIDERED
  }else{
    dat$MVR <- 1-(dat$MENDELIAN/dat$TRIOS_CONSIDERED)
  }
  dat$highMVR <- 0
  dat$highMVR[which(dat$MVR>cutoff)] <- 1
  dat$highMVR <- as.factor(dat$highMVR)
  return(dat)
}
#Filter data given a set of input conditions
filterDat <- function(dat,SVTYPE=NULL,FILTER=NULL,
                      SVLEN.category=NULL){
  dat.filt <- dat
  if(!is.null(SVTYPE)){
    dat.filt <- dat.filt[which(dat.filt$SVTYPE==SVTYPE),]
  }
  if(!is.null(FILTER)){
    dat.filt <- dat.filt[grep(FILTER,dat.filt$FILTER,fixed=T),]
  }
  if(!is.null(SVLEN.category)){
    dat.filt <- dat.filt[which(dat.filt$SVLEN.category==SVLEN.category),]
  }
  return(dat.filt)
}
#Calculate fraction of sites with high MVR
calcHighMVRFrac <- function(dat){
  if(nrow(dat)>0){
    return(length(which(dat$highMVR==1))/nrow(dat))
  }else{
    return(NA)
  }
}
#Calculate minimum median non-ref sample GQ to obtain target max MVR
optimizeNonrefGQ <- function(dat,step=5,target=0.05){
  m <- 0
  MVR <- calcHighMVRFrac(dat[which(dat$MEDIAN_NONREF_GQ>m),])
  while(!is.na(MVR) & MVR>target & m<999){
    m <- m+step
    MVR <- calcHighMVRFrac(dat[which(dat$MEDIAN_NONREF_GQ>m),])
  }
  if(m>999){
    m <- 999
  }
  return(m)
}


################
###RSCRIPT BLOCK
################
require(optparse,quietly=T)
###List of command-line options
option_list <- list(
  make_option(c("-p","--prefix"), type="character", default="lowQual_designation_script",
              help="Prefix to append to output files [default %default]"),
  make_option(c("-m","--min.variants"), type="integer", default=50,
              help="Minimum number of variants in order to consider category [default %default]"),
  make_option(c("-s","--step"), type="integer", default=5,
              help="Increments of median non-ref GQ to titrate [default %default]"),
  make_option(c("-d","--deNovo.only"), type="logical", default=FALSE,
              help="Restrict MVR analysis to consider only apparent de novo rate [default %default]"),
  make_option(c("--MVR"), type="numeric", default=0.1,
              help="Maximum Mendelian violation rate to tolerate for non-LOW_QUALITY variants [default %default]")
)

###Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog INFILE OUTDIR",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options


###Checks for appropriate positional arguments
if(length(args$args) != 2){
  stop("Incorrect number of required positional arguments\n")
}

###Writes args & opts to vars
INFILE <- args$args[1]
OUTDIR <- args$args[2]
prefix <- opts$prefix
min.variants <- opts$min.variants
step <- opts$step
target.MVR <- opts$MVR
deNovo.only <- opts$deNovo.only

#Dev parameters
INFILE <- "~/scratch/gnomAD_v2_SV_MASTER.merged_MVR_data.txt.gz"
OUTDIR <- "~/scratch/MVR_test/"
prefix <- "MVR_test"
min.variants <- 50
step <- 5
target.MVR <- 0.05
deNovo.only <- T

#Create output directory, if necessary
if(!dir.exists(OUTDIR)){
  dir.create(OUTDIR)
}

###Read & clean input data
dat <- read.table(INFILE,header=T,comment.char="",check.names=F)
dat <- cleanDat(dat)
dat <- classifyDatByMVR(dat,cutoff=target.MVR,deNovo.only=deNovo.only)

###Get vectors for titration
SVTYPES <- sort(unique(as.character(dat$SVTYPE)))
# FILTERS <- sort(as.character(unique(unlist(strsplit(dat$FILTER,split=",")))))
FILTERS <- unique(as.character(dat$FILTER))
SVLEN.categories <- c("SMALL","MEDIUM","LARGE")

###Titrate across all conditions and calculate MVR
base.MVRs <- as.data.frame(do.call("rbind", lapply(SVTYPES,function(SVTYPE){
  do.call("rbind", lapply(FILTERS,function(FILTER){
    do.call("rbind", lapply(SVLEN.categories,function(SVLEN.category){
      filt.dat <- filterDat(dat,SVTYPE=SVTYPE,FILTER=FILTER,
                            SVLEN.category=SVLEN.category)
      if(nrow(filt.dat)>0){
        MVR <- calcHighMVRFrac(filt.dat)
        res <- cbind("SVTYPE"=SVTYPE,
                     "FILTER"=FILTER,
                     "SVLEN.category"=SVLEN.category,
                     "SITES"=nrow(filt.dat),
                     "unfiltered.MVR"=MVR)
        return(res)
      }
    }))
  }))
})))

# 
# ###Restrict final filter panel to categories with:
# ### 1) MVR > target.MVR; and
# ### 2) SITES >= min.variants
# filter.panel <- MVR.results[which(MVR.results$SITES >= min.variants & 
#                                   MVR.results$MVR > target.MVR),
#                             1:6]
# 
# ###Write outfiles
# write.table(MVR.results,paste(OUTDIR,"/",prefix,".MVR_analysis_results.allCategores"))



# #Test:
# boxplot((1-(dat$APPARENT_DE_NOVO/dat$TRIOS_CONSIDERED))[which(dat$QUAL>=1 & dat$QUAL<200)],
#         (1-(dat$APPARENT_DE_NOVO/dat$TRIOS_CONSIDERED))[which(dat$QUAL>=200 & dat$QUAL<400)],
#         (1-(dat$APPARENT_DE_NOVO/dat$TRIOS_CONSIDERED))[which(dat$QUAL>=400 & dat$QUAL<600)],
#         (1-(dat$APPARENT_DE_NOVO/dat$TRIOS_CONSIDERED))[which(dat$QUAL>=600 & dat$QUAL<800)],
#         (1-(dat$APPARENT_DE_NOVO/dat$TRIOS_CONSIDERED))[which(dat$QUAL>=800 & dat$QUAL<1000)])
# plot(sapply(seq(0,999,5),function(i){
#   mean((dat$APPARENT_DE_NOVO/dat$TRIOS_CONSIDERED)[which(dat$QUAL>=i)],na.rm=T)
# }))
#Mean across variants
# plot(sapply(seq(0,999,5),function(i){
#   mean((1-(dat$MENDELIAN/dat$TRIOS_CONSIDERED))[which(dat$MEDIAN_NONREF_GQ>=i)],na.rm=T)
# }),ylab="")
# #Sum of all transmissions
# plot(sapply(seq(0,999,5),function(i){
#   sum(dat$APPARENT_DE_NOVO[which(dat$MEDIAN_NONREF_GQ>=i)])/sum(dat$TRIOS_CONSIDERED[which(dat$MEDIAN_NONREF_GQ>=i)])
# }),ylab="")
#Fraction of variants with de novo rate >5% (y)
# filt.dat <- filterDat(dat,SVTYPE="DEL",SVLEN.category="SMALL",FILTER="HIGH_SR_BACKGROUND")
# plot(sapply(seq(0,999,5),function(i){
#   calcHighMVRFrac(filt.dat[which(filt.dat$MEDIAN_NONREF_GQ>=i),])
# }),ylab="Mendelian Violation Rate",xlab="Median non-ref GQ (x5)",
# main="Small Deletions, High SR Background")
# plot(sapply(seq(0,999,5),function(i){
#   1-calcHighMVRFrac(dat[which(dat$QUAL>=i),],deNovo.only = T,cutoff = 0.05)
# }),ylab="")
#ROC: Fraction of sites with nonref GQ > i (x) vs Fraction of variants with de novo rate >5% (y)
# plot(sapply(seq(0,99,1),function(i){((nrow(dat)-length(which(dat$MEDIAN_NONREF_GQ>=i)))/nrow(dat))}),
#      sapply(seq(0,99,1),function(i){
#   1-calcHighMVRFrac(dat[which(dat$MEDIAN_NONREF_GQ>=i),],deNovo.only = T,cutoff = 0.01)
# }),ylab="")
# plot(c(mean((1-(dat$MENDELIAN/dat$TRIOS_CONSIDERED))[which(dat$QUAL>=1 & dat$QUAL<200)],na.rm=T),
#         mean((1-(dat$MENDELIAN/dat$TRIOS_CONSIDERED))[which(dat$QUAL>=200 & dat$QUAL<400)],na.rm=T),
#         mean((1-(dat$MENDELIAN/dat$TRIOS_CONSIDERED))[which(dat$QUAL>=400 & dat$QUAL<600)],na.rm=T),
#         mean((1-(dat$MENDELIAN/dat$TRIOS_CONSIDERED))[which(dat$QUAL>=600 & dat$QUAL<800)],na.rm=T),
#         mean((1-(dat$MENDELIAN/dat$TRIOS_CONSIDERED))[which(dat$QUAL>=800 & dat$QUAL<1000)]),na.rm=T),
#      ylab="")

#DEL PEAK ANALYSIS - JAN 8, 2019
#Isolate del peak
dels <- dat[which(dat$SVTYPE=="DEL"),]
dels$NCR <- dels$PCRMINUS_NULL_GTs/(dels$PCRMINUS_NULL_GTs+dels$PCRMINUS_REF_GTs+dels$PCRMINUS_NONREF_GTs)
hist(log10(dels$SVLEN),breaks=100,col="firebrick",
     main="Size (All Deletions)",xlab="log10(DEL size)",
     ylab="Variants (Count)")
abline(v=log10(c(350,1000)),lwd=2,lty=2)
dels.art <- dels[which(dels$SVLEN>350 & dels$SVLEN<1000),]
dels.noart <- dels[which(dels$SVLEN<=350 | dels$SVLEN>=1000),]
#Compare no-call rates
boxplot(100*dels$NCR,100*dels.art$NCR,100*dels.noart$NCR,outline=F,col="firebrick",lty=1,
        main="PCR- No-Call Rate",names=c("All DELs","Artifact\nDELs","Non-Artifact\nDELs"),
        ylab="No-Call Rate (%)")
boxplot(100*dels$NCR,100*dels.art$NCR,100*dels.noart$NCR,outline=T,col="firebrick",lty=1,
        main="PCR- No-Call Rate",names=c("All DELs","Artifact\nDELs","Non-Artifact\nDELs"),
        ylab="No-Call Rate (%)",cex=0.2)
boxplot(100*dels$NCR,100*dels.art$NCR,100*dels.noart$NCR,outline=F,col="firebrick",lty=1,
        main="PCR- No-Call Rate",names=c("All DELs","Artifact\nDELs","Non-Artifact\nDELs"),
        ylab="No-Call Rate (%)",ylim=c(0,5))
#Compare de novo rate by percent within the artifact peak
dnr.cdf <- sapply(seq(0,0.25,0.001),function(max.NCR){
  length(which(dels.art$DNR>0 & dels.art$NCR<=max.NCR))/length(which(dels.art$NCR<=max.NCR))
})
plot(100*seq(0,0.25,0.001),dnr.cdf,pch=19,cex=0.4,main="Fraction of Artifact DELs with >0% De Novo Rate vs No-Call Rate",
     xlab="Max No-Call Rate (Pct of PCR- Samples)",ylab="Fraction of Sites with any De Novos",
     ylim=c(0,max(dnr.cdf)),col="firebrick")
abline(v=1,lty=2,lwd=2)
plot(100*seq(0,0.25,0.001),dnr.cdf,pch=19,cex=0.4,main="Fraction of Artifact DELs with >0% De Novo Rate vs No-Call Rate",
     xlab="Max No-Call Rate (Pct of PCR- Samples)",ylab="Fraction of Sites with any De Novos",
     ylim=c(0,max(dnr.cdf)),col="firebrick",xlim=c(5,10))
dnr.cdf.inverse <- sapply(seq(0,0.25,0.001),function(max.NCR){
  length(which(dels.art$DNR>0 & dels.art$NCR>max.NCR))/length(which(dels.art$NCR>max.NCR))
})
plot(100*seq(0,0.25,0.001),dnr.cdf.inverse,pch=19,cex=0.4,main="% of Excluded Artifact DELs with >1% De Novo Genotypes by No-Call Rate",
     xlab="Max No-Call Rate (Pct of PCR- Samples)",ylab="Fraction of Sites with >1% De Novo Rate")
#Compute number of sites excluded as a function of max NCR
ncr.cdf <- sapply(seq(0,0.25,0.001),function(max.NCR){
  length(which(dels.art$NCR>max.NCR))
})
plot(100*seq(0,0.25,0.001),ncr.cdf,pch=19,cex=0.4,main="Number of Sites Excluded vs. Max NCR",
     xlab="Max No-Call Rate (Pct of PCR- Samples)",ylab="Sites Excluded (Count)",
     ylim=c(0,max(ncr.cdf)),col="firebrick")
abline(v=1,lty=2,lwd=2)
#Compute size distributions as a function of max NCR (applied only to the artifact peak)
ncr.steps <- seq(0,0.1,0.01)
par(mfrow=c(3,4))
sapply(ncr.steps,function(max.NCR){
  sizes <- c(dels.noart$SVLEN,dels.art$SVLEN[which(dels.art$NCR<=max.NCR)])
  sizes <- sizes[which(sizes<5000)]
  par(mar=c(2,2,2,0.5))
  hist(sizes,breaks=100,col="firebrick",
       main=max.NCR,xlab="",ylab="",ylim=c(0,10000),
       xlim=c(0,5000))
  abline(v=c(350,1000),lwd=2,lty=2)
})
par(mar=c(2,2,2,0.5))
hist(dels$SVLEN[which(dels$SVLEN<5000)],breaks=100,col="firebrick",
     main="Unfiltered",xlab="",ylab="",ylim=c(0,10000),
     xlim=c(0,5000))
abline(v=c(350,1000),lwd=2,lty=2)
#Compute rough estimate of # of sites excluded per sample vs various max NCRs
excl.persample <- sapply(seq(0,0.25,0.001),function(max.NCR){
  excl.idx <- which(dels.art$NCR>max.NCR)
  sum(dels.art$PCRMINUS_NONREF_GTs[excl.idx]/(dels.art$PCRMINUS_REF_GTs[excl.idx]+dels.art$PCRMINUS_NONREF_GTs[excl.idx]))
})
plot(100*seq(0,0.25,0.001),excl.persample,pch=19,cex=0.4,main="DEL Excluded per Sample vs. Max NCR",
     xlab="Max No-Call Rate (Pct of PCR- Samples)",ylab="Sites Excluded per Sample (Count)",
     ylim=c(0,max(excl.persample)),col="firebrick")
abline(v=c(1,2),lty=2,lwd=2)



