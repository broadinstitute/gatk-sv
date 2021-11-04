#!/usr/bin/env Rscript

# Generates a UMAP plot based on a GRM


###Set master parameters & load libraries
options(stringsAsFactors=F,scipen=1000)
require(umap)

#Inputs
grm.input <- "~/scratch/grm.txt.gz"
anc.input <- "~/Downloads/gnomAD_v2_SV.sample_population_assignments.rough_initial_simplified.txt"
seed <- 123456789

#Read input data
grm <- read.table(grm.input,header=T,check.names=F,comment.char="")
grm <- as.data.frame(t(grm[,-1]))
pop.train <- read.table(anc.input,header=F)
colnames(pop.train) <- c("ID","pop")
pop.train <- pop.train[which(pop.train$pop != "."),]
pops <- unique(pop.train$pop)
pop.cols <- c("firebrick2","darkorchid3","gray30","dodgerblue1","forestgreen","gray70")
pops.to.exclude <- c("OTH", "ASJ")

#UMAP
set.seed(seed)
umap.res <- umap(grm)
samp.cols <- sapply(rownames(grm),function(ID){
  pop <- pop.train[which(pop.train$ID==ID),]$pop
  if(length(pop)>0){
    return(as.character(pop.cols[which(pops==pop)]))
  }else{
    return("gray70")
  }
})
plot(umap.res$layout,col=samp.cols)
