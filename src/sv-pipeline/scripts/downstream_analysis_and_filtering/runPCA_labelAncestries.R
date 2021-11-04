#!/usr/bin/env Rscript

# Performs PCA based on an input VCF


###Set master parameters & load libraries
options(stringsAsFactors=F,scipen=1000)
require(e1071)


########################
###RSCRIPT FUNCTIONALITY
########################
require(optparse)

###List of command-line options
option_list <- list(
  make_option(c("-p", "--outprefix"), type="character", default="./SV_PCA",
              help="output file prefix to be appended to all outputs [default %default]",
              metavar="character"),
  make_option(c("--batchAssignments"), type="character", default=NULL,
              help="file listing batch assignments per sample [default %default]",
              metavar="character"),
  make_option(c("--PCRPLUSsamples"), type="character", default=NULL,
              help="list of PCR+ samples [default %default]",
              metavar="character"),
  make_option(c("--confidence"), type="numeric", default=0.8,
              help="minimum probability for a sample to be assigned to a population [default %default]",
              metavar="numeric")
)

###Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog GRM TRUTH_POP_ASSIGNMENTS POP_COLORS",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

###Checks for appropriate positional arguments
if(length(args$args) != 3){
  stop("Incorrect number of required positional arguments\n")
}

###Writes args & opts to vars
grm.in <- args$args[1]
pop.train.in <- args$args[2]
pop.colors.in <- args$args[3]
OUTPREFIX <- opts$outprefix
batch.assignments.in <- opts$batchAssignments
PCRPLUS.samples.in <- opts$PCRPLUSsamples
confidence <- opts$confidence

# #Dev inputs (local)
# # vcf.input <- "~/scratch/out.vcf.gz"
# # pop.train.input <- "~/Downloads/gnomAD_v2_SV.sample_population_assignments.rough_initial_simplified.txt"
# grm.in <- "~/scratch/gnomAD_v2_SV_MASTER.grm.txt.gz"
# pop.train.in <- "~/scratch/gnomAD_v2_SV.sample_population_assignments.cleaned_truth_subset_for_clustering.txt"
# pop.colors.in <- "~/Desktop/Collins/Talkowski/NGS/SV_Projects/gnomAD/gnomad_sv/FireCloud_formal_analysis/refs/gnomAD_population_colors.txt"
# OUTPREFIX <- "~/scratch/gnomAD_v2_SV_MASTER.PCA_test"
# batch.assignments.in <- "~/scratch/gnomAD_v2_SV.sample_batch_assignments.txt"
# PCRPLUS.samples.in <- "~/scratch/gnomAD_v2_SV_PCRPLUS.samples.list"
# confidence <- 0.8



###############
###PROCESS DATA
###############
#Read & clean GRM
grm <- as.data.frame(t(read.table(grm.in,header=T,comment.char="",sep="\t")[,-1]))
samples <- rownames(grm)
grm <- as.data.frame(apply(grm,2,function(values){
  values <- as.numeric(values)
  values[which(is.na(values))] <- mean(values,na.rm=T)
  return(values)
}))
rownames(grm) <- samples

#Read other input data
pop.train <- read.table(pop.train.in,header=F)
colnames(pop.train) <- c("ID","pop")
pop.train <- pop.train[which(pop.train$ID %in% rownames(grm)),]
pop.train$pop[which(pop.train$pop=="SAS")] <- "."
pop.train <- pop.train[which(pop.train$pop != "."),]
pop.train$pop[which(pop.train$pop %in% c("ASJ","FIN","NFE"))] <- "EUR"
pops <- unique(pop.train$pop)
pop.cols.df <- read.table(pop.colors.in,header=T,comment.char="",sep="\t")
pop.cols <- pop.cols.df$color
names(pop.cols) <- pop.cols.df$pop

# #Convert VCF to SNP GDS format
# snpgdsVCF2GDS(vcf.fn=vcf.input,
#               out.fn="input.gds",
#               method="copy.num.of.ref")
# genofile <- snpgdsOpen("input.gds")

#Run PCA
pca <- prcomp(grm)

# #Compute & plot variance explained
# eigs <- pca$sdev^2
# pct.explained <- eigs/sum(eigs)
# barplot(pct.explained[1:10])

#Format PC matrix & write out to file
pca.mat.all <- pca$x
colnames(pca.mat.all) <- paste("PC",1:ncol(pca.mat.all),sep="")
rownames(pca.mat.all) <- samples
pca.mat.out <- data.frame(samples,pca.mat.all[,1:100])
colnames(pca.mat.out)[1] <- "#sample"
write.table(pca.mat.out,paste(OUTPREFIX,".PCA_loadings.txt",sep=""),
            col.names=T,row.names=F,sep="\t",quote=F)
system(paste("gzip -f ",OUTPREFIX,".PCA_loadings.txt",sep=""),wait=T,intern=F)

#Subset to top 4 PCs for population clustering
pca.mat <- pca.mat.all[,1:4]

#Train SVM on samples with truth labels
train.mat <- pca.mat[which(rownames(pca.mat) %in% pop.train$ID),]
train.labels <- factor(sapply(rownames(train.mat),function(ID){as.character(pop.train$pop[which(pop.train$ID==ID)])}))
train.mat <- apply(train.mat,2,as.numeric)
train.dat <- data.frame("pop"=train.labels,train.mat)
set.seed(123456789)
svm.fit <- svm(pop ~ ., data=train.dat,
               type="C-classification",
               kernel="radial",
               cross=10,
               probability=T,
               na.action=na.omit)

#Apply SVM to classify all samples
classify.mat <- apply(pca.mat,2,as.numeric)
predictions <- predict(svm.fit, classify.mat, probability=T)
predicted.labels <- as.character(predictions)
names(predicted.labels) <- rownames(pca.mat)
prediction.pvals <- attr(predictions, "probabilities")
max.pred.pvals <- apply(prediction.pvals,1,max)
predicted.labels[which(max.pred.pvals<confidence)] <- "OTH"

#Write predictions to file
labels.out <- data.frame("sample"=names(predicted.labels),
                         "inferred_pop_final"=predicted.labels,
                         "pvalue"=max.pred.pvals,
                         "inferred_pop_raw"=as.character(predictions),
                         "groundTruth_pop"=sapply(names(predicted.labels),function(ID){
                           ID <- as.character(ID)
                           if(ID %in% pop.train$ID){
                             return(pop.train$pop[which(pop.train$ID==ID)])
                           }else{
                             return(NA)
                           }
                         }))
colnames(labels.out)[1] <- "#sample"
write.table(labels.out,paste(OUTPREFIX,".new_population_labels.txt",sep=""),
            col.names=T,row.names=F,sep="\t",quote=F)
system(paste("gzip -f ",OUTPREFIX,".new_population_labels.txt",sep=""),wait=T,intern=F)

# #Evaluate accuracy
# acc.table <- labels.out[which(!is.na(labels.out$groundTruth_pop) &
#                    labels.out$inferred_pop_final!="OTH"),
#            c(2,5)]
# accuracy <- length(which(acc.table[,1]==acc.table[,2]))/nrow(acc.table)
# fn.table <- labels.out[which(!is.na(labels.out$groundTruth_pop) &
#                                labels.out$inferred_pop_final=="OTH"),
#                        c(2,5)]


#Plot data with training labels
pca.cols.train <- unlist(sapply(rownames(pca.mat),function(ID){
  if(ID %in% pop.train$ID){
    pop <- pop.train[which(pop.train$ID==ID),]$pop
    return(as.character(pop.cols[which(names(pop.cols)==pop)]))
  }else{
    return(NA)
  }
}))
pca.cols.train.inverted <- pca.cols.train
pca.cols.train.inverted[which(is.na(pca.cols.train))] <- "gray30"
pca.cols.train.inverted[which(!is.na(pca.cols.train))] <- NA
png(paste(OUTPREFIX,".PCA.training_set.png",sep=""),height=2100,width=2100,res=300)
pairs(pca.mat,col=pca.cols.train,cex=0.3)
dev.off()
png(paste(OUTPREFIX,".PCA.unlabeled_samples.png",sep=""),height=2100,width=2100,res=300)
pairs(pca.mat,col=pca.cols.train.inverted,cex=0.3)
dev.off()


#Plot data with inferred labels
pca.cols.inferred <- unlist(sapply(rownames(pca.mat),function(ID){
  pop <- predicted.labels[which(names(predicted.labels)==as.character(ID))]
  return(as.character(pop.cols[which(names(pop.cols)==pop)]))
}))
png(paste(OUTPREFIX,".PCA.labeled_samples.png",sep=""),height=2100,width=2100,res=300)
pairs(pca.mat,col=pca.cols.inferred,cex=0.3)
dev.off()

#Plot data with PCR labels, if optioned
if(!is.null(PCRPLUS.samples.in)){
  plus.samples <- read.table(PCRPLUS.samples.in,header=F,sep="\t")[,1]
  pca.mat.pcr <- rbind(pca.mat[which(!(rownames(pca.mat) %in% plus.samples)),],
                       pca.mat[which(rownames(pca.mat) %in% plus.samples),])
  # pca.mat.pcr <- pca.mat[sample(1:nrow(pca.mat),size=nrow(pca.mat),replace=F),]
  pcr.cols <- unlist(sapply(rownames(pca.mat.pcr),function(ID){
    if(ID %in% plus.samples){
      return("#1F6AAB")
    }else{
      return("#E01925")
    }
  }))
  png(paste(OUTPREFIX,".PCA.PCR_status.png",sep=""),height=2100,width=2100,res=300)
  pairs(pca.mat.pcr,col=pcr.cols,cex=0.3)
  dev.off()
}

#Plot data with batch labels, if optioned
if(!is.null(PCRPLUS.samples.in)){
  batch.assignments <- read.table(batch.assignments.in,header=F,sep="\t")
  batches <- unique(batch.assignments[,2])
  batch.cols <- rainbow(length(unique(batch.assignments[,2])))
  pca.mat.batch <- pca.mat[sample(1:nrow(pca.mat),size=nrow(pca.mat),replace=F),]
  batch.cols <- unlist(sapply(rownames(pca.mat.batch),function(ID){
    batch <- batch.assignments[which(batch.assignments[,1]==ID),2]
    return(batch.cols[which(batches==batch)])
  }))
  png(paste(OUTPREFIX,".PCA.batch_assignment.png",sep=""),height=2100,width=2100,res=300)
  pairs(pca.mat.pcr,col=batch.cols,cex=0.3)
  dev.off()
}

