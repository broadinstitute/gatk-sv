#!/usr/bin/env Rscript

# Parse KING output to assign relatedness


###Set master parameters & load libraries
options(stringsAsFactors=F,scipen=1000)
require(e1071)


###################
###HELPER FUNCTIONS
###################
#Import fam file and split into relative classes
#Note: assumes parents of trios are not direct relatives
import.famfile <- function(TRUTH_TRIOS.in, supp.unrelateds.in=NULL){
  fam <- read.table(TRUTH_TRIOS.in,sep="\t")[,2:4]
  colnames(fam) <- c("pro","fa","mo")
  # fam <- fam[grep("SSC_",fam$pro,fixed=T),]
  #Get parent-child pairs
  parent_child <- rbind(data.frame("ID1"=fam$pro,"ID2"=fam$fa),
                        data.frame("ID1"=fam$pro,"ID2"=fam$mo),
                        data.frame("ID1"=fam$fa,"ID2"=fam$pro),
                        data.frame("ID1"=fam$mo,"ID2"=fam$pro))
  parent_child$collapsed_IDs <- apply(parent_child,1,paste,collapse="_")
  parent_child <- parent_child[which(!(duplicated(parent_child$collapsed_IDs))),]
  #Get unrelated pairs
  parent_parent <- rbind(data.frame("ID1"=fam$fa,"ID2"=fam$mo),
                         data.frame("ID1"=fam$mo,"ID2"=fam$fa))
  if(!is.null(supp.unrelateds.in)){
    supp.unrelateds <- read.table(supp.unrelateds.in, header=F, sep="\t")
    parent_parent <- rbind(parent_parent,
                           data.frame("ID1"=supp.unrelateds[, 1], "ID2"=supp.unrelateds[, 2]),
                           data.frame("ID1"=supp.unrelateds[, 2], "ID2"=supp.unrelateds[, 1]))
  }
  parent_parent$collapsed_IDs <- apply(parent_parent,1,paste,collapse="_")
  parent_parent <- parent_parent[which(!(duplicated(parent_parent$collapsed_IDs))),]
  #Get sibling pairs
  #Only use SSC siblings, and assume SSC nomenclature
  fam.ssc <- fam[grep("SSC_",fam$pro,fixed=T), ]
  sibling_sibling <- do.call("rbind", lapply(fam.ssc$pro,function(ID){
    base.famID <- unlist(strsplit(ID,split="_"))[2]
    pro.IDs <- paste("SSC",base.famID,"P1",sep="_")
    sib.IDs <- paste("SSC",base.famID,"S1",sep="_")
    if(all(sib.IDs %in% fam.ssc$pro)){
      return(rbind(data.frame("ID1"=pro.IDs,"ID2"=sib.IDs),
                   data.frame("ID1"=sib.IDs,"ID2"=pro.IDs)))
    }
  }))
  sibling_sibling$collapsed_IDs <- apply(sibling_sibling,1,paste,collapse="_")
  sibling_sibling <- sibling_sibling[which(!(duplicated(sibling_sibling$collapsed_IDs))),]
  return(list("parent_child"=parent_child,
              "parent_parent"=parent_parent,
              "sibling_sibling"=sibling_sibling))
}
#Import KING data and label based on truth trios
import.king <- function(KING.in,fam.labels){
  king <- read.table(KING.in,header=T,sep="\t",comment.char="")
  king$collapsed_IDs <- apply(king[,1:2],1,paste,collapse="_")
  king$label.truth <- NA
  king$label.truth[which(king$collapsed_IDs %in% fam.labels$parent_child$collapsed_IDs)] <- "parent_child"
  king$label.truth[which(king$collapsed_IDs %in% fam.labels$parent_parent$collapsed_IDs)] <- "unrelated"
  king$label.truth[which(king$collapsed_IDs %in% fam.labels$sibling_sibling$collapsed_IDs)] <- "siblings"
  return(king)
}
#Plot single KING metric based on truth labels
plot.king.metric <- function(king.train, metric){
  vals <- king.train[, which(colnames(king.train) == metric)]
  pc <- vals[which(king.train$label.truth=="parent_child")]
  sib <- vals[which(king.train$label.truth=="siblings")]
  unrel <- vals[which(king.train$label.truth=="unrelated")]
  par(mar=c(4,7,3,1))
  boxplot(pc, sib, unrel, 
          col=c("dodgerblue", "forestgreen", "gray50"),
          horizontal=T, las=2, outline=F, lty=1, staple.wex=0,
          names=c("Parent-Child", "Siblings", "Unrelated"),
          main=metric)
}
#Train & apply SVM to classify relationships
classify.relationships <- function(king,king.train){
  #Train SVM on pairs with truth labels
  svm.col.names <- c("HetHet","IBS0","HetConc","HomIBS0","Kinship","IBD1Seg","IBD2Seg","PropIBD")
  svm.col.idxs <- which(colnames(king) %in% svm.col.names)
  train.mat <- king.train[,svm.col.idxs]
  train.labels <- factor(king.train$label.truth)
  train.mat <- apply(train.mat,2,as.numeric)
  train.dat <- data.frame("relation"=train.labels,train.mat)
  set.seed(123456789)
  svm.fit <- svm(relation ~ ., data=train.dat,
                 type="C-classification",
                 kernel="linear",
                 cross=10,
                 probability=T,
                 na.action=na.omit)
  
  #Apply SVM to classify all samples
  classify.mat <- apply(king[,svm.col.idxs],2,as.numeric)
  predictions <- predict(svm.fit, classify.mat)
  predictions.out <- data.frame("ID1"=king$ID1,"ID2"=king$ID2,
                                "label.truth"=king$label.truth,
                                "label.inferred"=as.character(predictions),
                                "collapsed_IDs"=king$collapsed_IDs)
  return(predictions.out)
}
#Check accuracy of inferred relationship labels
check.label.accuracy <- function(relationships){
  acc.test.df <- relationships[which(!is.na(relationships$label.truth)),]
  acc.test.res <- as.data.frame(t(sapply(c("unrelated","all_related","parent_child","siblings"),function(label){
    if(label=="all_related"){
      related_labels <- c("parent_child", "siblings")
      sens <- length(which(acc.test.df$label.truth %in% related_labels & acc.test.df$label.inferred %in% related_labels))/length(which(acc.test.df$label.truth %in% related_labels))
      fdr <- length(which(!(acc.test.df$label.truth %in% related_labels) & acc.test.df$label.inferred %in% related_labels))/length(which(!(acc.test.df$label.truth %in% related_labels)))
      n.known.raw <- length(which(acc.test.df$label.truth %in% related_labels))
      n.known <- length(which(acc.test.df$label.truth %in% related_labels & acc.test.df$label.inferred %in% related_labels))
      n.new <- length(which(relationships$label.inferred %in% related_labels))-n.known
    }else{
      sens <- length(which(acc.test.df$label.truth==label & acc.test.df$label.inferred==label))/length(which(acc.test.df$label.truth==label))
      fdr <- length(which(acc.test.df$label.truth!=label & acc.test.df$label.inferred==label))/length(which(acc.test.df$label.truth!=label))
      n.known.raw <- length(which(acc.test.df$label.truth==label))
      n.known <- length(which(acc.test.df$label.truth==label & acc.test.df$label.inferred==label))
      n.new <- length(which(relationships$label.inferred==label))-n.known
    }
    return(c(sens,fdr,n.known.raw,n.known,n.new))
  })))
  colnames(acc.test.res) <- c("Sensitivity","FDR",
                              "N.known.possible","N.known.retained",
                              "N.new.predicted")
  acc.test.res <- cbind("Relationship"=rownames(acc.test.res),
                        acc.test.res)
  return(acc.test.res)
}
#Optimize list of samples to prune
optimize.prune.list <- function(relationships){
  #Get list of all related samples
  pairs <- relationships[which(relationships$label.inferred!="unrelated"),1:2]
  samples <- unique(as.character(c(pairs[,1],pairs[,2])))
  #Optimization loop
  samples.to.prune <- c()
  n.remaining.pairs <- length(which(complete.cases(pairs)))
  while(n.remaining.pairs>0){
    #Step 1: find sample with most remaining complete relationships
    pairs.per.sample <- table(as.character(c(pairs[which(complete.cases(pairs)),1],
                                             pairs[which(complete.cases(pairs)),2])),
                              useNA="no")
    next.sample.ID <- head(names(which(pairs.per.sample==max(pairs.per.sample))),1)
    #Step 2: prune that sample from all pairs
    pairs$ID1[which(pairs$ID1==next.sample.ID)] <- NA
    pairs$ID2[which(pairs$ID2==next.sample.ID)] <- NA
    samples.to.prune <- c(samples.to.prune,next.sample.ID)
    #Step 3: update number of pairs without at least one NA
    n.remaining.pairs <- length(which(complete.cases(pairs)))
  }
  return(samples.to.prune)
}


########################
###RSCRIPT FUNCTIONALITY
########################
require(optparse)

###List of command-line options
option_list <- list(
  make_option(c("-p", "--outprefix"), type="character", default="./KING_processed",
              help="output file prefix to be appended to all outputs [default %default]",
              metavar="character"),
  make_option(c("-u", "--unrelated"), type="character", default=NULL,
              help="supplementary pairs of samples to assume as unrelated [default %default]",
              metavar="character")
)

###Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog KING TRUTH_TRIOS",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

###Checks for appropriate positional arguments
if(length(args$args) != 2){
  stop("Incorrect number of required positional arguments\n")
}

###Writes args & opts to vars
KING.in <- args$args[1]
TRUTH_TRIOS.in <- args$args[2]
OUTPREFIX <- opts$outprefix
supp.unrelateds.in <- opts$unrelated

#Dev inputs (local)
# KING.in <- "~/scratch/gnomAD_v2_SV_MASTER.king_metrics.candidate_pairs_subsetted.txt.gz"
# KING.in <- "~/scratch/gnomAD_v2_SV_MASTER.king_metrics.train_truth_pairs.txt.gz"
# TRUTH_TRIOS.in <- "~/scratch/gnomAD_v2.trios.prepipeline_QC_pass.fam"
# OUTPREFIX <- "~/scratch/gnomAD_v2_SV_MASTER.KING_processed_test"
# OUTPREFIX <- "~/scratch/gnomAD_v2_SV_MASTER.KING_truth_test"
# supp.unrelateds.in <- "~/scratch/unrelated_sample_pairs.supplement.txt"


###############
###PROCESS DATA
###############
###Import data
fam.labels <- import.famfile(TRUTH_TRIOS.in, supp.unrelateds.in)
king <- import.king(KING.in,fam.labels)
king.train <- king[which(!is.na(king$label.truth)),]

###Summarize KING training data identified
cat(paste("\nKING training data contains", 
          prettyNum(length(which(king.train$label.truth=="parent_child")), big.mark=","), 
          "of", prettyNum(nrow(fam.labels$parent_child)/2, big.mark=","), 
          "known parent-child training pairs\n", sep=" "))
cat(paste("\nKING training data contains", 
          prettyNum(length(which(king.train$label.truth=="siblings")), big.mark=","), 
          "of", prettyNum(nrow(fam.labels$sibling_sibling)/2, big.mark=","), 
          "known sibling training pairs\n", sep=" "))
cat(paste("\nKING training data contains", 
          prettyNum(length(which(king.train$label.truth=="unrelated")), big.mark=","), 
          "of", prettyNum(nrow(fam.labels$parent_parent)/2, big.mark=","), 
          "known unrelated training pairs\n", sep=" "))

###Plot characteristics of training data
sapply(c("HetHet","IBS0","HetConc","HomIBS0","Kinship","IBD1Seg","IBD2Seg","PropIBD"), function(metric){
  pdf(paste(OUTPREFIX, metric, "training_metrics.pdf", sep="."), height=3, width=5)
  plot.king.metric(king.train, metric)
  dev.off()
})

###Assign relatedness labels with SVM
relationships <- classify.relationships(king, king.train)
#Check & report accuracy vs known training labels
acc.table <- check.label.accuracy(relationships)
write.table(acc.table,
            paste(OUTPREFIX,".relationship_inference_accuracy.txt",sep=""),
            col.names=T,row.names=F,sep="\t",quote=F)
#Write non-unrelated relationships to file
all.relationships.out <- relationships[which(relationships$label.inferred!="unrelated"),-c(3,5)]
colnames(all.relationships.out)[1] <- "#ID1"
write.table(all.relationships.out,
            paste(OUTPREFIX,".inferred_pairwise_relationships.txt",sep=""),
            col.names=T,row.names=F,sep="\t",quote=F)
system(paste("gzip -f ",OUTPREFIX,".inferred_pairwise_relationships.txt",sep=""),
       wait=T,intern=F)

###Compile list of samples to prune
samples.to.prune <- optimize.prune.list(relationships)
write.table(samples.to.prune,
            paste(OUTPREFIX,".related_samles_to_prune.txt",sep=""),
            col.names=F,row.names=F,sep="\t",quote=F)


