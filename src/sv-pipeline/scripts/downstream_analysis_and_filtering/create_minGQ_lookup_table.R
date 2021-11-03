#!/usr/bin/env Rscript

# Helper script to create final minGQ tree from
# combined outputs of ROC analysis at individual
# filter conditions


###Set master parameters
options(stringsAsFactors=F,scipen=1000)


###################
###HELPER FUNCTIONS
###################
#Read and merge all input data
process.inputs <- function(CONDITIONS,STATS,ROC_OPTS){
  #Read data
  cond <- read.table(CONDITIONS,comment.char="",header=T)
  colnames(cond)[1] <- "condition"
  stats <- read.table(STATS,comment.char="",header=F)
  colnames(stats) <- c("condition","perchild.median","perchild.q1","perchild.q3")
  opts <- read.table(ROC_OPTS,comment.char="",header=F)
  opts <- opts[which(opts[,1]!="condition"),]
  colnames(opts) <- c("condition","minGQ","inh.ret","inh.ret.rate","dn.ret","dn.rate")
  
  #Merge data
  merged <- merge(cond,stats,by="condition",all.x=T,all.y=F,sort=F)
  merged <- merge(merged,opts,by="condition",all.x=T,all.y=F,sort=F)
  merged <- merged[order(as.numeric(sapply(strsplit(as.character(merged$condition),split="_"),tail,1))),]
  rownames(merged) <- 1:nrow(merged)
  
  #Clean data
  numeric.cols.idx <- which(colnames(merged) %in% c("minSVLEN","maxSVLEN","minAF","maxAF",
                                                    "perchild.median","perchild.q1","perchild.q3",
                                                    "minGQ","inh.ret","inh.ret.rate","dn.ret","dn.rate"))
  merged[,numeric.cols.idx] <- apply(merged[,numeric.cols.idx],2,as.numeric)
  string.cols.idx <- setdiff(1:ncol(merged),numeric.cols.idx)
  merged[,string.cols.idx] <- apply(merged[,string.cols.idx],2,as.character)
  merged[which(merged$minSVLEN<=0),]$minSVLEN <- -2
  
  #Determine which filters were active for each condition
  merged$filt.SVLEN <- TRUE
  merged$filt.SVLEN[which(merged$maxSVLEN>300000000 & merged$minSVLEN<=1)] <- FALSE
  merged$filt.AF <- TRUE
  merged$filt.AF[which(merged$maxAF>=0.99 & merged$minAF<=0.001)] <- FALSE
  merged$filt.SVTYPE <- TRUE
  merged$filt.SVTYPE[which(merged$excludeSVTYPE=="None" & sapply(strsplit(merged$includeSVTYPE,split=","),length)>3)] <- FALSE
  merged$filt.FILTER <- TRUE
  merged$filt.FILTER[which(merged$includeFILTER=="None" & merged$excludeFILTER=="None")] <- FALSE
  merged$filt.EV <- TRUE
  merged$filt.EV[which(merged$includeEV=="None" & merged$excludeEV=="None")] <- FALSE
  merged$n.filt <- apply(merged[,which(colnames(merged) %in% paste("filt.",c("SVLEN","AF","SVTYPE","FILTER","EV"),sep=""))],
                         1,function(filts){
                           return(length(which(filts)))
                         })
  
  #Add filter category for each filter class
  merged$c.SVLEN <- paste(merged$minSVLEN,merged$maxSVLEN,sep="_")
  merged$c.AF <- paste(merged$minAF,merged$maxAF,sep="_")
  merged$c.SVTYPE <- paste(merged$includeSVTYPE,merged$excludeSVTYPE,sep="_")
  merged$c.FILTER <- paste(merged$includeFILTER,merged$excludeFILTER,sep="_")
  merged$c.EV <- paste(merged$includeEV,merged$excludeEV,sep="_")
  
  #Return processed data
  return(merged)
}
#Determine which additional filter maximizes inh.ret.rate
select.next.filter <- function(dat,filt.current=NULL){
  #Clean current filters and eligible filters
  all.filts <- paste("filt.",c("SVLEN","AF","SVTYPE","FILTER","EV"),sep="")
  n.filt.current <- length(filt.current)
  filt.unused <- all.filts[which(!(all.filts %in% filt.current))]
  n.filt.target <- n.filt.current + 1
  
  #Subset data
  dat.sub <- dat[which(dat$n.filt==n.filt.target),]
  dat.sub <- dat.sub[which(apply(as.data.frame(dat.sub[,which(colnames(dat.sub) %in% filt.current)]),1,all)),]
  
  #Iterate over eligible filters and compute total number of inherited variants retained for each
  m <- sapply(1:length(filt.unused),function(i){
    #Get filter indexes
    f.incl.idx <- which(colnames(dat.sub) %in% filt.unused[i])
    f.excl.idx <- which(colnames(dat.sub) %in% filt.unused[-i])
    
    #Get data subset
    f.sub <- dat.sub[which(dat.sub[,f.incl.idx] & !(apply(as.data.frame(dat.sub[,f.excl.idx]),1,all))),]
    
    #Calculate inh.ret sum
    m.idx <- which(!is.na(f.sub$inh.ret))
    if(length(m.idx)>0){
      m <- sum(f.sub$inh.ret[m.idx])
    }else{
      m <- NA
    }
    return(m)
  })
  
  #Select & return best filter
  #If no filter data are available, a filter is randomly chosen
  if(any(!is.na(m))){
    best <- max(m,na.rm=T)
    best.f <- filt.unused[which(m==max(m,na.rm=T))]
    #If two or more filters are equally good, a filter is randomly chosen
    if(length(best.f)>1){
      set.seed(123456789)
      best.f <- sample(best.f,size=1)
    }
  }else{
    set.seed(123456789)
    best.f <- sample(filt.unused,size=1)
  }
  return(best.f)
}
#Determine optimal ordering of filtering tree
build.tree <- function(dat,fixed.f1=NULL){
  #Create tree structure
  tree <- data.frame("filt.1"=NA,
                     "filt.1.category"=NA,
                     "filt.2"=NA,
                     "filt.2.category"=NA,
                     "filt.3"=NA,
                     "filt.3.category"=NA,
                     "filt.4"=NA,
                     "filt.4.category"=NA,
                     "filt.5"=NA,
                     "filt.5.category"=NA)
  tree <- tree[-1,]
  
  #Layer 1
  if(is.null(fixed.f1)){
    f1 <- select.next.filter(dat,filt.current=NULL)
  }else{
    f1 <- fixed.f1
  }
  f1.suf <- tail(unlist(strsplit(f1,split=".",fixed=T)),1)
  f1.f.idx <- which(colnames(dat) %in% f1)
  f1.c.idx <- which(colnames(dat) %in% paste("c",f1.suf,sep="."))
  
  #Layer 2
  for(f1.c in unique(dat[which(dat[,f1.f.idx]),f1.c.idx])){
    dat.f1 <- dat[which(dat[,f1.c.idx] == f1.c),]
    f2 <- select.next.filter(dat.f1,filt.current=f1)
    f2.suf <- tail(unlist(strsplit(f2,split=".",fixed=T)),1)
    f2.f.idx <- which(colnames(dat.f1) %in% f2)
    f2.c.idx <- which(colnames(dat.f1) %in% paste("c",f2.suf,sep="."))
    
    #Layer 3
    for(f2.c in unique(dat.f1[which(dat.f1[,f2.f.idx]),f2.c.idx])){
      dat.f2 <- dat.f1[which(dat.f1[,f2.c.idx] == f2.c),]
      f3 <- select.next.filter(dat.f2,filt.current=c(f1,f2))
      f3.suf <- tail(unlist(strsplit(f3,split=".",fixed=T)),1)
      f3.f.idx <- which(colnames(dat.f2) %in% f3)
      f3.c.idx <- which(colnames(dat.f2) %in% paste("c",f3.suf,sep="."))
      
      #Layer 4
      for(f3.c in unique(dat.f2[which(dat.f2[,f3.f.idx]),f3.c.idx])){
        dat.f3 <- dat.f2[which(dat.f2[,f3.c.idx] == f3.c),]
        f4 <- select.next.filter(dat.f3,filt.current=c(f1,f2,f3))
        f4.suf <- tail(unlist(strsplit(f4,split=".",fixed=T)),1)
        f4.f.idx <- which(colnames(dat.f3) %in% f4)
        f4.c.idx <- which(colnames(dat.f3) %in% paste("c",f4.suf,sep="."))
        
        #Layer 5
        for(f4.c in unique(dat.f3[which(dat.f3[,f4.f.idx]),f4.c.idx])){
          dat.f4 <- dat.f3[which(dat.f3[,f4.c.idx] == f4.c),]
          all.filts <- paste("filt.",c("SVLEN","AF","SVTYPE","FILTER","EV"),sep="")
          f5 <- all.filts[which(!(all.filts %in% c(f1,f2,f3,f4)))]
          f5.suf <- tail(unlist(strsplit(f5,split=".",fixed=T)),1)
          f5.f.idx <- which(colnames(dat.f4) %in% f5)
          f5.c.idx <- which(colnames(dat.f4) %in% paste("c",f5.suf,sep="."))
          
          #Add the entire branch of filter tree to the tree table
          for(f5.c in unique(dat.f4[which(dat.f4[,f5.f.idx]),f5.c.idx])){
            tree <- rbind(tree,c(f1,f1.c,f2,f2.c,f3,f3.c,f4,f4.c,f5,f5.c))
          }
        }
      }
    }
  }
  
  #Return tree
  colnames(tree) <- c("filt.1","filt.1.category",
                      "filt.2","filt.2.category",
                      "filt.3","filt.3.category",
                      "filt.4","filt.4.category",
                      "filt.5","filt.5.category")
  return(tree)
}
#Prune & build final tree
prune.tree <- function(dat,order.tree){
  #Simplify table
  res <- dat[which(dat$n.filt==max(dat$n.filt)),
             which(colnames(dat) %in% c("condition","minSVLEN","maxSVLEN","minAF","maxAF",
                                         "includeSVTYPE","excludeSVTYPE","includeFILTER","excludeFILTER",
                                         "includeEV","excludeEV","minGQ"))]
  res$source <- NA
  #Update minGQ
  for(i in 1:nrow(res)){
    #If minGQ is already set, leave as is & update source
    if(!is.na(res$minGQ[i])){
      res$source[i] <- "ROC"
      
      #If minGQ is not already set, determine closest parent node that is set
    }else{
      #Find matching data row
      r.idx <- which(dat$condition==res$condition[i])
      
      #Get filters for this category
      filts <- as.character(dat[r.idx,grep("c.",colnames(dat),fixed=T)])
      
      #Get branch to follow
      branch.order <- order.tree[which(order.tree$filt.1.category %in% filts & 
                                         order.tree$filt.2.category %in% filts & 
                                         order.tree$filt.3.category %in% filts & 
                                         order.tree$filt.4.category %in% filts & 
                                         order.tree$filt.5.category %in% filts),]
      f1.idx <- which(colnames(dat) == branch.order$filt.1)
      c1.idx <- grep(paste("c",unlist(strsplit(branch.order$filt.1,split=".",fixed=T))[2],sep="."),colnames(dat),fixed=T)
      f2.idx <- which(colnames(dat) == branch.order$filt.2)
      c2.idx <- grep(paste("c",unlist(strsplit(branch.order$filt.2,split=".",fixed=T))[2],sep="."),colnames(dat),fixed=T)
      f3.idx <- which(colnames(dat) == branch.order$filt.3)
      c3.idx <- grep(paste("c",unlist(strsplit(branch.order$filt.3,split=".",fixed=T))[2],sep="."),colnames(dat),fixed=T)
      f4.idx <- which(colnames(dat) == branch.order$filt.4)
      c4.idx <- grep(paste("c",unlist(strsplit(branch.order$filt.4,split=".",fixed=T))[2],sep="."),colnames(dat),fixed=T)
      f5.idx <- which(colnames(dat) == branch.order$filt.5)
      c5.idx <- grep(paste("c",unlist(strsplit(branch.order$filt.5,split=".",fixed=T))[2],sep="."),colnames(dat),fixed=T)
      
      #Walk up the branch until minGQ is set
      next.branch <- which(dat[,c1.idx]==branch.order$filt.1.category & 
                             dat[,c2.idx]==branch.order$filt.2.category & 
                             dat[,c3.idx]==branch.order$filt.3.category & 
                             dat[,c4.idx]==branch.order$filt.4.category & 
                             dat[,c5.idx]==branch.order$filt.5.category)
      n.minGQ <- dat$minGQ[next.branch]
      n.cond_id <- dat$condition[next.branch]
      
      if(is.na(n.minGQ)){
        next.branch <- which(dat[,c1.idx]==branch.order$filt.1.category & 
                               dat[,c2.idx]==branch.order$filt.2.category & 
                               dat[,c3.idx]==branch.order$filt.3.category & 
                               dat[,c4.idx]==branch.order$filt.4.category & 
                               dat[,f5.idx]==FALSE)
        n.minGQ <- dat$minGQ[next.branch]
        n.cond_id <- dat$condition[next.branch]
        
        if(is.na(n.minGQ)){
          next.branch <- which(dat[,c1.idx]==branch.order$filt.1.category & 
                                 dat[,c2.idx]==branch.order$filt.2.category & 
                                 dat[,c3.idx]==branch.order$filt.3.category & 
                                 dat[,f4.idx]==FALSE & 
                                 dat[,f5.idx]==FALSE)
          n.minGQ <- dat$minGQ[next.branch]
          n.cond_id <- dat$condition[next.branch]
          
          if(is.na(n.minGQ)){
            next.branch <- which(dat[,c1.idx]==branch.order$filt.1.category & 
                                   dat[,c2.idx]==branch.order$filt.2.category & 
                                   dat[,f3.idx]==FALSE & 
                                   dat[,f4.idx]==FALSE & 
                                   dat[,f5.idx]==FALSE)
            n.minGQ <- dat$minGQ[next.branch]
            n.cond_id <- dat$condition[next.branch]
            
            if(is.na(n.minGQ)){
              next.branch <- which(dat[,c1.idx]==branch.order$filt.1.category & 
                                     dat[,f2.idx]==FALSE & 
                                     dat[,f3.idx]==FALSE & 
                                     dat[,f4.idx]==FALSE & 
                                     dat[,f5.idx]==FALSE)
              n.minGQ <- dat$minGQ[next.branch]
              n.cond_id <- dat$condition[next.branch]
              
              if(is.na(n.minGQ)){
                next.branch <- which(dat[,f1.idx]==FALSE & 
                                       dat[,f2.idx]==FALSE & 
                                       dat[,f3.idx]==FALSE & 
                                       dat[,f4.idx]==FALSE & 
                                       dat[,f5.idx]==FALSE)
                n.minGQ <- dat$minGQ[next.branch]
                n.cond_id <- dat$condition[next.branch]
              }
            }
          }
        }
      }
      res$minGQ[i] <- n.minGQ
      res$source[i] <- n.cond_id
    }
  }
  
  #Return minGQ filtering table
  return(res)
}



################
###RSCRIPT BLOCK
################
require(optparse,quietly=T)

###List of command-line options
option_list <- list()

###Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog CONDITIONS STATS ROC_OPTS OUT_TREE OUT_TABLE",
                                option_list=option_list),
                   positional_arguments=TRUE)
options <- args$options


###Checks for appropriate positional arguments
if(length(args$args) != 5){
  stop("Incorrect number of required positional arguments\n")
}

###Writes args & opts to vars
CONDITIONS <- args$args[1]
STATS <- args$args[2]
ROC_OPTS <- args$args[3]
OUT_TREE <- args$args[4]
OUT_TABLE <- args$args[5]

###Process input data
dat <- process.inputs(CONDITIONS,STATS,ROC_OPTS)

###Build tree
order.tree <- build.tree(dat,fixed.f1="filt.SVTYPE")

###Create minGQ filtering lookup table
minGQ.table <- prune.tree(dat,order.tree)

###Write out final files
colnames(order.tree)[1] <- "#filt.1"
write.table(order.tree,OUT_TREE,col.names=T,row.names=F,quote=F,sep="\t")
colnames(minGQ.table)[1] <- "#condition"
write.table(minGQ.table,OUT_TABLE,col.names=T,row.names=F,quote=F,sep="\t")
