#!/usr/bin/env Rscript

#---------------------------
# Batch CNV Interval Genotyping
# Talkowski Laboratory
#
# Harrison Brand, Ryan Collins, and Joseph Glessner
# Update May 2017 (implementation for SV-pipeline)
# Update June 2016 (clean modules, smaller windows, multiallelic binning, denovo visual)
# Update October 2015 (load sample set prior)
# Update August 2015 (for incorporation into Holmes liWGS-SV 1.0)
#---------------------------

#Loads required packages; installs if necessary
RPackages <- c("optparse", "plyr", "MASS", "zoo","methods","metap", "e1071", "fpc", "BSDA", "DAAG", "pwr", "reshape", "perm", "hash")
for (i in RPackages)
{
  if (i %in% rownames(installed.packages()) == FALSE) {
    response <- readline("Install Required package (Y/N):")
    if (response == "Y") {
      install.packages((i), repos = "http://cran.rstudio.com")
      library(i, character.only = TRUE)
    } else {
      stop (paste("Unable to run script without package: ", i, sep = ""))
    }
  } else {
    library(i, character.only = TRUE)
  }
}

##build a list of command line options##
list <- structure(NA, class = "result")
"[<-.result" <- function(x, ..., value) {
  args <- as.list(match.call())
  args <- args[-c(1:2, length(args))]
  length(value) <- length(args)
  for (i in seq(along = args)) {
    a <- args[[i]]
    if (!missing(a))
      eval.parent(substitute(a <- v, list(a = a, v = value[[i]])))
  }
  x
}

#Command line options

option_list = list(
  make_option(c("-b", "--bed"), type="character", default=NULL,
              help="Bed file of CNVs to check. No header. Locus ID as fourth column. SampleIDs of interest comma delimited as fifth column. CNVtype (DEL,DUP) as the sixth column", metavar="character"),
  make_option(c("-c", "--coveragefile"), type="character", default=NULL,
              help="Full path to 1kb or 100bp binned coverage matrix for entire cohort", metavar="character"),
  make_option(c("-m", "--medianfile"), type="character", default=NULL,
              help="Full path to median intensity file with values for entire cohort", metavar="character"),
  make_option(c("-f", "--famfile"), type="character", default=NULL,
              help="Fam file FamID IndividualID(InCNVCallFile) FatherID MotherID Gender(1=male,2=female) Affected(1=unaffected,2=affected,-9=exclude)", metavar="character"),
  make_option(c("-o", "--outFolder"), type="character", default="./",
              help="Optional:Output folder", metavar="character"),
  make_option(c("-n", "--outputname"), type="character", default="out",
              help="Optional: Output file name for genotyping matrix.", metavar="character"),
  make_option(c("-r", "--refgeno"), type="character", default=NULL,
              help="Optional: File with precomputed genotype cutoffs; Requires -g TRUE", metavar="character"),
  make_option(c("-v", "--geno_call"), type="logical", default=FALSE,
              help="Optional:Only assign genotype for samples with SV, all others are assigned reference ; Requires -g TRUE. Default:FALSE", metavar="logical"),
  make_option(c("-g", "--rungenotype"), type="logical", default=FALSE,
              help="Optional:Peform genotyping on the cohort Default:FALSE", metavar="logical"),
  make_option(c("-d", "--denovo"), type="logical", default=FALSE,
              help="Optional:Call de novo per family (must only be single sample) Default:FALSE", metavar="logical"),
  make_option(c("-i", "--bins"), type="numeric", default=10,
              help="Optional:Number of bins", metavar="numeric"),
  make_option(c("-p", "--plot"), type="logical", default=FALSE,
              help="Optional:Plot JPG visualizations of CNV. Default:FALSE", metavar="logical"),
  make_option(c("-a", "--plotfamily"), type="logical", default=FALSE,
              help="Optional:Plot family based JPG visualizations; Requires -d TRUE. Default:FALSE", metavar="logical"),
  make_option(c("-j", "--runKmeans"), type="logical", default=FALSE,
              help="Optional: Run Kmeans", metavar="logical"),
  make_option(c("-e", "--Kintervalstart"), type="numeric", default=0.1,
              help="Optional:Lowest intesity diffrence between centers you want to test in kmeans Default:0.5", metavar="numeric"),
  make_option(c("-q", "--Kintervalend"), type="numeric", default=1,
              help="Optional:Highest intesity diffrence between centers you want to test in kmeans Default:0.5", metavar="numeric"),
  make_option(c("-t", "--Kinterval"), type="numeric", default=0.1,
              help="Optional:Intervals of intesity to test between interval start and end for example (start=0.1, end=0.5, interval=0.1) will test centers seperated by each of the following (0.1,0.2,0.3,0.4,0.5) default=0.5", metavar="numeric"),
  make_option(c("-k", "--plotK"), type="logical", default=FALSE,
              help="Optional:Plot JPG visualization of kmeans results (requires -j TRUE) . Default:FALSE", metavar="logical"),
  make_option(c("-s", "--sizefilter"), type="numeric", default=1000000,
              help="Optional:Restrict to large CNV to inner specified size Default:1000000", metavar="numeric"),
  make_option(c("-l", "--Blacklist"), type="character", default=NULL,
              help="Optional:Single column file with blacklist of samples to remove", metavar="character"),
  make_option(c("-w", "--Whitelist"), type="character", default=NULL,
              help="Optional:Single column file with whitelist of samples to include", metavar="character")
);

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

##QC check, see file inputs exist and are formated correctly and edit if neccessary##
if (is.null(opt$bed) || is.null(opt$coveragefile)) {
  print_help(opt_parser)
  stop("At least two arguments must be supplied (input bed and coverage file).",
       call. = FALSE)
}

##Add slash to output folder if necessary##
if (substr(opt$outFolder, nchar(opt$outFolder), nchar(opt$outFolder) + 1) !=
    "/") {
  opt$outFolder = paste(opt$outFolder, "/", sep = "")
} 

##check bed##
#Loads regions: chr start end locusID sampleID1,sampleID2,...
intervals <- read.table(opt$bed, sep = "\t", header = F)

#Make start and end numeric
intervals[, c(2:3)] <-
  apply(intervals[, c(2:3)], 2, function(vals) {
    return(as.numeric(as.character(vals)))
  })

##CNVtype check##
if (length(grep("del",intervals[,6],ignore.case = TRUE))+
    length(grep("dup",intervals[,6],ignore.case = TRUE))<length(intervals[,6]))
  {
    stop("WARNING: Incorrect CNV type specified")
  }

#Checks to make sure start less than end
if (length(which(intervals[,3]-intervals[,2]<0 ))>0) {
  stop("INPUT ERROR: Improper input coordinates. End must be greater than start.")
}

#Make sure cov bed is linked to gzipped file and tabix ready##
if (length(grep ("gz$",opt$coveragefile))<1)
{
  stop("Error cov file is not gzipped")
}

if (length(paste(opt$coveragefile,".tbi",sep=""))<1)
{
  stop("Error cov file is missing tabix index")
}


#Ensure Standard Decimal Notation
options(scipen = 1000)

##RdTest functions

#Rebinning helper function (df=dataframe,compression amount)
rebin <- function(df, compression) {
  Chr <- df[1, 1]
  Start <- df[1, 2]
  End <- df[compression, 3]
  for (i in 2:(floor(nrow(df) / compression))) {
    Chr <- c(Chr, as.character(df[((i - 1) * compression) + 1, 1]))
    Start <- c(Start, as.integer(df[((i - 1) * compression) + 1, 2]))
    End <- c(End, as.integer(df[i * compression, 3]))
  }
  newvals <- apply(df[, 4:ncol(df)], 2,
                   function(vals, compression) {
                     newcol <- sum(vals[1:compression])
                     for (i in 2:(floor(length(vals) / compression))) {
                       newcol <-
                         c(newcol, as.integer(sum(vals[(((i - 1) * compression) + 1):(i * compression)])))
                     }
                     return(newcol)
                   }, compression)
  return(as.data.frame(cbind(Chr, Start, End, newvals)))
}


#Reads coverage of specific queried region and compresses to a reasonable number of bins to create a region coverage matrix
#sampleIDs is comma specficed list of samples##

loadData <- function(chr, start, end, cnvID, sampleIDs,coveragefile,medianfile,bins)
  {
    #Take the coverage matrix header and tabix query the region in the .gz coverage matrix
    cov1 <-read.table(pipe(paste("tabix -h ",coveragefile," ", chr, ":", start, "-", end, " | sed 's/^#//'", sep = "")),sep = "\t", header = TRUE, check.names = FALSE)
    #Load plotting values if median coverage file generated by bincov##
    allnorm <- read.table(medianfile, header = TRUE, check.names = FALSE)
    ##remove when start or end pull in extra tabix line##
    cov1<-cov1[cov1$start!=end,]
    cov1<-cov1[cov1$end!=start,]
    #Check if no data
    if (nrow(cov1) < 1) {
      return("Failure")
    }
    #Find window bin size
    BinSize <- cov1$end[1] - cov1$start[1]
    
    ##Find variants with with some missing bins because bincov blacklist##
    if (nrow(cov1) < ((end - start) / BinSize)) {
      Rfinal = round_any(end, BinSize, floor)
      Rbeg = round_any(start, BinSize, ceiling)
      column_start = matrix(seq(Rbeg, Rfinal, by = BinSize), ncol = 1)
      column_end = matrix(seq(Rbeg + BinSize, Rfinal + BinSize, by = BinSize), ncol = 1)
      ncov_col = ncol(cov1)
      null_model <-
        cbind(chr, column_start, column_end, matrix(rep(0, times = nrow(column_start) *
                                                          (ncov_col - 3)), ncol = ncov_col - 3))
      colnames(null_model) <- colnames(cov1)
      covall <- rbind(cov1, null_model)
      cov1 <- covall[!duplicated(covall[, 2]), ]
      cov1 <- cov1[order(cov1[, 2]), ]
      ##Use sapply to convert files to numeric only more than one column in cov1 matrix. If not matrix will already be numeric##  
      if (nrow(cov1) > 1) {
        cov1 <- data.frame(sapply(cov1, as.numeric), check.names = FALSE)
      } else {cov1<-data.frame(t(sapply(cov1,as.numeric)),check.names=FALSE)}
    }
    
    
    #Round down the number of used bins events for smaller events (e.g at 100 bp bins can't have 10 bins if event is less than 1kb)
    if ((round_any(end, BinSize, floor) - round_any(start, BinSize, ceiling)) < bins * BinSize)
    {
      bins = (round_any(end, BinSize, floor) - round_any(start, BinSize, ceiling)) /
        BinSize
      if (bins <= 1)
      {
        Rstart <- round_any(start, BinSize, floor)
        Rend <- round_any(end, BinSize, ceiling)
        compression = 1
      }
    }
    
    #Round bins to ensure even compression (10 bins at 100bp should have a Rstart-Rend of 1kb)##
    if (!exists("compression"))
    {
      UnadjustedBins <-
        (round_any(end, BinSize, floor) - round_any(start, BinSize, ceiling)) /
        (bins * BinSize)
      RemainderForRemoval <-
        ##Need to account for round error by trunc so add the decimal####
      trunc(((
        UnadjustedBins - trunc(UnadjustedBins)
      ) * BinSize * bins / 2) + 0.000000001)
      RemainderFront <-
        round_any(RemainderForRemoval, BinSize, floor)
      RemainderBack <-
        round_any(RemainderForRemoval, BinSize, ceiling)
      Rstart <-
        round_any(start, BinSize, ceiling) + RemainderFront
      Rend <-
        round_any(end, BinSize, floor) - RemainderBack
      compression <- (Rend - Rstart) / (BinSize * bins)
    }
    #Cut bins down to those required for compressed clean size based on Rstart and Rend##
    cov1<-cov1[which(cov1[,3]>Rstart & cov1[,2]<Rend ),]
    
    #Samples Filter
    if ( !is.null(opt$Blacklist) ) {
      samplesBlacklist <- readLines(opt$Blacklist)
      IDsSamplesBlacklistButNotCoverage <-
        samplesBlacklist[!(samplesBlacklist %in% names(cov1))]
      if (length(IDsSamplesBlacklistButNotCoverage) > 0)
      {
        cat (
          " WARNING: IDs in samplesBlacklist but not coverage:",
          IDsSamplesBlacklistButNotCoverage
        )
      }
      ##Filter samples based of specified blacklist here##
      cov1 <- cov1[,!(names(cov1) %in% samplesBlacklist)]
      allnorm <-
        allnorm[,!(names(allnorm) %in% samplesBlacklist)]
    }
    
    ##Allow whitelist##
    if (!is.null(opt$Whitelist)) {
      samplesWhitelist <- readLines(opt$Whitelist)
      IDsSamplesWhitelistButNotCoverage <-
        samplesWhitelist[!(samplesWhitelist %in% names(cov1))]
      if (length(IDsSamplesWhitelistButNotCoverage) > 0)
      {
        cat (
          " WARNING: IDs in samplesWhitelist but not coverage:",
          IDsSamplesWhitelistButNotCoverage
        )
      }
      ##make sure to still include first three lines###
      cov1 <- cov1[,names(cov1) %in% c("chr","start","end",samplesWhitelist)]
      allnorm <-
      allnorm[, (names(allnorm) %in% samplesWhitelist)]
    }
    if (ncol(cov1) < 4)
    {
      stop (" WARNING: All samples excluded by filtering")
    } 
    #Approximates rebinned per-sample medians (approximated for speed & memory)
    allnorm[which(allnorm == 0)] <- 1
    allnorm <- compression * allnorm
    
    #Replace zero values with 1 for handling normalization
    cov1[cov1 == 0] <- 1
    
    #Rebins values
    if (compression > 1) {
      res <-
        rebin(cov1, compression)
      res <-
        apply(res[, 4:ncol(res)], 2, function(val) {
          as.numeric(as.matrix(val))
        })
    } else {
      res <- cov1[, 4:ncol(cov1)]
    }
  
    #Adds sample medians to df
    res0<-rbind((res), allnorm)
    
    #Scale each col within that sample
    res1<- apply(res0,2,
                 function(vals){
                   return(as.numeric(vals[1:(nrow(res0)-1)])/as.numeric(vals[nrow(res0)]))
                 })
  
    #need to transpose if more than one bin assessed
    if (ncol(as.matrix(res1)) > 1) {
      cnv_matrix <- t(res1)
    } else {
      cnv_matrix <- as.matrix(res1)
    }
    return(cnv_matrix)
  }

#Loads specified sample set in genotyping matrix based on the specified cnv type (del=1,dup=3) and unspecified samples as cn=2
#sampleIDs is comma specficed list of samples##
specified_cnv <- function(cnv_matrix, sampleIDs, cnvID, chr, start, end, cnvtype)
  {
    CNV <- matrix(c(cnvID, chr, start, end), nrow = 1)
    genotype_matrix <- cbind(CNV, t(matrix(seq(1, nrow(cnv_matrix)))))
    colnames(genotype_matrix) <- c("ID", "Chr", "Start", "End", rownames(cnv_matrix))
    samplenames <- colnames(as.matrix(genotype_matrix))
    columnswithsamp <- which(colnames(genotype_matrix) %in% unlist(strsplit(as.character(sampleIDs),split=","))) 
    if (length(columnswithsamp)==0) {
      ##"WARNING: No samples in coverage matrix for comparision check black/whitelist"##
      return ("No_Samples")
    }
    
    ##create genotype matrix##
    if (toupper(cnvtype) == "DEL")
    {
      genotype_matrix[1, columnswithsamp] = 1
      ##make sure first four columns are not modified##
      columnswithsamp <- c(columnswithsamp, 1, 2, 3, 4)
      genotype_matrix[1,-columnswithsamp] = 2
    } else if (toupper(cnvtype) == "DUP")
    {
      genotype_matrix[1,  columnswithsamp] = 3
      ##make sure first four columns are not modified##
      columnswithsamp <- c(columnswithsamp, 1, 2, 3, 4)
      genotype_matrix[1,-columnswithsamp] = 2
    } 
    return(genotype_matrix)
  }

###Kmeans multi-CNV Test##
##interval is measured by predicted copy state##
kMeans <-function(cnv_matrix,chr,start,end,cnvID,Kinterval,Kintervalstart,Kintervalend,outFolder,outputname)
  {
    samplenames <- rownames(cnv_matrix)
    #create Eucledian matrix###
    eucledianM <- dist(cnv_matrix, method = "euclidean")
    # counts clusters for different attempts at starting points for k-means k values
    ks = 0
    #avg. silhouette width to measuer success of each run
    avg.silwidth = 0
    count=0
    if (length(samplenames) > 100) {
      totalcopystate <- 100
    } else{
      totalcopystate <- length(samplenames) 
    }
    for (i in seq(Kintervalstart, Kintervalend, Kinterval)) {
      count=count+1  
      ##Need to make sure there are not more copy states than samples (max at a copy state of 100)##
      ##Run Kmeans##
      k <-kmeans(
        cnv_matrix,
        ##center assignment##    
        matrix(rep(seq(0, i*totalcopystate , by = i),ncol(cnv_matrix)),ncol = ncol(cnv_matrix)),
        iter.max = 100,
        algorithm = "Forgy"
      )
      # Number of clusters
      ##Count number of centers with values##
      a <- (count(k$centers > 0))
      ks[count] = a$freq[1]
      ##Get Cluster Stats and average silwidth#
      if (a$freq[1] > 1)
      {
        ClustSolStats = suppressWarnings(cluster.stats(eucledianM, k$cluster))
        avg.silwidth[count] = ClustSolStats$avg.silwidth
      } else  {
        avg.silwidth[count] = 0
      }
    }
    ##Select best cluster width###
    avg.silwidthINT = max(which(avg.silwidth == max(na.omit(avg.silwidth))))
    avg.silwidthClust = ks[avg.silwidthINT]
    finalinterval <-
      seq(Kintervalstart, Kintervalend, Kinterval)[avg.silwidthINT]
    finalk <- kmeans(cnv_matrix,
                     ##center assignment##
                     matrix(rep(
                       seq(0, finalinterval * totalcopystate , by = finalinterval),
                       ncol(cnv_matrix)
                     ), ncol = ncol(cnv_matrix)),
                     iter.max = 100,
                     algorithm = "Forgy")
    ##Make file to output K's##
    KclusterAsGenotype = cbind(t(as.data.frame(list(
      c(
        ID = cnvID,
        Chr = chr,
        Start = start,
        End = end
      )
    ))), t(as.matrix(finalk$cluster)))
    ##Count the number of copy states per SV##
    Kclustercount=matrix(c(KclusterAsGenotype[,1:4],length(unique(KclusterAsGenotype[,5:ncol(KclusterAsGenotype)]))),nrow=1)
    colnames(Kclustercount)<-c("ID","Chr","Start","End","N_CopyStates")
    if(file.exists(paste(outFolder,outputname,".clustercount",sep=""))) {
      #write.table(KclusterAsGenotype,paste(outFolder,outputname,".K",sep=""),quote=FALSE,append=TRUE,row.names=FALSE,col.names=FALSE,sep= "\t")  
      write.table(Kclustercount,paste(outFolder,outputname,".clustercount",sep=""),quote=FALSE,append=TRUE,row.names=FALSE,col.names=FALSE,sep= "\t")  
    } else { 
      #write.table(KclusterAsGenotype,paste(outFolder,outputname,".K",sep=""),quote=FALSE,row.names=FALSE,sep= "\t")  
      write.table(Kclustercount,paste(outFolder,outputname,".clustercount",sep=""),quote=FALSE,row.names=FALSE,sep= "\t")  
    }
  
    return(KclusterAsGenotype)
  }

#Seperate Samples into either Control or Treat group 
#Number of bins assessed is dependent on SV sample size
create_groups <- function(genotype_matrix, cnv_matrix)
{
  ##Remove outer two bins which tend to be noisy if CNV large enough##
  if (ncol(cnv_matrix) == 1) {
    ##a start bin, b end bin##
    Control <-cnv_matrix[which(genotype_matrix[, 5:ncol(genotype_matrix)] == 2), 1]
    Treat <- cnv_matrix[which(genotype_matrix[, 5:ncol(genotype_matrix)] != 2), 1]
    a<- 1
    b<- 1
  } else if (ncol(cnv_matrix) > 1 && ncol(cnv_matrix) < 4) {
    a <- 1
    b <- ncol(cnv_matrix)
    Control <- apply(cnv_matrix[which(genotype_matrix[, 5:ncol(genotype_matrix)] == 2), a:b, drop = F], 1, median)
    Treat <-apply(cnv_matrix[which(genotype_matrix[, 5:ncol(genotype_matrix)] != 2), a:b, drop =F], 1, median)
  } else {
    a <- 2
    b <- ncol(cnv_matrix) - 1
    Control <- apply(cnv_matrix[which(genotype_matrix[, 5:ncol(genotype_matrix)] == 2), a:b, drop = F], 1, median)
    Treat <- apply(cnv_matrix[which(genotype_matrix[, 5:ncol(genotype_matrix)] != 2), a:b, drop = F], 1, median)
  }
  output <- list(Control, Treat,a,b)
  names(output) <- c("Control", "Treat","a","b")
  return(output)
}

#Power
powerCalc <- function(genotype_matrix, cnv_matrix)
{
  #Call Treat (have SV) and Control Groups
  Control<-create_groups(genotype_matrix, cnv_matrix)$Control
  Treat<-create_groups(genotype_matrix, cnv_matrix)$Treat
  if (length(Control) > 1 && length(Treat) > 1)
  {
    #doesn't matter less or greater since absolute deviation, use 0.5 diffrence in mean between CNV to estimate effect size
    power <- pwr.t2n.test(n1 = length(Control), n2 = length(Treat), sig.level = 0.05,
                          alternative = "greater", d = (0.5 / sd(Control)))$power
  } else {
    power <- NA
  }
  return(power)
}

#OneSamplezscore to test single sample against everyone else; Can specifiy list of samples to exclude from analysis which may have CNV##
#samples exclude should be comma delimited list 
#singlesample being assessed must not be normal(cn=2) in the genotype matrix### 
onesamplezscore.median <- function(genotype_matrix,cnv_matrix,singlesample,cnvtype)
{
  #Call Treat (have SV) and Control Groups
  Control<-create_groups(genotype_matrix, cnv_matrix)$Control
  Treat<-create_groups(genotype_matrix, cnv_matrix)$Treat
  Treat<-Treat[singlesample]
  a<-create_groups(genotype_matrix, cnv_matrix)$a
  b<-create_groups(genotype_matrix, cnv_matrix)$b
  ##Calculate one-sided z score##
  if (toupper(cnvtype) == "DEL") {
    ztest.p <- pnorm((Treat - mean(Control)) / sd(Control))
  } else{
    ztest.p <- pnorm((mean(Control) - Treat) / sd(Control))
  }
  ##Find the secondest worst p-value and record as an assement metric## 
  plist <- c()
  i = 1
  for (column in a:b)
  {
    Control2 <-
      cnv_matrix[which(genotype_matrix[, 5:ncol(genotype_matrix)] == 2), column]
    Treat2 <-
      cnv_matrix[singlesample, column]
    if (toupper(cnvtype) == "DEL") {
      single.p <- pnorm((Treat2 - mean(Control2)) / sd(Control2))
    } else {
      single.p <- pnorm((mean(Control2) - Treat2) / sd(Control2))
    }
    #store diffrent z p-value by column##
    plist[i] <- single.p
    i = i + 1
  }
  if (length(plist) > 1)
  {
    mySecondMaxP <- sort(plist)[length(plist) - 1]
  } else {
    ##Note if only one bin, that bin is assigned as the second max P###
    mySecondMaxP <- plist[1]
  }
  output <- list(ztest.p, mySecondMaxP)
  names(output) <- c("singleZ_Pvalue", "Pmax_2nd")
  return(output)
}

#twosamplet t-test
twosamplezscore.median <- function(genotype_matrix,cnv_matrix,cnvtype)
{
  #Call Treat (have SV) and Control Groups
  Control<-create_groups(genotype_matrix, cnv_matrix)$Control
  Treat<-create_groups(genotype_matrix, cnv_matrix)$Treat
  a<-create_groups(genotype_matrix, cnv_matrix)$a
  b<-create_groups(genotype_matrix, cnv_matrix)$b
  if (toupper(cnvtype) == "DEL") {
    P_object <- permTS(Control, Treat, alternative = "greater", method = 'pclt')$p.value
  } else{ P_object <- permTS(Control, Treat, alternative = "less", method = 'pclt')$p.value }
  
  ##Find the secondest worst p-value and record as an assement metric#
  plist<-c()
  i=1
  for (column in a:b)
  {
    Control2 <- cnv_matrix[which(genotype_matrix[, 5:ncol(genotype_matrix)] == 2), column]
    Treat2 <- cnv_matrix[which(genotype_matrix[, 5:ncol(genotype_matrix)]!=2), column]
    if (toupper(cnvtype) == "DEL") {
      singlep <- permTS(Control2, Treat2, alternative = "greater", method = 'pclt')$p.value
    } else{
      singlep <- permTS(Control2, Treat2, alternative = "less", method = 'pclt')$p.value
    }
    #store diffrent z p-value by column##
    plist[i] <- singlep
    i=i+1
  }
  if (length(plist) > 1)
  {
    mySecondMaxP <- sort(plist)[length(plist) - 1]
  } else {
    ##Note if only one bin, that bin is assigned as the second max P###
    mySecondMaxP <- plist[1]
  }
  output<-list(P_object, mySecondMaxP)
  names(output)<-c("Pvalue","Pmax_2nd")
  return(output)  
}

##Provide a depth based rank of the sample##
##sample you want to pull out information, if NULL than will do treat vs control##
samprank_sep <- function(genotype_matrix,cnv_matrix,cnvtype,sample=NULL)
{
  #Call Treat (have SV) and Control Groups
  Control<-create_groups(genotype_matrix, cnv_matrix)$Control
  Treat<-create_groups(genotype_matrix, cnv_matrix)$Treat
  combined<-c(Treat,Control)
  if (toupper(cnvtype) == "DEL")
  {
    order.rank <- median(rank(combined)[1:length(Treat)])
    ##Seperation between Treatment and Control groups
    Sep = median(Control) - median(Treat)
  } else {
    order.rank <-
      length(combined) - median(rank(combined)[1:length(Treat)]) + 1
    Sep = median(Treat) - median(Control)
  }
  if (!is.null(sample)) {
    if (toupper(cnvtype) == "DEL")
    {
      order.rank = unname(rank(combined)[sample])
      Sep = median(Control) - median(combined[sample])
    } else {
      order.rank <- length(combined) - unname(rank(combined)[sample]) + 1
      Sep = median(combined[sample]) - median(Control)
    }
  }
  output <- list(order.rank, Sep)
  names(output) <- c("rank", "Sep")
  return(output)
}

##Plot of intensities across cohorts## 
plotJPG <- function(genotype_matrix,cnv_matrix,chr,start,end,cnvID,sampleIDs,outputname,cnvtype,plotK,plotfamily,famfile,outFolder)
{
  samplesPrior <- unlist(strsplit(as.character(sampleIDs),","))
  samplenames<-colnames(genotype_matrix)
  
  ##If only one bin##
  if(ncol(cnv_matrix)==1)
  {cnv_matrix<-cbind(cnv_matrix,cnv_matrix[,1])}
  
  ##File Output##
  jpeg(paste(outFolder,chr,"_",start,"_",end,"_",samplesPrior[1],"_",cnvID,"_",outputname,".jpg",sep=""),res=300, width=1800, height=1800)
  
  ##concatenate sample IDs if necessary##
  sampleIDs<-paste(sampleIDs,collapse=",")
  
  ##Limits number of sample Ids due to size limiations for readablity
  if(nchar(as.character(sampleIDs))>44){sampleIDsToDisplay<-paste(substr(sampleIDs,1,44),"...",sep="")}else{sampleIDsToDisplay<-sampleIDs}
  ##Title line 1##
  main1=paste("chr",chr,":",prettyNum(start,big.mark=","),"-",prettyNum(end,big.mark=",")," (hg19)",sep="")
  
  ###Add proper size abbr. for larger events
  size=end-start
  if(size<10000){mysize<-prettyNum(paste("(",size," bp)",sep=""), big.mark = ",")}
  if(size>=10000){mysize<-prettyNum(paste("(",signif(size/1000,3)," kb)",sep=""), big.mark = ",")}
  if(size>=1000000){mysize<-prettyNum(paste("(",signif(size/1000000,3)," Mb)",sep=""), big.mark = ",")}
  
  ##Formating##
  main2 = paste(sampleIDsToDisplay, " ", mysize, sep ="")
  mainText = paste(main1, "\n", main2, sep = "")
  maxcexXa <- par('pin')[1] / strwidth(main1, 'inches')
  maxcexXb <- par('pin')[1] / strwidth(main2, 'inches')
  maxcexXh <- min(maxcexXa, maxcexXb)
  if (maxcexXh < 0.5) {
    maxcexXh = 0.5
  }
  if (maxcexXh > 2.5) {
    maxcexXh = 2.5
  }
  par(mar = c(6.1, 6.1, 4.1, 2.1))
  
  ##Create matrix for plotting###
  columnstoshift <- which(rownames(cnv_matrix) %in% samplesPrior) 
  ##Place Samples with CNV on top###
  plot_cnvmatrix<-cbind(t(cnv_matrix)[,-columnstoshift],t(cnv_matrix)[,columnstoshift])
  ##column shift diffrent for genotype matrix because cnvID,chr,start,end
  columnstoshift <- which(colnames(genotype_matrix) %in% unlist(strsplit(as.character(samplesPrior),split=","))) 
  plot_colormatrix<-cbind(matrix(genotype_matrix[,-columnstoshift],nrow=1),matrix(genotype_matrix[,columnstoshift],nrow=1))
  endcolnormal<-ncol(plot_colormatrix)-(length(samplesPrior))
  plot_linematrix<-cbind(matrix(genotype_matrix[,-columnstoshift],nrow=1),matrix(genotype_matrix[,columnstoshift],nrow=1))

  ##Blue if Dup; Red if Del
  if ( plotK == TRUE ) {
    #keep plot_colormatrix
    main1=paste("chr",chr,":",prettyNum(start,big.mark=","),"-",prettyNum(end,big.mark=",")," (hg19)",sep="")
    mainText = paste(main1, "\n", "K-Test"," ", mysize, sep = "")  
  } else if (toupper(cnvtype) == "DEL") {
    plot_colormatrix[, (endcolnormal + 1):ncol(plot_colormatrix)] <- "red"
    plot_colormatrix[,5:endcolnormal]<-"grey"
    plot_linematrix[, (endcolnormal + 1):ncol(plot_colormatrix)] <- "3"
    plot_linematrix[,5:endcolnormal]<-"0.5"
  } else if (toupper(cnvtype) == "DUP") {
    plot_colormatrix[, (endcolnormal + 1):ncol(plot_colormatrix)] <- "blue"
    plot_colormatrix[,5:endcolnormal]<-"grey"
    plot_linematrix[, (endcolnormal + 1):ncol(plot_colormatrix)] <- "3"
    plot_linematrix[,5:endcolnormal]<-"0.5"
  } 
  
  ##Plotting Command##
  plot(as.zoo(plot_cnvmatrix),
    plot.type = "single",
    col = plot_colormatrix[1, 5:ncol(plot_colormatrix)],
    main = mainText,
    cex.main = maxcexXh,
    xlab = "Position (bp)",
    xaxt = 'n',
    ann = FALSE,
    ylab = "Intensity",
    lwd = plot_linematrix[1, 5:ncol(plot_linematrix)] 
  )
  mtext(
    side = 1,
    text = paste("chr", chr, " Position (bp)", sep = ""),
    line = 5
  )
  mtext(side = 2, text = "Normalized Read Depth Ratio", line = 3)
  myIntervalsXAxis <- round(seq(start,end,length.out=ncol(cnv_matrix)))
  axis(1,at=seq(1,ncol(cnv_matrix),by=1) ,labels = prettyNum(myIntervalsXAxis, big.mark = ","),las = 2,cex.axis = 0.8)
  ##Family-Based Plotting##
  if (plotfamily == TRUE ) {
    ##Call familes to plot###
    ##May have issues with multi-generation pedigress, Designed for Quad and Trio Families##
    family <- read.table(famfile)
    includedfams <-
      unique(family[which(family[, 2]  %in% samplesPrior), 1])
    proband_list <-as.character(
      family[which(family[, 1] %in% includedfams &
                     family[, 3] != 0 &
                     family[, 4] != 0 & family[, 6] == 2) , 2])
    sib_list <-as.character(
      family[which(family[, 1] %in% includedfams &
                     family[, 3] != 0 &
                     family[, 4] != 0 & family[, 6] == 1) , 2])
    father_list <-as.character(
      family[which(family[, 1] %in% includedfams &
                     family[, 3] == 0 &
                     family[, 4] == 0 &  family[, 5] == 1 & family[, 6] == 1) , 2])
    mother_list <-as.character(
      family[which(family[, 1] %in% includedfams &
                     family[, 3] == 0 &
                     family[, 4] == 0 &  family[, 5] == 2 & family[, 6] == 1) , 2])
    
    text(c(1:10), as.numeric(cnv_matrix[proband_list,]), "p", cex = 1)
    text(c(1:10), as.numeric(cnv_matrix[sib_list,]), "s", cex = 1)
    text(c(1:10), as.numeric(cnv_matrix[father_list,]), "fa", cex = 1)
    text(c(1:10), as.numeric(cnv_matrix[mother_list,]), "mo", cex = 1)
  }
  if (plotK == TRUE) {
    copy_states = order(as.numeric(unique(genotype_matrix[, 5:ncol(genotype_matrix)])), decreasing = TRUE)
    legend(
      ifelse(toupper(cnvtype) == "DEL", 'topright', 'bottomright'),
      paste("CN", copy_states) ,
      lty = 1,
      col = copy_states,
      cex = .3
    )
  } else if (toupper(cnvtype) == "DEL") {
    legend(
      'topright',
      c("Deletion", "Diploid"),
      lty = 1,
      col = c("red", "grey"),
      cex = .5
    )
  } else {
    legend(
      'topright',
      c("Diploid", "Duplication"),
      lty = 1,
      col = c("grey", "blue"),
      cex = .5
    )
  }
  dev.off()
}

##Provide genotype for VCF format##
genotype<- function(cnv_matrix,genotype_matrix,geno_call,refgeno,chr,start,end,cnvID,sampleIDs,cnvtype,outFolder,outputname)
{
 ##get depth intensities##  
 cnv_median <-c(create_groups(genotype_matrix, cnv_matrix)$Control,create_groups(genotype_matrix, cnv_matrix)$Treat) 
 ##order by names so same geno output for each variant##
 cnv_median<-cnv_median[order(names(cnv_median))]
 if (exists("refgeno") != TRUE) {
   ##use expected cutoffs for genotypeing (0.25,.75,1.25,1.75)
   cutoffs <- seq(0.25, 1.75, by = .5)
 } else {
   cutoffs <- unlist(c(read.table(refgeno, header = TRUE)[2]))
 }
 ##CN0##
 if (toupper(cnvtype) == "DEL") {
   cnv_median[which(cnv_median < cutoffs[1])] <- "1/1"
   ##CN1##
 }
 if (toupper(cnvtype) == "DEL") {
   cnv_median[which(cnv_median < cutoffs[2])] <- "0/1"
   ##CN3##
 }
 if (toupper(cnvtype) == "DUP") {
   cnv_median[which(cnv_median > cutoffs[3])] <- "0/1"
   ##CN4##
 }
 if (toupper(cnvtype) == "DUP") {
   cnv_median[which(cnv_median > cutoffs[4])] <- "1/1"
   ##CN2##
 }
 ##Assing everything else 0/0###
 if (length(grep("/1", cnv_median)) == 0) {
   cnv_median[which(cnv_median >= 0)] <- "0/0"
 } else {
   cnv_median[-grep("/1", cnv_median)] <- "0/0"
 }
 if (geno_call==TRUE)
 {
   samplesPrior<-unlist(strsplit(as.character(sampleIDs),","))
   cnv_median[-which(names(cnv_median) %in% samplesPrior)] <- "0/0"
   ##All predicted samples must have a genotype so assign 0/0 as 0/1###
   cnv_median[which(names(cnv_median) %in% samplesPrior & cnv_median=="0/0")] <- "0/1"
 }
 
 if(!file.exists(paste(outFolder,outputname,".geno",sep=""))) {
   ##write header##
   write.table(matrix(c("chr","start","end","cnvID",names(cnv_median)),nrow=1),paste(outFolder, outputname, ".geno", sep = ""),quote=FALSE,row.names=FALSE,col.names=FALSE,sep= "\t")  
 } 
 write.table(matrix(c(chr, start, end, cnvID,cnv_median),nrow=1),paste(outFolder, outputname, ".geno", sep = ""),
             quote = FALSE,col.names = FALSE, row.names = FALSE,append=TRUE,sep= "\t") 
}

runRdTest<-function(bed)
{ 
  chr<-as.character(bed[1])
  start<-as.numeric(bed[2])
  end<-as.numeric(bed[3])
  cnvID<-as.character(bed[4])
  sampleIDs<-as.character(bed[5])
  sampleOrigIDs<-as.character(bed[5])
  cnvtype<-as.character(bed[6])
  cnvtypeOrigIDs<-as.character(bed[6])
  
  ##Assign input values from opt list to variable##
  for (names in names(opt))
  {
    assign(names,unname(unlist(opt[names])))
  }
  #Speed up large cnvs by taking inner range of laregest desired size
  if (end - start  > sizefilter)
  {
    cat(paste(chr,":",start,"-",end,":Large size so subsampling in middle\n",sep=""))
    center=(start + end) / 2 
    start = round(center - (sizefilter/2))
    end = round(center + (sizefilter/2))
  } 
  
  if (end - start <= 0 )
  {
   end=start+1 
  }  
  ##Make sure region is in tabix##
  
  
  ##Get Intesity Data##
  cnv_matrix<-loadData(chr, start, end, cnvID, sampleIDs,coveragefile,medianfile,bins)
  
  if (cnv_matrix[1]=="Failure") {
    return(c(chr,start,end,cnvID,sampleOrigIDs,cnvtypeOrigIDs,"coverage_failure","coverage_failure","coverage_failure","coverage_failure","coverage_failure","coverage_failure"))
  }
  
  ##remove black or white list samples from sampleIDs###
  idsforsearch<-rownames(cnv_matrix)
  samplestokeep<-match(unlist(strsplit(sampleIDs,",")),idsforsearch)
  sampleIDs<-idsforsearch[na.omit(samplestokeep)]
  samplesPrior <-unlist(strsplit(as.character(sampleIDs),split=","))
  ##Run K Test if Specified##
  if (opt$runKmeans == TRUE) {
    k_matrix<-kMeans(cnv_matrix,chr,start,end,cnvID,Kinterval,Kintervalstart,Kintervalend,outFolder,outputname)
    if (opt$plotK==TRUE) {
      plotJPG(k_matrix,cnv_matrix,chr,start,end,cnvID,sampleIDs,outputname,cnvtype,plotK,plotfamily=FALSE,famfile,outFolder)
    }
  }
  ##Assign intial genotypes (del=1,dup=3,diploid=2)##
  genotype_matrix<-specified_cnv(cnv_matrix, sampleIDs, cnvID, chr, start, end, cnvtype)
  ##check if no samples are found in genotype matrix##
  if (as.matrix(genotype_matrix)[1,1]=="No_Samples") {
    return(c(chr,start,end,cnvID,sampleOrigIDs,cnvtypeOrigIDs,"No_samples_for_analysis","No_samples_for_analysis","No_samples_for_analysis","No_samples_for_analysis","No_samples_for_analysis","No_samples_for_analysis"))
  }
  
  ##QC on filtered sample counts##
  copystatecounts=table(genotype_matrix[1,5:ncol(genotype_matrix)])
  ##diploid count##
  dipcount=copystatecounts["2"]
  cnvcount=copystatecounts[ifelse(toupper(cnvtype)=="DEL","1","3")]
  if (is.na(dipcount)){
        return(c(chr,start,end,cnvID,sampleOrigIDs,cnvtypeOrigIDs,"All_samples_called_CNV_no_analysis","All_samples_called_CNV_no_analysis","All_samples_called_CNV_no_analysis","All_samples_called_CNV_no_analysis","All_samples_called_CNV_no_analysis","All_samples_called_CNV_no_analysis"))
  } 
  ##genotype and write to file##
  if (opt$rungenotype == TRUE) {
  genotype(cnv_matrix,genotype_matrix,geno_call,refgeno,chr,start,end,cnvID,sampleIDs,cnvtype,outFolder,outputname)
  }
  ##Plot JPG##
  if (opt$plot == TRUE){
    plotJPG(genotype_matrix,cnv_matrix,chr,start,end,cnvID,sampleIDs,outputname,cnvtype,plotK=FALSE,plotfamily=FALSE,famfile,outFolder)
  }
  ##De Novo Module##
  if (opt$denovo == TRUE) {
    ##Read in family file##
    family <- read.table(famfile)
    samplesPrior <- samplesPrior[which(grepl("s1|p1",samplesPrior)==TRUE)]
    ##If ID only has parents or all children removed by filtering##
    if (length(samplesPrior) == 0 ) {
      denovo_output <- cbind(chr,start,end,cnvID,cnvtype,"No_samples_for_analysis","No_samples_for_analysis","No_samples_for_analysis","No_samples_for_analysis","No_samples_for_analysis","No_samples_for_analysis","No_samples_for_analysis","No_samples_for_analysis","No_samples_for_analysis","No_samples_for_analysis",
                             "No_samples_for_analysis","No_samples_for_analysis","No_samples_for_analysis","No_samples_for_analysis",
                             "No_samples_for_analysis","No_samples_for_analysis","No_samples_for_analysis","No_samples_for_analysis")
      
      return(denovo_output)
    }
    includedfams <- unique(family[which(family[, 2]  %in% samplesPrior), 1])
    original_cnv_matrix<-cnv_matrix
    original_genotype_matrix<-genotype_matrix
    ##Find samples to include##
    proband_list <-as.character(
      family[which(family[, 1] %in% includedfams & family[, 3] != 0 & family[, 4] != 0 & family[, 6] == 2) , 2])
    sib_list <-as.character(
      family[which(family[, 1] %in% includedfams & family[, 3] != 0 & family[, 4] != 0 & family[, 6] == 1) , 2])
    father_list <-as.character(
      family[which(family[, 1] %in% includedfams & family[, 3] == 0 & family[, 4] == 0 &  family[, 5] == 1 & family[, 6] == 1) , 2])
    mother_list <-as.character(
      family[which(family[, 1] %in% includedfams & family[, 3] == 0 & family[, 4] == 0 &  family[, 5] == 2 & family[, 6] == 1) , 2])
    for (mem in c("mo","p1","s1","fa")) {
      eval(parse(text=paste(mem,".p.list<-c()",sep="")))
      eval(parse(text=paste(mem,".secmaxp.list<-c()",sep="")))
      eval(parse(text=paste(mem,".rankp.list<-c()",sep="")))
      eval(parse(text=paste(mem,".sepp.list<-c()",sep="")))
    } 
    affecteded_fam<-c()
    count=0
    for (i in includedfams) {
      count = count + 1
      for (singlesample in c(proband_list[count],sib_list[count],father_list[count],mother_list[count])) {
        mem=unlist(strsplit(singlesample,split=".",fixed=TRUE))[2]
        ##Add family members to sample include list##
        sampleID1s = unique(c(
            as.character(sampleIDs),
            father_list[count],
            mother_list[count],
            sib_list[count],
            proband_list[count]))
        ##gender restrict variants on X or Y###
        if ((chr == "X" ||
             chr == "Y")  &&
            (family[which(family[, 2] == singlesample), 5] == 2)) {
          cnv_matrix <-
            as.matrix(cnv_matrix[rownames(cnv_matrix)  %in%  family[which(family[, 5] == 2), 2], ])
        } else if ((chr == "X" ||
                    chr == "Y")  &&
                   (family[which(family[, 2] == singlesample), 5] == 1)) {
          cnv_matrix <-
            as.matrix(cnv_matrix[rownames(cnv_matrix)  %in%  family[which(family[, 5] == 1), 2], ])
        }
        ##remove sample of interest from sample exclude list and make new genotype matrix##
        genotype_matrix<-specified_cnv(cnv_matrix, sampleID1s, cnvID, chr, start, end, cnvtype)
        ##remove singlesample for exclusion list##
        p <-onesamplezscore.median(genotype_matrix,cnv_matrix,singlesample,cnvtype)
        ##write meteric for each family member
        eval(parse(text=paste(mem,".p.list[count]<-", p[1],sep="")))
        eval(parse(text=paste(mem,".secmaxp.list[count]<-", p[2],sep="")))
        rank_sep<-samprank_sep(genotype_matrix,cnv_matrix,cnvtype,singlesample)
        eval(parse(text=paste(mem,".rankp.list[count]<-", rank_sep[1],sep="")))
        eval(parse(text=paste(mem,".sepp.list[count]<-", rank_sep[2],sep="")))
        cnv_matrix<-original_cnv_matrix
      }
      affecteded_fam[count]<-paste(unique(grep(i,unlist(strsplit(as.character(sampleIDs),split=",")),value=TRUE)),collapse=",")
      if (opt$plotfamily==TRUE) {
        sampleID2s<-paste(unique(grep(i,unlist(strsplit(as.character(sampleIDs),split=",")),value=TRUE)),collapse=",")
        plotJPG(original_genotype_matrix,original_cnv_matrix,chr,start,end,cnvID,sampleIDs=sampleID2s,outputname=paste(outputname,"_",i,sep=""),cnvtype,plotK=FALSE,plotfamily,famfile,outFolder)
      }
    }
    denovo_output <- cbind(chr,start,end,cnvID,cnvtype,includedfams,affecteded_fam,p1.p.list,s1.p.list,fa.p.list,mo.p.list,
                          p1.secmaxp.list,s1.secmaxp.list,fa.secmaxp.list,mo.secmaxp.list,
                          p1.sepp.list,s1.sepp.list,fa.sepp.list,mo.sepp.list,
                          p1.rankp.list,s1.rankp.list,fa.rankp.list,mo.rankp.list)
    return(denovo_output)
  } 
  
  ##Flip samples and cnvtype to that with the lowest frequency##
  if(dipcount<cnvcount){
    cnvtype=ifelse(toupper(cnvtype)=="DEL","DUP","DEL")
    columnstoremove <- which(colnames(genotype_matrix) %in% unlist(strsplit(as.character(sampleIDs),split=","))) 
    ifelse(toupper(cnvtype)=="DEL",genotype_matrix[,-c(columnstoremove,1,2,3,4)]<-1,genotype_matrix[,-c(columnstoremove,1,2,3,4)]<-3)
    genotype_matrix[,columnstoremove]<-2
    ##Reclassify diplod as the acutal CNV##
    sampleIDs<-noquote(paste(names(genotype_matrix[,-c(columnstoremove,1,2,3,4)]),collapse=","))
    samplesPrior <-unlist(strsplit(as.character(sampleIDs),split=","))
  } 
  ##Power Calculation##
  power <- powerCalc(genotype_matrix, cnv_matrix)
  power<-ifelse(length(unlist(strsplit(as.character(sampleIDs), split = ","))) > 1,power,NA)
  if (!is.na(power) && power > 0.8) {
    p <- twosamplezscore.median(genotype_matrix, cnv_matrix, cnvtype)
    p[3]<-"twoSampPerm"
    names(p)<-c("Pvalue","Pmax_2nd","Test")
  } else {
    ##Need to break down underpowerd samples into multiple single z-tests##
    p.list<-c()
    p.2ndmax<-c()
    count=0
    for (i in samplesPrior) {
      count=count+1
      singlesample = i
      p <-onesamplezscore.median(genotype_matrix,cnv_matrix,singlesample,cnvtype)
      p.list[count]<-p[1]
      p.2ndmax[count]<-p[2]
    }
    ##Combine individual P-values with fisher.method##
    if (length(p.list) > 1) {
    ##Need to change 0 to 1e-300 for sumlog function##  
      p.list<-rapply(p.list,function(x) ifelse(x==0,1e-300,x), how = "replace")
      p.2ndmax<-rapply(p.2ndmax,function(x) ifelse(x==0,1e-300,x), how = "replace")
      p <- list(sumlog(unlist(p.list))$p, sumlog(unlist(p.2ndmax))$p)
    } else {
      p <- c(p.list[1], p.2ndmax[1])
    }
    p[3]<-"singlesampZ"
    names(p)<-c("Pvalue","Pmax_2nd","Test")
  } 
  rank_sep<-samprank_sep(genotype_matrix,cnv_matrix,cnvtype)
  output=matrix(unlist(c(chr,start,end,cnvID,sampleOrigIDs,cnvtypeOrigIDs,power,p[1],p[2],p[3],rank_sep[1],rank_sep[2])),nrow=1)
  return(output)
}

##Wrapper##
#Loads regions: chr start end locusID sampleID1,sampleID2,...
intervals <- read.table(opt$bed, sep = "\t", header = F)

#Make start and end numeric
intervals[, c(2:3)] <-
  apply(intervals[, c(2:3)], 2, function(vals) {
    return(as.numeric(as.character(vals)))
  })
intervals <- data.frame(lapply(intervals, as.character), stringsAsFactors=FALSE)
results<-alply(intervals,1,runRdTest)
if(class(results)=="list") {
  results<-do.call(rbind,results)
} else {
  results <- t(results)
}
#write ouputfile##
if (opt$denovo==FALSE) {
   if(!file.exists(paste(opt$outFolder,opt$outputname,".metrics",sep=""))) {
  ##write header##
  write.table(matrix(c("chr","Start","End","CNVID","SampleIDs","Type","Median_Power","P","2ndMaxP","Model","Median_Rank","Median_Separation"),nrow=1),paste(opt$outFolder, opt$outputname, ".metrics", sep = ""),quote=FALSE,row.names=FALSE,col.names=FALSE,sep= "\t")  
    } 
  write.table(results,paste(opt$outFolder, opt$outputname, ".metrics", sep = ""),
              quote = FALSE,col.names = FALSE, row.names = FALSE,append=TRUE,sep= "\t") 
} else {
  if(!file.exists(paste(opt$outFolder, opt$outputname,".denovo",sep=""))) {
    ##write header for de novo##
    write.table(matrix(c("chr","Start","End","CNVID","Type","Family","AffectedMember","Pro.P","Sib.P","Fa.P","Mo.P","Pro.secMaxP","Sib.secMaxP","Fa.secMaxP","Mo.secMaxP","Pro.Sep","Sib.Sep","Fa.Sep","Mo.Sep","Pro.rank","Sib.rank","Fa.rank","Mo.rank"),nrow=1),paste(opt$outFolder, opt$outputname, ".denovo", sep = ""),quote=FALSE,row.names=FALSE,col.names=FALSE,sep= "\t")  
  } 
  write.table(results,paste(opt$outFolder, opt$outputname, ".denovo", sep = ""),
              quote = FALSE,col.names = FALSE, row.names = FALSE,append=TRUE,sep= "\t") 
}

cat("FINISHED\n")
