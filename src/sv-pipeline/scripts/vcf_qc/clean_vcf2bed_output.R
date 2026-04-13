#!/usr/bin/env Rscript

# Helper script to clean VCF stats output from collectQC.sh

###Set master parameters
options(stringsAsFactors=F)

###Load optparse
require(optparse)

###List of command-line options
option_list <- list(
  make_option(c("-N", "--nsamp"), type="integer", default=NULL,
              help="number of samples to be used for allele frequency calculations [default %default]",
              metavar="integer"),
  make_option(c("-M", "--noLabelMCNV"), type="logical", action="store_false", default=TRUE,
              help="disable relabeling any multiallelic DEL or DUP records as MCNV [default: relabel MCNVs]",
              metavar="logical"),
  make_option(c("-R", "--noReclassify"), type="logical", action="store_false", default=TRUE,
              help="disable reclassifying any SV not matching an input SV class as 'OTH' [default: reclassify OTH]",
              metavar="logical"),
  make_option(c("--keepBNDsize"), type="logical", action="store_true", default=FALSE,
              help="keep BND size information [default: overwrite size as NA]",
              metavar="logical")
)

###Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog vcf2bed_cleaned.bed genotype_counts_per_SV.txt SVTYPES OUTFILE",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

###Checks for appropriate positional arguments
if(length(args$args) != 4){
  stop("Incorrect number of required positional arguments\n")
}

###Writes args & opts to vars
INFILE <- args$args[1]
GENOTYPES <- args$args[2]
svtypes.file <- args$args[3]
OUTFILE <- args$args[4]
nsamp <- opts$nsamp
labelMCNV <- opts$noLabelMCNV
reclassify <- opts$noReclassify
keepBNDsize <- opts$keepBNDsize

###Reads & cleans data
#Read data
dat <- read.table(INFILE,comment.char="",sep="\t",header=T,check.names=F)
#Fix leading column name
colnames(dat)[1] <- "chr"
drops <- c('AF','AN','AC')
dat=dat[,!colnames(dat)%in%drops]
#Read genotype counts
gt <- read.table(GENOTYPES,sep="\t",comment.char="",header=T,check.names=F)
#Merge data & genotypes
dat <- merge(dat,gt,sort=F,
                by.x=which(colnames(dat)=="name"),
                by.y=which(colnames(gt)=="VID"))
#Sets sv types
if(!is.null(svtypes.file)){
  svtypes <- read.table(svtypes.file,sep="\t",header=F,comment.char="",check.names=F)[,1]
}else{
  svtypes <- unique(dat$svtype)
}

###Processes VCF metrics
#Compute carrier frequency
if(is.null(nsamp)){
  nsamp <- length(unique(unlist(strsplit(as.character(dat$samples),split=","))))
}
dat$observations <- sapply(dat$samples,function(samples){
  length(unique(unlist(strsplit(as.character(samples),split=","))))
})
if (nrow(dat) > 0) {
  dat$carrierFreq <- dat$observations/nsamp
} else {
  dat$carrierFreq <- numeric()
}
if(any(dat$frequency > 1)){
  stop("Incorrect number of samples supplied; some carrier frequencies > 1")
}
#Subset to relevant columns
# Columns from bed (svtk vcf2bed -i ALL captures all INFO fields)
region.col <- if("REGION" %in% colnames(dat)) dat$REGION else NULL
ncr.col <- if("NCR" %in% colnames(dat)) dat$NCR else NULL
gmatch.col <- if("gnomAD_V4_match_ID" %in% colnames(dat)) dat$gnomAD_V4_match_ID else NULL
# Columns from genotype counts merge (QUAL/REF/ALT added to that file)
qual.col <- if("QUAL" %in% colnames(dat)) dat$QUAL else NULL
ref.col <- if("REF" %in% colnames(dat)) dat$REF else NULL
alt.col <- if("ALT" %in% colnames(dat)) dat$ALT else NULL
# Max motif length from MOTIFS INFO field (for TRV stratification)
motifs.raw <- if("MOTIFS" %in% colnames(dat)) as.character(dat$MOTIFS) else NULL
max.motif.len <- if(!is.null(motifs.raw)) sapply(motifs.raw, function(m){
  if(is.na(m) || m == ".") return(NA_integer_)
  mots <- strsplit(m, ",", fixed=TRUE)[[1]]
  max(nchar(mots), na.rm=TRUE)
}) else NULL
# TR_ENVELOPED flag (TRUE for SV sites fully enveloped by a TR/VNTR locus)
tr.env.col <- if("TR_ENVELOPED" %in% colnames(dat)) {
  v <- as.character(dat$TR_ENVELOPED)
  !is.na(v) & toupper(v) %in% c("TRUE","1","YES")
} else NULL
# TRV expansion ratio
trv.exp.col <- if("TRV_EXPANSION_RATIO" %in% colnames(dat)) as.numeric(as.character(dat$TRV_EXPANSION_RATIO)) else NULL
# VEP Consequence field (pipe-field index 1 per comma-chained annotation, semicolon-joined unique values)
vep.col <- if("vep" %in% colnames(dat)) {
  sapply(as.character(dat$vep), function(v){
    if(is.na(v) || v == "." || nchar(v) == 0) return(NA_character_)
    annots <- strsplit(v, ",", fixed=TRUE)[[1]]
    csqs <- unique(vapply(annots, function(a){
      parts <- strsplit(a, "|", fixed=TRUE)[[1]]
      if(length(parts) >= 2 && nchar(parts[2]) > 0) parts[2] else NA_character_
    }, character(1)))
    csqs <- csqs[!is.na(csqs)]
    if(length(csqs) == 0) return(NA_character_)
    paste(sort(csqs), collapse=";")
  })
} else NULL
# PREDICTED_* fields from SVAnnotate (TRUE if present and non-empty)
predicted.names <- grep("^PREDICTED_", colnames(dat), value=TRUE)
predicted.presence <- setNames(lapply(predicted.names, function(pname){
  v <- as.character(dat[[pname]])
  !is.na(v) & v != "." & v != "" & toupper(v) != "FALSE" & v != "0"
}), predicted.names)
dat <- data.frame("chr"=dat$chr,
                  "start"=dat$start,
                  "end"=dat$end,
                  "VID"=dat$name,
                  "svtype"=dat$SVTYPE,
                  "length"=dat$SVLEN,
                  "genotyped_samples"=dat$nsamp_gt,
                  "reference_gts"=dat$homref,
                  "het_gts"=dat$het,
                  "hom_gts"=dat$homalt,
                  "other_gts"=dat$other,
                  "missing_gts"=dat$unknown,
                  "carriers"=dat$observations,
                  "carrierFreq"=dat$carrierFreq,
                  "AN"=dat$AN,
                  "AC"=dat$AC,
                  "AF"=dat$AF)
if(!is.null(region.col)) dat$REGION <- region.col
if(!is.null(qual.col)) dat$QUAL <- as.numeric(as.character(qual.col))
if(!is.null(ref.col)) dat$REF <- as.character(ref.col)
if(!is.null(alt.col)) dat$ALT <- as.character(alt.col)
if(!is.null(ncr.col)) { ncr.v <- as.character(ncr.col); ncr.v[ncr.v=="."] <- NA; dat$NCR <- as.numeric(ncr.v) }
if(!is.null(gmatch.col)) dat$gnomAD_V4_match_ID <- as.character(gmatch.col)
if(!is.null(vep.col)) dat$VEP_consequences <- vep.col
for(pname in names(predicted.presence)) dat[[pname]] <- predicted.presence[[pname]]
if(!is.null(max.motif.len)) dat$max_motif_length <- as.integer(max.motif.len)
if(!is.null(trv.exp.col)) dat$TRV_EXPANSION_RATIO <- trv.exp.col
if(!is.null(tr.env.col)) dat$TR_ENVELOPED <- tr.env.col
zeroes <- which(dat$AC==0 | dat$carriers==0)
if(length(zeroes)>0){
  cat(paste("WARNING: ",prettyNum(length(zeroes),big.mark=","),"/",
              prettyNum(nrow(dat),big.mark=",")," (",
              round(100*length(zeroes)/nrow(dat),1),"%) of SV sites have no observed carriers or alleles.",
              " Excluding these sites from all downstream analyses.\n",sep=""))
  dat <- dat[-zeroes,]
}
#Relabel multiallelic DELs/DUPs and MCNVs, if optioned
if(labelMCNV==T){
  mCNV.to.label <- which(dat$svtype %in% c("DEL","DEL_SHORT","DEL_SV","DUP") & dat$other_gts>0)
  if(length(mCNV.to.label)>0){
    cat(paste("WARNING: ",prettyNum(length(mCNV.to.label),big.mark=","),"/",
              prettyNum(nrow(dat),big.mark=",")," (",
              round(100*length(zeroes)/nrow(dat),1),"%) of SV sites have non-Mendelian genotypes, but",
              " are not marked as multiallelic CNVs. Relableing these sites as 'MCNV'.\n",sep=""))
    
    dat[mCNV.to.label,]$svtype <- "MCNV"
  }
}
#Split TRV into TRV (small) and TRV_SV (large) by size threshold
trv.idx <- which(dat$svtype == "TRV" & !is.na(dat$length) & dat$length > 50)
if(length(trv.idx) > 0) dat$svtype[trv.idx] <- "TRV_SV"

#Set length of all BNDs to NA, if optioned
if(keepBNDsize==F){
  if(any(dat$svtype=="BND")){
    dat[which(dat$svtype=="BND"),]$length <- NA
  }
}
#Force all non-specified SV types to 'OTH', if optioned
if(reclassify==T){
  OTH.to.reclassify <- which(!(dat$svtype %in% svtypes))
  if(length(OTH.to.reclassify)>0){
    cat(paste("WARNING: ",prettyNum(length(OTH.to.reclassify),big.mark=","),"/",
              prettyNum(nrow(dat),big.mark=",")," (",
              round(100*length(zeroes)/nrow(dat),1),"%) of SV sites do not match a specified",
              " SV type. Relableing these sites as 'OTH'.\n",sep=""))
    dat[OTH.to.reclassify,]$svtype <- "OTH"
  }
}

#Write out results
write.table(dat,OUTFILE,sep="\t",quote=F,col.names=T,row.names=F)
