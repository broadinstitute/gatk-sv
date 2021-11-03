#!/usr/bin/env Rscript

# Talkowski SV pipeline downstream analysis helper script

# Identify variants with unstable frequencies before/after minGQ


### Set global parameters
options(stringsAsFactors=F, scipen=10)


### Import pre-minGQ frequency data

load.preMinGQ.freqs.PCRMINUS_only <- function(pre.freqs.PCRMINUS.in){
  minus.dat <- read.table(pre.freqs.PCRMINUS.in, header=F, sep="\t")
  colnames(minus.dat) <- c("VID", "AC.PCRMINUS", "AN.PCRMINUS")
  return(minus.dat)
  }


### Import batch effect-derived frequency data 
load.batcheffect.freqs <- function(freq.table.in){
  # Read data
  dat <- read.table(freq.table.in, header=T, sep="\t", comment.char="")
  colnames(dat)[1] <- "VID"
  # Drop multiallelics
  mcnv.idxs <- unique(unlist(apply(dat[, -c(1:3)], 2, function(vals){grep(",", vals, fixed=T)})))
  if(length(mcnv.idxs) > 0){
    dat <- dat[-mcnv.idxs, ]
  }
  # Convert all rows to numerics
  dat[, -c(1:3)] <- apply(dat[, -c(1:3)], 2, as.numeric)
  # Sum frequencies
  #PCRPLUS_AN.idxs <- intersect(grep("PCRPLUS", colnames(dat), fixed=T), grep("AN.", colnames(dat), fixed=T))
  #PCRPLUS_AN <- apply(dat[, PCRPLUS_AN.idxs], 1, sum)
  #PCRPLUS_AC.idxs <- intersect(grep("PCRPLUS", colnames(dat), fixed=T), grep("AC.", colnames(dat), fixed=T))
  #PCRPLUS_AC <- apply(dat[, PCRPLUS_AC.idxs], 1, sum)
  PCRMINUS_AN.idxs <-  grep("AN.", colnames(dat), fixed=T)
  PCRMINUS_AN <- apply(dat[, PCRMINUS_AN.idxs], 1, sum)
  PCRMINUS_AC.idxs <- grep("AC.", colnames(dat), fixed=T)
  PCRMINUS_AC <- apply(dat[, PCRMINUS_AC.idxs], 1, sum)
  # Prep cleaned data
  out.df <- data.frame("VID"=dat$VID,
                       "AC.PCRMINUS"=PCRMINUS_AC, "AN.PCRMINUS"=PCRMINUS_AN)
  return(out.df)
}


### Merge pre & post minGQ data, and convert ACs of all multiallelic sites
merge.clean.freqs <- function(pre.dat, post.dat){
  # Merge data
  merged.dat <- merge(pre.dat, post.dat, by="VID", sort=F, all.y=T, all.x=F, 
                      suffixes=c(".pre", ".post"))
  # Clean multiallelics
  # merged.dat[, -1] <- apply(merged.dat[, -1], 2, function(vals){
  #   midx <- grep(",", vals, fixed=T)
  #   vals[midx] <- sapply(midx, function(i){
  #     sum(as.numeric(unlist(strsplit(vals[i], split=",")))[-3])
  #   })
  #   vals <- as.numeric(vals)
  #   return(vals)
  # })
  mcnv.idxs <- grep(",", merged.dat$AC.PCRMINUS.pre, fixed=T)
  if(length(mcnv.idxs) > 0){
    merged.dat <- merged.dat[-mcnv.idxs, ]
  }
  # Convert all stats to numeric
  merged.dat[, -1] <- apply(merged.dat[, -1], 2, as.numeric)
  # Drop any records where any value is NA
  merged.dat <- merged.dat[which(complete.cases(merged.dat)), ]
  return(merged.dat)
}


### Run chi-squared test for all variants for a single PCR status
calc.pvals <- function(dat, PCR){
  # Get counts
  pre.AC <- dat[, which(colnames(dat) == paste("AC", PCR, "pre", sep="."))]
  pre.AN <- dat[, which(colnames(dat) == paste("AN", PCR, "pre", sep="."))]
  pre.ref <- pre.AN - pre.AC
  post.AC <- dat[, which(colnames(dat) == paste("AC", PCR, "post", sep="."))]
  post.AN <- dat[, which(colnames(dat) == paste("AN", PCR, "post", sep="."))]
  post.ref <- post.AN - post.AC
  # Run chisq tests
  sapply(1:nrow(dat), function(i){
    chisq.test(matrix(c(pre.ref[i], post.ref[i],
                        pre.AC[i], post.AC[i]),
                      nrow=2, byrow=T))$p.value
  })
}


### Estimate fraction of null GTs introduced by minGQ
estimate.null.gts <- function(dat, PCR){
  # Get counts
  pre.AN <- dat[, which(colnames(dat) == paste("AN", PCR, "pre", sep="."))]
  post.AN <- dat[, which(colnames(dat) == paste("AN", PCR, "post", sep="."))]
  # Compute normalized change in AN (mode takes into account excluded outlier samples)
  AN.delta <- pre.AN - post.AN
  delta.mode <- as.numeric(names(sort(-table(AN.delta))[1]))
  AN.delta.norm <- (AN.delta - delta.mode) / pre.AN
  return(AN.delta.norm)
}


### Plot frequencies before & after minGQ
plot.freqs <- function(dat, PCR){
  # Get data
  x.AC <- dat[, which(colnames(dat) == paste("AC", PCR, "pre", sep="."))]
  x.AN <- dat[, which(colnames(dat) == paste("AN", PCR, "pre", sep="."))]
  x <- x.AC / x.AN
  y.AC <- dat[, which(colnames(dat) == paste("AC", PCR, "post", sep="."))]
  y.AN <- dat[, which(colnames(dat) == paste("AN", PCR, "post", sep="."))]
  y <- y.AC / y.AN
  fail <- dat[, which(colnames(dat) == paste(PCR, "fail", sep="."))]
  # Prep plot
  par(mar=c(3, 3, 1.5, 0.5))
  plot(x=c(0, 1), y=c(0, 1), type="n",
       xaxt="n", yaxt="n", xlab="", ylab="")
  axis(1, at=seq(0, 1, 0.2), tick=F, line=-0.9)
  mtext(1, line=1.5, text="AF before minGQ")
  axis(2, at=seq(0, 1, 0.2), tick=F, line=-0.9, las=2)
  mtext(2, line=1.5, text="AF after minGQ")
  mtext(3, line=0.25, text=PCR, font=2)
  abline(h=seq(0, 1, 0.2), v=seq(0, 1, 0.2), col="gray90")
  # Add points
  pt.col <- rep("black", nrow(dat))
  pt.bg <- rep(NA, nrow(dat))
  pt.col[which(fail)] <- "firebrick"
  pt.bg[which(fail)] <- "red"
  points(x=x, y=y, cex=0.1, lwd=0.2, pch=21, col=pt.col, bg=pt.bg)
}


#################
### RScript block
#################

### Read command-line arguments
args <- commandArgs(trailingOnly=T)
#pre.freqs.PCRPLUS.in <- as.character(args[1])
pre.freqs.PCRMINUS.in <- as.character(args[1])
freq.table.in <- as.character(args[2])
OUTDIR <- as.character(args[3])
prefix <- as.character(args[4])

# #DEV:
# # pre.freqs.PCRPLUS.in <- "~/scratch/gnomAD_v2_SV_MASTER.PCRPLUS.AF_preMinGQ.txt"
# pre.freqs.PCRMINUS.in <- "~/scratch/minGQ_test/PD_1perc_fdr.PCRMINUS.AF_preMinGQ.txt"
# freq.table.in <- "~/scratch/minGQ_test/mod07_mcnv_test.merged_AF_table.txt.gz"
# OUTDIR <- "~/scratch/minGQ_test/"
# prefix <- "mod07_mcnv_test"

# Read frequency data
pre.dat <- load.preMinGQ.freqs.PCRMINUS_only(pre.freqs.PCRMINUS.in)
post.dat <- load.batcheffect.freqs(freq.table.in)

# Merge & clean frequency data
merged.dat <- merge.clean.freqs(pre.dat, post.dat)

# Compute p-values
#p.PCRPLUS <- suppressWarnings(calc.pvals(merged.dat, "PCRPLUS"))
p.PCRMINUS <- suppressWarnings(calc.pvals(merged.dat, "PCRMINUS"))

# Estimate pct of null-GT samples
#PCRPLUS.nullGTs <- estimate.null.gts(merged.dat, "PCRPLUS")
PCRMINUS.nullGTs <- estimate.null.gts(merged.dat, "PCRMINUS")

# Write lists of failing variants
PCRMINUS.fail.idx <- which(p.PCRMINUS < 0.05/length(p.PCRMINUS)
                           & PCRMINUS.nullGTs >= 0.02)
write.table(merged.dat$VID[PCRMINUS.fail.idx],
            paste(OUTDIR, "/", prefix, "PCRMINUS_minGQ_AF_prePost_fails.VIDs.list", sep=""),
            col.names=F, row.names=F, quote=F)

#generate an empty list for PCRPLUS
#PCRPLUS.fail.idx <- which(p.PCRMINUS < -0.05/length(p.PCRMINUS))
#write.table(merged.dat$VID[PCRPLUS.fail.idx],
#           paste(OUTDIR, "/", prefix, "PCRPLUS_minGQ_AF_prePost_fails.VIDs.list", sep=""),
#            col.names=F, row.names=F, quote=F)


# Write table of all data, for reference
#merged.dat$p.PCRPLUS <- p.PCRPLUS
#merged.dat$PCRPLUS.nullGTs <- PCRPLUS.nullGTs
#merged.dat$PCRPLUS.fail <- FALSE
#merged.dat$PCRPLUS.fail[PCRPLUS.fail.idx] <- TRUE
merged.dat$p.PCRMINUS <- p.PCRMINUS
merged.dat$PCRMINUS.nullGTs <- PCRMINUS.nullGTs
merged.dat$PCRMINUS.fail <- FALSE
merged.dat$PCRMINUS.fail[PCRMINUS.fail.idx] <- TRUE
colnames(merged.dat)[1] <- "#VID"
out.data.path <- paste(OUTDIR, "/", prefix, "minGQ_AF_prePost_comparison.data.txt", sep="")
write.table(merged.dat, out.data.path,
            col.names=T, row.names=F, quote=F, sep="\t")
system(paste("gzip -f ", out.data.path), wait=T)

# Generate plots of AF before/after with fail labels
png(paste(OUTDIR, "/", prefix, "minGQ_AF_prePost_comparison.plot.png", sep=""),
    height=500, width=1000)
par(mfrow=c(1, 2))
#plot.freqs(merged.dat, "PCRPLUS")
plot.freqs(merged.dat, "PCRMINUS")
dev.off()

