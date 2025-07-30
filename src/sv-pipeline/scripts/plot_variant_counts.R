#!/usr/bin/env Rscript

library(vioplot)

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
  stop("Usage: plot_variant_counts.R <variant_count_file> <caller>")
}

variant_count_file <- args[1]
caller <- args[2]

if (!file.exists(variant_count_file)) {
  stop(paste("File does not exist:", variant_count_file))
}

df <- read.table(variant_count_file, header=TRUE, sep="\t", 
                 check.names=FALSE, stringsAsFactors=FALSE)

y_vals <- c("BND","DEL","DUP","INS","INV")
all_contigs <- c(paste0("chr",1:22),"chrX","chrY")

for (Y in y_vals) {
  present_cols <- grep(paste0("^",caller,"_",Y,"_"), names(df), value=TRUE)
  if (length(present_cols)==0) next
  contigs_present <- sub(paste0("^",caller,"_",Y,"_"), "", present_cols)
  contigs_present <- intersect(all_contigs, contigs_present)
  
  if (length(contigs_present) > 0) {
    ncol <- min(4, length(contigs_present))
    nrow <- ceiling(length(contigs_present) / ncol)
    
    cell_w <- 3.0
    cell_h <- 2.0
    fig_w <- ncol * cell_w
    fig_h <- nrow * cell_h
    
    out_file <- paste0(caller,"_",Y,".png")
    
    png(out_file, width=fig_w, height=fig_h, units="in", res=150)
    
    par(mfrow=c(nrow, ncol), mar=c(4,4,2,1))
    
    for (i in seq_along(contigs_present)) {
      contig <- contigs_present[i]
      vals <- df[[paste(caller,Y,contig,sep="_")]]
      vals <- vals[vals > 0]
      
      if (length(vals) > 0) {
        vioplot(vals, main=contig, col="skyblue", 
                cex.main=1.2, border="black", xaxt="n")
        points(jitter(rep(1, length(vals)), amount=0.1), vals, 
                col="darkblue", pch=16, cex=0.8)
      } else {
        plot.new()
        title(main=contig, cex.main=1.2)
        text(0.5, 0.5, "No data", cex=1.5, col="gray")
      }
    }
    
    dev.off()
    cat("Generated plot:", out_file, "\n")
  }
}

cat("Finished processing variant counts for caller:", caller, "\n") 