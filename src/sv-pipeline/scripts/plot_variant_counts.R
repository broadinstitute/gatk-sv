#!/usr/bin/env Rscript

library(ggplot2)
library(gridExtra)

violin_single_group <- function(y, n=NA, log_scale=FALSE,
                                title="Violin Plot with Jittered Points",
                                ylab="Value") {
  y <- y[y > 0]
  data <- data.frame(group=factor("Group"), value=y)
  p <- ggplot(data, aes(x=group, y=value)) +
    geom_violin(trim=FALSE, fill="skyblue", color="black", alpha=0.7) +
    geom_jitter(width=0.1, size=1.5, alpha=0.6, color="darkblue") +
    labs(title=title, x=NULL, y=ylab) +
    theme_minimal(base_size=14) +
    theme(
      axis.title=element_text(size=18),
      axis.text=element_text(size=18),
      plot.title=element_text(size=20, face="bold", hjust=0.5),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank()
    )
  if (log_scale) p <- p + scale_y_log10()
  if (!is.na(n)) p <- p + geom_hline(yintercept=n, color="red",
                                     linetype="dashed", size=1)
  return(p)
}

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
  plot_list <- lapply(contigs_present, function(Z) {
    vals <- df[[paste(caller,Y,Z,sep="_")]]
    violin_single_group(vals, NA, FALSE,
                        title=Z, ylab=Z)
  })
  nplots <- length(plot_list)
  if (nplots > 0) {
    ncol   <- min(4, nplots)
    nrow   <- ceiling(nplots / ncol)
    grid   <- arrangeGrob(grobs=plot_list, nrow=nrow, ncol=ncol)
    out_file <- paste0(caller,"_",Y,".jpg")
    cell_w <- 4.0
    cell_h <- 2.5
    fig_w  <- ncol * cell_w
    fig_h  <- nrow * cell_h
    ggsave(out_file, grid,
           width  = fig_w, height = fig_h,
           units  = "in", dpi = 300)
    cat("Generated plot:", out_file, "\n")
  }
}

cat("Finished processing variant counts for caller:", caller, "\n") 