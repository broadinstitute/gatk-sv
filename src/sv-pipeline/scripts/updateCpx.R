#!/usr/bin/env Rscript

# Load libraries
library(tidyr)
library(data.table)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
    stop("Usage: updateCpx.R <input_bed>")
}

input_bed <- args[1]

# Read data
denovo <- fread(input_bed)
out_path <- "final.denovo.merged.cpx_split.bed"

# Split dataframe into complex and non-complex calls
denovo_non_cpx <- subset(denovo, svtype != "CPX", select = c("chrom", "start", "end", "name", "svtype", "sample"))
denovo_cpx <- subset(denovo, svtype == "CPX")

# Check if there are CPX variants to process
if (nrow(denovo_cpx) > 0 && "CPX_INTERVALS" %in% colnames(denovo_cpx)) {
    # Reformat cpx calls
    denovo_cpx %>% 
        separate_rows(CPX_INTERVALS, sep = ",") %>%
        select(c("chrom", "start", "end", "name", "svtype", "sample", "CPX_INTERVALS")) %>%
        separate(CPX_INTERVALS, sep = "_", into = c("cpx_type", "cpx_coord")) %>%
        separate(cpx_coord, sep = ":", into = c("cpx_chr", "cpx_pos")) %>%
        separate(cpx_pos, sep = "-", into = c("cpx_start", "cpx_end")) -> denovo_cpx
    
    # Update name of cpx calls and change column names
    denovo_cpx$name <- paste0(denovo_cpx$name, "_", denovo_cpx$cpx_type)
    denovo_cpx <- subset(denovo_cpx, select = c("cpx_chr", "cpx_start", "cpx_end", "name", "cpx_type", "sample"))
    names(denovo_cpx) <- names(denovo_non_cpx)
    
    # Merge all de novo calls
    denovo_split <- rbind(denovo_cpx, denovo_non_cpx)
} else {
    # No CPX variants to split, just use non-CPX
    denovo_split <- denovo_non_cpx
}

# Write output
write.table(denovo_split, out_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE) 