# This script reformats the output bed files from module 01 with raw call evidence
# to match desired input for de novo SV filtering

# Define args and load libraries
args = commandArgs(trailingOnly=TRUE)
require(plyr)
library(data.table)
library(tidyverse)

input_bed <- args[1]
input_ped <- args[2]
out_probands <- args[3]
out_parents <- args[4]

# Read input bed
bed <- fread(input_bed)
ped <- fread(input_ped)

# adding back column names
colnames(bed) <- c("CHROM", "start", "end", "name", "svtype", "samples", "SVTYPE")

# Split into one sample per row
bed %>%
  subset(select = c(CHROM, start, end, SVTYPE, samples)) %>%
  separate_rows(samples, sep = ",", convert = T) -> bed_split

# split into parents and probands
bed_probands <- bed_split
bed_probands <- subset(bed_probands, !(samples %in% ped$MotherID) & !(samples %in% ped$FatherID))

bed_parents <- bed_split
bed_parents <- subset(bed_parents, samples %in% ped$MotherID | samples %in% ped$FatherID)

# add in a family ID column to bed_parents
ped_subset <- subset(ped, select = c('FamID', 'IndividualID'))
colnames(ped_subset) <- c('FamID', 'samples')

bed_parents <- merge(bed_parents,ped_subset, by="samples", all.x = TRUE)

# Reformat chromosome column
bed_probands$CHROM <- paste0(bed_probands$CHROM, "_", bed_probands$SVTYPE, "_", bed_probands$samples)
bed_parents$CHROM <- paste0(bed_parents$CHROM, "_", bed_parents$SVTYPE, "_", bed_parents$FamID)

# reorder the columns of bed_parents
bed_parents2 <- bed_parents[, c(2,3,4,5,1)]

# Write output file
write.table(bed_probands, out_probands, sep = "\t", quote = F, row.names = F, col.names = F)
write.table(bed_parents2, out_parents, sep = "\t", quote = F, row.names = F, col.names = F)
