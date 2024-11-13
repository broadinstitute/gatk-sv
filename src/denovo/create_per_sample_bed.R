require(plyr)
library(data.table)
library(tidyverse)


args = commandArgs(trailingOnly=TRUE)
input_gd <- args[1]
proband_gd <- args[2]
parents_gd <- args[3]
ped <- args[4]
chromosome <- args[5]

df_ori <- as.data.frame(read.table(input_gd, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
colnames(df_ori) <- c("chrom", "start", "end", "name", "annotation", "svtype")
df <- subset(df_ori, chrom == chromosome)

df_DUP <- subset(df, svtype == 'DUP' | svtype == 'DEL;DUP')
if (nrow(df_DUP) > 0) {
  df_DUP$svtype <- 'DUP'
}

df_DEL <- subset(df, svtype == 'DEL' | svtype == 'DEL;DUP')
if (nrow(df_DEL) > 0) {
  df_DEL$svtype <- 'DEL'
}

if (nrow(df_DEL) > 0 & nrow(df_DUP) > 0) {
  df2 = rbind(df_DUP,df_DEL)
} else if (nrow(df_DEL) > 0) {
  df2 = df_DEL
} else if (nrow(df_DUP) > 0) {
  df2 = df_DUP
}

if (exists('df2')) {
  ped <- as.data.frame(read.table(ped, header = TRUE, sep="\t",stringsAsFactors=FALSE, quote=""))
  samples <- ped$IndividualID
  df2$samples <- gsub(" ", "",toString(samples))

  df2 %>%
    subset(select = c(chrom, start, end, name, svtype, samples)) %>%
    separate_rows(samples, sep = ",", convert = T) -> df3

  df3$chrom <- paste0(df3$chrom, "_", df3$svtype, "_", df3$samples)
  df3 <- df3[, c(1,2,3,4)]

  ped_subset <- subset(ped, select = c('FamID', 'IndividualID'))
  colnames(ped_subset) <- c('FamID', 'samples')

  df2 %>%
    subset(select = c(chrom, start, end, name, svtype, samples)) %>%
    separate_rows(samples, sep = ",", convert = T) -> new_df

  new_df <- merge(new_df,ped_subset, by="samples", all.x = TRUE)
  new_df$chrom <- paste0(new_df$chrom, "_", new_df$svtype, "_", new_df$FamID)
  new_df <- new_df[, c(2,3,4,5)]

  #df4 = rbind(df3,new_df)

  write.table(df3, proband_gd, sep = "\t", quote = F, row.names = F, col.names = F)
  write.table(new_df, parents_gd, sep = "\t", quote = F, row.names = F, col.names = F)
} else{
  df3 <- data.frame()
  write.table(df3, proband_gd, sep = "\t", quote = F, row.names = F, col.names = F)
  write.table(df3, parents_gd, sep = "\t", quote = F, row.names = F, col.names = F)
}