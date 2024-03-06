# Description: Post-hoc QC analysis of filtered de novo SVs


# Load libraries
library(tidyr)
library(ggplot2)
library(UpSetR)
library(data.table)
library(dplyr)
library(gridExtra)
library(grid)
library(tidyverse)
library("grid")
library("ggplotify")
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

input_bed <- args[1]
input_ped <- args[2]
out_file <- args[3]

#input_outliers <-args [2]
#input_ped <- args[3]
#out_file <- args[4]

ped <- fread(input_ped)

denovo <- as.data.frame(read.table(input_bed,header = TRUE, sep="\t",stringsAsFactors=FALSE, quote=""))
#outliers <- as.data.frame(read.table(input_outliers,header = TRUE, sep="\t",stringsAsFactors=FALSE, quote=""))

denovo$SVTYPE <- factor(denovo$SVTYPE, levels = rev(c("DEL", "DUP", "INS", "INV", "CPX", "CTX")))

denovo_ins <- subset(denovo, SVTYPE == "INS")
denovo_del <- subset(denovo, SVTYPE == "DEL")
denovo_dup <- subset(denovo, SVTYPE == "DUP")
denovo_inv <- subset(denovo, SVTYPE == "INV")
denovo_cpx <- subset(denovo, SVTYPE == "CPX")
denovo_ctx <- subset(denovo, SVTYPE == "CTX")

# Define colors
del_col <- "#D43925"
dup_col <- "#2376B2"
ins_col <- "#D474E0"
inv_col <- "#FA9627"
ctx_col <- "#638E6C"
cpx_col <- "#4DA1A9"

# Get some numbers
length(unique(denovo$name)) #of unique denovos
length(unique(denovo$sample)) #number of samples with denovos

#length(unique(outliers$name))
#length(unique(outliers$sample))
#outliers_in_gd <- subset(outliers, select= c(name, in_gd))

denovo %>%
  # subset(chrom != "chrX") %>%
  select(sample, name) %>%
  group_by(sample) %>%
  tally() -> sample_count

denovo %>% 
  count(sample, SVTYPE, name= "svtype_per_sample") -> type

type_boxplot <- ggplot(type, aes(x=SVTYPE, fill = SVTYPE, y=svtype_per_sample)) + 
  geom_jitter(position = position_jitter(seed = 1, width = 0.2), color = 'grey', size = 4.0) + geom_boxplot(outlier.shape=NA) + labs(title = "Number of De Novo SVs per Type", y = "Number of de novo SVs", x = "SV Type") + scale_fill_manual(values = rev(c("DEL"=del_col, "DUP"=dup_col, "INS"=ins_col, "INV"=inv_col, "CPX"=cpx_col, "CTX"=ctx_col))) + 
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 70),
    axis.title = element_text(size = 70),
    plot.title = element_text(size=70),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.title = element_text(size=30),
    legend.text = element_text(size=25))
ggsave("type_boxplot.png", type_boxplot, width = 30, height = 40, limitsize = FALSE) 

sample_boxplot <- ggplot(sample_count, aes(x=factor(0), y=n)) + 
  geom_jitter(position = position_jitter(seed = 1, width = 0.2), color = 'gray47', size = 4.0) + geom_boxplot(outlier.shape=NA) + labs(title = "Number of de Novo SVs per Sample", y = "Number of de novo SVs", x = "Samples") + scale_fill_manual(values = rev(c("DEL"=del_col, "DUP"=dup_col, "INS"=ins_col, "INV"=inv_col, "CPX"=cpx_col, "CTX"=ctx_col))) + 
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 70),
    axis.title = element_text(size = 70),
    plot.title = element_text(size=70),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.title = element_text(size=30),
    legend.text = element_text(size=25))
ggsave("sample_boxplot.png", sample_boxplot, width = 30, height = 35, limitsize = FALSE) 

sample_count %>%
  ggplot(aes(y = n)) +
  geom_boxplot() -> boxplot

# Make plots
denovo$chrom <- factor(denovo$chrom, levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"))

# Count by SV Type
denovo %>%
  select(SVTYPE, name, SVTYPE) %>%
  unique() %>%
  ggplot(aes(x = SVTYPE, fill = SVTYPE)) +
  geom_bar() +
  labs(x="SV Type", y="Number de novo SVs", title= "Number of De Novo SVs per SV Type") +
  scale_fill_manual(values = rev(c("DEL"=del_col, "DUP"=dup_col, "INS"=ins_col, "INV"=inv_col, "CPX"=cpx_col, "CTX"=ctx_col))) + 
  theme_classic() +
  scale_x_discrete(limits = rev(levels(denovo$SVTYPE))) +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 70),
    axis.title = element_text(size = 70),
    plot.title = element_text(size=70),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.title = element_text(size=30),
    legend.text = element_text(size=25)) -> p_type_count
p_type_count

ggsave("per_type.png", p_type_count, width = 30, height = 35, limitsize = FALSE)

# Count by chromosome
denovo %>%
  select(chrom, name, SVTYPE) %>%
  unique() %>%
  ggplot(aes(x = chrom, fill = SVTYPE)) +
  geom_bar() +
  labs(x="Chromosome", y="Number de novo SVs", title= "Number of De Novo SVs per Chromosome") +
  scale_fill_manual(values = rev(c("DEL"=del_col, "DUP"=dup_col, "INS"=ins_col, "INV"=inv_col, "CPX"=cpx_col, "CTX"=ctx_col))) + 
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 70),
    axis.title = element_text(size = 70),
    plot.title = element_text(size=70),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.title = element_text(size=30),
    legend.text = element_text(size=25)) -> p_chr_count
p_chr_count

ggsave("per_chrom.png", p_chr_count, width = 30, height = 45, limitsize = FALSE)

# Count by sample
df <- as.data.frame(table(denovo$sample))
df$samples <- as.character(df$Var1)
df$num_of_denovos <- as.numeric(as.character(df$Freq))
df$samples <- as.factor(df$samples)

p <- ggplot(df, aes(y=num_of_denovos)) + geom_boxplot() + labs(title = "Number of de novo SVs per sample", y = "Number of de novos", x = "Samples") + theme_classic() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 70),
    axis.title = element_text(size = 70),
    plot.title = element_text(size=70),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.title = element_text(size=30),
    legend.text = element_text(size=25))
ggsave("per_sample.png", p, width = 30, height = 35, limitsize = FALSE)

# De novo count by allele frequency
denovo$AF <- as.numeric(denovo$AF)  
denovo$AC <- as.numeric(denovo$AC) 
ac_1_freq <- min(subset(denovo, AC == 1)$AF)
max_AF_str <- toString(round(max(denovo$AF), digits=3))
# denovo$AF_interv <- cut(denovo$AF, breaks = c(0, ac_1_freq, 0.001, 0.01, 0.03, max(denovo$AF)),
#                         labels = c("AC=1", "AC 1 - AF<=0.001", "AF 0.001-0.01", "AF 0.01-0.03", "AF>0.03"))
denovo$AF_interv <- cut(denovo$AF, breaks = c(0, ac_1_freq, 0.001, max(denovo$AF)),
                        labels = c("AC=1", "AC 1-0.001", "AF>0.001"))

denovo %>%
  select(AF_interv, name, SVTYPE) %>%
  unique() %>%
  ggplot(aes(x = AF_interv, fill = SVTYPE)) +
  geom_bar() +
  # scale_y_log10()+
  labs(x="Allele Frequency", y="Number de novo SVs", title="Number of De Novo\n SVs by Frequency") +
  scale_fill_manual(values = rev(c("DEL"=del_col, "DUP"=dup_col, "INS"=ins_col, "INV"=inv_col, "CPX"=cpx_col, "CTX"=ctx_col))) + 
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 70),
    axis.title = element_text(size = 70),
    plot.title = element_text(size=70),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.title = element_text(size=30),
    legend.text = element_text(size=25)) -> p_af_count
ggsave("per_freq.png", p_af_count, width = 30, height = 35, limitsize = FALSE)

denovo_in_gd <- subset(denovo, in_gd == "True")
denovo_not_in_gd <- subset(denovo, in_gd == "False")

denovo_in_gd %>%
  select(AF_interv, name, SVTYPE) %>%
  unique() %>%
  ggplot(aes(x = AF_interv, fill = SVTYPE)) +
  geom_bar() +
  # scale_y_log10()+
  labs(x="Allele Frequency", y="Number de novo SVs", title= "Number of De Novo SVs \nin Genomic Disorder Regions by Frequency") +
  scale_fill_manual(values = rev(c("DEL"=del_col, "DUP"=dup_col, "INS"=ins_col, "INV"=inv_col, "CPX"=cpx_col, "CTX"=ctx_col))) + 
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 70),
    axis.title = element_text(size = 70),
    plot.title = element_text(size=70),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.title = element_text(size=30),
    legend.text = element_text(size=25)) -> p_af_count_in_gd
ggsave("per_freq_gd.png", p_af_count_in_gd, width = 30, height = 35, limitsize = FALSE)

denovo_not_in_gd %>%
  select(AF_interv, name, SVTYPE) %>%
  unique() %>%
  ggplot(aes(x = AF_interv, fill = SVTYPE)) +
  geom_bar() +
  # scale_y_log10()+
  labs(x="Allele Frequency", y="Number de novo SVs", title= "Number of De Novo SVs not in \n Genomic Disorder Regions by Frequency") +
  scale_fill_manual(values = rev(c("DEL"=del_col, "DUP"=dup_col, "INS"=ins_col, "INV"=inv_col, "CPX"=cpx_col, "CTX"=ctx_col))) + 
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 70),
    axis.title = element_text(size = 70),
    plot.title = element_text(size=70),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.title = element_text(size=30),
    legend.text = element_text(size=25)) -> p_af_count_not_in_gd
ggsave("per_freq_not_gd.png", p_af_count_not_in_gd, width = 30, height = 35, limitsize = FALSE) 

# De novo count by size
denovo$SVLEN <- as.numeric(denovo$SVLEN)  
#denovo[denovo$SVLEN < 1,]$SVLEN <- 1
denovo$SVLEN_interv <- cut(denovo$SVLEN, breaks = c(0, 100, 500, 2500, 10000, 50000, max(denovo$SVLEN)),
                           labels = c("<100bp", "101-500bp", "501bp-2.5kb", "2.5kb-10kb", "10kb-50kb", ">50kb"))
denovo %>%
  select(SVLEN_interv, name, SVTYPE) %>%
  unique() %>%
  ggplot(aes(x = SVLEN_interv, fill = SVTYPE)) +
  geom_bar() +
  # scale_y_log10()+
  labs(x="Size of SV", y="Number de novo SVs", title= "De Novo SVs by Size") +
  scale_fill_manual(values = rev(c("DEL"=del_col, "DUP"=dup_col, "INS"=ins_col, "INV"=inv_col, "CPX"=cpx_col, "CTX"=ctx_col))) + 
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 70),
    axis.title = element_text(size = 70),
    plot.title = element_text(size=70),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.title = element_text(size=30),
    legend.text = element_text(size=25)) -> p_size_count
ggsave("size.png", p_size_count, width = 25, height = 20) 

##Evidence plot
SR_denovo <- subset(denovo, EVIDENCE_FIX == "SR", select = c(ALGORITHMS, GQ))
caller_quality_denovos <- subset(denovo, select = c(name,ALGORITHMS, GQ, EVIDENCE))

denovo %>%
  select(name, SVTYPE, EVIDENCE_FIX) %>%
  unique() %>%
  ggplot(aes(x = forcats::fct_infreq(EVIDENCE_FIX), fill = SVTYPE)) +
  geom_bar() +
  # scale_y_log10()+
  labs(x="Evidence", y="Number de novo SVs", title= "De Novo SVs by Evidence") +
  scale_fill_manual(values = rev(c("DEL"=del_col, "DUP"=dup_col, "INS"=ins_col, "INV"=inv_col, "CPX"=cpx_col, "CTX"=ctx_col))) + 
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 70),
    axis.title = element_text(size = 70),
    plot.title = element_text(size=70),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.title = element_text(size=30),
    legend.text = element_text(size=25)) -> p_evidence
ggsave("evidence.png", p_evidence, width = 25, height = 20) 

# Annotation plot
annot <- grep("_", names(denovo), value = T)
annot <- grep("LINCRNA", annot, value = T, invert = T)
# annot <- annot[!annot %in% c("PROTEIN_CODING__NEAREST_TSS", "PROTEIN_CODING__INTERGENIC")]

lof <- grep("LOF", annot, value = T)
cg <- grep("COPY_GAIN", annot, value = T)
msv <- grep("MSV_EXON_OVR", annot, value = T)
inv_span <- grep("INV_SPAN", annot, value = T)
dup_partial <- grep("DUP_PARTIAL", annot, value = T)
intronic <- grep("INTRONIC", annot, value = T)
utr <- grep("UTR", annot, value = T)
promoter<- grep("PROMOTER", annot, value = T)
nearest_tss<- grep("NEAREST_TSS", annot, value = T)
intergenic<- grep("INTERGENIC", annot, value = T)
coding <- grep("CODING", annot, value = T)
coding <- coding[!coding %in% c('PROTEIN_CODING__NEAREST_TSS', "PROTEIN_CODING__INTERGENIC")]

denovo$is_lof <- apply(denovo[,lof, drop = FALSE], 1, function(r) !all(is.na(r) | r == 'NA' | r == ''))
denovo$is_cg <- apply(denovo[,cg, drop = FALSE], 1, function(r) !all(is.na(r) | r == 'NA' | r == ''))
denovo$is_msv <- apply(denovo[,msv, drop = FALSE], 1, function(r) !all(is.na(r) | r == 'NA' | r == ''))
denovo$is_inv_span <- apply(denovo[,inv_span, drop = FALSE], 1, function(r) !all(is.na(r) | r == 'NA' | r == ''))
denovo$is_dup_partial <- apply(denovo[,dup_partial, drop = FALSE], 1, function(r) !all(is.na(r) | r == 'NA' | r == ''))
denovo$is_intronic <- apply(denovo[,intronic, drop = FALSE], 1, function(r) !all(is.na(r) | r == 'NA' | r == ''))
denovo$is_utr <- apply(denovo[,utr, drop = FALSE], 1, function(r) !all(is.na(r) | r == 'NA' | r == ''))
denovo$is_promoter <- apply(denovo[,promoter, drop = FALSE], 1, function(r) !all(is.na(r) | r == 'NA' | r == ''))
denovo$is_nearest_tss <- apply(denovo[,nearest_tss, drop = FALSE], 1, function(r) !all(is.na(r) | r == 'NA' | r == ''))
denovo$is_intergenic <- apply(denovo[,intergenic, drop = FALSE], 1, function(r) !all(is.na(r) | r == 'NA' | r == '' | r == FALSE))
denovo$is_coding <- apply(denovo[,coding, drop = FALSE], 1, function(r) !all(is.na(r) | r == 'NA' | r == ''))

keep_cols <- c(grep("SVTYPE", names(denovo)), (ncol(denovo)-9-1):ncol(denovo))
denovo_cons <- denovo[,keep_cols]

ref_cols <- grep("is_", names(denovo_cons))
denovo_cons_type <- cbind(denovo_cons[,'SVTYPE', drop = FALSE], (denovo_cons[,ref_cols])*1)

sv_cols <- list(
  list(query = elements, 
       params = list("SVTYPE", c("DEL","DUP", "INS", "INV", "CPX", "CTX")), color = del_col, active = T
       , query.name = "DEL"
  ), #red
  list(query = elements, 
       params = list("SVTYPE", c("DUP", "INS", "INV", "CPX", "CTX")), color = dup_col, active = T
       , query.name = "DUP"
  ), #blue
  list(query = elements, 
       params = list("SVTYPE", c("INS", "INV", "CPX", "CTX")), color = ins_col, active = T
       , query.name = "INS"
  ), #orange
  list(query = elements, 
       params = list("SVTYPE", c("INV", "CPX", "CTX")), color = inv_col, active = T
       , query.name = "INV"
  ), #green
  list(query = elements, 
       params = list("SVTYPE", "CPX", "CTX"), color = cpx_col, active = T
       , query.name = "CPX" )
  # ), #purple
  #list(query = elements, 
  #params = list("SVTYPE", "CTX"), color = ctx_col, active = T,
  #query.name = "CTX"
  #)
)

denovo_cons_type %>% 
  upset(sets = grep("is_", names(denovo_cons_type), value = T),
        #queries = sv_cols,
        order.by = "freq", 
        set_size.show = TRUE,
        # set_size.scale_max = 15000, 
        text.scale = 4, 
        query.legend = "top",
        ) -> p_denovo_upset_all

grob_annotation_upset_plot <- as.grob(p_denovo_upset_all)
ggsave("annotation.png", grob_annotation_upset_plot, width = 30, height = 20)

# Creating a panel of plots
lay <- rbind(c(1,1,1,2,2,2,3,3,3),
             c(1,1,1,2,2,2,3,3,3),
             c(4,4,4,4,4,4,4,4,4),
             c(4,4,4,4,4,4,4,4,4),
             c(5,5,NA,6,6,NA,7,7,NA),
             c(5,5,NA,6,6,NA,7,7,NA),
             c(8,8,8,8,NA,9,9,9,9),
             c(8,8,8,8,NA,9,9,9,9),
             c(10,10,10,10,10,10,NA,NA,NA),
             c(10,10,10,10,10,10,NA,NA,NA))

ml <- grid.arrange(sample_boxplot, type_boxplot, p_type_count, p_chr_count, p_af_count, p_af_count_in_gd, p_af_count_not_in_gd, p_size_count, p_evidence, grob_annotation_upset_plot, layout_matrix = lay, top=textGrob("De Novo SV Data", gp=gpar(fontsize=75)))
ggsave(out_file, ml, width = 80, height = 150, limitsize = FALSE)
