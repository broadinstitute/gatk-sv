# ...existing code...

#!R
library("optparse")

# ...existing code...
option_list <- list(
  make_option(c("-p", "--permutate"), type = "character", default = NULL,
              help = "permutation round label, e.g. permi_1", metavar = "character"),
  make_option(c("-s", "--sv-info-file"), type = "character", default = NULL,
              help = "SV info table (tsv/tsv.gz)", metavar = "character"),
  make_option(c("-g", "--gene-info"), type = "character", default = NULL,
              help = "Gene info table (tsv/tsv.gz)", metavar = "character"),
  make_option(c("-r", "--reanno-file"), type = "character", default = NULL,
              help = "Re-annotated SV-vs-gene file", metavar = "character"),
  make_option(c("-o", "--output-rdata"), type = "character", default = NULL,
              help = "Output .RData path", metavar = "character")
)

# ...existing code...

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

required_args <- c("sv_info_file", "gene_info", "reanno_file", "output_rdata", "permutate")
missing_args <- required_args[sapply(required_args, function(x) is.null(opt[[x]]) || opt[[x]] == "")]
if (length(missing_args) > 0) {
  stop(paste("Missing required arguments:", paste(missing_args, collapse = ", ")))
}

print('reading in sv info ...')
sv_info_file <- opt$sv_info_file
sv_info <- read.table(sv_info_file, header = TRUE, comment.char = "", sep = '\t')

print('reading in gene info ...')
permu <- opt$permutate
gene_info <- opt$gene_info
snv.data.permu <- read.table(gene_info, header = TRUE, sep = '\t', comment.char = "")

print('reading in sv vs. genes ...')
re_anno_file <- opt$reanno_file
re_anno_svtype.permutate <- readin.re_anno_svtype(re_anno_file, sv_info)

# ...existing code...

save(gene.data.reanno.permu, file = opt$output_rdata)

# ...existing code...