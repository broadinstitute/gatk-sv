#!/usr/bin/env Rscript
# integrate_rmc_status.R
#
# Adds an rmc_status column to the gene_information tsv.gz produced by
# Rplot.integrate_loeuf_v4_and_gene_blacklist.R.
#
# Matching logic:
#   - transcript_id in the tsv.gz (e.g. ENST00000456328.2) is stripped of the
#     version suffix (everything after the last '.') before joining.
#   - rmc transcript IDs are also stripped the same way.
#   - Precedence: "outlier" overrides "PASS" (if a transcript appears in both).
#
# rmc_status values:
#   "PASS"    - transcript found in rmc_pass_transcripts.tsv
#   "outlier" - transcript found in rmc_outlier_transcripts.tsv
#   NA        - transcript not found in either file

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-i", "--input"),   type="character", default=NULL,
              help="Input tsv.gz: gene_information_loeufv4_pHaplo_pTriplo.with_GATKSV_blacklist.txt.gz",
              metavar="FILE"),
  make_option(c("--rmc-pass"),      type="character", default=NULL,
              help="rmc_pass_transcripts.tsv (transcript IDs, one per row; first column used)",
              metavar="FILE"),
  make_option(c("--rmc-outlier"),   type="character", default=NULL,
              help="rmc_outlier_transcripts.tsv (transcript IDs, one per row; first column used)",
              metavar="FILE"),
  make_option(c("--transcript-col"),type="character", default="transcript_id",
              help="Column name for transcript ID in --input [default: transcript_id]",
              metavar="STRING"),
  make_option(c("-o", "--output"),  type="character", default=NULL,
              help="Output tsv.gz path [default: <input_basename>.rmc.txt.gz]",
              metavar="FILE")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$input))       stop("--input is required")
if (is.null(opt$rmc_pass))    stop("--rmc-pass is required")
if (is.null(opt$rmc_outlier)) stop("--rmc-outlier is required")

# ── output path ───────────────────────────────────────────────────────────────
if (is.null(opt$output)) {
  bn  <- sub("[.]txt[.]gz$", "", sub("[.]tsv[.]gz$", "", basename(opt$input)))
  opt$output <- file.path(dirname(opt$input), paste0(bn, ".rmc.txt.gz"))
}

# ── helper: strip version suffix after last '.' ───────────────────────────────
strip_version <- function(x) sub("[.][^.]+$", "", x)

# ── helper: read first column of a TSV (with or without header) ───────────────
read_transcript_list <- function(path) {
  dat <- read.table(path, header=TRUE, sep="\t", comment.char="",
                    check.names=FALSE, quote="", stringsAsFactors=FALSE)
  strip_version(as.character(dat[[1]]))
}

# ── load rmc lists ────────────────────────────────────────────────────────────
cat("Reading rmc PASS list:", opt$rmc_pass, "\n")
pass_ids    <- unique(read_transcript_list(opt$rmc_pass))
cat("  ->", length(pass_ids), "unique transcript base IDs\n")

cat("Reading rmc outlier list:", opt$rmc_outlier, "\n")
outlier_ids <- unique(read_transcript_list(opt$rmc_outlier))
cat("  ->", length(outlier_ids), "unique transcript base IDs\n")

# ── load main table ───────────────────────────────────────────────────────────
cat("Reading input:", opt$input, "\n")
dat <- read.table(gzfile(opt$input), header=TRUE, sep="\t",
                  comment.char="", check.names=FALSE, quote="",
                  stringsAsFactors=FALSE)
cat("  ->", nrow(dat), "rows,", ncol(dat), "columns\n")

tcol <- opt$transcript_col
if (!tcol %in% colnames(dat)) {
  stop("Column '", tcol, "' not found in input. Available columns:\n  ",
       paste(colnames(dat), collapse=", "))
}

# ── join ──────────────────────────────────────────────────────────────────────
base_ids <- strip_version(as.character(dat[[tcol]]))

rmc_status <- rep(NA_character_, nrow(dat))
rmc_status[base_ids %in% pass_ids]    <- "PASS"
rmc_status[base_ids %in% outlier_ids] <- "outlier"   # outlier overrides PASS

n_pass    <- sum(rmc_status == "PASS",    na.rm=TRUE)
n_outlier <- sum(rmc_status == "outlier", na.rm=TRUE)
n_na      <- sum(is.na(rmc_status))

cat(sprintf("rmc_status assigned: PASS=%d, outlier=%d, NA=%d\n",
            n_pass, n_outlier, n_na))

dat$rmc_status <- rmc_status

# ── write output ──────────────────────────────────────────────────────────────
cat("Writing:", opt$output, "\n")
write.table(dat, gzfile(opt$output), sep="\t", quote=FALSE, row.names=FALSE)
cat("Done.\n")
