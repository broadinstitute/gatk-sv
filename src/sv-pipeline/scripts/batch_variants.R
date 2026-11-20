# Usage: Rscript batch_variants.R <cnvs> <sample_table> <outdir> <pad>
#
# Creates manifests that can be used to subset the bincov matrices. A manifest
# is created for each batch that would need to be queried.
# Each manifest is an R environment containing a list for each variant that
# needs to be queried from that batch. The lists are assigned to the variant
# names and each list has three components: the ranges to query from the bincov
# matrix, the carrier samples and the number of background samples. Each manifest
# is then serialized as an RDS file with the name "{batch}.rdx".
#
#      batch                +---------------------------+
# +=============+           |  $ranges: data.frame      |
# |  variant_1 -|---------> |  $carriers: character()   |
# |  variant_2  |           |  $bg_count: integer(1)    |
# |     ...     |           +---------------------------+
# +=============+                       list
#   environment

# Constants -------------------------------------------------------------------

# The intervals in the bincov matrices are 100 bp so the matrices of large CNVs
# are noisy and take up a lot of memory. For these CNVs, we sample the median
# coverage at equally spaced windows across the SV.

# number of samples to take
LARGE_CNV_SAMPLE_COUNT <- 500L
# size of the sampling window
LARGE_CNV_SAMPLE_WINDOW_SIZE <- 2000L
# minimum CNV size to use sampling strategy
LARGE_CNV_SIZE <- LARGE_CNV_SAMPLE_COUNT * LARGE_CNV_SAMPLE_WINDOW_SIZE

# maximum size of region that tabix can index
TABIX_MAX_SEQLEN <- 536870912L

MAX_BACKGROUND_SAMPLES <- 500L

# Functions ------------------------------------------------------------------

usage <- function(con) {
    cat("usage: Rscript batch_variants.R <cnvs> <sample_table> <outdir> <pad>\n",
        file = con)
}

#' Read the CNVs to shard.
#'
#' File should be six tab-separated columns without a header
#' 1. chromosome
#' 2. start (1-based, inclusive)
#' 3. end (1-based, inclusive)
#' 4. variant ID
#' 5. SV type
#' 6. sample ID
#'
#' @param path `character(1)` Path to file.
#' @returns `data.table` The table in `path`.
read_cnvs <- function(path) {
    stopifnot("CNVs path must be a string" = is.character(path) && length(path) == 1L)
    cnvs <- fread(path, sep = "\t", header = FALSE,
                  col.names = c("chr", "start", "end", "vid", "svtype", "samples"),
                  colClasses = c("character", "integer", "integer", "character",
                                 "character", "character"))
    if (nrow(cnvs) == 0) {
        stop("no variants found")
    }

    if (any(is.na(c(cnvs$start, cnvs$end)))) {
        stop("CNV coordinates must not be `NA`")
    }

    if (!all(cnvs$start <= cnvs$end)) {
        stop("CNV start must be less than or equal to end")
    }

    if (!all(grepl("DEL|DUP", cnvs$svtype))) {
        stop("only DEL and DUP SV types are allowed")
    }

    if (!all(nzchar(cnvs$vid))) {
        stop("all variant IDs must be non-empty")
    }

    if (anyDuplicated(cnvs$vid) != 0) {
        stop("variant IDs must be unique")
    }

    if (any(grepl(.Platform$file.sep, cnvs$vid, fixed = TRUE))) {
        stop("variant IDs must not contain path separators")
    }

    if (!all(nzchar(cnvs$samples))) {
        stop("all sample IDs must be non-empty")
    }

    cnvs
}

#' Read samples and their batches.
#'
#' The file must contain all samples in the cohort. The format is a TSV with
#' two columns:
#' 1. sample ID
#' 2. batch ID
#'
#' @param path `character(1)` Path to file.
#' @returns `data.table` Sample table in `path`.
read_sample_table <- function(path) {
    stopifnot("CNVs path must be a string" = is.character(path) && length(path) == 1L)
    s <- fread(path, sep = "\t", header = FALSE,
               col.names = c("sample", "batch"),
               colClasses = c("character", "character"))

    if (nrow(s) == 0) {
        stop("no samples found")
    }

    if (!all(nzchar(s$sample))) {
        stop("all sample IDs must be non-empty")
    }

    if (anyDuplicated(s$sample) != 0) {
        stop("sample IDs must be unique")
    }

    if (!all(nzchar(s$batch))) {
        stop("all batch IDs must be non-empty")
    }

    s
}

make_hashmap <- function(keys, vals) {
    h <- hashtab(type = "identical", length(keys))
    invisible(mapply(\(k, v) sethash(h, k, v), keys, vals))

    h
}

expand_cnvs <- function(x, pad_factor) {
    sizes <- x$end - x$start + 1L
    pad <- ceiling(sizes * pad_factor)
    x[, start := pmax(1L, start - pad)]
    x[, end := pmin(end + pad, TABIX_MAX_SEQLEN)]

    x
}

tile_cnv <- function(cnv) {
    gaps <- LARGE_CNV_SAMPLE_COUNT - 1L
    total_gap_size <- cnv$end - cnv$start + 1L - LARGE_CNV_SIZE
    pad <- floor(total_gap_size / gaps)
    pads <- rep(pad, gaps)
    remaining_pad <- total_gap_size - pad * gaps
    if (remaining_pad > 0) {
        i <- (gaps - remaining_pad + 1):gaps
        pads[i] <- pads[i] + 1L
    }
    steps <- c(cnv$start, pads + LARGE_CNV_SAMPLE_WINDOW_SIZE)
    starts <- cumsum(steps)
    ends <- starts + LARGE_CNV_SAMPLE_WINDOW_SIZE - 1L

    if (any(is.na(starts))) {
        dump.frames(to.file = TRUE)
        stop("tiling CNVs resulted in `NA` starts")
    }

    if (any(is.na(ends))) {
        stop("tiling CNVs resulted in `NA` ends")
    }

    data.frame(chr = cnv$chr, start = starts, end = ends)
}

make_manifest <- function(cnv, sample_table, sample_map, store) {
    cnv_size <- cnv$end - cnv$start + 1L
    if (cnv_size >= LARGE_CNV_SIZE) {
        ranges <- tile_cnv(cnv)
    } else {
        ranges <- data.frame(chr = cnv$chr, start = cnv$start, end = cnv$end)
    }

    carriers <- strsplit(cnv$samples, split = ",", fixed = TRUE)[[1]]
    batches <- lapply(carriers, \(s) gethash(sample_map, s))
    missing_batches_idx <- which(vapply(batches, is.null, logical(1)))
    if (length(missing_batches_idx) != 0) {
        stop(sprintf("sample '%s' does not have a batch", carriers[missing_batches_idx[[1]]]))
    }

    batches <- unlist(batches)
    batches_uniq <- unique(batches)
    batch_samples <- sample_table[batches_uniq, ]
    bg_samples <- batch_samples[!sample %in% carriers, ]
    if (nrow(bg_samples) == 0) {
        bg_counts <- vector(mode = "list", length = length(batches_uniq))
        names(bg_counts) <- batches_uniq
    } else {
        bg_counts <- sample(bg_samples$batch, min(MAX_BACKGROUND_SAMPLES, nrow(bg_samples))) |>
            table() |>
            as.list()
    }
    carriers_grouped <- split(carriers, batches)

    for (b in batches_uniq) {
        batch_store <- get0(b, envir = store, mode = "environment", inherits = FALSE)
        if (is.null(batch_store)) {
            batch_store <- new.env(parent = emptyenv())
            assign(b, batch_store, pos = store)
        }
        bc <- bg_counts[[b]]
        m <- list(ranges = ranges,
                  carriers = carriers_grouped[[b]],
                  bg_count = if (is.null(bc)) 0 else bc)
        assign(cnv$vid, m, pos = batch_store)
    }
}

write_manifests <- function(store, outdir) {
    batches <- ls(store)
    batch_stores <- mget(batches, envir = store, mode = "environment")
    for (i in seq_along(batches)) {
        saveRDS(batch_stores[[i]], file = file.path(outdir, paste0(batches[[i]], ".rdx")))
    }
}

main <- function() {
    argv <- commandArgs(trailingOnly = TRUE)
    if (length(argv) != 4) {
        usage(stderr())
        quit(save = "no", status = 2)
    }

    suppressPackageStartupMessages(library(data.table))

    padding <- suppressWarnings(as.double(argv[[4]]))
    if (is.na(padding) || is.nan(padding) || !is.finite(padding) || padding < 0) {
        stop("padding must be a finite float >= 0")
    }
    outdir <- argv[[3]]
    dir.create(outdir)
    message("reading CNVs")
    cnvs <- read_cnvs(argv[[1]])
    message("reading sample table")
    sample_table <- read_sample_table(argv[[2]])
    sample_map <- make_hashmap(sample_table$sample, sample_table$batch)
    setkey(sample_table, batch)
    message("expanding CNVs")
    cnvs <- expand_cnvs(cnvs, padding)
    store <- new.env(parent = emptyenv())
    for (i in seq_len(nrow(cnvs))) {
        tmp <- as.list(cnvs[i, ])
        message(paste0("making manifest for ", tmp$vid))
        make_manifest(tmp, sample_table, sample_map, store)
    }
    message("writing manifests")
    write_manifests(store, outdir)
    message("done")
}

main()
