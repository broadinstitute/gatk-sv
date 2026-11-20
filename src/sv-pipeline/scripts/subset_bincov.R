# Usage: Rscript subset_bincov.R <manifest> <bincov> <medians> <outdir>
#
# Generate the subset and normalized coverage matrix for each variant in
# <manifest>. Each subset is then serialized as an RDS file.

# Functions ------------------------------------------------------------------

usage <- function(con) {
    cat("usage: Rscript subset_bincov.R <manifest> <bincov> <medians> <outdir>\n",
        file = con)
}

#' Retrieve the header of a bincov matrix.
#'
#' @param con `TabixFile` Connection to matrix.
#' @returns `character` Header of the matrix.
bincov_header <- function(con) {
    header <- headerTabix(con)[["header"]]
    if (is.null(header)) {
        stop(sprintf("bincov matrix is missing a header", con$ath))
    }
    header <- strsplit(header, split = "\t", fixed = TRUE)[[1]]
    if (length(header) <= 4 || !all(header[1:3] == c("#Chr", "Start", "End"))) {
        stop(sprintf("bincov matrix has an invalid header", con$path))
    }
    header[1:3] <- c("chr", "start", "end")

    header
}

#' Read coverage medians file into a named vector.
#'
#' @param path `character(1)` Path to the coverage medians file.
#' @returns `double` The coverage medians. The names of the vector are the
#'   sample IDs.
read_medians <- function(path) {
    lines <- readLines(path, n = 2L, ok = FALSE)
    ids <- strsplit(lines[[1]], split = "\t", fixed = TRUE)[[1]]
    if (anyDuplicated(ids) != 0) {
        stop("duplicate samples in medians file")
    }
    medians <- as.double(strsplit(lines[[2]], split = "\t", fixed = TRUE)[[1]])

    names(medians) <- ids

    medians
}

#' Query a bincov matrix.
#'
#' @param con `TabixFile` Connection to bincov matrix.
#' @param header `character` The parsed header line from the bincov matrix.
#'   Assumes the first three elements are "chr", "start", "end" and the
#'   remaining elements are sample IDs.
#' @param select `integer` Which columns to keep.
#' @param ranges `GRanges` Genomic ranges to check for overlap.
#' @returns `data.table` The bincov matrix. If `ranges` contains a single
#'   range, the overlapping bincov rows will be returned as is. Otherwise,
#'   the median coverage of each sample for the bins overlapping each
#'   range will be returned. Ranges without any overlapping bins will have
#'   `NA` coverage values unless all ranges do not have any overlapping bins
#'   in which case an error is signaled.
query_bincov <- function(con, header, select, ranges) {
    tabix <- scanTabix(con, param = ranges)
    if (all(lengths(tabix) == 0)) {
        stop("no overlapping intervals in coverage matrix")
    }

    coltypes <- list("character" = 1L,
                     "integer" = c(2L, 3L),
                     "double" = seq.int(4, length(header)))
    qcoords <- data.table(
        qchr = as.character(seqnames(ranges)),
        qstart = start(ranges),
        qend = end(ranges)
    )
    sub_header <- header[select]

    parse <- function(qchr, qstart, qend, lines) {
        if (length(lines) == 0) {
            return(NULL)
        }
        # select determines order of columns
        d <- fread(text = lines, sep = "\t", header = FALSE,
                   colClasses = coltypes, select = select)
        colnames(d) <- sub_header
        # convert from 0-start, exclusive to 1-start, inclusive
        d[, start := start + 1L]
        d[, c("qchr", "qstart", "qend") := list(qchr, qstart, qend)]
        d
    }

    mat <- mapply(parse, qcoords$qchr, qcoords$qstart, qcoords$qend, tabix,
                  SIMPLIFY = FALSE, USE.NAMES = FALSE) |>
        rbindlist(use.names = TRUE)

    # a single range indicates the query region is the span of the CNV while multiple regions
    # indicates interval sampling
    if (length(ranges) > 1) {
        samples <- setdiff(sub_header, c("chr", "start", "end"))
        mat <- mat[qcoords, on = c("qchr", "qstart", "qend")]
        mat <- mat[, lapply(.SD, median, na.rm = TRUE), .SDcols = samples, by = c("qchr", "qstart", "qend")]
        setnames(mat, colnames(qcoords), c("chr", "start", "end"))
    } else {
        mat[, c("qchr", "qstart", "qend") := list(NULL)]
    }
    setkey(mat, chr, start, end)

    mat
}

make_bincov_subset <- function(query, header, bincov_con, outdir) {
    ranges_df <- query[["ranges"]]
    ranges_gr <- GRanges(ranges_df$chr, IRanges(ranges_df$start, ranges_df$end))

    carriers <- query[["carriers"]]
    carrier_idx <- which(header %in% carriers)
    if (length(carriers) != length(carrier_idx)) {
        stop("some carrier samples not found in bincov")
    }

    bg_count <- query[["bg_count"]]
    if (bg_count > 0L) {
        bg_idx <- which(!header %in% c("chr", "start", "end", carriers))
        if (length(bg_idx) < bg_count) {
            warning("fewer background samples in bincov than expected", immediate. = TRUE)
        }
        keep_idx <- c(1:3, carrier_idx, sample(bg_idx, bg_count))
    } else {
        keep_idx <- c(1:3, carrier_idx)
    }

    keep_idx <- sort(keep_idx)

    query_bincov(bincov_con, header, keep_idx, ranges_gr)
}

#' Normalize coverage values.
#'
#' Normalization is done by dividing each sample's coverage values by the
#' genome-wide median coverage of the sample.
#'
#' @param x `data.table` Coverage matrix.
#' @param medians `double` Genome-wide median coverage for each sample in `x`.
#' @param `data.table` Normalized coverage matrix.
normalize_cov <- function(x, medians) {
    cov_cols <- colnames(x)
    cov_samples <- cov_cols[!cov_cols %in% c("chr", "start", "end")]
    medians_samples <- names(medians)
    common_samples <- intersect(cov_samples, medians_samples)

    medians <- medians[common_samples]
    x[, names(.SD) := mapply(`/`, .SD, medians, SIMPLIFY = FALSE, USE.NAMES = FALSE), .SDcols = names(medians)]

    x
}

main <- function() {
    argv <- commandArgs(trailingOnly = TRUE)
    if (length(argv) != 4) {
        usage(stderr())
        quit(save = "no", status = 2)
    }

    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(GenomicRanges))
    suppressPackageStartupMessages(library(Rsamtools))

    outdir <- argv[[4]]
    dir.create(outdir)
    manifest <- readRDS(argv[[1]])
    bincov <- argv[[2]]
    bincov_con <- TabixFile(bincov)
    header <- bincov_header(bincov_con)
    medians <- read_medians(argv[[3]])
    if (!all(header[-(1:3)] %in% names(medians))) {
        stop("some bincov samples are not in the medians file")
    }

    variants <- ls(manifest)
    for (v in variants) {
        message(sprintf("query variant '%s'", v))
        query <- get(v, pos = manifest, mode = "list", inherits = FALSE)
        subset <- make_bincov_subset(query, header, bincov_con, outdir)
        norms <- normalize_cov(subset, medians)

        saveRDS(norms, file.path(outdir, paste0(v, ".rdx")))
    }
}

main()