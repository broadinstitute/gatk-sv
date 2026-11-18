# Usage: Rscript visualize_cnvs.R <cnvs> <sample_table> <batch_dir_paths> \
#   <outdir> <one_sample_per_plot>
#
# Plot CNVs

# Constants -------------------------------------------------------------------

# must be odd
SMOOTH_WINDOW <- 21

# number of intervals to plot
INTERVAL_PLOT_COUNT <- 26L

# Functions ------------------------------------------------------------------

usage <- function(con) {
    cat("usage: Rscript visualize_cnvs.R <cnvs> <sample_table> <batch_tars_table> <outdir> <one_sample_per_plot>\n",
        file = con)
}

#' Read the CNVs to plot.
#'
#' File should be six tab-separated columns without a header
#' 1. chromosome
#' 2. start (1-based, inclusive)
#' 3. end (1-based, inclusive)
#' 4. variant ID
#' 5. SV type
#' 6. sample ID
#'
#' @param path `character(1)` Path to the file.
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

#' Read a two-column, tab-separated file into a hash table.
#'
#' The first column contains the keys and the second column contains the
#' values.
#'
#' @param path `character(1)` Path to the file.
#' @returns `utils::hashtab` Hash table.
read_keyed_tsv <- function(path) {
    d <- fread(path, sep = "\t", header = FALSE,
                    col.names = c("V1", "V2"),
                    colClasses = c("character", "character"))
    if (nrow(d) == 0) {
        stop("at least one row must be given in mapping")
    }

    if (anyDuplicated(d$V1) != 0) {
        stop("TSV keys must be unique")
    }

    h <- hashtab(type = "identical", nrow(d))
    mapply(\(k, v) sethash(h, k, v), d$V1, d$V2)

    h
}

select_spaced_intervals <- function(n) {
    if (n <= INTERVAL_PLOT_COUNT) {
        return(1:n)
    }
    gaps <- INTERVAL_PLOT_COUNT - 1L
    total_gap_size <- n - INTERVAL_PLOT_COUNT
    remaining_pad <- total_gap_size %% gaps
    pads <- rep(floor(total_gap_size / gaps), gaps)
    i <- seq_len(remaining_pad)
    pads[i] <- pads[i] + 1L
    steps <- c(1, pads + 1L)

    cumsum(steps)
}

rollmean <- function(x) {
    n <- length(x)
    if (n < SMOOTH_WINDOW) {
        return(frollmean(x, 1:n, adaptive = TRUE))
    }
    partial <- ceiling(SMOOTH_WINDOW / 2)
    a <- vapply(seq(partial, SMOOTH_WINDOW - 1),
                \(y) mean(x[1:y], na.rm = TRUE),
                double(1))
    b <- vapply(seq(SMOOTH_WINDOW - 1, partial),
                \(y) mean(x[(n - y + 1):n], na.rm = TRUE),
                double(1))
    # using algo = "fast" makes NA propagate even with na.rm = TRUE
    m <- frollmean(x, SMOOTH_WINDOW, algo = "exact", align = "center", na.rm = TRUE)
    m[1:(partial - 1)] <- a
    m[(n - partial + 2):n] <- b

    m
}

#' Prettify the CNV size.
#'
#' @param cnv `list` Information of the CNV.
#' @returns `character(1)` A CNV size formatted as a pretty string.
pretty_cnv_size <- function(cnv) {
    cnv_size <- cnv$end - cnv$start + 1L
    if (cnv_size <= 1000) {
        size_pretty <- sprintf("%d b", cnv_size)
    } else if (cnv_size <= 100000) {
        size_pretty <- sprintf("%.2f Kb", cnv_size / 1000)
    } else {
        size_pretty <- sprintf("%.2f Mb", cnv_size / 1000000)
    }

    size_pretty
}

plot_singleton_intervals <- function(x, mids, col, cex) {
    idx <- which(!is.na(x))
    na_idx <- which(is.na(x))
    na_after <- idx[(idx + 1) %in% na_idx]
    na_before <- idx[(idx - 1) %in% na_idx]
    singleton_idx <- sort(intersect(na_after, na_before))
    points(mids[singleton_idx], x[singleton_idx], cex = cex, pch = 19, col = col)

    # Handle intervals at the start and end of vector
    if (length(x) < 2) {
        return()
    }
    if (!is.na(x[[1]]) && is.na(x[[2]])) {
        points(mids[[1]], x[[1]], cex = cex, pch = 19, col = col)
    }
    if (!is.na(x[[length(x)]]) && is.na(x[[length(x) - 1]])) {
        points(mids[[length(mids)]], x[[length(x)]], cex = cex, pch = 19, col = col)
    }
}

#' Plot coverage over a CNV region.
#'
#' @param cnv `list` Information of the CNV.
#' @param norm_cov `data.table` Normalized coverage matrix for the CNV region. Should be smoothed
#'   sampled.
#' @param carriers `character` Carrier samples.
#' @param outfile The path to the plot.
plot_cnv <- function(cnv, norm_cov, carriers, outfile) {
    # make plot main
    size_pretty <- pretty_cnv_size(cnv)
    start_pretty <- formatC(cnv$start, big.mark = ",", format = "d")
    end_pretty <- formatC(cnv$end, big.mark = ",", format = "d")
    # if we are plotting one sample per plot, the carrier sample ID could be
    # different from the comma-separated list of carrier sample IDs
    if (length(carriers) == 1) {
        sample_string <- carriers
    } else {
        sample_string <- cnv$samples
    }

    if (nchar(sample_string) > 30) {
        samples_pretty <- paste0(substr(sample_string, 1, 30), "...")
    } else {
        samples_pretty <- sample_string
    }
    main <- sprintf("%s:%s-%s (hg38)\n%s (%s)",
                    cnv$chr, start_pretty, end_pretty, samples_pretty, size_pretty)

    mids <- norm_cov$start + floor((norm_cov$end - norm_cov$start + 1) / 2)
    xlim <- c(min(norm_cov$start), max(norm_cov$end))
    jpeg(outfile, res = 300, width = 1800, height = 1800)
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    par(mar = c(6.1, 4.1, 4.1, 2.1))
    cnv_col <- if (cnv$svtype == "DEL") "red" else "blue"
    # plot background
    matplot(mids, as.matrix(norm_cov[, .SD, .SDcols = !patterns("chr|start|end")]),
            type = "l", lty = 1, col = "grey", lwd = 0.5,
            xlim = xlim,
            main = main,
            ylab = "Normalized Read Depth Ratio",
            xlab = "",
            xaxs = "i",
            xaxt = "n")
    # plot carriers
    matlines(mids, as.matrix(norm_cov[, .SD, .SDcols = carriers]),
             lty = 1, lwd = 2,
             col = cnv_col)
    # lonely intervals - background
    bg_samples <- setdiff(colnames(norm_cov), c("chr", "start", "end", carriers))
    for (bg in norm_cov[, .SD, .SDcols = bg_samples]) {
        plot_singleton_intervals(bg, mids, "grey", 0.5)
    }
    # lonely intervals - carriers
    for (fg in norm_cov[, .SD, .SDcols = carriers]) {
        plot_singleton_intervals(fg, mids, cnv_col, 2)
    }
    # add padding indicators
    ylim <- par("usr")[3:4]
    rect(xlim[[1]], ylim[[1]], cnv$start, ylim[[2]], col = "#FFAF0044", border = NA)
    rect(cnv$end, ylim[[1]], xlim[[2]], ylim[[2]], col = "#FFAF0044", border = NA)
    # add x-axis
    xticks <- axisTicks(xlim, log = FALSE, nint = 10)
    xlabs <- formatC(xticks, big.mark = ",", format = "d")
    axis(1, at = xticks, labels = xlabs, las = 2, cex.axis = 0.8)
    title(xlab = sprintf("%s Position (bp)", cnv$chr), line = 5)
    dev.off()
}

get_bincov_from_dir <- function(path, variant, all_carriers, batch_carriers) {
    rdx <- file.path(path, sprintf("%s.rdx", variant))
    mat <- readRDS(rdx)
    mat_samples <- setdiff(colnames(mat), c("chr", "start", "end"))
    # find and remove carrier samples that are not carrier samples for this batch
    # this addresses the case where samples are mistakenly assigned to multiple
    # batches and become a carrier sample in one batch, but a background sample
    # in another for the same variant
    dups <- setdiff(intersect(mat_samples, all_carriers), batch_carriers)
    if (length(dups) > 0) {
        mat[, (dups) := list(NULL)]
    }

    mat
}

make_plot <- function(cnv, sample_map, batch_dir_map, one_sample_per_plot, outdir) {
    samples <- strsplit(cnv$sample, split = ",", fixed = TRUE)[[1]]
    batches <- vapply(samples, \(x) gethash(sample_map, x), character(1), USE.NAMES = FALSE)
    grouped_samples <- split(samples, batches)
    batch_dirs <- vapply(names(grouped_samples), \(x) gethash(batch_dir_map, x), character(1),
                         USE.NAMES = FALSE)
    bincovs <- mapply(\(x, y) get_bincov_from_dir(x, cnv$vid, samples, y),
                      batch_dirs, grouped_samples, SIMPLIFY = FALSE, USE.NAMES = FALSE)
    # assume bincovs are keyed on chr, start, end
    merged <- Reduce(merge, bincovs)

    # smooth
    merged[, names(.SD) := lapply(.SD, rollmean), .SDcols = !patterns("chr|start|end")]

    # sample intervals
    merged <- merged[select_spaced_intervals(.N), ]

    if (one_sample_per_plot) {
        bg_samples <- setdiff(colnames(merged), c("chr", "start", "end", samples))
        for (s in samples) {
            plot_name <- sprintf("%s_%d-%d_%s_%s_%s.jpg",
                                 cnv$chr, cnv$start, cnv$end, cnv$vid,
                                 cnv$svtype, s)
            plot_path <- file.path(outdir, plot_name)
            keep <- c("chr", "start", "end", s, bg_samples)
            plot_cnv(cnv, merged[, .SD, .SDcols = keep], s, plot_path)
        }
    } else {
        plot_name <- sprintf("%s_%d-%d_%s_%s.jpg",
                             cnv$chr, cnv$start, cnv$end, cnv$vid, cnv$svtype)
        plot_path <- file.path(outdir, plot_name)
        plot_cnv(cnv, merged, samples, plot_path)
    }
}

main <- function() {
    argv <- commandArgs(trailingOnly = TRUE)
    if (length(argv) != 5) {
        usage(stderr())
        quit(save = "no", status = 2)
    }

    suppressPackageStartupMessages(library(data.table))

    one_sample_per_plot <- if (argv[[5]] == "0") FALSE else TRUE
    outdir <- argv[[4]]
    dir.create(outdir)
    message("reading CNVs")
    cnvs <- read_cnvs(argv[[1]])
    message("reading sample table")
    sample_map <- read_keyed_tsv(argv[[2]])
    message("reading batch directory table")
    batch_dir_map <- read_keyed_tsv(argv[[3]])

    for (i in seq_len(nrow(cnvs))) {
        tmp <- as.list(cnvs[i, ])
        message(paste0("plotting ", tmp$vid))
        make_plot(tmp, sample_map, batch_dir_map, one_sample_per_plot, outdir)
    }
    message("done")
}

main()
