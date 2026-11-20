# Visualize the PE, SR, and RD evidence of CNV carriers.
#
# Adapted from visualize_denovos.R in
# https://github.com/talkowski-lab/gatk-sv-utils for non-trio CNV plots. One
# plot is made per variant and carrier sample. The PE and SR tracks show the
# evidence of the carrier and the RD track shows the carrier against a random
# set of non-carrier samples from the same batch.
#
# The evidence files can be local paths or gs:// URIs. Remote files are
# streamed with the system tabix, which must be built with GCS support.
#
# Usage:
# Rscript visualize_cnv.R <variants> <evidence> <samples> <outdir> <exclusions>
#
# <variants>
# TSV of variants to visualize with a header line. The required columns are
# chr, start, end, svlen, vid, svtype, and sample_id. Other columns are ignored.
#
# <evidence>
# Headerless TSV with paths to the SV evidence files with columns in the order
# batch ID, PE, SR, RD, median coverage.
#
# <samples>
# Headerless TSV mapping samples to batches with columns in the order
# sample ID, batch ID.
#
# <outdir>
# Directory to write the plots.
#
# <exclusions>
# File to write the variants that were skipped.

REQ_VARIANT_FIELDS <- c("chr", "start", "end", "svlen", "vid", "svtype", "sample_id")

# maximum number of non-carrier samples to plot in the RD track
MAX_BG_SAMPLES <- 100L

# colors for the carrier evidence
CARRIER_COLOR <- "#bc5d41"

# relative heights of the PE (multiplied by the PE track scale), SR, RD,
# repeats, and genes tracks
TRACK_HEIGHTS <- c(0.6, 0.3, 1, 0.1, 0.6)

usage <- function(con) {
    cat(
        "usage: Rscript visualize_cnv.R <variants> <evidence> <samples> <outdir> <exclusions>\n",
        file = con
    )
}

write_exclusion <- function(con, sample_id, vid, reason) {
    writeLines(sprintf("%s\t%s\t%s", sample_id, vid, reason), con)
}

plot_carrier <- function(sve, carrier, bg_samples, vid, plot_path) {
    carrier_sve <- subset_samples(sve, carrier)
    rd_sve <- subset_samples(sve, c(carrier, bg_samples))

    pe <- pe_plotter(carrier_sve$pe)

    # make the SR track plot over the same interval as the other tracks
    synced_sr <- carrier_sve$sr
    synced_sr$region <- carrier_sve$pe$region
    sr <- sr_plotter(synced_sr)

    rd <- rd_plotter(rd_sve$rd, median_region = sve$region)

    region <- sve$region
    genes <- gorilla:::new_genes_plotter(region$contig, region$qstart, region$qend)
    repeats <- gorilla:::new_repeats_plotter(region$contig, region$qstart, region$qend)

    min_lanes <- gorilla:::MIN_PE_LANES_PER_SAMPLE
    max_lanes <- gorilla:::MAX_PE_LANES_PER_SAMPLE
    pe_lanes <- if (any(pe$lanes_per_sample$lanes > min_lanes)) max_lanes else min_lanes
    heights <- TRACK_HEIGHTS * c(pe_lanes / min_lanes, 1, 1, 1, 1)

    main <- sprintf(
        "%s %s\n%s (%s)",
        sve$svtype,
        vid,
        gorilla:::pretty_sample_id(carrier),
        gorilla:::pretty_size(region$end - region$start + 1)
    )

    height_scale <- sum(heights) / sum(TRACK_HEIGHTS)
    png(plot_path, width = 3840, height = 2160 * height_scale, res = 300)
    on.exit(dev.off(), add = TRUE)
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par), add = TRUE, after = FALSE)

    layout(matrix(1:5, nrow = 5), heights = heights)

    par(mar = c(0, 4.1, 4.1, 2.1), mgp = c(1, 1, 0))
    plot(pe, samples = carrier, main = main, col = CARRIER_COLOR)
    gorilla:::draw_padding(region$start, region$end)

    par(mar = c(0, 4.1, 0, 2.1), mgp = c(1, 1, 0))
    plot(sr, samples = carrier, col = CARRIER_COLOR)
    gorilla:::draw_padding(region$start, region$end)

    par(mgp = c(3, 1, 0))
    plot(
        rd,
        samples = carrier,
        bg_samples = if (length(bg_samples) > 0) bg_samples else NULL,
        col = CARRIER_COLOR
    )
    gorilla:::draw_padding(region$start, region$end)

    plot(repeats)

    par(mar = c(3.1, 4.1, 0, 2.1))
    plot(genes)
    gorilla:::draw_genomic_coords_axis(region$contig)
}

main <- function() {
    argv <- commandArgs(trailingOnly = TRUE)
    if (length(argv) != 5) {
        usage(stderr())
        quit(save = "no", status = 2)
    }

    suppressPackageStartupMessages(library(gorilla))

    outdir <- argv[[4]]
    dir.create(outdir, showWarnings = FALSE)

    variants <- fread(
        argv[[1]], header = TRUE, sep = "\t",
        colClasses = list(character = c("chr", "vid", "svtype", "sample_id"))
    )
    if (!all(REQ_VARIANT_FIELDS %in% colnames(variants))) {
        stop(
            paste0(
                "all of these fields must be in the variant file: ",
                paste0(REQ_VARIANT_FIELDS, collapse = ", ")
            ),
            call. = FALSE
        )
    }

    evidence_paths <- fread(
        argv[[2]], header = FALSE, sep = "\t",
        colClasses = rep("character", 5),
        col.names = c("batch_id", "pe_path", "sr_path", "rd_path", "medians_path"),
        key = "batch_id"
    )

    samples <- fread(
        argv[[3]], header = FALSE, sep = "\t",
        colClasses = c("character", "character"),
        col.names = c("sample_id", "batch_id"),
        key = "sample_id"
    )

    exclusions_fp <- file(argv[[5]], open = "wt")
    on.exit(close(exclusions_fp), add = TRUE)
    writeLines("sample_id\tvid\texclusion_reason", exclusions_fp)

    tmpdir <- tempfile("temp", ".")
    if (!dir.create(tmpdir, showWarnings = FALSE)) {
        stop("failed to create temporary directory in working directory")
    }

    if (anyDuplicated(samples$sample_id) != 0) {
        stop("sample IDs in the sample table must be unique")
    }

    set.seed(42)
    batch_files <- list()
    variants[, batch_id := samples[J(variants$sample_id), batch_id]]
    for (i in which(is.na(variants$batch_id))) {
        write_exclusion(exclusions_fp, variants$sample_id[[i]], variants$vid[[i]], "sample has no batch")
    }

    # query the evidence once for each variant and batch
    groups <- unique(variants[!is.na(batch_id), list(vid, batch_id)])
    for (g in seq_len(nrow(groups))) {
        vid_g <- groups$vid[[g]]
        batch_g <- groups$batch_id[[g]]
        rows <- variants[vid == vid_g & batch_id == batch_g, ]
        v <- rows[1, ]
        message(sprintf("plotting %s in batch %s", vid_g, batch_g))

        if (!v$svtype %in% c("DEL", "DUP")) {
            for (s in rows$sample_id) {
                write_exclusion(exclusions_fp, s, vid_g, "SV type is not DEL or DUP")
            }
            next
        }

        if (is.null(batch_files[[batch_g]])) {
            paths <- evidence_paths[batch_g, nomatch = NULL]
            if (nrow(paths) == 0) {
                stop(sprintf("batch '%s' has no evidence files", batch_g))
            }
            cachedir <- file.path(tmpdir, batch_g)
            batch_files[[batch_g]] <- list(
                pe = pe_file(paths$pe_path, cachedir),
                sr = sr_file(paths$sr_path, cachedir),
                rd = rd_file(paths$rd_path, paths$medians_path, cachedir)
            )
        }
        bf <- batch_files[[batch_g]]

        pad <- max(1000 / (v$end - v$start + 1), 0.3)
        sve <- svevidence(v$chr, v$start, v$end, bf$pe, bf$sr, bf$rd, v$svtype, pad = pad)

        if (nrow(sve$rd$mat) == 0) {
            for (s in rows$sample_id) {
                write_exclusion(exclusions_fp, s, vid_g, "no RD bins in region")
            }
            next
        }

        rd_samples <- colnames(sve$rd$mat)
        carriers <- intersect(rows$sample_id, rd_samples)
        for (s in setdiff(rows$sample_id, rd_samples)) {
            write_exclusion(exclusions_fp, s, vid_g, "sample missing from RD matrix")
        }
        if (length(carriers) == 0) {
            next
        }

        non_carriers <- setdiff(rd_samples, variants[vid == vid_g, sample_id])
        bg_samples <- if (length(non_carriers) > MAX_BG_SAMPLES) {
            sample(non_carriers, MAX_BG_SAMPLES)
        } else {
            non_carriers
        }

        for (s in carriers) {
            plot_path <- file.path(
                outdir,
                sprintf("%s_%d-%d_%s_%s_%s.png", v$chr, v$start, v$end, vid_g, v$svtype, s)
            )
            plot_carrier(sve, s, bg_samples, vid_g, plot_path)
        }
    }
}

main()
