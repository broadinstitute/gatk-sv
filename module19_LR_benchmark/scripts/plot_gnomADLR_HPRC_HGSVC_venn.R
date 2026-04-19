#!/usr/bin/env Rscript
# plot_gnomADLR_HPRC_HGSVC_venn.R
# Parse pairwise-overlap stats and draw 3-way Venn diagrams
# (one PDF per combination of variant_type x Genomic_Location)

suppressPackageStartupMessages({
  library(VennDiagram)
  library(grid)
})

# ── 1. Locate input file ──────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  input_file <- args[1]
} else {
  input_file <- file.path(
    path.expand("~"),
    "Library/CloudStorage/GoogleDrive-zxf0419@gmail.com",
    "My Drive/Talkowski_Lab/gnomAD_LR/final_callsets/benchmark",
    "gnomADLR_HPRC_HGSVC.venn.stat"
  )
}
output_dir <- if (length(args) >= 2) args[2] else dirname(input_file)
cat("Input :", input_file, "\n")
cat("Output:", output_dir, "\n\n")

# ── 2. Parse the stats tables ────────────────────────────────────────────────
# The file is whitespace-separated with either:
#   - 3 stacked pairwise blocks (old format), or
#   - 1 all-3-overlap block followed by 3 pairwise blocks (new format).

raw_lines <- readLines(input_file)
raw_lines <- trimws(raw_lines)

parse_block <- function(lines) {
  # Read lines with a header and data rows separated by blank lines.
  data.frame(
    t(sapply(lines, function(l) strsplit(l, "\\s+")[[1]])),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

# Identify header row indices and data ranges
header_idx <- which(grepl("^variant_type", raw_lines))
block_list <- list()
for (i in seq_along(header_idx)) {
  start <- header_idx[i]
  end   <- if (i < length(header_idx)) header_idx[i + 1] - 1 else length(raw_lines)
  block_lines <- raw_lines[start:end]
  block_lines <- block_lines[nzchar(block_lines)]   # drop blank lines
  block_list[[i]] <- block_lines
}

read_block <- function(bl) {
  hdr  <- strsplit(bl[1], "\\s+")[[1]]
  if (identical(hdr, c("variant_type", "Genomic_Location", "all", "3", "overlap"))) {
    hdr <- c("variant_type", "Genomic_Location", "all_3_overlap")
  }
  rows <- lapply(bl[-1], function(l) strsplit(l, "\\s+")[[1]])
  df   <- as.data.frame(
    do.call(rbind, rows), stringsAsFactors = FALSE
  )
  colnames(df) <- hdr
  for (col in colnames(df)[-c(1, 2)]) df[[col]] <- as.numeric(df[[col]])
  df
}

if (length(block_list) == 4) {
  tbl_all3        <- read_block(block_list[[1]])  # explicit all-3 overlap
  tbl_hprc_hgsvc  <- read_block(block_list[[2]])  # HPRC vs HGSVC
  tbl_gADLR_hprc  <- read_block(block_list[[3]])  # gnomADLR vs HPRC
  tbl_gADLR_hgsvc <- read_block(block_list[[4]])  # gnomADLR vs HGSVC
} else {
  tbl_all3        <- NULL
  tbl_hprc_hgsvc  <- read_block(block_list[[1]])  # HPRC vs HGSVC
  tbl_gADLR_hprc  <- read_block(block_list[[2]])  # gnomADLR vs HPRC
  tbl_gADLR_hgsvc <- read_block(block_list[[3]])  # gnomADLR vs HGSVC
}

# keep keys for merging
if (!is.null(tbl_all3)) {
  tbl_all3[["key"]] <- paste(tbl_all3$variant_type, tbl_all3$Genomic_Location)
}
tbl_hprc_hgsvc[["key"]]  <- paste(tbl_hprc_hgsvc$variant_type, tbl_hprc_hgsvc$Genomic_Location)
tbl_gADLR_hprc[["key"]]  <- paste(tbl_gADLR_hprc$variant_type, tbl_gADLR_hprc$Genomic_Location)
tbl_gADLR_hgsvc[["key"]] <- paste(tbl_gADLR_hgsvc$variant_type, tbl_gADLR_hgsvc$Genomic_Location)

# ── 3. Derive 3-way set sizes ─────────────────────────────────────────────────
# For each key we have three pairwise overlap tables.  Notation:
#   AB = gnomADLR ∩ HPRC  (symmetric average of both reported directions)
#   AC = gnomADLR ∩ HGSVC
#   BC = HPRC ∩ HGSVC
#
# Pairwise totals (two estimates each – take the average):
#   A = gnomADLR_Uniq(tbl2) + AB   =  gnomADLR_Uniq(tbl3) + AC
#   B = HPRC_Uniq(tbl2)    + AB   =  HPRC_Uniq(tbl1)    + BC
#   C = HGSVC_Uniq(tbl3)   + AC   =  HGSVC_Uniq(tbl1)   + BC
#
# 3-way overlap (ABC):
#   - if the new all_3_overlap block is present, use it directly
#   - otherwise estimate via Bonferroni lower bound using HPRC as pivot:
#       ABC_est = max(0, AB + BC − B)
#     where B is the HPRC total (average of two estimates).

keys  <- tbl_hprc_hgsvc$key
stats <- lapply(keys, function(k) {
  r0 <- if (!is.null(tbl_all3)) tbl_all3[tbl_all3$key == k, ] else NULL
  r1 <- tbl_hprc_hgsvc[tbl_hprc_hgsvc$key == k, ]
  r2 <- tbl_gADLR_hprc[tbl_gADLR_hprc$key == k, ]
  r3 <- tbl_gADLR_hgsvc[tbl_gADLR_hgsvc$key == k, ]

  AB <- (r2$gnomADLR_HPRC + r2$HPRC_gnomADLR) / 2
  AC <- (r3$gnomADLR_HGSVC + r3$HGSVC_gnomADLR) / 2
  BC <- (r1$HPRC_HGSVC + r1$HGSVC_HPRC) / 2

  # Total size of each set (two independent estimates – average them)
  A <- ((r2$gnomADLR_Uniq + AB) + (r3$gnomADLR_Uniq + AC)) / 2
  B <- ((r2$HPRC_Uniq     + AB) + (r1$HPRC_Uniq     + BC)) / 2
  C <- ((r3$HGSVC_Uniq    + AC) + (r1$HGSVC_Uniq    + BC)) / 2

  if (!is.null(r0) && nrow(r0) == 1) {
    ABC <- r0$all_3_overlap
  } else {
    ABC <- max(0, AB + BC - B)
    ABC <- min(ABC, AB, AC, BC)
  }

  list(
    variant_type     = r1$variant_type,
    Genomic_Location = r1$Genomic_Location,
    A = A, B = B, C = C,
    AB = AB, AC = AC, BC = BC, ABC = ABC
  )
})
stats_df <- do.call(rbind, lapply(stats, as.data.frame))

# ── 4. Formatting helpers ─────────────────────────────────────────────────────
fmt_number <- function(x) {
  x <- round(x)
  if (x >= 1e6)       sprintf("%.2fM", x / 1e6)
  else if (x >= 1e3)  sprintf("%.1fK", x / 1e3)
  else                as.character(x)
}

nice_type <- c(snv = "SNV", indel = "Indel", sv = "SV")
nice_loc  <- c(sd_sr = "SD/SR", us_rm = "Unique/RM")

label_colors <- c(
  gnomADLR = "#E64B35",   # red-ish
  HPRC     = "#4DBBD5",   # blue
  HGSVC    = "#00A087"    # green
)

# ── 5. One PDF per (variant_type x Genomic_Location) ─────────────────────────
variant_types <- c("snv", "indel", "sv")
genomic_locs  <- c("sd_sr", "us_rm")

for (vtype in variant_types) {
  for (gloc in genomic_locs) {

    r <- stats_df[stats_df$variant_type == vtype &
                  stats_df$Genomic_Location == gloc, ]
    if (nrow(r) == 0) next

    A   <- r$A;  B  <- r$B;  C  <- r$C
    AB  <- r$AB; AC <- r$AC; BC <- r$BC; ABC <- r$ABC

    A <- round(A)
    B <- round(B)
    C <- round(C)
    AB <- round(AB)
    AC <- round(AC)
    BC <- round(BC)
    ABC <- round(ABC)

    area1_only <- max(0, A - AB - AC + ABC)
    area2_only <- max(0, B - AB - BC + ABC)
    area3_only <- max(0, C - AC - BC + ABC)
    area12_only <- max(0, AB - ABC)
    area13_only <- max(0, AC - ABC)
    area23_only <- max(0, BC - ABC)

    # —— Annotation labels ———————————————————————————————————————————————
    title_str <- paste0(
      nice_type[vtype], "  |  ", nice_loc[gloc], "\n",
      "gnomADLR=", fmt_number(A),
      "  HPRC=",  fmt_number(B),
      "  HGSVC=", fmt_number(C)
    )
    # —— Render to PDF ————————————————————————————————————————————————————
    pdf_name <- file.path(
      output_dir,
      sprintf("venn_%s_%s_gnomADLR_HPRC_HGSVC.pdf", vtype, gloc)
    )
    pdf(pdf_name, width = 7, height = 6.5)

    grid.newpage()
    grid.text(
      title_str,
      x = 0.5,
      y = 0.96,
      gp = gpar(fontsize = 11, fontface = "bold")
    )

    venn_grob <- draw.triple.venn(
      area1 = A,
      area2 = B,
      area3 = C,
      n12 = AB,
      n23 = BC,
      n13 = AC,
      n123 = ABC,
      category = c("gnomADLR", "HPRC", "HGSVC"),
      fill = unname(label_colors),
      alpha = rep(0.25, 3),
      cex = 1.0,
      cat.cex = 1.0,
      cat.col = unname(label_colors),
      cat.pos = c(-20, 20, 180),
      cat.dist = c(0.06, 0.06, 0.05),
      lwd = 2,
      scaled = FALSE,
      euler.d = FALSE,
      direct.area = TRUE,
      area.vector = c(
        area1_only,
        area12_only,
        area2_only,
        area13_only,
        ABC,
        area23_only,
        area3_only
      ),
      ind = FALSE,
      margin = 0.12
    )
    grid.draw(venn_grob)

    abc_label <- if (!is.null(tbl_all3)) "ABC_all3" else "ABC_est"

    # Add a legend-like subtitle with raw overlap numbers
    grid.text(
      label = sprintf(
        "AB_overlap=%.0f  AC_overlap=%.0f  BC_overlap=%.0f  %s=%.0f",
        AB, AC, BC, abc_label, ABC
      ),
      x = 0.5, y = 0.03,
      gp = gpar(fontsize = 7.5, col = "grey40")
    )

    dev.off()
    cat(sprintf("Written: %s\n", pdf_name))
  }
}

cat("\nDone. 6 PDFs created.\n")
