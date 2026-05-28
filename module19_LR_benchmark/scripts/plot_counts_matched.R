#!/usr/bin/env Rscript
# Bar plot: 5 or 7 variant classes, overlaid with dbSNP Matched slash strip
library(ggplot2)
library(scales)

ANNO_DIR <- "/Users/xzhao/Library/CloudStorage/GoogleDrive-zxf0419@gmail.com/My Drive/Talkowski_Lab/gnomAD_LR/annotation"

VAR_COLS_7 <- c("SNV", "INS 1-49bp", "DEL 1-49bp",
                "INS 50-499bp", "DEL 50-499bp",
                "INS >499bp", "DEL >499bp")

BAR_COLORS_5 <- c(
  "SNV"        = "#7B2D8B",
  "INS 1-49bp" = "#FFB6C1",
  "DEL 1-49bp" = "#CC0000",
  "INS >=50bp" = "#FF69B4",
  "DEL >=50bp" = "#FF4444"
)
VAR_COLS_5 <- names(BAR_COLORS_5)

BAR_COLORS_7 <- c(
  "SNV"          = "#7B2D8B",
  "INS 1-49bp"   = "#FFB6C1",
  "DEL 1-49bp"   = "#CC0000",
  "INS 50-499bp" = "#FF69B4",
  "DEL 50-499bp" = "#FF4444",
  "INS >499bp"   = "#FF1493",
  "DEL >499bp"   = "#FF6666"
)

# Overlay colour
DBSNP_COLOR <- "#4477CC"   # blue (slash line colour)

# ── Helpers ──────────────────────────────────────────────────────────────────
fmt_count <- function(x) {
  ifelse(x >= 1e9, sprintf("%.1fB", x / 1e9),
  ifelse(x >= 1e6, sprintf("%.1fM", x / 1e6),
  ifelse(x >= 1e3, sprintf("%.1fK", x / 1e3),
         sprintf("%.0f", x))))
}

merge5 <- function(v7) {
  c(v7["SNV"], v7["INS 1-49bp"], v7["DEL 1-49bp"],
    v7["INS 50-499bp"] + v7["INS >499bp"],
    v7["DEL 50-499bp"] + v7["DEL >499bp"])
}

extract_counts <- function(df, cat1, cat2 = "All", cat3 = "All", cat4 = "All") {
  row <- df[df$category     == cat1 &
            df$sub_category  == cat2 &
            df$tr_status     == cat3 &
            df$region        == cat4, ]
  row <- row[!grepl("prop", row$sub_category, ignore.case = TRUE), ]
  if (nrow(row) == 0) return(setNames(rep(0, length(VAR_COLS_7)), VAR_COLS_7))
  setNames(as.numeric(row[1, VAR_COLS_7]), VAR_COLS_7)
}

# Diagonal slash segments within a horizontal band [y_bots[i], y_tops[i]]
make_slash_strip <- function(x_pos, half_w, y_bots, y_tops, n_lines = 10) {
  rows <- lapply(seq_along(x_pos), function(i) {
    yb <- y_bots[i];  yt <- y_tops[i]
    if (is.na(yt) || is.na(yb) || yt <= yb) return(NULL)
    hb <- log10(max(yb, 1));  ht <- log10(max(yt, 1.001))
    if (ht <= hb) return(NULL)
    xl <- x_pos[i] - half_w;  xr <- x_pos[i] + half_w
    w  <- xr - xl;  dh <- ht - hb;  slope <- dh / w
    offsets <- seq(xl - w, xr, length.out = n_lines + 2)
    offsets <- offsets[-c(1, length(offsets))]
    segs <- lapply(offsets, function(x0) {
      xs <- max(xl, x0);  xe <- min(xr, x0 + w)
      lys <- (xs - x0) * slope;  lye <- (xe - x0) * slope
      if (lys < 0)  { xs <- x0;              lys <- 0   }
      if (lye > dh) { xe <- x0 + dh / slope; lye <- dh  }
      xs <- max(xs, xl);  xe <- min(xe, xr)
      if (is.na(xs) || is.na(xe) || xs >= xe) return(NULL)
      data.frame(x = xs, xend = xe,
                 y    = 10^(hb + max(lys, 0)),
                 yend = 10^(hb + min(lye, dh)),
                 stringsAsFactors = FALSE)
    })
    do.call(rbind, Filter(Negate(is.null), segs))
  })
  do.call(rbind, Filter(Negate(is.null), rows))
}

# ── Core plotting function ────────────────────────────────────────────────────
# var_labels    : character vector of bar names
# bar_colors_vec: fill colour per bar (same length)
# all_n, dbsnp_n: numeric vectors (same length)
build_match_plot <- function(var_labels, bar_colors_vec,
                             all_n, dbsnp_n,
                             title_str, y_label, out_file, plot_width) {

  all_n      <- pmax(all_n, 1)
  prop_dbsnp <- dbsnp_n / all_n
  prop_dbsnp[is.nan(prop_dbsnp)] <- 0

  log_all   <- log10(all_n)
  top_dbsnp <- 10^(log_all * prop_dbsnp)

  n      <- length(var_labels)
  x_pos  <- seq_len(n)
  half_w <- 0.4

  bar_df <- data.frame(
    var_class   = factor(var_labels, levels = var_labels),
    fill_col    = bar_colors_vec,
    count_all   = all_n,
    count_label = fmt_count(all_n),
    label_y     = all_n * 2.5,
    stringsAsFactors = FALSE
  )

  # Slash strip for dbSNP Matched
  db_slash <- make_slash_strip(x_pos, half_w, rep(1, n), top_dbsnp)

  # % label just above strip top
  lbl_df <- data.frame(
    x     = x_pos,
    y     = top_dbsnp * 1.25,
    label = paste0(round(prop_dbsnp * 100, 1), "%"),
    stringsAsFactors = FALSE
  )
  lbl_df <- lbl_df[top_dbsnp > 1.15, ]

  # Legend
  mx  <- max(all_n)
  lx0 <- n + 0.65;  lx1 <- n + 0.90
  leg_yb <- mx * 0.25;  leg_yt <- mx * 0.55

  p <- ggplot() +
    geom_col(data = bar_df, aes(x = var_class, y = count_all),
             fill = bar_df$fill_col, color = NA, width = 2 * half_w)

  if (!is.null(db_slash) && nrow(db_slash) > 0)
    p <- p + geom_segment(data = db_slash,
                          aes(x = x, xend = xend, y = y, yend = yend),
                          color = DBSNP_COLOR, linewidth = 0.55, alpha = 0.9,
                          inherit.aes = FALSE)

  if (nrow(lbl_df) > 0)
    p <- p + geom_text(data = lbl_df,
                       aes(x = x, y = y, label = label),
                       size = 2.5, color = DBSNP_COLOR, fontface = "bold",
                       vjust = 0, inherit.aes = FALSE)

  p <- p +
    geom_text(data = bar_df,
              aes(x = var_class, y = label_y, label = count_label),
              size = 3.2, fontface = "bold", color = "#333333") +
    scale_y_log10(
      labels = label_number(scale_cut = cut_short_scale()),
      expand = expansion(mult = c(0, 0)),
      limits = c(1, max(all_n) * 10)
    ) +
    labs(
      x        = "Variant class",
      y        = y_label,
      title    = title_str,
      subtitle = "Slash strip = dbSNP Matched; height in log10-space = log10(N) x proportion"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x        = element_text(angle = 25, hjust = 1, size = 11),
      axis.text.y        = element_text(size = 10),
      axis.title         = element_text(size = 12),
      plot.title         = element_text(face = "bold", size = 14),
      plot.subtitle      = element_text(size = 8, color = "grey40"),
      panel.grid.major.y = element_line(linetype = "dashed", color = "grey80",
                                        linewidth = 0.4)
    )

  # Legend: mini slash box
  leg_slash <- make_slash_strip(c((lx0 + lx1) / 2), (lx1 - lx0) / 2,
                                c(leg_yb), c(leg_yt), n_lines = 5)
  p <- p +
    annotate("rect", xmin = lx0, xmax = lx1, ymin = leg_yb, ymax = leg_yt,
             fill = NA, color = DBSNP_COLOR, linewidth = 0.5)
  if (!is.null(leg_slash) && nrow(leg_slash) > 0)
    p <- p + annotate("segment",
                      x = leg_slash$x, xend = leg_slash$xend,
                      y = leg_slash$y, yend = leg_slash$yend,
                      color = DBSNP_COLOR, linewidth = 0.5)
  p <- p +
    annotate("text", x = lx1 + 0.08, y = sqrt(leg_yb * leg_yt),
             label = "dbSNP Matched", hjust = 0, size = 3.0, fontface = "bold",
             color = DBSNP_COLOR) +
    coord_cartesian(xlim = c(0.5, n + 2.8), clip = "off")

  ggsave(out_file, plot = p, width = plot_width, height = 6.5, device = "pdf")
  cat("Saved:", out_file, "\n")
}

# ── Dataset loop ──────────────────────────────────────────────────────────────
datasets <- list(
  list(file   = "aou_I.counts_sites.tsv",
       prefix = "aou_I.counts_sites.barplot_matched",
       label  = "AoU variant site counts",
       ylabel = "Number of variant sites"),
  list(file   = "aou_I.counts_samples.tsv",
       prefix = "aou_I.counts_samples.barplot_matched",
       label  = "AoU per-sample variant counts",
       ylabel = "Mean variant count per sample"),
  list(file   = "hgsv_hprc.counts_sites.tsv",
       prefix = "hgsv_hprc.counts_sites.barplot_matched",
       label  = "HGSV/HPRC variant site counts",
       ylabel = "Number of variant sites"),
  list(file   = "hgsv_hprc.counts_samples.tsv",
       prefix = "hgsv_hprc.counts_samples.barplot_matched",
       label  = "HGSV/HPRC per-sample variant counts",
       ylabel = "Mean variant count per sample")
)

for (ds in datasets) {
  cat("\n── Processing:", ds$file, "\n")
  df <- read.table(file.path(ANNO_DIR, ds$file),
                   header = TRUE, sep = "\t", check.names = FALSE)

  ca_all   <- extract_counts(df, "All")[VAR_COLS_7]
  ca_dbsnp <- extract_counts(df, "Matched", cat2 = "dbSNP Matched")[VAR_COLS_7]

  total_N <- sum(ca_all)
  title_str <- sprintf("%s (N = %s)", ds$label, fmt_count(total_N))

  # ── 5-bar version
  m5 <- function(v) unname(merge5(v))
  build_match_plot(
    var_labels     = VAR_COLS_5,
    bar_colors_vec = unname(BAR_COLORS_5),
    all_n          = m5(ca_all),
    dbsnp_n        = m5(ca_dbsnp),
    title_str      = title_str,
    y_label        = ds$ylabel,
    out_file       = file.path(ANNO_DIR, paste0(ds$prefix, "_5bars.pdf")),
    plot_width     = 9
  )

  # ── 7-bar version
  build_match_plot(
    var_labels     = VAR_COLS_7,
    bar_colors_vec = unname(BAR_COLORS_7),
    all_n          = unname(ca_all),
    dbsnp_n        = unname(ca_dbsnp),
    title_str      = title_str,
    y_label        = ds$ylabel,
    out_file       = file.path(ANNO_DIR, paste0(ds$prefix, "_7bars.pdf")),
    plot_width     = 10.5
  )
}
