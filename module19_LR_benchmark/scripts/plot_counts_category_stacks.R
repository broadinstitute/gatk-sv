#!/usr/bin/env Rscript
# Bar plot: 5 variant classes, bars split into base / ME / Duplication stacks
library(ggplot2)
library(scales)

ANNO_DIR <- "/Users/xzhao/Library/CloudStorage/GoogleDrive-zxf0419@gmail.com/My Drive/Talkowski_Lab/gnomAD_LR/annotation"

VAR_COLS_7 <- c("SNV", "INS 1-49bp", "DEL 1-49bp",
                "INS 50-499bp", "DEL 50-499bp",
                "INS >499bp", "DEL >499bp")

# Bar colours (same as existing palette)
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

# Overlay fill colours (semi-transparent)
ME_FILL    <- "#AADDFF"   # light blue  (ME/All)
DUP_FILL   <- "#FFDD88"   # light amber (Duplication/All)
TR_FILL    <- "#99DD99"   # light green (All/TR)
METR_FILL  <- "#2266BB"   # mid blue    (ME/TR)
DUPTR_FILL <- "#CC7700"   # dark amber  (Duplication/TR)

ME_ALPHA    <- 0.50
DUP_ALPHA   <- 0.50
TR_ALPHA    <- 0.55
METR_ALPHA  <- 0.65
DUPTR_ALPHA <- 0.65

# Helper: format numbers
fmt_count <- function(x) {
  ifelse(x >= 1e9, sprintf("%.1fB", x / 1e9),
  ifelse(x >= 1e6, sprintf("%.1fM", x / 1e6),
  ifelse(x >= 1e3, sprintf("%.1fK", x / 1e3),
         sprintf("%.0f", x))))
}

# Helper: merge INS/DEL 50+ into >=50 category
merge5 <- function(ca7) {
  c(
    ca7["SNV"],
    ca7["INS 1-49bp"],
    ca7["DEL 1-49bp"],
    ca7["INS 50-499bp"] + ca7["INS >499bp"],
    ca7["DEL 50-499bp"] + ca7["DEL >499bp"]
  )
}

# Extract a named count vector for a given category row
extract_counts <- function(df, cat1, cat2 = "All", cat3 = "All", cat4 = "All") {
  row <- df[df$category    == cat1 &
            df$sub_category == cat2 &
            df$tr_status    == cat3 &
            df$region       == cat4, ]
  # exclude "(prop)" rows
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
      data.frame(x    = xs, xend = xe,
                 y    = 10^(hb + max(lys, 0)),
                 yend = 10^(hb + min(lye, dh)),
                 stringsAsFactors = FALSE)
    })
    do.call(rbind, Filter(Negate(is.null), segs))
  })
  do.call(rbind, Filter(Negate(is.null), rows))
}

# Shared helper: build stacked overlay data frames and plot layers
# Counts vectors must be the same length; var_labels names the bars.
build_stack_plot <- function(var_labels, bar_colors_vec,
                             all_n, me_n, dup_n,
                             all_tr_n, me_tr_n, dup_tr_n,
                             title_str, y_label,
                             out_file, plot_width) {

  all_n  <- pmax(all_n,  1)

  prop_me     <- me_n     / all_n
  prop_dup    <- dup_n    / all_n
  prop_tr     <- all_tr_n / all_n
  prop_me_tr  <- me_tr_n  / all_n
  prop_dup_tr <- dup_tr_n / all_n
  for (v in c("prop_me","prop_dup","prop_tr","prop_me_tr","prop_dup_tr")) {
    x <- get(v); x[is.nan(x)] <- 0; assign(v, x)
  }

  cum1 <- prop_me
  cum2 <- cum1 + prop_dup
  cum3 <- cum2 + prop_tr
  cum4 <- cum3 + prop_me_tr
  cum5 <- cum4 + prop_dup_tr

  log_all <- log10(all_n)
  t0  <- rep(1, length(all_n))
  t1  <- 10^(log_all * cum1)
  t2  <- 10^(log_all * cum2)
  t3  <- 10^(log_all * cum3)
  t4  <- 10^(log_all * cum4)
  t5  <- 10^(log_all * cum5)

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

  # Solid overlay data frames for ME and Dup
  solid_ovl <- function(ymin, ymax, prop) {
    data.frame(
      xmin    = x_pos - half_w, xmax = x_pos + half_w,
      ymin    = ymin,           ymax = ymax,
      pct_lbl = paste0(round(prop * 100, 1), "%"),
      lbl_y   = ifelse(ymax > ymin * 1.2, sqrt(ymin * ymax), NA_real_),
      stringsAsFactors = FALSE
    )
  }

  me_df  <- solid_ovl(t0, t1, prop_me)
  dup_df <- solid_ovl(t1, t2, prop_dup)

  # Slash strips for TR layers
  tr_slash    <- make_slash_strip(x_pos, half_w, t2, t3)
  me_tr_slash <- make_slash_strip(x_pos, half_w, t3, t4)
  dtr_slash   <- make_slash_strip(x_pos, half_w, t4, t5)

  # Label for All-TR: % of total, shown just above the strip top
  tr_lbl_df <- data.frame(
    x     = x_pos,
    y     = t3 * 1.25,
    label = paste0(round(prop_tr * 100, 1), "%"),
    stringsAsFactors = FALSE
  )
  tr_lbl_df <- tr_lbl_df[t3 > t2 * 1.15, ]

  # Legend positions
  leg_x0 <- n + 0.65;  leg_x1 <- n + 0.9
  mx <- max(all_n)
  leg_tops <- c(mx*0.60, mx*0.30, mx*0.15, mx*0.075, mx*0.035)
  leg_bots <- c(mx*0.30, mx*0.15, mx*0.075, mx*0.035, mx*0.015)
  leg_colors <- c(ME_FILL, DUP_FILL, TR_FILL, METR_FILL, DUPTR_FILL)
  leg_labels <- c("ME", "Duplication", "All TR", "ME (TR)", "Dup (TR)")

  p <- ggplot() +
    geom_col(data = bar_df, aes(x = var_class, y = count_all),
             fill = bar_df$fill_col, color = NA, width = 2 * half_w) +
    # ME solid overlay
    geom_rect(data = me_df,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = ME_FILL, alpha = ME_ALPHA, color = NA, inherit.aes = FALSE) +
    geom_text(data = me_df[!is.na(me_df$lbl_y), ],
              aes(x = (xmin + xmax) / 2, y = lbl_y, label = pct_lbl),
              size = 2.7, color = "black", fontface = "bold", inherit.aes = FALSE) +
    # Dup solid overlay
    geom_rect(data = dup_df,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = DUP_FILL, alpha = DUP_ALPHA, color = NA, inherit.aes = FALSE) +
    geom_text(data = dup_df[!is.na(dup_df$lbl_y), ],
              aes(x = (xmin + xmax) / 2, y = lbl_y, label = pct_lbl),
              size = 2.7, color = "black", fontface = "bold", inherit.aes = FALSE) +
    # Count labels above bars
    geom_text(data = bar_df, aes(x = var_class, y = label_y, label = count_label),
              size = 3.2, fontface = "bold", color = "#333333")

  # TR slash strips
  if (!is.null(tr_slash) && nrow(tr_slash) > 0)
    p <- p + geom_segment(data = tr_slash,
                          aes(x = x, xend = xend, y = y, yend = yend),
                          color = TR_FILL, linewidth = 0.55, alpha = 0.9,
                          inherit.aes = FALSE)
  if (!is.null(me_tr_slash) && nrow(me_tr_slash) > 0)
    p <- p + geom_segment(data = me_tr_slash,
                          aes(x = x, xend = xend, y = y, yend = yend),
                          color = METR_FILL, linewidth = 0.55, alpha = 0.9,
                          inherit.aes = FALSE)
  if (!is.null(dtr_slash) && nrow(dtr_slash) > 0)
    p <- p + geom_segment(data = dtr_slash,
                          aes(x = x, xend = xend, y = y, yend = yend),
                          color = DUPTR_FILL, linewidth = 0.55, alpha = 0.9,
                          inherit.aes = FALSE)

  # All-TR % label just above strip
  if (nrow(tr_lbl_df) > 0)
    p <- p + geom_text(data = tr_lbl_df,
                       aes(x = x, y = y, label = label),
                       size = 2.5, color = TR_FILL, fontface = "bold",
                       vjust = 0, inherit.aes = FALSE)

  p <- p +
    scale_y_log10(
      labels = label_number(scale_cut = cut_short_scale()),
      expand = expansion(mult = c(0, 0)),
      limits = c(1, max(all_n) * 10)
    ) +
    labs(x = "Variant class", y = y_label, title = title_str,
         subtitle = paste0(
           "Solid: ME (blue) | Dup (amber)    Slash strips: All-TR | ME-TR | Dup-TR\n",
           "Heights in log10-space = log10(N) x cumulative proportion")) +
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

  # Manual legend: solid rects for ME/Dup; outlined box + slashes for TR
  p <- p +
    annotate("rect", xmin = leg_x0, xmax = leg_x1,
             ymin = leg_bots[1], ymax = leg_tops[1],
             fill = ME_FILL,  alpha = ME_ALPHA,  color = NA) +
    annotate("text", x = leg_x1 + 0.08, y = sqrt(leg_bots[1] * leg_tops[1]),
             label = "ME",          hjust = 0, size = 3.2, fontface = "bold") +
    annotate("rect", xmin = leg_x0, xmax = leg_x1,
             ymin = leg_bots[2], ymax = leg_tops[2],
             fill = DUP_FILL, alpha = DUP_ALPHA, color = NA) +
    annotate("text", x = leg_x1 + 0.08, y = sqrt(leg_bots[2] * leg_tops[2]),
             label = "Duplication", hjust = 0, size = 3.2, fontface = "bold")

  for (ii in 1:3) {
    lcol <- c(TR_FILL, METR_FILL, DUPTR_FILL)[ii]
    llbl <- leg_labels[ii + 2]
    yb   <- leg_bots[ii + 2];  yt <- leg_tops[ii + 2]
    # outlined box
    p <- p +
      annotate("rect", xmin = leg_x0, xmax = leg_x1, ymin = yb, ymax = yt,
               fill = NA, color = lcol, linewidth = 0.5) +
      annotate("text", x = leg_x1 + 0.08, y = sqrt(yb * yt),
               label = llbl, hjust = 0, size = 3.2, fontface = "bold",
               color = lcol)
    # draw mini slash segments inside the legend box
    sldf <- make_slash_strip(c((leg_x0 + leg_x1) / 2), (leg_x1 - leg_x0) / 2,
                             c(yb), c(yt), n_lines = 5)
    if (!is.null(sldf) && nrow(sldf) > 0)
      p <- p + annotate("segment",
                        x = sldf$x, xend = sldf$xend,
                        y = sldf$y, yend = sldf$yend,
                        color = lcol, linewidth = 0.5)
  }

  p <- p + coord_cartesian(xlim = c(0.5, n + 2.5), clip = "off")

  ggsave(out_file, plot = p, width = plot_width, height = 6.5, device = "pdf")
  cat("Saved:", out_file, "\n")
}

# ── 5-bar version ────────────────────────────────────────────────────────────
make_category_stack_plot <- function(data_file, out_file, title_str, y_label) {
  df <- read.table(data_file, header = TRUE, sep = "\t", check.names = FALSE)
  m  <- function(x) unname(merge5(x))
  build_stack_plot(
    var_labels    = VAR_COLS_5,
    bar_colors_vec= unname(BAR_COLORS_5),
    all_n         = m(extract_counts(df, "All")),
    me_n          = m(extract_counts(df, "ME")),
    dup_n         = m(extract_counts(df, "Duplication")),
    all_tr_n      = m(extract_counts(df, "All",         cat3 = "TR")),
    me_tr_n       = m(extract_counts(df, "ME",          cat3 = "TR")),
    dup_tr_n      = m(extract_counts(df, "Duplication", cat3 = "TR")),
    title_str     = title_str,
    y_label       = y_label,
    out_file      = out_file,
    plot_width    = 9
  )
}

# ── 7-bar version ────────────────────────────────────────────────────────────
make_category_stack_plot_7 <- function(data_file, out_file, title_str, y_label) {
  df <- read.table(data_file, header = TRUE, sep = "\t", check.names = FALSE)
  g  <- function(cat1, cat3 = "All") unname(extract_counts(df, cat1, cat3 = cat3)[VAR_COLS_7])
  build_stack_plot(
    var_labels    = VAR_COLS_7,
    bar_colors_vec= unname(BAR_COLORS_7),
    all_n         = g("All"),
    me_n          = g("ME"),
    dup_n         = g("Duplication"),
    all_tr_n      = g("All",         cat3 = "TR"),
    me_tr_n       = g("ME",          cat3 = "TR"),
    dup_tr_n      = g("Duplication", cat3 = "TR"),
    title_str     = title_str,
    y_label       = y_label,
    out_file      = out_file,
    plot_width    = 10
  )
}

# ── Run for each dataset ──────────────────────────────────────────────────────
datasets <- list(
  list(file   = "aou_I.counts_sites.tsv",
       prefix = "aou_I.counts_sites.barplot_categories",
       label  = "AoU variant site counts",
       ylabel = "Number of variant sites"),
  list(file   = "aou_I.counts_samples.tsv",
       prefix = "aou_I.counts_samples.barplot_categories",
       label  = "AoU per-sample variant counts",
       ylabel = "Mean variant count per sample"),
  list(file   = "hgsv_hprc.counts_sites.tsv",
       prefix = "hgsv_hprc.counts_sites.barplot_categories",
       label  = "HGSV/HPRC variant site counts",
       ylabel = "Number of variant sites"),
  list(file   = "hgsv_hprc.counts_samples.tsv",
       prefix = "hgsv_hprc.counts_samples.barplot_categories",
       label  = "HGSV/HPRC per-sample variant counts",
       ylabel = "Mean variant count per sample")
)

for (ds in datasets) {
  cat("\n── Processing:", ds$file, "\n")
  df_tmp <- read.table(file.path(ANNO_DIR, ds$file),
                       header = TRUE, sep = "\t", check.names = FALSE)
  ca_all <- extract_counts(df_tmp, "All")
  total_N <- sum(unname(merge5(ca_all)))
  title_str <- sprintf("%s (N = %s)", ds$label, fmt_count(total_N))

  make_category_stack_plot(
    data_file = file.path(ANNO_DIR, ds$file),
    out_file  = file.path(ANNO_DIR, paste0(ds$prefix, "_5bars.pdf")),
    title_str = title_str,
    y_label   = ds$ylabel
  )

  total_N7 <- sum(unname(ca_all[VAR_COLS_7]))
  title_str7 <- sprintf("%s (N = %s)", ds$label, fmt_count(total_N7))
  make_category_stack_plot_7(
    data_file = file.path(ANNO_DIR, ds$file),
    out_file  = file.path(ANNO_DIR, paste0(ds$prefix, "_7bars.pdf")),
    title_str = title_str7,
    y_label   = ds$ylabel
  )
}
