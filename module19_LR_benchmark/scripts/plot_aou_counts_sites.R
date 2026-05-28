#!/usr/bin/env Rscript
library(ggplot2)
library(scales)

ANNO_DIR <- "/Users/xzhao/Library/CloudStorage/GoogleDrive-zxf0419@gmail.com/My Drive/Talkowski_Lab/gnomAD_LR/annotation"

VAR_COLS_7 <- c("SNV", "INS 1-49bp", "DEL 1-49bp",
                "INS 50-499bp", "DEL 50-499bp",
                "INS >499bp", "DEL >499bp")

COLORS_7 <- c(
  "SNV"          = "#7B2D8B",
  "INS 1-49bp"   = "#FFB6C1",
  "DEL 1-49bp"   = "#CC0000",
  "INS 50-499bp" = "#FF69B4",
  "DEL 50-499bp" = "#FF4444",
  "INS >499bp"   = "#FF1493",
  "DEL >499bp"   = "#FF6666"
)

COLORS_5 <- c(
  "SNV"        = "#7B2D8B",
  "INS 1-49bp" = "#FFB6C1",
  "DEL 1-49bp" = "#CC0000",
  "INS >=50bp" = "#FF69B4",
  "DEL >=50bp" = "#FF4444"
)

# ── Helpers ──────────────────────────────────────────────────────────────────
fmt_count <- function(x) {
  ifelse(x >= 1e9, sprintf("%.1fB", x / 1e9),
  ifelse(x >= 1e6, sprintf("%.1fM", x / 1e6),
  ifelse(x >= 1e3, sprintf("%.1fK", x / 1e3),
         sprintf("%.0f", x))))
}

make_slash_df <- function(x_pos, half_w, y_tops, n_lines = 10) {
  rows <- lapply(seq_along(x_pos), function(i) {
    yt <- y_tops[i]
    if (is.na(yt) || yt <= 1) return(NULL)
    xl <- x_pos[i] - half_w;  xr <- x_pos[i] + half_w
    w  <- xr - xl;            hy <- log10(yt);  slope <- hy / w
    offsets <- seq(xl - w, xr, length.out = n_lines + 2)
    offsets <- offsets[-c(1, length(offsets))]
    segs <- lapply(offsets, function(x0) {
      xs <- max(xl, x0);  xe <- min(xr, x0 + w)
      lys <- (xs - x0) * slope;  lye <- (xe - x0) * slope
      if (lys < 0)  { xs <- x0;              lys <- 0  }
      if (lye > hy) { xe <- x0 + hy / slope; lye <- hy }
      xs <- max(xs, xl);  xe <- min(xe, xr)
      if (is.na(xs) || is.na(xe) || xs >= xe) return(NULL)
      data.frame(x = xs, xend = xe,
                 y = 10^max(lys, 0), yend = 10^min(lye, hy),
                 stringsAsFactors = FALSE)
    })
    do.call(rbind, Filter(Negate(is.null), segs))
  })
  do.call(rbind, Filter(Negate(is.null), rows))
}

# ── V1: 7-bar plot ──────────────────────────────────────────────────────────
save_v1 <- function(counts_all, counts_tr, title_str, y_label, out_file) {
  var_cols    <- VAR_COLS_7
  bar_colors  <- COLORS_7
  proportions <- counts_tr / counts_all
  strip_top   <- 10^(log10(counts_all) * proportions)
  x_pos       <- seq_len(length(var_cols))
  half_w      <- 0.4

  bar_df <- data.frame(
    var_class   = factor(var_cols, levels = var_cols),
    count_all   = counts_all,
    count_label = fmt_count(counts_all),
    label_y     = counts_all * 2.5,
    stringsAsFactors = FALSE
  )
  strip_df <- data.frame(
    xmin      = x_pos - half_w,
    xmax      = x_pos + half_w,
    ymin      = 1,
    ymax      = strip_top,
    pct_label = paste0(round(proportions * 100, 1), "%"),
    label_y   = strip_top * 1.8,
    stringsAsFactors = FALSE
  )
  slash_df <- make_slash_df(x_pos, half_w, strip_top)

  p <- ggplot() +
    geom_col(data = bar_df,
             aes(x = var_class, y = count_all, fill = var_class),
             color = "white", linewidth = 0.5, width = 2 * half_w) +
    geom_rect(data = strip_df,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = NA, color = "black", linewidth = 0.6, inherit.aes = FALSE) +
    geom_segment(data = slash_df,
                 aes(x = x, xend = xend, y = y, yend = yend),
                 color = "black", linewidth = 0.4, inherit.aes = FALSE) +
    geom_text(data = strip_df,
              aes(x = (xmin + xmax) / 2, y = label_y, label = pct_label),
              size = 3.2, color = "black", fontface = "bold", inherit.aes = FALSE) +
    geom_text(data = bar_df,
              aes(x = var_class, y = label_y, label = count_label),
              size = 3.5, fontface = "bold", color = "#333333") +
    scale_fill_manual(values = bar_colors, guide = "none") +
    scale_y_log10(
      labels = label_number(scale_cut = cut_short_scale()),
      expand = expansion(mult = c(0, 0)),
      limits = c(1, max(counts_all) * 10)
    ) +
    labs(x = "Variant class", y = y_label, title = title_str,
         subtitle = "Hatched strip = TR fraction  (strip height in log10-space = log10(N) x TR/N)") +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x        = element_text(angle = 25, hjust = 1, size = 11),
      axis.text.y        = element_text(size = 10),
      axis.title         = element_text(size = 12),
      plot.title         = element_text(face = "bold", size = 14),
      plot.subtitle      = element_text(size = 9, color = "grey40"),
      panel.grid.major.y = element_line(linetype = "dashed", color = "grey80",
                                        linewidth = 0.4)
    )

  ggsave(out_file, plot = p, width = 10, height = 7, device = "pdf")
  cat("Saved:", out_file, "\n")
}

# ── V2: 5-bar grouped layout ─────────────────────────────────────────────────
save_v2 <- function(counts_all_7, counts_tr_7, title_str, y_label,
                    out_file_a, out_file_b) {
  VAR_COLS_5 <- c("SNV", "INS 1-49bp", "DEL 1-49bp", "INS >=50bp", "DEL >=50bp")
  counts_all <- unname(c(
    counts_all_7["SNV"],
    counts_all_7["INS 1-49bp"],
    counts_all_7["DEL 1-49bp"],
    counts_all_7["INS 50-499bp"] + counts_all_7["INS >499bp"],
    counts_all_7["DEL 50-499bp"] + counts_all_7["DEL >499bp"]
  ))
  counts_tr <- unname(c(
    counts_tr_7["SNV"],
    counts_tr_7["INS 1-49bp"],
    counts_tr_7["DEL 1-49bp"],
    counts_tr_7["INS 50-499bp"] + counts_tr_7["INS >499bp"],
    counts_tr_7["DEL 50-499bp"] + counts_tr_7["DEL >499bp"]
  ))

  half_w5  <- 0.28
  x_pos5   <- c(1.0, 2.5, 3.1, 4.3, 4.9)
  x_snv    <- x_pos5[1]
  x_ins1   <- x_pos5[2];  x_del1 <- x_pos5[3]
  x_ins2   <- x_pos5[4];  x_del2 <- x_pos5[5]

  proportions <- counts_tr / counts_all
  strip_top   <- 10^(log10(counts_all) * proportions)

  bar_df <- data.frame(
    x           = x_pos5,
    fill_col    = unname(COLORS_5[VAR_COLS_5]),
    count_all   = counts_all,
    count_label = fmt_count(counts_all),
    label_x     = x_pos5 - half_w5 + 0.02,
    label_y     = counts_all * 1.9,
    stringsAsFactors = FALSE
  )
  strip_df <- data.frame(
    xmin      = x_pos5 - half_w5,
    xmax      = x_pos5 + half_w5,
    ymin      = 1,
    ymax      = strip_top,
    pct_label = paste0(round(proportions * 100, 1), "%"),
    label_x   = x_pos5 - half_w5 + 0.03,
    label_y   = sqrt(strip_top),
    stringsAsFactors = FALSE
  )
  slash_df <- make_slash_df(x_pos5, half_w5, strip_top)

  groups_df <- data.frame(
    label   = c("1-49bp", ">=50bp"),
    x_ctr   = c((x_ins1 + x_del1) / 2, (x_ins2 + x_del2) / 2),
    x_left  = c(x_ins1 - half_w5,       x_ins2 - half_w5),
    x_right = c(x_del1 + half_w5,       x_del2 + half_w5),
    stringsAsFactors = FALSE
  )

  v2_theme <- theme_classic(base_size = 13) +
    theme(
      axis.text.x        = element_text(size = 12, face = "bold"),
      axis.text.y        = element_text(size = 11),
      axis.title.x       = element_blank(),
      axis.title.y       = element_text(size = 13, face = "bold"),
      plot.title         = element_text(face = "bold", size = 15, hjust = 0.5,
                                        lineheight = 1.2),
      panel.grid.major.y = element_line(linetype = "dotted", color = "grey65",
                                        linewidth = 0.45),
      plot.margin        = margin(t = 10, r = 20, b = 55, l = 10)
    )

  build <- function(add_strips) {
    p <- ggplot() +
      geom_col(data = bar_df, aes(x = x, y = count_all),
               fill = bar_df$fill_col, color = NA, width = 2 * half_w5)
    if (add_strips) {
      p <- p +
        geom_rect(data = strip_df,
                  aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                  fill = NA, color = "black", linewidth = 0.6, inherit.aes = FALSE) +
        geom_segment(data = slash_df,
                     aes(x = x, xend = xend, y = y, yend = yend),
                     color = "black", linewidth = 0.35, inherit.aes = FALSE) +
        geom_text(data = strip_df,
                  aes(x = label_x, y = label_y, label = pct_label),
                  hjust = 0, vjust = 0.5, size = 3.3, fontface = "bold",
                  color = "black", inherit.aes = FALSE)
    }
    p +
      geom_text(data = bar_df,
                aes(x = label_x, y = label_y, label = count_label),
                hjust = 0, vjust = 0, size = 4, fontface = "bold",
                color = "grey25") +
      geom_segment(data = groups_df,
                   aes(x = x_left, xend = x_right, y = 0.42, yend = 0.42),
                   color = "black", linewidth = 0.55, inherit.aes = FALSE) +
      geom_text(data = groups_df,
                aes(x = x_ctr, y = 0.15, label = label),
                hjust = 0.5, vjust = 1, size = 3.8, color = "black",
                inherit.aes = FALSE) +
      scale_x_continuous(limits = c(0.5, 5.45),
                         breaks = x_pos5,
                         labels = c("SNV", "INS", "DEL", "INS", "DEL"),
                         expand = expansion(mult = 0)) +
      scale_y_log10(labels = label_number(scale_cut = cut_short_scale()),
                    expand = expansion(mult = c(0, 0)),
                    limits = c(1, max(counts_all) * 8)) +
      coord_cartesian(clip = "off") +
      labs(y = y_label, title = title_str) +
      v2_theme
  }

  ggsave(out_file_a, plot = build(FALSE), width = 6, height = 6, device = "pdf")
  cat("Saved:", out_file_a, "\n")
  ggsave(out_file_b, plot = build(TRUE),  width = 6, height = 6, device = "pdf")
  cat("Saved:", out_file_b, "\n")
}

# ── Main: loop over datasets ─────────────────────────────────────────────────
datasets <- list(
  list(
    file     = "aou_I.counts_sites.tsv",
    prefix   = "aou_I.counts_sites.barplot",
    dataset  = "AoU",
    is_sites = TRUE
  ),
  list(
    file     = "aou_I.counts_samples.tsv",
    prefix   = "aou_I.counts_samples.barplot",
    dataset  = "AoU",
    is_sites = FALSE
  ),
  list(
    file     = "hgsv_hprc.counts_sites.tsv",
    prefix   = "hgsv_hprc.counts_sites.barplot",
    dataset  = "HGSV/HPRC",
    is_sites = TRUE
  ),
  list(
    file     = "hgsv_hprc.counts_samples.tsv",
    prefix   = "hgsv_hprc.counts_samples.barplot",
    dataset  = "HGSV/HPRC",
    is_sites = FALSE
  )
)

for (ds in datasets) {
  cat("\n── Processing:", ds$file, "\n")
  df <- read.table(file.path(ANNO_DIR, ds$file),
                   header = TRUE, sep = "\t", check.names = FALSE)

  row_all <- df[df$category    == "All" &
                df$sub_category == "All" &
                df$tr_status    == "All" &
                df$region       == "All", ][1, ]
  row_tr  <- df[df$category    == "All" &
                df$sub_category == "All" &
                df$tr_status    == "TR"  &
                df$region       == "All", ][1, ]

  ca7 <- setNames(as.numeric(row_all[, VAR_COLS_7]), VAR_COLS_7)
  ct7 <- setNames(as.numeric(row_tr[,  VAR_COLS_7]), VAR_COLS_7)

  total_N <- sum(ca7)
  y_label <- if (ds$is_sites) "Count of variant sites" else "Mean variant count per sample"
  kind    <- if (ds$is_sites) "variant site counts" else "per-sample variant counts"
  title_7 <- sprintf("%s %s (N = %s)", ds$dataset, kind, fmt_count(total_N))
  title_5 <- sprintf("%s %s\n(N = %s)",  ds$dataset, kind, fmt_count(total_N))

  save_v1(unname(ca7), unname(ct7), title_7, y_label,
          file.path(ANNO_DIR, paste0(ds$prefix, "_v1_7bars.pdf")))

  save_v2(ca7, ct7, title_5, y_label,
          file.path(ANNO_DIR, paste0(ds$prefix, "_v2a_5bars.pdf")),
          file.path(ANNO_DIR, paste0(ds$prefix, "_v2b_5bars.pdf")))
}
