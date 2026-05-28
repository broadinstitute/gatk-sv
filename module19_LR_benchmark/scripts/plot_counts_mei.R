#!/usr/bin/env Rscript
# 5-bar plot with MEI (ME/All/All/All) proportion shade overlay
# Aspect ratio 6:4 (h:w), narrow bars
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

MEI_FILL  <- "#AADDFF"   # light blue
MEI_ALPHA <- 0.65

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

# ── Plot function ─────────────────────────────────────────────────────────────
build_mei_plot <- function(var_labels, bar_colors_vec,
                           all_n, me_n,
                           title_str, y_label, out_file) {

  all_n    <- pmax(all_n, 1)
  prop_me  <- me_n / all_n
  prop_me[is.nan(prop_me)] <- 0

  log_all  <- log10(all_n)
  top_me   <- 10^(log_all * prop_me)

  n      <- length(var_labels)
  x_pos  <- seq_len(n)
  half_w <- 0.28   # narrower bars

  bar_df <- data.frame(
    var_class   = factor(var_labels, levels = var_labels),
    fill_col    = bar_colors_vec,
    count_all   = all_n,
    count_label = fmt_count(all_n),
    label_y     = all_n * 2.2,
    stringsAsFactors = FALSE
  )

  me_df <- data.frame(
    xmin    = x_pos - half_w,
    xmax    = x_pos + half_w,
    ymin    = 1,
    ymax    = top_me,
    pct_lbl = paste0(round(prop_me * 100, 1), "%"),
    lbl_y   = ifelse(top_me > 1.3, sqrt(top_me), NA_real_),
    stringsAsFactors = FALSE
  )

  # legend: relative to data max
  mx  <- max(all_n)
  lx0 <- n + 0.55;  lx1 <- n + 0.75
  leg_yb <- mx * 0.30;  leg_yt <- mx * 0.60

  p <- ggplot() +
    # Main bars
    geom_col(data = bar_df, aes(x = var_class, y = count_all),
             fill = bar_df$fill_col, color = NA, width = 2 * half_w) +
    # MEI shade overlay
    geom_rect(data = me_df,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = MEI_FILL, alpha = MEI_ALPHA, color = NA,
              inherit.aes = FALSE) +
    # % label inside the shaded region
    geom_text(data = me_df[!is.na(me_df$lbl_y), ],
              aes(x = (xmin + xmax) / 2, y = lbl_y, label = pct_lbl),
              size = 2.6, color = "#114477", fontface = "bold",
              inherit.aes = FALSE) +
    # Count labels above each bar
    geom_text(data = bar_df,
              aes(x = var_class, y = label_y, label = count_label),
              size = 3.0, fontface = "bold", color = "#333333") +
    scale_y_log10(
      labels = label_number(scale_cut = cut_short_scale()),
      expand = expansion(mult = c(0, 0)),
      limits = c(1, max(all_n) * 9)
    ) +
    labs(
      x        = "Variant class",
      y        = y_label,
      title    = title_str,
      subtitle = "Shaded region = MEI proportion; height in log10-space = log10(N) x proportion"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x        = element_text(angle = 30, hjust = 1, size = 11),
      axis.text.y        = element_text(size = 10),
      axis.title         = element_text(size = 12),
      plot.title         = element_text(face = "bold", size = 13),
      plot.subtitle      = element_text(size = 8, color = "grey40"),
      panel.grid.major.y = element_line(linetype = "dashed", color = "grey80",
                                        linewidth = 0.4)
    ) +
    # Manual legend
    annotate("rect", xmin = lx0, xmax = lx1, ymin = leg_yb, ymax = leg_yt,
             fill = MEI_FILL, alpha = MEI_ALPHA, color = NA) +
    annotate("text", x = lx1 + 0.08, y = sqrt(leg_yb * leg_yt),
             label = "MEI", hjust = 0, size = 3.2, fontface = "bold",
             color = "#114477") +
    coord_cartesian(xlim = c(0.5, n + 1.8), clip = "off")

  # 6:4 h:w ratio → height = 1.5 × width; use width = 6.5
  ggsave(out_file, plot = p, width = 6.5, height = 9.75, device = "pdf")
  cat("Saved:", out_file, "\n")
}

# ── Dataset loop ──────────────────────────────────────────────────────────────
datasets <- list(
  list(file   = "aou_I.counts_sites.tsv",
       prefix = "aou_I.counts_sites.barplot_mei",
       label  = "AoU variant site counts",
       ylabel = "Number of variant sites"),
  list(file   = "aou_I.counts_samples.tsv",
       prefix = "aou_I.counts_samples.barplot_mei",
       label  = "AoU per-sample variant counts",
       ylabel = "Mean variant count per sample"),
  list(file   = "hgsv_hprc.counts_sites.tsv",
       prefix = "hgsv_hprc.counts_sites.barplot_mei",
       label  = "HGSV/HPRC variant site counts",
       ylabel = "Number of variant sites"),
  list(file   = "hgsv_hprc.counts_samples.tsv",
       prefix = "hgsv_hprc.counts_samples.barplot_mei",
       label  = "HGSV/HPRC per-sample variant counts",
       ylabel = "Mean variant count per sample")
)

for (ds in datasets) {
  cat("\n── Processing:", ds$file, "\n")
  df <- read.table(file.path(ANNO_DIR, ds$file),
                   header = TRUE, sep = "\t", check.names = FALSE)

  ca_all <- unname(merge5(extract_counts(df, "All")))
  ca_me  <- unname(merge5(extract_counts(df, "ME")))

  total_N   <- sum(ca_all)
  title_str <- sprintf("%s (N = %s)", ds$label, fmt_count(total_N))

  build_mei_plot(
    var_labels     = VAR_COLS_5,
    bar_colors_vec = unname(BAR_COLORS_5),
    all_n          = ca_all,
    me_n           = ca_me,
    title_str      = title_str,
    y_label        = ds$ylabel,
    out_file       = file.path(ANNO_DIR, paste0(ds$prefix, "_5bars.pdf"))
  )
}
