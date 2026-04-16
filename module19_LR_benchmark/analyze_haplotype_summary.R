library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# ── Paths ─────────────────────────────────────────────────────────────────────
base_dir <- "/Users/xzhao/zxf0419@gmail.com - Google Drive/My Drive/Talkowski_Lab/gnomAD_LR/final_callsets/haplotypes"
out_dir  <- file.path(base_dir, "figures")
dir.create(out_dir, showWarnings = FALSE)

files <- list(
  "All variants"  = file.path(base_dir, "hgsv_hprc.with_gene_anno.haplotypes.summary.tsv"),
  "AF > 1%"       = file.path(base_dir, "hgsv_hprc.with_gene_anno.haplotypes_common_AF_1%.summary.tsv"),
  "AF > 10%"      = file.path(base_dir, "hgsv_hprc.with_gene_anno.haplotypes_common_AF_10%.summary.tsv")
)

dat <- bind_rows(lapply(names(files), function(label) {
  read.table(files[[label]], header = TRUE, sep = "\t", stringsAsFactors = FALSE) |>
    mutate(af_filter = label)
})) |>
  mutate(af_filter = factor(af_filter, levels = c("All variants", "AF > 1%", "AF > 10%")))

cat("Rows per group:\n")
print(table(dat$af_filter))

palette <- c("All variants" = "#4C72B0", "AF > 1%" = "#DD8452", "AF > 10%" = "#55A868")

save_plot <- function(p, name, w = 8, h = 5) {
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 150)
  message("Saved: ", name)
}

# ── 1. Variants per gene distribution (density) ───────────────────────────────
p1 <- ggplot(dat, aes(x = n_variants, colour = af_filter, fill = af_filter)) +
  geom_density(alpha = 0.25) +
  scale_x_log10(labels = comma) +
  scale_colour_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(title = "Distribution of variants per gene",
       x = "Number of variants (log10)", y = "Density",
       colour = "AF filter", fill = "AF filter") +
  theme_bw(base_size = 13)
save_plot(p1, "1_variants_per_gene_density")

# ── 2. Haplotype patterns per gene distribution ───────────────────────────────
p2 <- ggplot(dat, aes(x = n_patterns, colour = af_filter, fill = af_filter)) +
  geom_density(alpha = 0.25) +
  scale_x_log10(labels = comma) +
  scale_colour_manual(values = palette) +
  scale_fill_manual(values = palette) +
  labs(title = "Distribution of unique haplotype patterns per gene",
       x = "Number of patterns (log10)", y = "Density",
       colour = "AF filter", fill = "AF filter") +
  theme_bw(base_size = 13)
save_plot(p2, "2_patterns_per_gene_density")

# ── 3. Variants vs Patterns scatter ──────────────────────────────────────────
p3 <- dat |>
  filter(n_variants > 0, n_patterns > 0) |>
  ggplot(aes(x = n_variants, y = n_patterns, colour = af_filter)) +
  geom_point(alpha = 0.15, size = 0.6) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  scale_x_log10(labels = comma) +
  scale_y_log10(labels = comma) +
  scale_colour_manual(values = palette) +
  facet_wrap(~af_filter) +
  labs(title = "Variants vs haplotype patterns per gene",
       x = "n_variants (log10)", y = "n_patterns (log10)",
       colour = "AF filter") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")
save_plot(p3, "3_variants_vs_patterns_scatter", w = 12, h = 4)

# ── 4. Pattern class composition (stacked bar, mean per group) ─────────────────
pattern_class <- dat |>
  group_by(af_filter) |>
  summarise(
    Singletons  = mean(singletons,        na.rm = TRUE),
    Doubletons  = mean(doubletons,        na.rm = TRUE),
    Common      = mean(common_patterns,   na.rm = TRUE),
    .groups = "drop"
  ) |>
  pivot_longer(cols = c(Singletons, Doubletons, Common),
               names_to = "class", values_to = "mean_count") |>
  mutate(class = factor(class, levels = c("Common", "Doubletons", "Singletons")))

p4 <- ggplot(pattern_class, aes(x = af_filter, y = mean_count, fill = class)) +
  geom_col() +
  scale_fill_manual(values = c(Singletons = "#E15759", Doubletons = "#F28E2B", Common = "#76B7B2")) +
  labs(title = "Mean pattern class composition per gene",
       x = "AF filter", y = "Mean count per gene", fill = "Pattern class") +
  theme_bw(base_size = 13)
save_plot(p4, "4_pattern_class_composition")

# ── 5. Max pattern AC boxplot ─────────────────────────────────────────────────
p5 <- ggplot(dat, aes(x = af_filter, y = max_pattern_ac, fill = af_filter)) +
  geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.3) +
  scale_y_log10(labels = comma) +
  scale_fill_manual(values = palette) +
  labs(title = "Max haplotype pattern allele count per gene",
       x = "AF filter", y = "Max pattern AC (log10)") +
  theme_bw(base_size = 13) +
  theme(legend.position = "none")
save_plot(p5, "5_max_pattern_AC_boxplot")

# ── 6. Gene size vs n_variants ────────────────────────────────────────────────
p6 <- dat |>
  filter(gene_size > 0, n_variants > 0) |>
  ggplot(aes(x = gene_size, y = n_variants, colour = af_filter)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  scale_x_log10(labels = comma) +
  scale_y_log10(labels = comma) +
  scale_colour_manual(values = palette) +
  labs(title = "Gene size vs number of variants",
       x = "Gene size (bp, log10)", y = "n_variants (log10)",
       colour = "AF filter") +
  theme_bw(base_size = 13)
save_plot(p6, "6_gene_size_vs_variants", w = 8, h = 5)

# ── 7. Summary stats table ────────────────────────────────────────────────────
summary_tbl <- dat |>
  group_by(af_filter) |>
  summarise(
    n_genes        = n(),
    median_variants  = median(n_variants,    na.rm = TRUE),
    median_patterns  = median(n_patterns,    na.rm = TRUE),
    median_singletons = median(singletons,   na.rm = TRUE),
    pct_common_gt0 = mean(common_patterns > 0, na.rm = TRUE) * 100,
    .groups = "drop"
  )
cat("\n── Summary table ──\n")
print(summary_tbl)
write.table(summary_tbl, file.path(out_dir, "summary_stats.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("\nAll figures saved to: ", out_dir)
