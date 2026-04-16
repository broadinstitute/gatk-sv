library(ggplot2)

plot_data <- data.frame(
  label = c(
    "HiPhase",
    "SHAPEIT",
    "BackbonePhase",
    "HiPhase",
    "SHAPEIT",
    "BackbonePhase",
    "BackbonePhase, Normalized",
    "BackbonePhase, len(REF) < 100",
    "BackbonePhase, abs(len(REF) - len(ALT)) < 3",
    "BackbonePhase, len(MOTIF) > 2",
    "BackbonePhase, len(MOTIF) > 1"
  ),
  concordance = c(56.67, 61.88, 96.69, 53.83, 56.65, 84.85, 84.90, 86.73, 94.82, 95.39, 96.58),
  group = c(rep("No TRGT", 3), rep("TRGT", 8)),
  x = c(1, 2, 3, 5.3, 6.3, 7.5, 8.5, 9.7, 10.9, 12.1, 13.3)
)

p <- ggplot(plot_data, aes(x = x, y = concordance, fill = group)) +
  geom_hline(
    yintercept = c(20, 40, 60, 80, 100),
    color = "grey82",
    linewidth = 0.7
  ) +
  geom_hline(
    yintercept = c(10, 30, 50, 70, 90),
    color = "grey70",
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_col(width = 0.8) +
  geom_text(
    aes(label = sprintf("%.2f%%", concordance)),
    vjust = -0.35,
    fontface = "bold",
    size = 7
  ) +
  scale_fill_manual(values = c("No TRGT" = "#d92525", "TRGT" = "#2e7db6")) +
  scale_x_continuous(breaks = plot_data$x, labels = plot_data$label) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 20), expand = c(0, 0)) +
  labs(x = NULL, y = "Concordance Rate (%)", fill = NULL) +
  theme_classic(base_size = 22) +
  theme(
    legend.position = c(0.02, 0.92),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = grDevices::adjustcolor("white", alpha.f = 0.92), color = "grey80"),
    axis.text.x = element_text(angle = 40, hjust = 1, size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.y = element_text(size = 22),
    legend.text = element_text(size = 18),
    panel.grid = element_blank()
  )

ggsave("/Users/xzhao/Downloads/plot.png", p, width = 16, height = 10, dpi = 150)