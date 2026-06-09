library(ggplot2)

setwd("/home/carlos/VizMut-lolliplot")

bench_df <- read.csv("benchmarks/benchmark_results.csv")

p <- ggplot(bench_df, aes(x = n_variants, y = time_mean)) +
  geom_ribbon(aes(ymin = time_mean - time_sd,
                  ymax = time_mean + time_sd),
              fill = "#4A90D9", alpha = 0.2) +
  geom_line(color = "#2C5F8A", linewidth = 1) +
  geom_point(color = "#2C5F8A", size = 3) +
  geom_text(aes(label = paste0(round(time_mean, 1), "s")),
            vjust = -1.2, size = 3.5, color = "gray30") +
  scale_x_continuous(breaks = bench_df$n_variants) +
  scale_y_continuous(limits = c(0, max(bench_df$time_mean + bench_df$time_sd) * 1.2),
                     labels = function(x) paste0(x, "s")) +
  labs(
    title    = "VizMut-lolliplot — Enrichment mode performance",
    subtitle = "Time includes ClinVar, NCBI Variation Services and gnomAD v4 API calls\nMean ± SD of 3 replicates",
    x        = "Number of variants",
    y        = "Elapsed time (seconds)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(color = "gray50", size = 10),
    panel.grid.minor = element_blank()
  )

ggsave("benchmarks/benchmark_plot.png", p,
       width = 8, height = 5, dpi = 150)

message("Plot guardado en benchmarks/benchmark_plot.png")


