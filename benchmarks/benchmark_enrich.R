#' Benchmark del modo --enrich de VizMut-lolliplot
#' Mide el tiempo de enriquecimiento para diferentes tamaños de input

library(dplyr)

setwd("/home/carlos/VizMut-lolliplot")
source("R/enrich_variants.R")

# cargar variantes base
variants_base <- read.csv("data/variants.csv")
cat("Total variantes disponibles:", nrow(variants_base), "\n\n")

# tamaños a evaluar
sizes <- c(5, 10, 25, 50, 68)

results <- list()

for (n in sizes) {
  cat("Benchmarking n =", n, "variantes...\n")

  # subset de n variantes unicas
  variants_sub <- variants_base %>%
    distinct(hgvs_c, .keep_all = TRUE) %>%
    head(n)

  # medir tiempo
  times <- numeric(3)  # 3 repeticiones
  for (rep in 1:3) {
    t_start <- proc.time()
    tryCatch(
      enrich_variants(variants_sub),
      error = function(e) NULL
    )
    t_end <- proc.time()
    times[rep] <- (t_end - t_start)["elapsed"]
    Sys.sleep(2)  # evitar rate limiting
  }

  results[[as.character(n)]] <- data.frame(
    n_variants  = n,
    time_mean   = mean(times),
    time_sd     = sd(times),
    time_min    = min(times),
    time_max    = max(times)
  )

  cat("  Tiempo medio:", round(mean(times), 1), "s\n\n")
}

bench_df <- do.call(rbind, results)
cat("\n=== Resultados del benchmark ===\n")
print(bench_df)

# guardar resultados
write.csv(bench_df, "benchmarks/benchmark_results.csv", row.names=FALSE)
cat("\nResultados guardados en benchmarks/benchmark_results.csv\n")
