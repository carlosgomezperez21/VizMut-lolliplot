#' Deduplicate variants and add count column
#'
#' @param variants data.frame de variantes
#' @param by columna para agrupar (default: "hgvs_c")
#' @return data.frame deduplicado con columna count

dedup_variants <- function(variants, by = "hgvs_c") {

  library(dplyr)

  if (!by %in% names(variants)) {
    warning("Columna '", by, "' no encontrada. No se deduplicara.")
    variants$count <- 1
    return(variants)
  }

  n_before <- nrow(variants)

  # columnas a conservar (primera ocurrencia)
  variants_dedup <- variants %>%
    group_by(across(all_of(by))) %>%
    summarise(
      across(everything(), function(x) x[1]),
      count = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(count))

  n_after <- nrow(variants_dedup)
  n_removed <- n_before - n_after

  message("Deduplicacion completada:")
  message("  Variantes totales:  ", n_before)
  message("  Variantes unicas:   ", n_after)
  message("  Duplicados removidos: ", n_removed)

  if (n_removed > 0) {
    message("  Variantes mas frecuentes:")
    top <- variants_dedup %>%
      filter(count > 1) %>%
      select(all_of(by), count) %>%
      head(5)
    for (i in seq_len(nrow(top))) {
      message("    - ", top[[by]][i], ": n=", top$count[i])
    }
  }

  return(variants_dedup)
}
