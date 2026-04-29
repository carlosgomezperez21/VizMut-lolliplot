#' Plot chromosome ideogram with gene location highlighted
#'
#' @param cytobands data.frame de fetch_cytobands()
#' @param gene_start posicion genomica de inicio del gen
#' @param gene_end posicion genomica de fin del gen
#' @param chrom nombre del cromosoma
#' @param gene_name nombre del gen para la etiqueta
#' @return ggplot object

plot_ideogram <- function(cytobands,
                          gene_start,
                          gene_end,
                          chrom,
                          gene_name = "") {

  library(ggplot2)
  library(dplyr)

  colors <- giemsa_colors()

  # asignar colores a cada banda
  cytobands <- cytobands %>%
    mutate(
      fill_color = ifelse(gieStain %in% names(colors),
                          colors[gieStain],
                          "#FFFFFF"),
      is_centro  = gieStain == "acen"
    )

  chr_end <- max(cytobands$chromEnd)

  p <- ggplot() +

    # cuerpo del cromosoma
    geom_rect(data = cytobands,
              aes(xmin = chromStart, xmax = chromEnd,
                  ymin = 0.2, ymax = 0.8),
              fill  = cytobands$fill_color,
              color = "gray40",
              linewidth = 0.15) +

    # centromero en forma de triángulo (adelgazamiento)
    geom_rect(data = cytobands %>% filter(is_centro),
              aes(xmin = chromStart, xmax = chromEnd,
                  ymin = 0.3, ymax = 0.7),
              fill  = "#D92F27",
              color = "gray40",
              linewidth = 0.15) +

    # resaltar región del gen
    annotate("rect",
             xmin = gene_start, xmax = gene_end,
             ymin = 0.1, ymax = 0.9,
             fill  = "#FF6B35",
             color = "#CC3300",
             alpha = 0.8,
             linewidth = 0.5) +

    # líneas de delimitación del gen
    annotate("segment",
             x = gene_start, xend = gene_start,
             y = 0, yend = 1,
             color = "#CC3300",
             linewidth = 0.6,
             linetype  = "dashed") +

    annotate("segment",
             x = gene_end, xend = gene_end,
             y = 0, yend = 1,
             color = "#CC3300",
             linewidth = 0.6,
             linetype  = "dashed") +

    # etiqueta del gen
    annotate("text",
             x     = (gene_start + gene_end) / 2,
             y     = 1.3,
             label = gene_name,
             size  = 3,
             fontface = "bold",
             color = "#CC3300") +

    # posiciones genomicas
    annotate("text",
             x = gene_start, y = -0.2,
             label = scales::comma(gene_start),
             size  = 2.5,
             color = "gray30",
             hjust = 0.5) +

    annotate("text",
             x = gene_end, y = -0.2,
             label = scales::comma(gene_end),
             size  = 2.5,
             color = "gray30",
             hjust = 0.5) +

    # nombre del cromosoma
    annotate("text",
             x = -chr_end * 0.02, y = 0.5,
             label = chrom,
             size  = 3.5,
             fontface = "bold",
             hjust = 1) +

    # limites del eje x
    scale_x_continuous(
      limits = c(-chr_end * 0.05, chr_end * 1.02),
      labels = function(x) paste0(round(x/1e6, 0), "Mb")
    ) +

    scale_y_continuous(limits = c(-0.4, 1.6)) +

theme_void() +
    theme(
      axis.text.x  = element_text(size = 8, color = "gray40"),
      axis.ticks.x = element_line(color = "gray40"),
      plot.margin  = margin(8, 15, 0, 15)
    )

  return(p)
}
