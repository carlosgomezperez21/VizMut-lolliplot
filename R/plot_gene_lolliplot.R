#' Plot gene structure lolliplot
#'
#' @param variants data.frame resultado de parse_variants()
#' @param transcript_structure data.frame resultado de fetch_transcript_structure()
#' @param gene_name nombre del gen para el título
#' @param label_type "protein" (p.) o "cds" (c.) para etiquetas
#' @param grid logical, dividir por variant_type
#' @return ggplot object

plot_gene_lolliplot <- function(variants,
                                transcript_structure,
                                gene_name = "Gene",
                                label_type = "protein",
                                grid = FALSE) {

  library(ggplot2)
  library(dplyr)
  library(ggrepel)

  #---------------------------
  # ACMG colors
  #---------------------------
  acmg_colors <- c(
    "P"   = "#D73027",
    "LP"  = "#FC8D59",
    "VUS" = "#BDBDBD",
    "LB"  = "#91BFDB",
    "B"   = "#1A9850"
  )

  #---------------------------
  # Preparar estructura del transcrito
  #---------------------------
  exons   <- transcript_structure %>% filter(type == "exon")
  introns <- transcript_structure %>% filter(type == "intron")
  utr5    <- transcript_structure %>% filter(type == "utr5")
  utr3    <- transcript_structure %>% filter(type == "utr3")

  gene_start <- min(transcript_structure$start)
  gene_end   <- max(transcript_structure$end)

  #---------------------------
  # Preparar variantes
  #---------------------------
  variants$variant_type <- ifelse(is.na(variants$variant_type),
                                  "unknown", variants$variant_type)
  variants$ACMG <- as.character(variants$ACMG)
  variants <- variants[!is.na(variants$pos), ]

  # etiqueta según disponibilidad de columnas
  if (label_type == "cds" && "hgvs_c" %in% names(variants)) {
    variants$label <- ifelse(!is.na(variants$hgvs_c),
                             variants$hgvs_c, variants$variant_id)
  } else if ("protein_change" %in% names(variants)) {
    variants$label <- ifelse(!is.na(variants$protein_change),
                             variants$protein_change, variants$variant_id)
  } else {
    variants$label <- variants$variant_id
  }

  # etiqueta según disponibilidad de columnas
  if (label_type == "cds" && "hgvs_c" %in% names(variants)) {
    variants <- variants %>%
      mutate(label = ifelse(!is.na(hgvs_c), hgvs_c, variant_id))
  } else if (!is.na(match("protein_change", names(variants)))) {
    variants <- variants %>%
      mutate(label = ifelse(!is.na(protein_change), protein_change, variant_id))
  } else {
    variants <- variants %>%
      mutate(label = variant_id)
  }

  #--------------------------
  # Base plot
  #---------------------------
  p <- ggplot() +

    # línea base del gen
    annotate("segment",
             x = gene_start, xend = gene_end,
             y = 0, yend = 0,
             linewidth = 0.5,
             color = "gray40")

  #---------------------------
  # UTR5
  #---------------------------
  if (nrow(utr5) > 0) {
    p <- p +
      geom_rect(data = utr5,
                aes(xmin = start, xmax = end,
                    ymin = -0.15, ymax = 0.15),
                fill = "#B0BEC5",
                color = "gray40",
                linewidth = 0.3)
  }

  #---------------------------
  # UTR3
  #---------------------------
  if (nrow(utr3) > 0) {
    p <- p +
      geom_rect(data = utr3,
                aes(xmin = start, xmax = end,
                    ymin = -0.15, ymax = 0.15),
                fill = "#B0BEC5",
                color = "gray40",
                linewidth = 0.3)
  }

  #---------------------------
  # Intrones (línea delgada con dirección)
  #---------------------------
  if (nrow(introns) > 0) {

    introns <- introns %>%
      mutate(mid = (start + end) / 2)

    p <- p +
      geom_segment(data = introns,
                   aes(x = start, xend = mid,
                       y = 0, yend = 0.08),
                   color = "gray50",
                   linewidth = 0.4) +
      geom_segment(data = introns,
                   aes(x = mid, xend = end,
                       y = 0.08, yend = 0),
                   color = "gray50",
                   linewidth = 0.4)
  }

  #---------------------------
  # Exones
  #---------------------------
  if (nrow(exons) > 0) {
    p <- p +
      geom_rect(data = exons,
                aes(xmin = start, xmax = end,
                    ymin = -0.25, ymax = 0.25),
                fill = "#4A90D9",
                color = "#2C5F8A",
                linewidth = 0.3)
  }

  #---------------------------
  # Lollipops de variantes
  #---------------------------
  p <- p +
    geom_segment(data = variants,
                 aes(x = pos, xend = pos,
                     y = 0.25, yend = 1,
                     color = ACMG),
                 linewidth = 0.8) +

    geom_point(data = variants,
               aes(x = pos, y = 1,
                   color = ACMG),
               size = 3) +

    scale_color_manual(values = acmg_colors, name = "ACMG")

  #---------------------------
  # Etiquetas P y LP
  #---------------------------
  p <- p +
    geom_text_repel(
      data = variants %>% filter(ACMG %in% c("P", "LP")),
      aes(x = pos, y = 1, label = label),
      size = 2,
      min.segment.length = 0
    )

  #---------------------------
  # Grid por variant_type
  #---------------------------
  if (grid) {
    p <- p + facet_wrap(~variant_type, ncol = 1)
  }

  #---------------------------
  # Theme
  #---------------------------
  p <- p +
    theme_minimal() +
    labs(
      title = paste(gene_name, "- Gene Structure & Variants"),
      x = "Genomic position",
      y = ""
    ) +
    scale_x_continuous(labels = scales::comma) +
    theme(
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(size = 7, face = "bold"),
      legend.text  = element_text(size = 6),
      legend.key.size = unit(0.25, "cm")
    )

  return(p)
}
