#' Generate protein lolliplot
#'
#' @param variants data.frame resultado de parse_variants()
#' @param features data.frame resultado de parse_features()
#' @param meta lista con metadata del gen y proteína
#' @param grid logical, dividir por variant_type
#' @return ggplot object

plot_protein_lolliplot <- function(variants,
                                   features = NULL,
                                   meta,
                                   grid = FALSE) {

  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(RColorBrewer)
  library(ggnewscale)

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

  protein_length <- meta$protein$length

  #---------------------------
  # features
  #---------------------------
  features <- features %>%
    mutate(type = tolower(type)) %>%
    mutate(type = gsub(" ", "_", type))

  domains <- features %>% filter(type == "domain")
  motifs  <- features %>% filter(type == "motif")
  ptm     <- features %>% filter(type == "ptm")
  zf      <- features %>% filter(type %in% c("zf", "zinc_finger"))

  if (nrow(domains) > 0) {
    domains <- domains %>%
      mutate(mid = (start + end) / 2)
  }

  #---------------------------
  # palettes dinamicas
  #---------------------------
  domain_colors <- setNames(
    brewer.pal(max(3, min(8, length(unique(domains$name)))), "Set3"),
    unique(domains$name)
  )

  motif_colors <- setNames(
    brewer.pal(max(3, min(8, length(unique(motifs$name)))), "Pastel1"),
    unique(motifs$name)
  )

  ptm_colors <- setNames(
    brewer.pal(max(3, min(8, length(unique(ptm$name)))), "Set2"),
    unique(ptm$name)
  )

  zf_colors <- setNames(
    brewer.pal(max(3, min(8, length(unique(zf$name)))), "Dark2"),
    unique(zf$name)
  )

  #---------------------------
  # variants
  #---------------------------

  # filtrar variantes sin posicion proteica
  n_intronic <- sum(is.na(variants$protein_pos))
  if (n_intronic > 0) {
    message("Nota: ", n_intronic,
            " variantes intronicas/UTR excluidas del plot proteico")
    variants <- variants %>% filter(!is.na(protein_pos))
  }

  variants <- variants %>%
    mutate(
      variant_type = ifelse(is.na(variant_type), "unknown", variant_type),
      ACMG  = as.character(ACMG),
      label = ifelse(!is.na(protein_change),
                     protein_change,
                     variant_id)
    )

  #---------------------------
  # base plot
  #---------------------------
  p <- ggplot(variants) +
    annotate("segment",
             x = 0, xend = protein_length,
             y = 0, yend = 0,
             linewidth = 2.5)

  #---------------------------
  # DOMAINS
  #---------------------------
  if (nrow(domains) > 0) {
    p <- p +
      geom_rect(data = domains,
                aes(xmin = start, xmax = end,
                    ymin = -0.15, ymax = -0.05,
                    fill = name),
                alpha = 0.9) +
      scale_fill_manual(values = domain_colors, name = "Domain")
  }

  #---------------------------
  # MOTIFS
  #---------------------------
  if (nrow(motifs) > 0) {
    p <- p +
      new_scale_fill() +
      geom_rect(data = motifs,
                aes(xmin = start, xmax = end,
                    ymin = -0.30, ymax = -0.18,
                    fill = name),
                alpha = 0.9) +
      scale_fill_manual(values = motif_colors, name = "Motif")
  }

  #---------------------------
  # PTM
  #---------------------------
  if (nrow(ptm) > 0) {
    p <- p +
      new_scale_fill() +
      geom_segment(data = ptm,
                   aes(x = start, xend = start,
                       y = -0.05, yend = -0.35),
                   color = "gray50",
                   linewidth = 0.4) +
      geom_point(data = ptm,
                 aes(x = start, y = -0.42, fill = name),
                 shape = 21,
                 size  = 3,
                 color = "black") +
      scale_fill_manual(values = ptm_colors, name = "PTM")
  }

  #---------------------------
  # ZINC FINGERS
  #---------------------------
  if (nrow(zf) > 0) {
    p <- p +
      new_scale_fill() +
      geom_rect(data = zf,
                aes(xmin = start, xmax = end,
                    ymin = -0.58, ymax = -0.46,
                    fill = name),
                alpha = 0.9) +
      scale_fill_manual(values = zf_colors, name = "Zinc Finger")
  }

  #---------------------------
  # VARIANTS
  #---------------------------
  p <- p +
    geom_segment(data = variants,
                 aes(x = protein_pos, xend = protein_pos,
                     y = 0.05, yend = 1,
                     color = ACMG),
                 linewidth = 0.8) +
    geom_point(data = variants,
               aes(x = protein_pos, y = 1,
                   color = ACMG),
               size = 3) +
    scale_color_manual(values = acmg_colors, name = "ACMG")

  # labels P y LP
  p <- p +
    geom_text_repel(
      data = variants %>% filter(ACMG %in% c("P", "LP")),
      aes(x = protein_pos, y = 1, label = label),
      size = 2
    )

  #---------------------------
  # Grid por variant_type
  #---------------------------
  if (grid) {
    p <- p + facet_wrap(~variant_type, ncol = 1, scales = "free_y")
  }

  #---------------------------
  # Escala Y fija cuando no hay grid
  #---------------------------
  if (!grid) {
    p <- p + coord_cartesian(ylim = c(-0.70, 1.1), expand = FALSE)
  }

#---------------------------
  # Theme
  #---------------------------
  p <- p +
    theme_minimal() +
    labs(
      title = paste(meta$gene$name, "- Protein Variants"),
      x     = "Amino acid position",
      y     = ""
    ) +
    scale_x_continuous(
      breaks = c(0, 500, 1000, 1500, 2000, protein_length),
      labels = c("0", "500", "1000", "1500", "2000",
                 as.character(protein_length)),
      expand = c(0.1, 0)
    ) +
    theme(
      legend.position  = "right",
      legend.title     = element_text(size = 7, face = "bold"),
      legend.text      = element_text(size = 6),
      legend.key.size  = unit(0.25, "cm"),
      legend.box       = "vertical",
      axis.title.y     = element_blank(),
      axis.text.y      = element_blank(),
      axis.ticks.y     = element_blank()
    ) +
    guides(
      fill  = guide_legend(ncol = 1, title.position = "top"),
      color = guide_legend(ncol = 1, title.position = "top")
    )

  return(p)
}
