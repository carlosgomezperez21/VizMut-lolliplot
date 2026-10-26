#' Generate protein lolliplot
#'
#' @param variants data.frame resultado de parse_variants()
#' @param features data.frame resultado de parse_features()
#' @param meta lista con metadata del gen y proteína
#' @param grid logical, dividir por variant_type
#' @param label_all logical, etiquetar todas las variantes
#' @return ggplot object

plot_protein_lolliplot <- function(variants,
                                   features = NULL,
                                   meta,
                                   grid      = FALSE,
                                   label_all = FALSE) {

  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(RColorBrewer)
  library(ggnewscale)

  # ----------------------------------------------------------------
  # ACMG colors
  # ----------------------------------------------------------------
  acmg_colors <- c(
    "P"           = "#D73027",
    "LP"          = "#FC8D59",
    "P/LP"        = "#F46D43",
    "VUS"         = "#BDBDBD",
    "LB"          = "#91BFDB",
    "B"           = "#1A9850",
    "Conflicting" = "#9E66AB",
    "RF"          = "#D9B600",
    "DR"          = "#3288BD",
    "NP"          = "#969696",
    "Other"       = "#636363",
    "NA"          = "#525252"
  )

  protein_length <- meta$protein$length

  # ----------------------------------------------------------------
  # Preparar features
  # ----------------------------------------------------------------
  if (!is.null(features) && nrow(features) > 0) {
    # Normalizar nombre de columna: fetch_uniprot devuelve "feature_type" y "feature_name"
    if ("feature_type" %in% names(features) && !"type" %in% names(features))
      features <- features %>% rename(type = feature_type)
    if ("feature_name" %in% names(features) && !"name" %in% names(features))
      features <- features %>% rename(name = feature_name)
    features <- features %>%
      mutate(type = tolower(gsub(" ", "_", type)))
  } else {
    features <- data.frame(type=character(), name=character(),
                           start=integer(), end=integer(),
                           stringsAsFactors=FALSE)
  }

  domains     <- features %>% filter(type == "domain")
  motifs      <- features %>% filter(type == "motif")
  ptm         <- features %>% filter(type == "ptm")
  zf          <- features %>% filter(type %in% c("zf", "zinc_finger"))
  ss_helix    <- features %>% filter(type == "helix")
  ss_beta     <- features %>% filter(type == "beta_strand")
  ss_turn     <- features %>% filter(type == "turn")
  mutagenesis <- features %>% filter(type == "mutagenesis")

  if (nrow(domains) > 0) domains <- domains %>% mutate(mid = (start + end) / 2)

  # ----------------------------------------------------------------
  # Paletas dinámicas
  # ----------------------------------------------------------------
  make_colors <- function(names_vec, palette) {
    n <- length(unique(names_vec))
    if (n == 0) return(character(0))
    cols <- colorRampPalette(brewer.pal(min(8, max(3, n)), palette))(n)
    setNames(cols, unique(names_vec))
  }

  domain_colors <- make_colors(domains$name, "Set3")
  motif_colors  <- make_colors(motifs$name,  "Pastel1")
  ptm_colors    <- make_colors(ptm$name,     "Set2")
  zf_colors     <- make_colors(zf$name,      "Dark2")

  helix_color <- "#E8A598"
  beta_color  <- "#98C4E8"
  turn_color  <- "#B8E8A0"
  mut_color   <- "#F5D76E"

  # ----------------------------------------------------------------
  # Y layout  (backbone en y=0, lollipops hacia arriba)
  # Capas debajo del backbone (negativo):
  #   SS:          -0.06 a -0.14
  #   Dominios:    -0.18 a -0.30
  #   Motivos:     -0.34 a -0.44
  #   Mutagénesis: -0.47 a -0.49 (línea fina)
  #   PTM:         -0.54 (círculos)
  #   ZF:          -0.60 a -0.70
  # ----------------------------------------------------------------

  # ----------------------------------------------------------------
  # Preparar variantes
  # ----------------------------------------------------------------
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
      label = ifelse(!is.na(protein_change), protein_change, variant_id)
    )

  # ----------------------------------------------------------------
  # Base plot (backbone)
  # ----------------------------------------------------------------
  p <- ggplot(variants) +
    annotate("segment",
             x = 0, xend = protein_length,
             y = 0, yend = 0,
             linewidth = 2.5,
             color = "black") +
    scale_x_continuous(
      limits = c(-protein_length * 0.04, protein_length * 1.04),
      labels = scales::comma
    )

  # ----------------------------------------------------------------
  # Layer 1: Estructura secundaria — franja SOBRE el backbone
  # (entre backbone y los palitos de lollipop)
  # helix/beta:  y = 0.02 a 0.07
  # turn:        y = 0.025 a 0.065 (algo más delgada)
  # ----------------------------------------------------------------
  if (nrow(ss_helix) > 0) {
    p <- p +
      geom_rect(data = ss_helix %>%
                  mutate(end = ifelse(end == start, start + 1, end)),
                aes(xmin = start, xmax = end,
                    ymin = 0.02, ymax = 0.07),
                fill = helix_color, color = NA,
                inherit.aes = FALSE)
  }
  if (nrow(ss_beta) > 0) {
    p <- p +
      geom_rect(data = ss_beta %>%
                  mutate(end = ifelse(end == start, start + 1, end)),
                aes(xmin = start, xmax = end,
                    ymin = 0.02, ymax = 0.07),
                fill = beta_color, color = NA,
                inherit.aes = FALSE)
  }
  if (nrow(ss_turn) > 0) {
    p <- p +
      geom_rect(data = ss_turn %>%
                  mutate(end = ifelse(end == start, start + 1, end)),
                aes(xmin = start, xmax = end,
                    ymin = 0.025, ymax = 0.065),
                fill = turn_color, color = NA,
                inherit.aes = FALSE)
  }

  # ----------------------------------------------------------------
  # Layer 2: Dominios
  # ----------------------------------------------------------------
  if (nrow(domains) > 0) {
    p <- p +
      geom_rect(data = domains,
                aes(xmin = start, xmax = end,
                    ymin = -0.18, ymax = -0.30,
                    fill = name),
                alpha = 0.9,
                inherit.aes = FALSE) +
      scale_fill_manual(values = domain_colors, name = "Domain",
                        guide = guide_legend(ncol = 1,
                                             title.position = "top"))
  }

  # ----------------------------------------------------------------
  # Layer 3: Motivos
  # ----------------------------------------------------------------
  if (nrow(motifs) > 0) {
    p <- p +
      new_scale_fill() +
      geom_rect(data = motifs %>%
                  mutate(end = ifelse(end == start, start + 2, end)),
                aes(xmin = start, xmax = end,
                    ymin = -0.34, ymax = -0.44,
                    fill = name),
                alpha = 0.9,
                inherit.aes = FALSE) +
      scale_fill_manual(values = motif_colors, name = "Motif",
                        guide = guide_legend(ncol = 1,
                                             title.position = "top"))
  }

  # ----------------------------------------------------------------
  # Layer 4: Mutagénesis (banda delgada)
  # ----------------------------------------------------------------
  if (nrow(mutagenesis) > 0) {
    p <- p +
      geom_rect(data = mutagenesis %>%
                  mutate(end = ifelse(end == start, start + 1, end)),
                aes(xmin = start - 0.5, xmax = end + 0.5,
                    ymin = -0.47, ymax = -0.49),
                fill = mut_color, color = NA,
                inherit.aes = FALSE)
  }

  # ----------------------------------------------------------------
  # Layer 5: PTMs (círculos)
  # ----------------------------------------------------------------
  if (nrow(ptm) > 0) {
    p <- p +
      new_scale_fill() +
      geom_segment(data = ptm,
                   aes(x = start, xend = start,
                       y = -0.01, yend = -0.50),
                   color = "gray70", linewidth = 0.3,
                   inherit.aes = FALSE) +
      geom_point(data = ptm,
                 aes(x = start, y = -0.54, fill = name),
                 shape = 21, size = 2.5,
                 color = "black", stroke = 0.3,
                 inherit.aes = FALSE) +
      scale_fill_manual(values = ptm_colors, name = "PTM",
                        guide = guide_legend(ncol = 1,
                                             title.position = "top"))
  }

  # ----------------------------------------------------------------
  # Layer 6: Zinc Fingers
  # ----------------------------------------------------------------
  if (nrow(zf) > 0) {
    p <- p +
      new_scale_fill() +
      geom_rect(data = zf,
                aes(xmin = start, xmax = end,
                    ymin = -0.60, ymax = -0.70,
                    fill = name),
                alpha = 0.9,
                inherit.aes = FALSE) +
      scale_fill_manual(values = zf_colors, name = "Zinc Finger",
                        guide = guide_legend(ncol = 1,
                                             title.position = "top"))
  }

  # ----------------------------------------------------------------
  # Leyenda de estructura secundaria (anotación manual como guide)
  # Solo se agrega si hay al menos un tipo de SS
  # ----------------------------------------------------------------
  has_ss <- nrow(ss_helix) > 0 || nrow(ss_beta) > 0 || nrow(ss_turn) > 0
  if (has_ss) {
    ss_legend_data <- data.frame(
      x    = numeric(0),
      y    = numeric(0),
      fill = character(0),
      stringsAsFactors = FALSE
    )
    if (nrow(ss_helix) > 0)
      ss_legend_data <- rbind(ss_legend_data,
        data.frame(x=NA_real_, y=NA_real_, fill="Helix", stringsAsFactors=FALSE))
    if (nrow(ss_beta) > 0)
      ss_legend_data <- rbind(ss_legend_data,
        data.frame(x=NA_real_, y=NA_real_, fill="Beta strand", stringsAsFactors=FALSE))
    if (nrow(ss_turn) > 0)
      ss_legend_data <- rbind(ss_legend_data,
        data.frame(x=NA_real_, y=NA_real_, fill="Turn", stringsAsFactors=FALSE))

    ss_fill_vals <- c("Helix" = helix_color,
                      "Beta strand" = beta_color,
                      "Turn" = turn_color)

    p <- p +
      new_scale_fill() +
      geom_tile(data = ss_legend_data,
                aes(x = x, y = y, fill = fill),
                inherit.aes = FALSE) +
      scale_fill_manual(values = ss_fill_vals,
                        name  = "2° Structure",
                        guide = guide_legend(ncol = 1,
                                             title.position = "top"))
  }

  # ----------------------------------------------------------------
  # Lollipops (variantes) — encima del backbone
  # ----------------------------------------------------------------
  p <- p +
    new_scale_fill() +
    geom_segment(data = variants,
                 aes(x = protein_pos, xend = protein_pos,
                     y = 0.09, yend = 1.0,
                     color = ACMG),
                 linewidth = 0.8) +
    geom_point(data = variants,
               aes(x = protein_pos, y = 1.0,
                   color = ACMG),
               size = 3) +
    scale_color_manual(values = acmg_colors, name = "ACMG",
                       guide = guide_legend(ncol = 1,
                                            title.position = "top"))

  # ----------------------------------------------------------------
  # Etiquetas (geom_text_repel sin fill para evitar el error de dimensión)
  # ----------------------------------------------------------------
  label_data <- if (label_all) {
    variants
  } else {
    variants %>% filter(ACMG %in% c("P", "LP", "P/LP", "Conflicting"))
  }

  if (nrow(label_data) > 0) {
    # Color manual para labels — no usa fill, solo color
    label_colors <- acmg_colors[as.character(label_data$ACMG)]
    label_colors[is.na(label_colors)] <- "#525252"

    p <- p +
      geom_text_repel(
        data               = label_data,
        aes(x = protein_pos, y = 1.0,
            label = label, color = ACMG),
        size               = 2.2,
        fontface           = "bold",
        max.overlaps       = 25,
        min.segment.length = 0,
        box.padding        = 0.35,
        nudge_y            = 0.22,
        direction          = "x",
        segment.size       = 0.3,
        segment.color      = "gray60",
        show.legend        = FALSE,
        seed               = 42
      )
  }

  # ----------------------------------------------------------------
  # Grid
  # ----------------------------------------------------------------
  if (grid) {
    p <- p + facet_wrap(~variant_type, ncol = 1, scales = "free_y")
  }

  # ----------------------------------------------------------------
  # Límites Y dinámicos
  # ----------------------------------------------------------------
  y_min <- -0.40   # sin features debajo del backbone
  if (nrow(domains) > 0 || nrow(motifs) > 0) y_min <- -0.52
  if (nrow(ptm) > 0)         y_min <- -0.62
  if (nrow(zf) > 0)          y_min <- -0.78

  if (!grid) {
    p <- p + coord_cartesian(ylim = c(y_min, 1.42), expand = FALSE)
  }

  # ----------------------------------------------------------------
  # Theme
  # ----------------------------------------------------------------
  p <- p +
    theme_bw() +
    labs(
      title = paste(meta$gene$name, "- Protein Variants"),
      x     = "Amino acid position",
      y     = ""
    ) +
    theme(
      panel.background   = element_rect(fill = "white", color = NA),
      plot.background    = element_rect(fill = "white", color = NA),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_line(color = "gray92", linewidth = 0.3),
      panel.grid.minor.x = element_blank(),
      panel.border       = element_blank(),
      axis.title.y       = element_blank(),
      axis.text.y        = element_blank(),
      axis.ticks.y       = element_blank(),
      axis.text.x        = element_text(size = 9),
      axis.title.x       = element_text(size = 11, face = "bold"),
      plot.title         = element_text(size = 13, face = "bold"),
      legend.position    = "right",
      legend.direction   = "vertical",
      legend.box         = "vertical",
      legend.box.just    = "left",
      legend.title       = element_text(size = 7, face = "bold"),
      legend.text        = element_text(size = 6),
      legend.key.size    = unit(0.3, "cm"),
      legend.spacing.y   = unit(0.2, "cm"),
      legend.margin      = margin(0, 0, 0, 5),
      legend.background  = element_rect(fill = "white", color = "gray80",
                                        linewidth = 0.3),
      legend.box.margin  = margin(5, 5, 5, 5)
    ) +
    guides(
      color = guide_legend(ncol = 1, title.position = "top",
                           override.aes = list(size = 3))
    )

  return(p)
}
