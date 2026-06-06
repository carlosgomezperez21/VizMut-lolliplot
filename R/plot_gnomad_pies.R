add_gnomad_pies <- function(p, variants, acmg_colors) {

  library(dplyr)

  pop_colors <- c(
    "AFR" = "#E41A1C",
    "AMR" = "#FF7F00",
    "EAS" = "#4DAF4A",
    "EUR" = "#377EB8",
    "SAS" = "#984EA3"
  )

  gene_range <- max(variants$pos, na.rm=TRUE) - min(variants$pos, na.rm=TRUE)
  r_y        <- 0.06
  r_x        <- r_y * (gene_range / 1.5) * (7 / 14)
  n_points   <- 50
  pie_data   <- list()

  for (i in seq_len(nrow(variants))) {

    chr_i <- as.character(variants$chr[i])
    pos_i <- as.numeric(variants$pos[i])
    vid_i <- as.character(variants$variant_id[i])
    y_i   <- as.numeric(variants$y_top[i])

    if (is.na(pos_i) || is.na(chr_i)) next
      if (length(vid_i) == 0 || is.na(vid_i) || vid_i == "") next
      if (length(y_i) == 0 || is.na(y_i)) next

    af_raw <- tryCatch(
      fetch_gnomad_af(chr_i, pos_i),
      error = function(e) NULL
    )

    if (is.null(af_raw) || nrow(af_raw) == 0) next

    # forzar tipos y verificar
    af_data <- data.frame(
      pop = as.character(af_raw$pop),
      af  = as.numeric(af_raw$af),
      stringsAsFactors = FALSE
    )

    # verificar que pop no tiene vacios
    af_data <- af_data[nchar(af_data$pop) > 0 & !is.na(af_data$pop), ]
    if (nrow(af_data) == 0) next

   af_sum <- sum(af_data$af)
    if (af_sum == 0 || !is.finite(af_sum)) next

    current_angle <- -pi / 2

    for (j in seq_len(nrow(af_data))) {

      pop_j  <- af_data$pop[j]
      af_j   <- af_data$af[j]
      prop_j <- af_j / af_sum
      arc_j  <- prop_j * 2 * pi

      if (!is.finite(arc_j) || arc_j < 0.01) {
        current_angle <- current_angle + arc_j
        next
      }

      n_seg    <- max(3L, as.integer(round(n_points * prop_j)))
      angles   <- seq(current_angle, current_angle + arc_j, length.out = n_seg)
      sx       <- c(pos_i, pos_i + r_x * cos(angles), pos_i)
      sy       <- c(y_i,   y_i   + r_y * sin(angles), y_i)
      n_pts    <- length(sx)

      pie_data[[length(pie_data) + 1L]] <- data.frame(
        x          = sx,
        y          = sy,
        pop        = rep(pop_j,                      n_pts),
        variant_id = rep(vid_i,                      n_pts),
        group      = rep(paste(vid_i, pop_j, sep="_"), n_pts),
        stringsAsFactors = FALSE
      )

      current_angle <- current_angle + arc_j
    }
  }

  if (length(pie_data) == 0) {
    message("No se encontraron frecuencias en gnomAD para las variantes.")
    return(p)
  }

  pie_df  <- do.call(rbind, pie_data)

  n_total <- sum(!is.na(variants$pos))
  n_pies  <- length(unique(pie_df$variant_id))
  n_no_af <- n_total - n_pies

  message("\nResumen gnomAD v4:")
  message("  Variantes con frecuencia poblacional: ", n_pies)
  message("  Variantes sin datos en gnomAD: ", n_no_af)

  if (n_no_af > 0 && "hgvs_c" %in% names(variants)) {
    no_af <- variants$hgvs_c[!variants$variant_id %in%
                               unique(pie_df$variant_id) &
                               !is.na(variants$pos)]
    message("  Lista:")
    for (v in unique(no_af)) message("    - ", v)
  }

  p <- p +
    geom_polygon(
      data        = pie_df,
      aes(x=x, y=y, group=group, fill=pop),
      color       = "white",
      linewidth   = 0.2,
      alpha       = 0.9,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(
      values = pop_colors,
      name   = "Population (gnomAD v4)",
      guide  = guide_legend(nrow=1, direction="horizontal")
    )

  return(p)
}
