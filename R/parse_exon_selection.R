#' Parse exon selection string into a vector of exon numbers
#'
#' @param exon_str string con formato "1,5-7,10"
#' @return vector numerico con los numeros de exon seleccionados

parse_exon_selection <- function(exon_str) {

  if (is.null(exon_str) || exon_str == "all") return(NULL)

  parts  <- trimws(strsplit(exon_str, ",")[[1]])
  result <- c()

  for (part in parts) {
    if (grepl("-", part)) {
      # rango: "5-7" -> c(5,6,7)
      bounds <- as.integer(strsplit(part, "-")[[1]])
      if (length(bounds) != 2 || any(is.na(bounds))) {
        stop("Rango invalido en --exons: '", part, "'")
      }
      result <- c(result, seq(bounds[1], bounds[2]))
    } else {
      # numero simple: "1" -> c(1)
      n <- suppressWarnings(as.integer(part))
      if (is.na(n)) stop("Valor invalido en --exons: '", part, "'")
      result <- c(result, n)
    }
  }

  result <- sort(unique(result))
  message("Exones seleccionados: ", paste(result, collapse=", "))
  return(result)
}

#' Filter transcript structure to selected exons
#' and add '//' break markers between non-consecutive exons
#'
#' @param struct data.frame de fetch_transcript_structure()
#' @param exon_numbers vector de numeros de exon seleccionados
#' @return list con $structure (filtrada) y $breaks (posiciones de //)

filter_exons <- function(struct, exon_numbers) {

  if (is.null(exon_numbers)) {
    return(list(structure = struct, breaks = NULL))
  }

  # extraer exones con su numero de orden
  exons <- struct[struct$type == "exon", ]
  exons <- exons[order(exons$start), ]
  exons$exon_number <- seq_len(nrow(exons))

  # validar que los numeros pedidos existen
  max_exon <- nrow(exons)
  invalid  <- exon_numbers[exon_numbers > max_exon | exon_numbers < 1]
  if (length(invalid) > 0) {
    stop("Numeros de exon invalidos: ", paste(invalid, collapse=", "),
         ". El transcrito tiene ", max_exon, " exones.")
  }

  # filtrar exones seleccionados
  exons_sel <- exons[exons$exon_number %in% exon_numbers, ]

  # calcular posiciones de breaks (//) entre exones no consecutivos
  breaks <- c()
  if (nrow(exons_sel) > 1) {
    for (i in seq_len(nrow(exons_sel) - 1)) {
      curr_num <- exons_sel$exon_number[i]
      next_num <- exons_sel$exon_number[i + 1]
      if (next_num - curr_num > 1) {
        # break entre el fin del exon actual y el inicio del siguiente
        break_pos <- (exons_sel$end[i] + exons_sel$start[i + 1]) / 2
        breaks    <- c(breaks, break_pos)
      }
    }
  }

  # reconstruir estructura con solo exones seleccionados
  # e intrones SOLO entre exones consecutivos seleccionados
  keep_starts <- exons_sel$start
  keep_ends   <- exons_sel$end

  struct_filtered <- struct[struct$type %in% c("utr5", "utr3") |
    (struct$type == "exon" & struct$start %in% keep_starts), ]

  # agregar solo intrones entre exones seleccionados consecutivos
  if (nrow(exons_sel) > 1) {
    for (i in seq_len(nrow(exons_sel) - 1)) {
      curr_num <- exons_sel$exon_number[i]
      next_num <- exons_sel$exon_number[i + 1]

      # solo agregar intron si los exones son consecutivos
      if (next_num - curr_num == 1) {
        intron_between <- struct[
          struct$type == "intron" &
          struct$start >= exons_sel$end[i] &
          struct$end   <= exons_sel$start[i + 1], ]
        struct_filtered <- bind_rows(struct_filtered, intron_between)
      }
    }
  }

  struct_filtered <- struct_filtered %>% arrange(start)


message("Exones en el plot: ", nrow(exons_sel),
          " de ", max_exon)
  if (length(breaks) > 0) {
    message("Breaks (//) en posiciones: ",
            paste(round(breaks), collapse=", "))
  }

  return(list(
    structure    = struct_filtered,
    breaks       = breaks,
    exon_numbers = exons_sel$exon_number
  ))
}
