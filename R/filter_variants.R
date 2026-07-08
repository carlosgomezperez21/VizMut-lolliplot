#' Filtrar variantes para graficar segun criterios del usuario
#'
#' @param variants data.frame de variantes
#' @param show_only string con filtros: "acmg=P,LP type=SNV phenotype=ALS"
#' @return data.frame filtrado

filter_variants <- function(variants, show_only) {

  if (is.null(show_only) || show_only == "") return(variants)

  # parsear los filtros — formato: "key=val1,val2 key2=val3"
  tokens <- strsplit(trimws(show_only), "\\s+")[[1]]

  for (token in tokens) {
    parts <- strsplit(token, "=")[[1]]
    if (length(parts) != 2) {
      warning("Filtro ignorado (formato invalido): ", token)
      next
    }

    key    <- trimws(parts[1])
    values <- trimws(strsplit(parts[2], ",")[[1]])

    if (key == "acmg") {
      if ("ACMG" %in% names(variants)) {
        variants <- variants[!is.na(variants$ACMG) &
                             variants$ACMG %in% values, ]
        message("  Filtro ACMG=", paste(values, collapse=","),
                ": ", nrow(variants), " variantes")
      } else {
        warning("Columna ACMG no encontrada — filtro acmg ignorado")
      }
    } else if (key == "type") {
      if ("variant_type" %in% names(variants)) {
        # normalizar aliases comunes
        type_aliases <- list(
          "SNV"       = c("SNV", "single nucleotide variant",
                          "single_nucleotide_variant"),
          "DEL"       = c("DEL", "Deletion", "deletion"),
          "INS"       = c("INS", "Insertion", "insertion"),
          "DUP"       = c("DUP", "Duplication", "duplication"),
          "INDEL"     = c("INDEL", "Indel", "indel"),
          "Microsatellite" = c("Microsatellite", "microsatellite")
        )
        # expandir los valores del usuario con sus aliases
        expanded <- unlist(lapply(values, function(v) {
          matched <- type_aliases[[toupper(v)]]
          if (is.null(matched)) {
            # buscar por alias inverso
            for (canonical in names(type_aliases)) {
              if (v %in% type_aliases[[canonical]]) return(type_aliases[[canonical]])
            }
            return(v)
          }
          return(matched)
        }))
        variants <- variants[!is.na(variants$variant_type) &
                             variants$variant_type %in% expanded, ]
        message("  Filtro type=", paste(values, collapse=","),
                ": ", nrow(variants), " variantes")
      } else {
        warning("Columna variant_type no encontrada — filtro type ignorado")
      }
    } else if (key == "phenotype") {
      if ("phenotype" %in% names(variants)) {
        pattern <- paste(values, collapse="|")
        variants <- variants[!is.na(variants$phenotype) &
                             grepl(pattern, variants$phenotype,
                                   ignore.case=TRUE), ]
        message("  Filtro phenotype=", paste(values, collapse=","),
                ": ", nrow(variants), " variantes")
      } else {
        warning("Columna phenotype no encontrada — filtro phenotype ignorado")
      }
    } else {
      warning("Filtro desconocido ignorado: ", key)
    }
  }

  if (nrow(variants) == 0) {
    stop("--show_only no retorno variantes. Verifica los criterios de filtro.")
  }

  return(variants)
}
