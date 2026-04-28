#' Fetch chromosome cytobands from UCSC API
#'
#' @param chrom cromosoma (ej. "chr19")
#' @param genome ensamble (default: "hg38")
#' @return data.frame con columnas: chrom, chromStart, chromEnd, name, gieStain

fetch_cytobands <- function(chrom, genome = "hg38") {

  library(httr)
  library(jsonlite)

  message("Obteniendo citobandas de ", chrom, " (", genome, ")...")

  resp <- GET(paste0(
    "https://api.genome.ucsc.edu/getData/track",
    "?genome=", genome,
    ";track=cytoBand",
    ";chrom=", chrom
  ))

  if (http_error(resp)) {
    stop("Error al obtener citobandas para: ", chrom)
  }

  data <- fromJSON(content(resp, as = "text", encoding = "UTF-8"),
                   simplifyVector = TRUE)

  bands <- data$cytoBand
  message("Citobandas obtenidas: ", nrow(bands), " bandas")

  return(bands)
}

#' Get standard colors for Giemsa staining patterns
#'
#' @return named vector con colores por tipo de tincion

giemsa_colors <- function() {
  c(
    "gneg"   = "#FFFFFF",  # blanco - negativo
    "gpos25" = "#C8C8C8",  # gris claro
    "gpos50" = "#969696",  # gris medio
    "gpos75" = "#646464",  # gris oscuro
    "gpos100"= "#000000",  # negro - muy denso
    "acen"   = "#D92F27",  # rojo - centromero
    "gvar"   = "#DCDCDC",  # gris muy claro - variable
    "stalk"  = "#647FA4"   # azul - stalk
  )
}
