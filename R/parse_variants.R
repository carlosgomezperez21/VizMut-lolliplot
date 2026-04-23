#' Parse variant file
#'
#' @param df data.frame leído desde el CSV de variantes
#' @return data.frame con columnas chr, pos, ref, alt, protein_pos

parse_variants <- function(df) {

  library(dplyr)
  library(stringr)

  df_parsed <- df %>%

    mutate(
      chr = str_extract(gen_position, "chr[0-9XY]+"),
      pos = as.numeric(str_extract(gen_position, "(?<=-)[0-9]+")),
      ref = str_extract(gen_position, "(?<= )[A-Z]"),
      alt = str_extract(gen_position, "(?<=>)[A-Z]")
    ) %>%

    mutate(
      protein_pos = as.numeric(str_extract(protein_change, "[0-9]+"))
    )

  return(df_parsed)
}
