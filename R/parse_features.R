#' Parse features file
#'
#' @param df data.frame leído desde el CSV de features
#' @return data.frame con columnas type, name, start, end

parse_features <- function(df) {

  library(dplyr)

  df <- df %>%
    rename(
      type = feature_type,
      name = feature_name,
      start = start,
      end = end
    )

  return(df)
}

#' Split features by type
#'
#' @param features data.frame resultado de parse_features()
#' @return lista con dominios, motifs, ptm y zinc fingers

split_features <- function(features) {

  list(
    domains = features %>% dplyr::filter(type == "domain"),
    motifs  = features %>% dplyr::filter(type == "motif"),
    ptm     = features %>% dplyr::filter(type == "ptm"),
    zf      = features %>% dplyr::filter(type == "Zinc finger")
  )
}
