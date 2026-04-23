#' Validar que las variantes caen dentro de los límites del gen/proteía
#'
#' @param variants data.frame resultado de parse_variants()
#' @param meta lista con metadata del gen y proteína

validate_variants <- function(variants, meta) {

  if (any(variants$pos < meta$gene$start |
          variants$pos > meta$gene$end)) {

    warning("Some variants fall outside gene boundaries")
  }
}
