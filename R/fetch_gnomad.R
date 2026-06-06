#' Fetch allele frequencies from gnomAD v4 API
#'
#' @param chr cromosoma (ej. "chrX")
#' @param pos posicion genomica hg38
#' @return data.frame con columnas: pop, af, ac, an, variant_id

fetch_gnomad_af <- function(chr, pos) {

  library(httr)
  library(jsonlite)

  chr_num <- sub("chr", "", chr)

  query <- paste0(
    '{region(chrom:"', chr_num, '",start:', pos-1, ',stop:', pos+1,
    ',reference_genome:GRCh38){',
    'variants(dataset:gnomad_r4){',
    'variant_id,',
    'exome{populations{id,ac,an}}',
    'genome{populations{id,ac,an}}',
    '}}}'
  )

  resp <- tryCatch(
    POST("https://gnomad.broadinstitute.org/api",
         body   = list(query = query),
         encode = "json",
         add_headers("Content-Type" = "application/json")),
    error = function(e) NULL
  )

  if (is.null(resp) || http_error(resp)) return(NULL)

  data <- fromJSON(content(resp, as="text", encoding="UTF-8"),
                   simplifyVector=FALSE)

  variants <- tryCatch(
    data$data$region$variants,
    error = function(e) NULL
  )

  if (is.null(variants) || length(variants) == 0) return(NULL)

  target_pops <- c("afr", "amr", "eas", "nfe", "sas")
  pop_labels  <- c(
    "afr"="AFR", "amr"="AMR", "eas"="EAS",
    "nfe"="EUR", "sas"="SAS"
  )

  # tomar la primera variante en esa posicion
  # (la mas comun o la unica)
  v      <- variants[[1]]
  vid    <- v$variant_id

  result_list <- list()

  for (pop in target_pops) {
    ac_total <- 0
    an_total <- 0

    for (src in c("exome", "genome")) {
      src_data <- tryCatch(v[[src]]$populations,
                           error=function(e) NULL)
      if (is.null(src_data)) next
      for (p in src_data) {
        if (p$id == pop) {
          ac_total <- ac_total + p$ac
          an_total <- an_total + p$an
          break
        }
      }
    }

    af <- ifelse(an_total > 0, ac_total / an_total, 0)
    result_list[[pop]] <- data.frame(
      pop        = pop_labels[pop],
      af         = af,
      ac         = ac_total,
      an         = an_total,
      variant_id = vid,
      stringsAsFactors = FALSE
    )
  }

  result <- do.call(rbind, result_list)
  return(result)
}
