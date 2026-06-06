#' Fetch allele frequencies from gnomAD v4 API
#'
#' @param chr cromosoma (ej. "chr19")
#' @param pos posicion genomica hg38
#' @param ref alelo referencia
#' @param alt alelo alternativo
#' @return data.frame con columnas: pop, af, ac, an, ref, alt

fetch_gnomad_af <- function(chr, pos, ref=NULL, alt=NULL) {

  library(httr)
  library(jsonlite)

  chr_num <- sub("chr", "", chr)

  target_pops <- c("afr", "amr", "eas", "nfe", "sas")
  pop_labels  <- c(
    "afr"="AFR", "amr"="AMR", "eas"="EAS",
    "nfe"="EUR", "sas"="SAS"
  )

  #---------------------------
  # Si tenemos ref y alt: buscar por variantId exacto
  #---------------------------
  if (!is.null(ref) && !is.null(alt) &&
      !is.na(ref) && !is.na(alt) &&
      ref != "na" && alt != "na" &&
      nchar(ref) <= 10 && nchar(alt) <= 10) {

    variant_id <- paste(chr_num, pos, ref, alt, sep="-")

    query <- paste0(
      '{variant(variantId:"', variant_id, '",dataset:gnomad_r4){',
      'exome{populations{id,ac,an}}',
      'genome{populations{id,ac,an}}',
      '}}'
    )

    resp <- tryCatch(
      POST("https://gnomad.broadinstitute.org/api",
           body   = list(query = query),
           encode = "json",
           add_headers("Content-Type" = "application/json")),
      error = function(e) NULL
    )

    if (!is.null(resp) && !http_error(resp)) {
      data    <- fromJSON(content(resp, as="text", encoding="UTF-8"),
                          simplifyVector=FALSE)
      variant <- tryCatch(data$data$variant, error=function(e) NULL)

      if (!is.null(variant)) {
        result_list <- list()
        for (pop in target_pops) {
          ac_total <- 0
          an_total <- 0
          for (src in c("exome", "genome")) {
            src_data <- tryCatch(variant[[src]]$populations,
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
            ref        = ref,
            alt        = alt,
            variant_id = variant_id,
            stringsAsFactors = FALSE
          )
        }
        result <- do.call(rbind, result_list)
        return(result)
      }
    }
  }

  #---------------------------
  # Fallback: buscar por region
  #---------------------------
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

  data     <- fromJSON(content(resp, as="text", encoding="UTF-8"),
                       simplifyVector=FALSE)
  variants <- tryCatch(data$data$region$variants,
                       error=function(e) NULL)

  if (is.null(variants) || length(variants) == 0) return(NULL)

  v   <- variants[[1]]
  vid <- v$variant_id

  # extraer ref y alt del variant_id
  vid_parts  <- strsplit(vid, "-")[[1]]
  ref_gnomad <- ifelse(length(vid_parts) >= 3, vid_parts[3], NA)
  alt_gnomad <- ifelse(length(vid_parts) >= 4, vid_parts[4], NA)

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
      ref        = ref_gnomad,
      alt        = alt_gnomad,
      variant_id = vid,
      stringsAsFactors = FALSE
    )
  }

  return(do.call(rbind, result_list))
}
