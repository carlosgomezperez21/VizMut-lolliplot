#' Fetch canonical RefSeq transcript ID for a gene symbol
#'
#' @param gene_symbol simbolo del gen (ej. "KMT2B")
#' @return string con el transcript ID canonico (ej. "NM_014727.3")

fetch_canonical_transcript <- function(gene_symbol) {

  library(httr)
  library(jsonlite)

  base_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

  #---------------------------
  # 1. Buscar Gene ID por simbolo
  #---------------------------
  message("Buscando transcrito canonico para: ", gene_symbol)

  resp <- GET(paste0(base_url, "/esearch.fcgi",
                     "?db=gene",
                     "&term=", gene_symbol, "[gene]+AND+homo+sapiens[organism]",
                     "&retmode=json"))

  ids <- fromJSON(content(resp, as="text", encoding="UTF-8"))$esearchresult$idlist

  if (length(ids) == 0) {
    stop("No se encontro el gen: ", gene_symbol)
  }

  gene_id <- ids[1]
  message("Gene ID NCBI: ", gene_id)

  #---------------------------
  # 2. Obtener product report
  #---------------------------
  resp_prod <- GET(
    paste0("https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id/",
           gene_id, "/product_report"),
    add_headers("Accept" = "application/json")
  )

  data <- fromJSON(content(resp_prod, as="text", encoding="UTF-8"),
                   simplifyVector = FALSE)

  #---------------------------
  # 3. Buscar transcrito canonico
  #---------------------------
  
  #---------------------------
  # 3. Buscar transcrito MANE Select primero,
  #    fallback a primer NM_ con coordenadas genomicas
  #---------------------------
  mane_select  <- NULL
  first_nm     <- NULL

  for (report in data$reports) {
    transcripts <- tryCatch(report$product$transcripts,
                            error = function(e) list())
    for (tx in transcripts) {
      acc      <- tryCatch(tx$accession_version, error=function(e) "")
      category <- tryCatch(tx$select_category,   error=function(e) "")
      locs     <- tryCatch(tx$genomic_locations,  error=function(e) NULL)

      if (!grepl("^NM_", acc)) next
      if (is.null(locs) || length(locs) == 0) next

      # priorizar MANE Select
      if (!is.null(category) && length(category) > 0 &&
          category == "MANE_SELECT") {
        mane_select <- acc
        break
      }

      # guardar primer NM_ como fallback
      if (is.null(first_nm)) first_nm <- acc
    }
    if (!is.null(mane_select)) break
  }

  canonical <- if (!is.null(mane_select)) mane_select else first_nm

  if (is.null(canonical)) {
    stop("No se encontro transcrito canonico NM_ para: ", gene_symbol)
  }

  if (!is.null(mane_select)) {
    message("Transcrito MANE Select: ", canonical)
  } else {
    message("Transcrito canonico (fallback): ", canonical)
    message("Nota: no se encontro MANE Select para ", gene_symbol)
  }

  return(canonical)

 }
