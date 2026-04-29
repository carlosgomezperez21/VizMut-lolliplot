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
  canonical <- NULL

  for (report in data$reports) {
    transcripts <- tryCatch(report$product$transcripts,
                            error = function(e) list())
    for (tx in transcripts) {
      acc <- tryCatch(tx$accession_version, error = function(e) "")
      # transcrito canonico empieza con NM_ y tiene genomic_locations
      if (grepl("^NM_", acc)) {
        locs <- tryCatch(tx$genomic_locations, error = function(e) NULL)
        if (!is.null(locs) && length(locs) > 0) {
          canonical <- acc
          break
        }
      }
    }
    if (!is.null(canonical)) break
  }

  if (is.null(canonical)) {
    stop("No se encontro transcrito canonico NM_ para: ", gene_symbol)
  }

  message("Transcrito canonico: ", canonical)
  return(canonical)
}
