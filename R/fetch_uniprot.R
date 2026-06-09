#' Fetch protein features from UniProt REST API
#'
#' @param gene_symbol simbolo del gen (ej. "KMT2B")
#' @return data.frame compatible con parse_features()

fetch_uniprot_features <- function(gene_symbol) {

  library(httr)
  library(jsonlite)
  library(dplyr)

  message("Buscando features en UniProt para: ", gene_symbol)

  #---------------------------
  # 1. Obtener accession UniProt
  #---------------------------
  resp_search <- GET(
    "https://rest.uniprot.org/uniprotkb/search",
    query = list(
      query  = paste0("gene:", gene_symbol,
                      " AND organism_id:9606 AND reviewed:true"),
      fields = "accession,protein_name,length",
      format = "json",
      size   = 1
    )
  )

  if (http_error(resp_search)) {
    stop("Error al consultar UniProt para: ", gene_symbol)
  }

  search_data <- fromJSON(content(resp_search, as="text", encoding="UTF-8"),
                          simplifyVector=FALSE)

  if (length(search_data$results) == 0) {
    stop("No se encontro entrada UniProt para: ", gene_symbol)
  }

  accession    <- search_data$results[[1]]$primaryAccession
  protein_name <- tryCatch(
    search_data$results[[1]]$proteinDescription$recommendedName$fullName$value,
    error = function(e) gene_symbol
  )

  message("  Accession UniProt: ", accession, " (", protein_name, ")")

  #---------------------------
  # 2. Obtener features
  #---------------------------
  resp_entry <- GET(
    paste0("https://rest.uniprot.org/uniprotkb/", accession),
    query = list(format = "json")
  )

  if (http_error(resp_entry)) {
    stop("Error al obtener entrada UniProt: ", accession)
  }

  entry <- fromJSON(content(resp_entry, as="text", encoding="UTF-8"),
                    simplifyVector=FALSE)

  features_raw <- entry$features
  protein_len  <- entry$sequence$length

  message("  Longitud proteica: ", protein_len, " aa")
  message("  Features totales en UniProt: ", length(features_raw))

  #---------------------------
  # 3. Mapear tipos UniProt a VizMut
  #---------------------------
  type_map <- list(
    "Domain"            = "domain",
    "Region"            = "domain",
    "Transmembrane"     = "domain",
    "Topological domain"= "domain",
    "DNA binding"       = "motif",
    "Binding site"      = "motif",
    "Active site"       = "motif",
    "Site"              = "motif",
    "Zinc finger"       = "Zinc finger",
    "Modified residue"  = "ptm",
    "Cross-link"        = "ptm",
    "Lipidation"        = "ptm",
    "Glycosylation"     = "ptm"
  )

  # tipos a incluir
  include_types <- names(type_map)

  result <- list()

  for (feat in features_raw) {
    ftype <- tryCatch(feat$type, error=function(e) NULL)
    if (is.null(ftype)) next
    if (!ftype %in% include_types) next

    fname <- tryCatch(feat$description, error=function(e) ftype)
    if (is.null(fname) || fname == "") fname <- ftype
    # truncar nombres muy largos
    if (nchar(fname) > 30) fname <- paste0(substr(fname, 1, 28), "..")



    fstart <- tryCatch(as.integer(feat$location$start$value),
                       error=function(e) NA)
    fend   <- tryCatch(as.integer(feat$location$end$value),
                       error=function(e) NA)

    if (is.na(fstart) || is.na(fend)) next
    if (fstart > fend) next

    result[[length(result) + 1]] <- data.frame(
      feature_type = type_map[[ftype]],
      feature_name = fname,
      start        = fstart,
      end          = fend,
      stringsAsFactors = FALSE
    )
  }

  if (length(result) == 0) {
    warning("No se encontraron features relevantes en UniProt para: ", gene_symbol)
    return(NULL)
  }

  features_df <- do.call(rbind, result)

  message("  Features obtenidos: ", nrow(features_df),
          " (domains: ", sum(features_df$feature_type == "domain"),
          ", motifs: ", sum(features_df$feature_type == "motif"),
          ", ptm: ", sum(features_df$feature_type == "ptm"),
          ", zinc fingers: ", sum(features_df$feature_type == "Zinc finger"), ")")

  attr(features_df, "accession")   <- accession
  attr(features_df, "protein_len") <- protein_len

  return(features_df)
}
