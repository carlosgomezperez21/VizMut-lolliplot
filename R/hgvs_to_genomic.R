#' Fetch variant coordinates from ClinicalTables NLM API + liftover hg19->hg38
#'
#' @param gene_symbol simbolo del gen (ej. "KMT2B")
#' @param transcript_id RefSeq transcript ID (ej. "NM_014727.3")
#' @param max_results numero maximo de variantes a obtener (default: 500)
#' @param chain_file ruta al chain file hg19ToHg38 (se descarga si no existe)
#' @return data.frame con columnas: hgvs_c, hgvs_p, chr, pos, ref, alt, variant_type

fetch_variants_clinvar <- function(gene_symbol,
                                   transcript_id,
                                   max_results  = 500,
                                   chain_file   = "/tmp/hg19ToHg38.over.chain") {

  library(httr)
  library(jsonlite)
  library(dplyr)
  library(rtracklayer)
  library(GenomicRanges)

  #---------------------------
  # 1. Descargar chain file si no existe
  #---------------------------
  chain_gz <- paste0(chain_file, ".gz")

  if (!file.exists(chain_file)) {
    message("Descargando chain file hg19->hg38...")
    chain_url <- paste0(
      "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/",
      "hg19ToHg38.over.chain.gz"
    )
    download.file(chain_url, chain_gz, quiet = TRUE)
    R.utils::gunzip(chain_gz, destname = chain_file, remove = FALSE)
  }

  chain <- import.chain(chain_file)

  #---------------------------
  # 2. Obtener variantes de ClinicalTables (hg19)
  #---------------------------
  message("Buscando variantes de ", gene_symbol,
          " (", transcript_id, ") en ClinVar...")

  url <- paste0(
    "https://clinicaltables.nlm.nih.gov/api/variants/v4/search",
    "?terms=", utils::URLencode(transcript_id, reserved = TRUE),
    "&maxList=", max_results,
    "&ef=HGVS_c,HGVS_p,Chromosome,Start,Stop,",
    "ReferenceAllele,AlternateAllele,Type,AminoAcidChange,",
    "ClinicalSignificance,PhenotypeList,dbSNP"
  )

  resp <- GET(url)
  if (http_error(resp)) stop("Error al conectar con ClinicalTables API")

  data   <- fromJSON(content(resp, as = "text", encoding = "UTF-8"),
                     simplifyVector = FALSE)
  total  <- data[[1]]
  extras <- data[[3]]
  n      <- length(extras$HGVS_c)

  message("Total variantes en ClinVar: ", total,
          " | obtenidas: ", n)

  if (n == 0) stop("No se encontraron variantes para: ", transcript_id)

  safe_get <- function(lst, i) {
    val <- tryCatch(lst[[i]], error = function(e) NA)
    if (is.null(val) || length(val) == 0) NA else val
  }

  #---------------------------
  # 3. Construir data.frame hg19
  #---------------------------
  result <- tibble(
    hgvs_c         = sapply(seq_len(n), function(i) safe_get(extras$HGVS_c, i)),
    hgvs_p         = sapply(seq_len(n), function(i) safe_get(extras$HGVS_p, i)),
    chr_hg19       = paste0("chr", sapply(seq_len(n), function(i)
                       safe_get(extras$Chromosome, i))),
    pos_hg19       = as.numeric(sapply(seq_len(n), function(i)
                       safe_get(extras$Start, i))),
    ref            = sapply(seq_len(n), function(i)
                       safe_get(extras$ReferenceAllele, i)),
    alt            = sapply(seq_len(n), function(i)
                       safe_get(extras$AlternateAllele, i)),
    variant_type   = sapply(seq_len(n), function(i)
                       safe_get(extras$Type, i)),
    protein_change = sapply(seq_len(n), function(i)
                       safe_get(extras$AminoAcidChange, i)),
    clinsig        = sapply(seq_len(n), function(i)
                       safe_get(extras$ClinicalSignificance, i)),
    phenotype      = sapply(seq_len(n), function(i)
                       safe_get(extras$PhenotypeList, i)),
    dbsnp          = sapply(seq_len(n), function(i)
                       safe_get(extras$dbSNP, i))
  ) %>%
    mutate(
      ref = ifelse(ref == "na" | ref == "", NA, ref),
      alt = ifelse(alt == "na" | alt == "", NA, alt),
      gene = gene_symbol
    )

#---------------------------
  # Mapeo ClinicalSignificance -> ACMG
  #---------------------------
  result <- result %>%
    mutate(
      ACMG = case_when(
        grepl("^Pathogenic$", clinsig, ignore.case=TRUE)                    ~ "P",
        grepl("Pathogenic/Likely pathogenic", clinsig, ignore.case=TRUE)    ~ "LP",
        grepl("^Likely pathogenic$", clinsig, ignore.case=TRUE)             ~ "LP",
        grepl("^Uncertain significance$", clinsig, ignore.case=TRUE)        ~ "VUS",
        grepl("^Likely benign$", clinsig, ignore.case=TRUE)                 ~ "LB",
        grepl("Benign/Likely benign", clinsig, ignore.case=TRUE)            ~ "LB",
        grepl("^Benign$", clinsig, ignore.case=TRUE)                        ~ "B",
        TRUE                                                                 ~ "VUS"
      )
    )
  #---------------------------
  # 4. Liftover hg19 -> hg38
  #---------------------------
  message("Realizando liftover hg19 -> hg38...")

  valid <- !is.na(result$pos_hg19)

  gr_hg19 <- GRanges(
    seqnames = result$chr_hg19[valid],
    ranges   = IRanges(
      start = result$pos_hg19[valid],
      end   = result$pos_hg19[valid]
    )
  )

  gr_hg38  <- liftOver(gr_hg19, chain)
  lifted   <- sapply(gr_hg38, function(x) {
    if (length(x) == 1) start(x) else NA
  })

  result$pos <- NA_real_
  result$chr <- NA_character_
  result$pos[valid] <- lifted
  result$chr[valid] <- result$chr_hg19[valid]

  # reporte
  n_lifted  <- sum(!is.na(result$pos))
  n_failed  <- sum(is.na(result$pos))
  message("Liftover exitoso: ", n_lifted, " | fallido: ", n_failed)

  result <- result %>%
    select(hgvs_c, hgvs_p, chr, pos, ref, alt,
           variant_type, protein_change, clinsig,
           ACMG, phenotype, dbsnp, gene) %>%
    filter(!is.na(pos))

  return(result)
}
