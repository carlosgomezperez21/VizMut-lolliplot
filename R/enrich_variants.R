#' Enrich a minimal variant CSV with genomic coordinates,
#' ACMG classification and protein notation from ClinVar
#'
#' @param variants data.frame con columnas minimas: variant_id, hgvs_c
#'                 donde hgvs_c tiene formato "GENE:c.XXXX"
#' @param genome ensamble objetivo: "hg38" | "hg19"
#' @return data.frame enriquecido compatible con el pipeline

enrich_variants <- function(variants, genome = "hg38") {

  library(dplyr)
  library(httr)
  library(jsonlite)

  source("R/fetch_canonical_transcript.R")
  source("R/hgvs_to_genomic.R")

  #---------------------------
  # funcion auxiliar: coordenadas desde NCBI
  #---------------------------
  get_coords_from_hgvs <- function(hgvs_full) {
    tryCatch({
      url1 <- paste0(
        "https://api.ncbi.nlm.nih.gov/variation/v0/hgvs/",
        utils::URLencode(hgvs_full, reserved=TRUE),
        "/contextuals"
      )
      resp1 <- GET(url1)
      if (http_error(resp1)) return(NULL)

      data1    <- fromJSON(content(resp1, as="text", encoding="UTF-8"),
                           simplifyVector=FALSE)
      spdi     <- data1$data$spdis[[1]]
      spdi_str <- paste0(spdi$seq_id, ":", spdi$position, ":",
                         spdi$deleted_sequence, ":", spdi$inserted_sequence)

      url2  <- paste0(
        "https://api.ncbi.nlm.nih.gov/variation/v0/spdi/",
        utils::URLencode(spdi_str, reserved=TRUE), "/rsids"
      )
      resp2 <- GET(url2)
      if (http_error(resp2)) return(NULL)

      data2 <- fromJSON(content(resp2, as="text", encoding="UTF-8"),
                        simplifyVector=FALSE)
      rsids <- data2$data$rsids
      if (length(rsids) == 0) return(NULL)
      rsid  <- rsids[[1]]

      url3  <- paste0("https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/", rsid)
      resp3 <- GET(url3)
      if (http_error(resp3)) return(NULL)

      data3     <- fromJSON(content(resp3, as="text", encoding="UTF-8"),
                            simplifyVector=FALSE)
      movements <- data3$present_obs_movements

      acc_to_chr <- c(
        "NC_000001.11"="chr1",  "NC_000002.12"="chr2",
        "NC_000003.12"="chr3",  "NC_000004.12"="chr4",
        "NC_000005.10"="chr5",  "NC_000006.12"="chr6",
        "NC_000007.14"="chr7",  "NC_000008.11"="chr8",
        "NC_000009.12"="chr9",  "NC_000010.11"="chr10",
        "NC_000011.10"="chr11", "NC_000012.12"="chr12",
        "NC_000013.11"="chr13", "NC_000014.9"="chr14",
        "NC_000015.10"="chr15", "NC_000016.10"="chr16",
        "NC_000017.11"="chr17", "NC_000018.10"="chr18",
        "NC_000019.10"="chr19", "NC_000020.11"="chr20",
        "NC_000021.9"="chr21",  "NC_000022.11"="chr22",
        "NC_000023.11"="chrX",  "NC_000024.10"="chrY"
      )

      for (mov in movements) {
        allele <- tryCatch(mov$allele_in_cur_release, error=function(e) NULL)
        if (is.null(allele)) next
        seq_id <- allele$seq_id
        if (grepl("^NC_", seq_id)) {
          chr <- ifelse(seq_id %in% names(acc_to_chr),
                        acc_to_chr[seq_id], seq_id)
          return(list(
            chr  = chr,
            pos  = as.numeric(allele$position) + 1,
            ref  = allele$deleted_sequence,
            alt  = allele$inserted_sequence,
            rsid = paste0("rs", rsid)
          ))
        }
      }
      return(NULL)
    }, error = function(e) NULL)
  }

  #---------------------------
  # 1. Validar columnas minimas
  #---------------------------
  required <- c("variant_id", "hgvs_c")
  missing  <- setdiff(required, names(variants))
  if (length(missing) > 0) {
    stop("Columnas requeridas faltantes: ", paste(missing, collapse=", "))
  }

  #---------------------------
  # 2. Extraer gen y construir hgvs_c completo
  #---------------------------
  message("Extrayendo simbolos de genes desde hgvs_c...")

  variants <- variants %>%
    mutate(
      gene  = sub(":.*", "", hgvs_c),
      c_dot = sub(".*:", "", hgvs_c)
    )

  genes <- unique(variants$gene)
  message("Genes identificados: ", paste(genes, collapse=", "))

  #---------------------------
  # 3. Obtener transcrito canonico por gen
  #---------------------------
  message("Obteniendo transcritos canonicos...")

  transcript_map <- list()
  for (g in genes) {
    tx <- tryCatch(
      fetch_canonical_transcript(g),
      error = function(e) {
        message("Warning: no se encontro transcrito para ", g)
        NULL
      }
    )
    transcript_map[[g]] <- tx
  }

  variants <- variants %>%
    rowwise() %>%
    mutate(
      transcript  = transcript_map[[gene]],
      hgvs_c_full = ifelse(!is.null(transcript),
                            paste0(transcript, ":", c_dot), NA)
    ) %>%
    ungroup()

  #---------------------------
  # 4. Inicializar columnas resultado
  #---------------------------
  result                  <- variants
  result$chr              <- NA_character_
  result$pos              <- NA_real_
  result$ref              <- NA_character_
  result$alt              <- NA_character_
  result$protein_change   <- NA_character_
  result$variant_type     <- "unknown"
  result$ACMG             <- NA_character_
  result$clinsig          <- "Not found in ClinVar"
  result$phenotype        <- NA_character_
  result$dbsnp            <- NA_character_
  result$found_clinvar    <- FALSE
  result$coords_from_ncbi <- FALSE

  #---------------------------
  # 5. Buscar cada variante en ClinVar
  #---------------------------
  message("Buscando variantes en ClinVar...")

  for (i in seq_len(nrow(result))) {
    hgvs_full <- result$hgvs_c_full[i]
    if (is.na(hgvs_full)) next

    message("  Buscando: ", hgvs_full)

    url <- paste0(
      "https://clinicaltables.nlm.nih.gov/api/variants/v4/search",
      "?terms=", utils::URLencode(hgvs_full, reserved=TRUE),
      "&maxList=1",
      "&ef=HGVS_c,HGVS_p,Chromosome,Start,ClinicalSignificance,",
      "Type,AminoAcidChange,ReferenceAllele,AlternateAllele,PhenotypeList,dbSNP"
    )

    resp  <- tryCatch(GET(url), error = function(e) NULL)
    if (is.null(resp)) next

    data   <- fromJSON(content(resp, as="text", encoding="UTF-8"),
                       simplifyVector=FALSE)
    total  <- data[[1]]
    extras <- data[[3]]
    if (total == 0) next

    safe <- function(lst, key) {
      val <- tryCatch(lst[[key]][[1]], error=function(e) NA)
      if (is.null(val) || length(val) == 0) NA else val
    }

    clinsig_raw <- safe(extras, "ClinicalSignificance")
    acmg <- dplyr::case_when(
      grepl("^Pathogenic$", clinsig_raw, ignore.case=TRUE)                 ~ "P",
      grepl("Pathogenic/Likely pathogenic", clinsig_raw, ignore.case=TRUE) ~ "LP",
      grepl("^Likely pathogenic$", clinsig_raw, ignore.case=TRUE)          ~ "LP",
      grepl("^Uncertain significance$", clinsig_raw, ignore.case=TRUE)     ~ "VUS",
      grepl("^Likely benign$", clinsig_raw, ignore.case=TRUE)              ~ "LB",
      grepl("Benign/Likely benign", clinsig_raw, ignore.case=TRUE)         ~ "LB",
      grepl("^Benign$", clinsig_raw, ignore.case=TRUE)                     ~ "B",
      TRUE                                                                  ~ "VUS"
    )

    result$chr[i]            <- paste0("chr", safe(extras, "Chromosome"))
    result$pos[i]            <- as.numeric(safe(extras, "Start"))
    result$ref[i]            <- safe(extras, "ReferenceAllele")
    result$alt[i]            <- safe(extras, "AlternateAllele")
    result$protein_change[i] <- safe(extras, "AminoAcidChange")
    result$variant_type[i]   <- safe(extras, "Type")
    result$ACMG[i]           <- acmg
    result$clinsig[i]        <- clinsig_raw
    result$phenotype[i]      <- safe(extras, "PhenotypeList")
    result$dbsnp[i]          <- safe(extras, "dbSNP")
    result$found_clinvar[i]  <- TRUE
  }

  #---------------------------
  # 6. Coordenadas desde NCBI para variantes no en ClinVar
  #---------------------------
  not_found_idx <- which(result$found_clinvar == FALSE &
                         !is.na(result$hgvs_c_full))

  if (length(not_found_idx) > 0) {
    message("\nNota: ", length(not_found_idx),
            " variantes no encontradas en ClinVar.")
    message("Buscando coordenadas en NCBI Variation Services...")

    for (i in not_found_idx) {
      hgvs_full <- result$hgvs_c_full[i]
      message("  Consultando NCBI para: ", hgvs_full)

      coords <- get_coords_from_hgvs(hgvs_full)

      if (!is.null(coords)) {
        result$chr[i]            <- coords$chr
        result$pos[i]            <- coords$pos
        result$ref[i]            <- coords$ref
        result$alt[i]            <- coords$alt
        result$dbsnp[i]          <- coords$rsid
        result$variant_type[i]   <- "SNV"
        result$clinsig[i]        <- "Not in ClinVar - coords from NCBI"
        result$coords_from_ncbi[i] <- TRUE
        message("    Coordenadas obtenidas: ", coords$chr, ":", coords$pos)
      } else {
        message("    No se encontraron coordenadas para: ", hgvs_full)
      }
    }
  }

  #---------------------------
  # 7. Liftover hg19 -> hg38 solo para coordenadas de ClinVar
  #---------------------------
  clinvar_idx <- which(result$found_clinvar == TRUE)

  if (genome == "hg38" && length(clinvar_idx) > 0) {
    message("Realizando liftover hg19 -> hg38 para coordenadas ClinVar...")

    library(rtracklayer)
    library(GenomicRanges)

    chain_file <- "/tmp/hg19ToHg38.over.chain"
    chain_gz   <- paste0(chain_file, ".gz")

    if (!file.exists(chain_file)) {
      message("Descargando chain file...")
      download.file(
        "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
        chain_gz, quiet = TRUE
      )
      R.utils::gunzip(chain_gz, destname = chain_file, remove = FALSE)
    }

    chain   <- import.chain(chain_file)
    gr_hg19 <- GRanges(
      seqnames = result$chr[clinvar_idx],
      ranges   = IRanges(start = result$pos[clinvar_idx],
                         end   = result$pos[clinvar_idx])
    )
    gr_hg38 <- liftOver(gr_hg19, chain)
    lifted  <- sapply(gr_hg38, function(x) if (length(x) == 1) start(x) else NA)

    result$pos[clinvar_idx] <- lifted
    message("Liftover completado: ", sum(!is.na(lifted)), " coordenadas convertidas")
  }

  #---------------------------
  # 8. Construir gen_position
  #---------------------------
  result <- result %>%
    mutate(
      gen_position = ifelse(
        !is.na(pos),
        paste0(chr, "-", pos, " ",
               ifelse(is.na(ref), "N", ref), ">",
               ifelse(is.na(alt), "N", alt)),
        NA
      )
    )

  #---------------------------
  # 9. Resumen y output final
  #---------------------------
  n_found    <- sum(result$found_clinvar == TRUE)
  n_notfound <- sum(result$found_clinvar == FALSE)

  message("\nResumen del enriquecimiento:")
  message("  Encontradas en ClinVar: ", n_found)
  message("  No encontradas en ClinVar: ", n_notfound)

  if (n_notfound > 0) {
    not_found_vars <- result[result$found_clinvar == FALSE,
                             c("variant_id", "hgvs_c")]
    message("\nNota: las siguientes variantes no se identificaron en ClinVar")
    message("y se graficaran sin clasificacion ACMG con etiqueta c.:")
    for (i in seq_len(nrow(not_found_vars))) {
      message("  - ", not_found_vars$variant_id[i],
              ": ", not_found_vars$hgvs_c[i])
    }
  }

  result <- result %>%
    select(variant_id, gene, gen_position,
           protein_change, hgvs_c, variant_type,
           ACMG, clinsig, phenotype, dbsnp,
           chr, pos)

  return(result)
}
