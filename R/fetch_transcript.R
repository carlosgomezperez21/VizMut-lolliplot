#' Fetch transcript structure from NCBI Datasets API
#'
#' @param transcript_id RefSeq transcript ID (ej. NM_014727.3)
#' @return data.frame con columnas: type, start, end, strand, chr, name, length

fetch_transcript_structure <- function(transcript_id) {

  library(httr)
  library(jsonlite)
  library(dplyr)

  #---------------------------
  # 1. Obtener Gene ID via NCBI Datasets API
  #---------------------------
  message("Buscando transcrito en NCBI: ", transcript_id)

  # extraer simbolo del gen desde el transcript_id no es posible
  # usamos esearch solo como fallback, primero intentamos Datasets
  # buscando el transcrito directamente en Datasets por accession

  resp_acc <- GET(
    paste0("https://api.ncbi.nlm.nih.gov/datasets/v2/gene/accession/",
           transcript_id, "/product_report"),
    add_headers("Accept" = "application/json")
  )

  if (!http_error(resp_acc)) {
    data_acc <- fromJSON(content(resp_acc, as="text", encoding="UTF-8"),
                         simplifyVector=FALSE)
    gene_id <- tryCatch(
      data_acc$reports[[1]]$product$gene_id,
      error = function(e) NULL
    )
  } else {
    gene_id <- NULL
  }

  # fallback a eutils si Datasets no funciona
  if (is.null(gene_id) || gene_id == "") {
    message("  Intentando via eutils...")
    get_with_retry <- function(url, max_attempts=3, wait=5) {
      for (attempt in seq_len(max_attempts)) {
        Sys.sleep(wait)
        resp <- GET(url)
        txt  <- content(resp, as="text", encoding="UTF-8")
        if (!grepl("<!DOCTYPE|Search Backend failed|ERROR", txt)) return(txt)
        message("  Reintentando (", attempt, "/", max_attempts, ")...")
        wait <- wait * 1.5
      }
      stop("NCBI no respondió después de ", max_attempts, " intentos")
    }

    txt_search <- get_with_retry(paste0(
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
      "?db=nuccore&term=", transcript_id, "&retmode=json"
    ))
    ncbi_id <- fromJSON(txt_search)$esearchresult$idlist[1]

    txt_link <- get_with_retry(paste0(
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi",
      "?dbfrom=nuccore&db=gene&id=", ncbi_id, "&retmode=json"
    ))
    link_data <- fromJSON(txt_link, simplifyVector=FALSE)
    gene_id <- tryCatch(
      link_data$linksets$linksetdbs[[1]]$links[[1]][1],
      error = function(e) stop("No se pudo obtener Gene ID para: ", transcript_id)
    )
  }

  message("Gene ID NCBI: ", gene_id)
  
  #---------------------------
  # 2. Obtener product report
  #---------------------------
  message("Obteniendo estructura del transcrito...")

  Sys.sleep(1)
  resp_prod <- GET(
    paste0("https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id/",
           gene_id, "/product_report"),
    add_headers("Accept" = "application/json")
  )

  if (http_error(resp_prod)) {
    stop("Error al obtener product report para gene_id: ", gene_id)
  }

  data <- fromJSON(content(resp_prod, as="text", encoding="UTF-8"),
                   simplifyVector=FALSE)

  #---------------------------
  # 3. Encontrar el transcrito especifico
  #---------------------------
  target_transcript <- NULL

  for (report in data$reports) {
    transcripts <- tryCatch(report$product$transcripts,
                            error = function(e) list())
    for (tx in transcripts) {
      acc      <- tryCatch(tx$accession_version, error = function(e) "")
      acc_base <- sub("\\..*", "", acc)
      tid_base <- sub("\\..*", "", transcript_id)
      if (acc == transcript_id || acc_base == tid_base) {
        target_transcript <- tx
        break
      }
    }
    if (!is.null(target_transcript)) break
  }

  if (is.null(target_transcript)) {
    stop("No se encontró el transcrito ", transcript_id,
         " en el product report.")
  }

  #---------------------------
  # 4. Extraer exones
  #---------------------------
  genomic_loc <- target_transcript$genomic_locations[[1]]

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

  acc_version <- genomic_loc$genomic_accession_version
  chr    <- ifelse(acc_version %in% names(acc_to_chr),
                   acc_to_chr[acc_version],
                   paste0("chr_", acc_version))
  strand <- ifelse(genomic_loc$genomic_range$orientation == "plus", 1, -1)

  exons_raw <- genomic_loc$exons
  exon_df   <- tibble(
    type   = "exon",
    start  = as.numeric(sapply(exons_raw, function(e) e$begin)) + 1,
    end    = as.numeric(sapply(exons_raw, function(e) e$end)),
    strand = strand,
    chr    = chr,
    name   = paste0("exon_", sapply(exons_raw, function(e) e$order))
  ) %>% arrange(start)

  #---------------------------
  # 5. Calcular UTRs desde CDS en coordenadas mRNA
  #---------------------------
  utr5_df <- tibble(type=character(), start=numeric(), end=numeric(),
                    strand=integer(), chr=character(), name=character())
  utr3_df <- utr5_df

  # obtener rango CDS en mRNA
  cds_mrna <- tryCatch(
    target_transcript$cds$range[[1]],
    error = function(e) NULL
  )

  if (!is.null(cds_mrna)) {
    cds_mrna_start <- as.numeric(cds_mrna$begin) + 1  # 1-based
    cds_mrna_end   <- as.numeric(cds_mrna$end)

    message("CDS en mRNA: ", cds_mrna_start, " - ", cds_mrna_end)

    # mapear coordenadas mRNA a genomico recorriendo exones en orden
    
   strand_val    <- exon_df$strand[1]
    exons_ordered <- if (strand_val == 1) {
      exon_df %>% arrange(start)
    } else {
      exon_df %>% arrange(desc(start))
    }

    mrna_pos <- 0  # posicion acumulada en mRNA

    utr5_regions <- list()
    utr3_regions <- list()

    for (k in seq_len(nrow(exons_ordered))) {
      ex_start  <- exons_ordered$start[k]
      ex_end    <- exons_ordered$end[k]
      ex_len    <- ex_end - ex_start + 1

      mrna_ex_start <- mrna_pos + 1
      mrna_ex_end   <- mrna_pos + ex_len

      # region de este exon en mRNA
      # overlap con 5'UTR (antes del CDS)
      if (mrna_ex_start < cds_mrna_start) {
        utr5_end_mrna <- min(mrna_ex_end, cds_mrna_start - 1)
        utr5_len      <- utr5_end_mrna - mrna_ex_start + 1

        if (strand == 1) {
          g_start <- ex_start
          g_end   <- ex_start + utr5_len - 1
        } else {
          g_end   <- ex_end
          g_start <- ex_end - utr5_len + 1
        }

        utr5_regions[[length(utr5_regions)+1]] <- tibble(
          type   = "utr5",
          start  = g_start,
          end    = g_end,
          strand = strand_val,
          chr    = chr,
          name   = "5'UTR"
        )
      }

      # overlap con 3'UTR (despues del CDS)
      if (mrna_ex_end > cds_mrna_end) {
        utr3_start_mrna <- max(mrna_ex_start, cds_mrna_end + 1)
        utr3_len        <- mrna_ex_end - utr3_start_mrna + 1
        offset          <- utr3_start_mrna - mrna_ex_start

        if (strand == 1) {
          g_start <- ex_start + offset
          g_end   <- ex_end
        } else {
          g_end   <- ex_end - offset
          g_start <- ex_start
        }

        utr3_regions[[length(utr3_regions)+1]] <- tibble(
          type   = "utr3",
          start  = g_start,
          end    = g_end,
          strand = strand,
          chr    = chr,
          name   = "3'UTR"
        )
      }

      mrna_pos <- mrna_ex_end
    }

    if (length(utr5_regions) > 0)
      utr5_df <- bind_rows(utr5_regions)
    if (length(utr3_regions) > 0)
      utr3_df <- bind_rows(utr3_regions)

    message("UTRs calculados: ",
            nrow(utr5_df), " regiones 5'UTR, ",
            nrow(utr3_df), " regiones 3'UTR")
  } else {
    message("No se encontraron coordenadas del CDS en mRNA")
  }

  #---------------------------
  # 6. Intrones
  #---------------------------
  intron_df <- tibble(type=character(), start=numeric(), end=numeric(),
                      strand=integer(), chr=character(), name=character())

  if (nrow(exon_df) > 1) {
    intron_df <- tibble(
      type   = "intron",
      start  = exon_df$end[-nrow(exon_df)] + 1,
      end    = exon_df$start[-1] - 1,
      strand = strand,
      chr    = chr,
      name   = paste0("intron_", seq_len(nrow(exon_df)-1))
    )
  }

  #---------------------------
  # 7. Combinar
  #---------------------------
  structure_df <- bind_rows(exon_df, utr5_df, utr3_df, intron_df) %>%
    mutate(length = end - start + 1) %>%
    arrange(start)

  message("Transcrito obtenido: ",
          nrow(exon_df), " exones, ",
          nrow(intron_df), " intrones, ",
          nrow(utr5_df) + nrow(utr3_df), " UTRs")

  return(structure_df)
}
