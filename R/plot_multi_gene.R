#' Plot multi-gene lolliplot with genomic structure
#'
#' @param variants data.frame con columna 'gene' que separa los genes
#' @param gene_list vector de genes a plotear o "all"
#' @param genome ensamble: "hg38" | "hg19"
#' @param grid logical, dividir por variant_type
#' @param label_type "protein" | "cds"
#' @return patchwork plot object

plot_multi_gene <- function(variants,
                            gene_list  = "all",
                            genome     = "hg38",
                            grid       = FALSE,
                            label_type = "protein") {

  library(dplyr)
  library(patchwork)

  source("R/fetch_canonical_transcript.R")
  source("R/fetch_transcript.R")
  source("R/fetch_cytobands.R")
  source("R/plot_ideogram.R")
  source("R/plot_gene_lolliplot.R")

  #---------------------------
  # 1. Determinar genes a plotear
  #---------------------------
  if (!("gene" %in% names(variants))) {
    stop("El CSV debe tener una columna 'gene'")
  }

  all_genes <- unique(variants$gene)

  if (identical(gene_list, "all")) {
    genes_to_plot <- all_genes
  } else {
    genes_to_plot <- intersect(gene_list, all_genes)
    missing <- setdiff(gene_list, all_genes)
    if (length(missing) > 0) {
      warning("Genes no encontrados en el CSV: ",
              paste(missing, collapse=", "))
    }
  }

  if (length(genes_to_plot) == 0) {
    stop("No hay genes validos para plotear")
  }

  message("Genes a plotear: ", paste(genes_to_plot, collapse=", "))

  #---------------------------
  # 2. Generar plot por gen
  #---------------------------
  plots <- lapply(genes_to_plot, function(gene) {

    message("\n=== Procesando: ", gene, " ===")

    # filtrar variantes del gen
    vars_gene <- variants %>% filter(.data$gene == !!gene)

    # obtener transcrito canonico
    transcript_id <- tryCatch(
      fetch_canonical_transcript(gene),
      error = function(e) {
        message("Error obteniendo transcrito para ", gene, ": ", e$message)
        return(NULL)
      }
    )

    if (is.null(transcript_id)) return(NULL)

    # obtener estructura del transcrito
    struct <- tryCatch(
      fetch_transcript_structure(transcript_id),
      error = function(e) {
        message("Error obteniendo estructura para ", gene, ": ", e$message)
        return(NULL)
      }
    )

    if (is.null(struct)) return(NULL)

    # filtrar variantes dentro del transcrito
    tx_start    <- min(struct$start)
    tx_end      <- max(struct$end)
    vars_in     <- vars_gene[!is.na(vars_gene$pos) &
                              vars_gene$pos >= tx_start &
                              vars_gene$pos <= tx_end, ]

    n_out <- nrow(vars_gene) - nrow(vars_in)
    if (n_out > 0) {
      message("  ", n_out, " variantes fuera del transcrito excluidas")
    }

    # obtener citobandas
    chr   <- unique(struct$chr)[1]
    bands <- tryCatch(
      fetch_cytobands(chr),
      error = function(e) NULL
    )

    # generar plot
    p <- tryCatch(
      plot_gene_lolliplot(
        variants             = vars_in,
        transcript_structure = struct,
        gene_name            = gene,
        label_type           = label_type,
        grid                 = grid,
        cytobands            = bands
      ),
      error = function(e) {
        message("Error generando plot para ", gene, ": ", e$message)
        return(NULL)
      }
    )

    return(p)
  })

  # eliminar NULLs
  plots <- Filter(Negate(is.null), plots)

  if (length(plots) == 0) {
    stop("No se pudo generar ningun plot")
  }

  #---------------------------
  # 3. Combinar plots
  #---------------------------
  message("\nCombinando ", length(plots), " plots...")

  final_plot <- wrap_plots(plots, ncol = 1)

  return(final_plot)
}
