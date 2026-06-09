#' Write a run log file documenting the pipeline execution
#'
#' @param log_data list con toda la informacion del run
#' @param log_path ruta del archivo log

write_run_log <- function(log_data, log_path) {

  dir.create(dirname(log_path), showWarnings=FALSE, recursive=TRUE)

  lines <- c(
    "================================================================",
    "VizMut-lolliplot — Run Log",
    paste("Date:   ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste("Version:", "1.3.0"),
    "================================================================",
    "",
    "=== Command ===",
    log_data$command,
    "",
    "=== Input ===",
    paste("Variants file: ", log_data$variants_file),
    paste("Plot type:     ", log_data$plot_type),
    paste("Gene name:     ", log_data$gene_name),
    ""
  )

  # enriquecimiento
  if (!is.null(log_data$enrich)) {
    lines <- c(lines,
      "=== Enrichment summary ===",
      paste("Total variants:            ", log_data$enrich$total),
      paste("Found in ClinVar:          ", log_data$enrich$found_clinvar),
      paste("Not found in ClinVar:      ", log_data$enrich$not_found_clinvar)
    )
    if (length(log_data$enrich$not_found_list) > 0) {
      lines <- c(lines, "  Variants not in ClinVar:")
      for (v in log_data$enrich$not_found_list) {
        lines <- c(lines, paste("    -", v))
      }
    }
    lines <- c(lines, "")
  }

  # gnomAD
  if (!is.null(log_data$gnomad)) {
    lines <- c(lines,
      "=== gnomAD v4 summary ===",
      paste("Variants with frequency data:      ", log_data$gnomad$found),
      paste("Variants with AF=0 (rare/P/LP):    ", log_data$gnomad$af_zero),
      paste("Variants not in gnomAD:            ", log_data$gnomad$not_found)
    )
    if (length(log_data$gnomad$not_found_list) > 0) {
      lines <- c(lines, "  Variants not in gnomAD:")
      for (v in log_data$gnomad$not_found_list) {
        lines <- c(lines, paste("    -", v))
      }
    }
    lines <- c(lines, "")
  }

  # liftover
  if (!is.null(log_data$liftover)) {
    lines <- c(lines,
      "=== Liftover summary ===",
      paste("Coordinates converted hg19->hg38:  ", log_data$liftover$success),
      paste("Failed liftover:                   ", log_data$liftover$failed),
      ""
    )
  }

  # exon selection
  if (!is.null(log_data$exons)) {
    lines <- c(lines,
      "=== Exon selection ===",
      paste("Selected exons:            ", log_data$exons$selected),
      paste("Total exons in transcript: ", log_data$exons$total),
      paste("Variants excluded:         ", log_data$exons$excluded),
      ""
    )
  }

  # variantes excluidas
  if (!is.null(log_data$excluded)) {
    lines <- c(lines,
      "=== Variants excluded ===",
      paste("Outside transcript:        ", log_data$excluded$outside_transcript),
      ""
    )
  }

  # output files
  lines <- c(lines,
    "=== Output files ===",
    paste("Plot:     ", log_data$output_plot)
  )
  if (!is.null(log_data$output_enrich)) {
    lines <- c(lines, paste("Enriched: ", log_data$output_enrich))
  }
  lines <- c(lines,
    paste("Log:      ", log_path),
    "",
    "================================================================"
  )

  writeLines(lines, log_path)
  message("Log guardado en: ", log_path)
}
