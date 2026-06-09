library(testthat)
source("/home/carlos/VizMut-lolliplot/R/write_run_log.R")

test_that("write_run_log crea el archivo correctamente", {
  log_path <- tempfile(fileext = ".log")
  log_data <- list(
    command       = "Rscript main.R --test",
    variants_file = "data/test.csv",
    plot_type     = "single_gene",
    gene_name     = "GJB1",
    output_plot   = "output/test.png",
    output_enrich = NULL
  )
  write_run_log(log_data, log_path)
  expect_true(file.exists(log_path))
  on.exit(unlink(log_path))
})

test_that("write_run_log incluye informacion basica", {
  log_path <- tempfile(fileext = ".log")
  log_data <- list(
    command       = "Rscript main.R --test",
    variants_file = "data/test.csv",
    plot_type     = "single_gene",
    gene_name     = "GJB1",
    output_plot   = "output/test.png",
    output_enrich = NULL
  )
  write_run_log(log_data, log_path)
  content <- readLines(log_path)
  expect_true(any(grepl("VizMut-lolliplot", content)))
  expect_true(any(grepl("GJB1", content)))
  expect_true(any(grepl("single_gene", content)))
  on.exit(unlink(log_path))
})

test_that("write_run_log incluye resumen de enriquecimiento", {
  log_path <- tempfile(fileext = ".log")
  log_data <- list(
    command       = "Rscript main.R --test",
    variants_file = "data/test.csv",
    plot_type     = "single_gene",
    gene_name     = "GJB1",
    output_plot   = "output/test.png",
    output_enrich = NULL,
    enrich = list(
      total             = 10,
      found_clinvar     = 8,
      not_found_clinvar = 2,
      not_found_list    = c("GJB1:c.285C>T", "GJB1:c.477G>A")
    )
  )
  write_run_log(log_data, log_path)
  content <- readLines(log_path)
  expect_true(any(grepl("Enrichment summary", content)))
  expect_true(any(grepl("GJB1:c.285C>T", content)))
  on.exit(unlink(log_path))
})
