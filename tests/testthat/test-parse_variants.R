library(testthat)
source("/home/carlos/VizMut-lolliplot/R/parse_variants.R")

test_that("parse_variants acepta CSV valido", {
  df <- data.frame(
    variant_id     = c("cv_1", "cv_2"),
    gen_position   = c("chr19-35718270 G>A", "chr19-35720473 C>T"),
    protein_change = c("p.Trp84Ter", "p.Lys376del"),
    hgvs_c         = c("NM_014727.3:c.252G>A", "NM_014727.3:c.1126del"),
    variant_type   = c("SNV", "DEL"),
    ACMG           = c("P", "B"),
    stringsAsFactors = FALSE
  )
  result <- parse_variants(df)
  expect_s3_class(result, "data.frame")
  expect_true("pos" %in% names(result))
  expect_true("chr" %in% names(result))
  expect_equal(nrow(result), 2)
})

test_that("parse_variants extrae chr y pos correctamente", {
  df <- data.frame(
    variant_id     = "cv_1",
    gen_position   = "chr19-35718270 G>A",
    protein_change = "p.Trp84Ter",
    hgvs_c         = "NM_014727.3:c.252G>A",
    variant_type   = "SNV",
    ACMG           = "P",
    stringsAsFactors = FALSE
  )
  result <- parse_variants(df)
  expect_equal(result$chr, "chr19")
  expect_equal(result$pos, 35718270)
})

test_that("parse_variants maneja NA en gen_position", {
  df <- data.frame(
    variant_id     = "cv_1",
    gen_position   = NA,
    protein_change = NA,
    hgvs_c         = "KMT2B:c.252G>A",
    variant_type   = "SNV",
    ACMG           = NA,
    stringsAsFactors = FALSE
  )
  result <- parse_variants(df)
  expect_true(is.na(result$pos))
})
