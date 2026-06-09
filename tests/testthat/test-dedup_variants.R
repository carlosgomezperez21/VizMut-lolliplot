library(testthat)
source("/home/carlos/VizMut-lolliplot/R/dedup_variants.R")

test_that("dedup_variants deduplica correctamente", {
  df <- data.frame(
    variant_id = c("cv_1", "cv_2", "cv_3"),
    hgvs_c     = c("KMT2B:c.252G>A", "KMT2B:c.252G>A", "KMT2B:c.500A>T"),
    ACMG       = c("P", "P", "VUS"),
    stringsAsFactors = FALSE
  )
  result <- suppressMessages(dedup_variants(df))
  expect_equal(nrow(result), 2)
  expect_true("count" %in% names(result))
})

test_that("dedup_variants cuenta frecuencias correctamente", {
  df <- data.frame(
    variant_id = c("cv_1", "cv_2", "cv_3", "cv_4", "cv_5"),
    hgvs_c     = c("A", "A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )
  result <- suppressMessages(dedup_variants(df))
  expect_equal(result$count[result$hgvs_c == "A"], 3)
  expect_equal(result$count[result$hgvs_c == "B"], 2)
})

test_that("dedup_variants sin duplicados conserva todos", {
  df <- data.frame(
    variant_id = c("cv_1", "cv_2"),
    hgvs_c     = c("A", "B"),
    stringsAsFactors = FALSE
  )
  result <- suppressMessages(dedup_variants(df))
  expect_equal(nrow(result), 2)
  expect_true(all(result$count == 1))
})
