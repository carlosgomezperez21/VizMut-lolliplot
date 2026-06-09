library(testthat)
source("/home/carlos/VizMut-lolliplot/R/fetch_uniprot.R")

test_that("fetch_uniprot_features devuelve data.frame para gen conocido", {
  skip_on_ci()
  result <- fetch_uniprot_features("GJB1")
  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) > 0)
  expect_true(all(c("feature_type", "feature_name", "start", "end") %in% names(result)))
})

test_that("fetch_uniprot_features devuelve longitud proteica correcta para GJB1", {
  skip_on_ci()
  result <- fetch_uniprot_features("GJB1")
  expect_equal(attr(result, "protein_len"), 283)
})

test_that("fetch_uniprot_features devuelve longitud proteica correcta para KMT2B", {
  skip_on_ci()
  result <- fetch_uniprot_features("KMT2B")
  expect_equal(attr(result, "protein_len"), 2715)
})

test_that("fetch_uniprot_features falla para gen inexistente", {
  skip_on_ci()
  expect_error(fetch_uniprot_features("GENEINEXISTENTE99999"))
})

test_that("fetch_uniprot_features devuelve tipos validos", {
  skip_on_ci()
  result <- fetch_uniprot_features("KMT2B")
  valid_types <- c("domain", "motif", "ptm", "Zinc finger")
  expect_true(all(result$feature_type %in% valid_types))
})
