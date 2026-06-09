library(testthat)
source("/home/carlos/VizMut-lolliplot/R/parse_exon_selection.R")

test_that("parse_exon_selection parsea numeros simples", {
  result <- suppressMessages(parse_exon_selection("1,3,5"))
  expect_equal(result, c(1, 3, 5))
})

test_that("parse_exon_selection parsea rangos", {
  result <- suppressMessages(parse_exon_selection("1-3"))
  expect_equal(result, c(1, 2, 3))
})

test_that("parse_exon_selection parsea combinacion", {
  result <- suppressMessages(parse_exon_selection("1,5-7,10"))
  expect_equal(result, c(1, 5, 6, 7, 10))
})

test_that("parse_exon_selection devuelve NULL para all", {
  result <- parse_exon_selection("all")
  expect_null(result)
})

test_that("parse_exon_selection devuelve NULL para NULL", {
  result <- parse_exon_selection(NULL)
  expect_null(result)
})

test_that("parse_exon_selection elimina duplicados y ordena", {
  result <- suppressMessages(parse_exon_selection("3,1,2,1"))
  expect_equal(result, c(1, 2, 3))
})
