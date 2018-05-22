context("run supertest")

test_that("summarizeOmicsResByMinPvalue", {
  expect_equal(summarizeOmicsResByMinPvalue(2:3, mat=matrix(c(1,2,4,1,2,5), nrow=2)),
               c(2,1))
})

test_that("minOrNA", {
  expect_equal(minOrNA(c(1,5,0.1,NA)),0.1)
  expect_true(is.na(minOrNA(c(NA,NA,NA))))
})

test_that("createSignificanceMask", {
  l <- list(a=c(0.03, 0.06, 1), b=c(0.6, 0.01, 0.5))
  expected <- data.frame(a=c(1,0,0), b=c(0,1,0), stringsAsFactors = F)
  expect_identical(createSignificanceMask(l), expected)
})
