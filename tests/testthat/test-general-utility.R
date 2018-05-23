context("General Utility")

test_that("mapFromDict", {
  res = mapFromDict("a", c(a = 5))
  expect_equal(res, c(a=5))
})

require("org.Hs.eg.db")
test_that("convertIds", {
  s <- convertIds("ENTREZID:7157", graphiteStyle=T)
  expect_equal(s, paste0("SYMBOL:TP53"))
})

test_that("reverseDict", {
  res <- lapply(reverseDict(list(a = c("tp53", "tp53a"), b="tp53b")), identity)
  expected <- list(tp53="a", tp53a="a", tp53b="b")
  expect_identical(res, expected)
})

test_that("randomExpression", {
  expect_equal(dim(randomExpression(4,10)),c(10,4))
})

