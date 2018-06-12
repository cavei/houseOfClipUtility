context("PCA")
require(graph)
set.seed(1234)
exprs <- randomExpression(30,20)
nodes <- letters[1:20]
degrees <- makeDegreesEven(sample(c(0,1,2,3,4), 20, replace=T))
names(degrees) <- nodes
row.names(exprs) <- nodes
colnames(exprs) <- paste0("P",seq_len(ncol(exprs)))
graph <- graph::randomNodeGraph(degrees)
cliques <- extractCliquesFromDag(graph)

resamples <- lapply(seq_len(10), function(x) {
  sample(colnames(exprs), ncol(exprs)-2)
})

test_that("PCA regular", {
  pcs <- computePCs(t(exprs), maxPCs = 2)
  expect_equal(ncol(pcs$x), 2)
})

test_that("PCA regular Resampled", {
  pcs <- computePCsWithResampling(resampling = resamples, t(exprs), maxPCs = 2)
  expect_equal(ncol(pcs$x), 2)
})

test_that("PCA sparse", {
  pcs <- computePCs(t(exprs), method = "sparse", maxPCs = 2)
  expect_equal(ncol(pcs$x), 2)
})

test_that("PCA sparse Resampled", {
  pcs <- computePCsWithResampling(resampling = resamples, t(exprs), method = "sparse", maxPCs = 2)
  expect_equal(ncol(pcs$x), 2)
})

test_that("PCA topological", {
  pcs <- computePCs(t(exprs), shrink = TRUE, method = "topological", cliques = cliques, maxPCs = 2)
  expect_equal(ncol(pcs$x), 2)
})


test_that("PCA topological Resampled", {
  pcs <- computePCsWithResampling(resampling = resamples, t(exprs), shrink = TRUE, method = "topological", cliques = cliques, maxPCs = 2)
  expect_equal(ncol(pcs$x), 2)
})

exprs <- exprs[,1:10]
test_that("PCA topological auto shrink", {
  pcs <- computePCs(t(exprs), shrink = FALSE, method = "topological", cliques = cliques, maxPCs = 2)
  expect_equal(ncol(pcs$x), 2)
})

resamples <- lapply(seq_len(10), function(x) {
  sample(colnames(exprs), ncol(exprs)-1)
})

test_that("PCA topological Resampled auto shrink", {
  pcs <- computePCsWithResampling(resampling = resamples, t(exprs), shrink = FALSE, method = "topological", cliques = cliques, maxPCs = 2)
  expect_equal(ncol(pcs$x), 2)
})
