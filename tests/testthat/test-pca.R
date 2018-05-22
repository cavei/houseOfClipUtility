context("PCA")
require(graph)
set.seed(1234)
exprs <- randomExpression(3, 2, dimension=c(20,20))
nodes <- letters[1:20]
degrees <- sample(c(0,1,2,3), 20, replace=T)
names(degrees) <- nodes
row.names(exprs) <- nodes
graph <- graph::randomNodeGraph(degrees)
cliques <- extractCliquesFromDag(graph)

test_that("PCA regular", {
  pcs <- computePCs(t(exprs), shrink = FALSE, method = "regular", maxPCs = 2)
  expect_equal(ncol(pcs$x), 2)
})

test_that("PCA sparse", {
  pcs <- computePCs(t(exprs), shrink = FALSE, method = "sparse", maxPCs = 2)
  expect_equal(ncol(pcs$x), 2)
})

test_that("PCA topological", {
  pcs <- computePCs(t(exprs), shrink = TRUE, method = "topological", cliques = cliques, maxPCs = 2)
  expect_equal(ncol(pcs$x), 2)
})
