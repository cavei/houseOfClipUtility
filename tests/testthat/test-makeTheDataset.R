context("makeTheDataset")
require(graphite)
require(igraph)
pathway <- pathways("hsapiens", "kegg")[["EGFR tyrosine kinase inhibitor resistance"]]
survAnn <- simulateData(beta=list(0.6), x=list(c("bern", 0.4)))
dataset <- makeTheDataset(pathway, pathBoundaries = list(start="ENTREZID:3082",end="ENTREZID:6774"),
                          deInChain = "ENTREZID:3716", deInPath = NULL,
                          w=5, ann = survAnn)

test_that("survAnnotation",{
  expect_equal(nrow(survAnn), 300)
})

test_that("theDataset", {
  expect_equal(ncol(dataset$exprs), 300)
  expect_equal(nrow(dataset$exprs), length(graphite::nodes(dataset$graph)))
})

