context("makeTheDataset")
require(graphite)
require(igraph)
pathway <- pathways("hsapiens", "kegg")[["EGFR tyrosine kinase inhibitor resistance"]]
pathway <- filterDirected(pathway)
survAnn <- simulateData(beta=list(0.6), x=list(c("bern", 0.4)))

dataset <- makeTheExpressionDataset(pathway,
                                    pathBoundaries = list(start="ENTREZID:3082",end="ENTREZID:6774"),
                                    deInChain = "ENTREZID:3716", deInPath = NULL,
                                    w=5, ann = survAnn)

test_that("survAnnotation",{
  expect_equal(nrow(survAnn), 300)
})

test_that("expression", {
  expect_equal(ncol(dataset$exprs), 300)
  expect_equal(nrow(dataset$exprs), length(graphite::nodes(dataset$graph)))
  expect_identical(colnames(dataset$exprs), row.names(survAnn))
})

i_means <- list("ENTREZID:3480" = c(0.8, 0.4))
dataset_met <- makeTheMethylationDataset(pathway, names(i_means), i_means, ann = survAnn)
sel_cl1 <- survAnn$x==1
m_cl1 <- mean(dataset_met$exprs[names(i_means), sel_cl1])
m_cl2 <- mean(dataset_met$exprs[names(i_means), !sel_cl1])

test_that("methylation", {
  expect_equal(ncol(dataset_met$exprs), 300)
  expect_equal(nrow(dataset_met$exprs), length(graphite::nodes(dataset_met$graph)))
  expect_identical(colnames(dataset_met$exprs), row.names(survAnn))
  expect_true(m_cl1 - m_cl2 > 0.2)
})

genes = c("ENTREZID:558","ENTREZID:1956", "ENTREZID:2064", "ENTREZID:5335")
patients_freqs <- c(0.05, 0.1, 0.2)
dataset_mut <- makeTheMutationDataset(pathway, genes, patients_freqs, ann = survAnn)

test_that("mutation", {
  expect_equal(ncol(dataset_mut$exprs), 300)
  expect_equal(nrow(dataset_mut$exprs), length(graphite::nodes(dataset_mut$graph)))
  expect_identical(colnames(dataset_mut$exprs), row.names(survAnn))
})

