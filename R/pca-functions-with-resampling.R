#' Compute PCs.
#'
#' For internal use only. Performs Principal Componenent analysis.
#'
#' Three methods are implemented:
#'   * regular: a regular PCA ('prcomp')
#'   * topological: PCA using a pathway topology.
#'   * sparse: sparse PCA analysis implemented by 'elasticnet'
#'
#' @inheritParams topoCompPCsWithResampling
#' @param method one of 'regular', 'topological' and 'sparse'
#' @param maxPCs the maximum number of PCs to consider
#'
#' @return a list with the following elements:
#'   \item{x}{the computed PCs}
#'   \item{sdev}{the standard deviation captured by the PCs}
#'   \item{loadings}{the loadings}
#'
#' @examples
#'   fakeExp <- randomExpression(4)
#'   computePCs(t(fakeExp))
#'
#' @importFrom FactoMineR estim_ncp
#'
#' @export
#'
computePCsWithResampling <- function(resampling, exp, shrink=FALSE, method=c("regular", "topological", "sparse"),
                       cliques=NULL, maxPCs=3) {
  k<- min(FactoMineR::estim_ncp(exp,scale=FALSE,ncp.min=1)$ncp, maxPCs)
  switch(method[1],
         "regular"     = compPCsWithResampling(resampling, exp=exp, shrink=shrink, k=k),
         "topological" = topoCompPCsWithResampling(resampling, exp=exp, shrink=shrink, cliques=cliques, k=k),
         "sparse"      = sparseCompPCsWithResampling(resampling, exp=exp, shrink=shrink, k=k))
}

#' Topological PCA
#'
#' @param resampling list of resampled columns
#'
#' @inheritParams topoCompPCs
#'
#' @rdname computePCs
#'
#' @importFrom qpgraph qpIPF
#' @export
#'
topoCompPCsWithResampling <- function(resampling, exp, shrink, cliques, k) {
  if (is.null(cliques))
    stop("Cliques argument is needed")
  nms <- colnames(exp)

  resamplesEigenMats <- lapply(resampling, function(sample) {
    expSampled <- exp[sample, ,drop=T]
    covmat <- estimateExprCov(expSampled, shrink) ## Consider collapse with the following line!
    covmat <- makePositiveDefinite(covmat)$m1
    cliquesIdx <- lapply(cliques, function(c) match(c, row.names(covmat)))
    if (any(sapply(cliquesIdx, function(x) {any(is.na(x))})))
      stop("Some genes in the cliques are not present as expression.")
    scovmat <- qpgraph::qpIPF(covmat, cliquesIdx)
    pcCov <- base::eigen(scovmat)
    eigenvector <- pcCov$vectors[, seq_len(k), drop=F]
    eigenvector
  })
  eigenvector <- averageReplicateMatrices(resamplesEigenMats, na.rm=T)
  scalee <- scale(exp, scale=FALSE)
  scores <- scalee%*%eigenvector
  colnames(scores) <- paste0("PC", seq_len(k))
  colnames(eigenvector) <- paste0("PC", seq_len(k))
  row.names(eigenvector) <- nms
  sd<-apply(scores, 2, sd)
  return(list(x=scores, sdev=sd, loadings=eigenvector))
}

#' Sparse PCA
#'
#' @inheritParams topoCompPCsWithResampling
#'
#' @rdname computePCs
#' @importFrom elasticnet spca
#' @export
#'
sparseCompPCsWithResampling <- function(resampling, exp, shrink, k) {
  nms <- colnames(exp)

  resamplesEigenMats <- lapply(resampling, function(sample) {
    expSampled <- exp[sample, ,drop=T]
    covmat <- estimateExprCov(expSampled, shrink) ## Consider collapse with the following line!
    covmat <- makePositiveDefinite(covmat)$m1
    paraSingle <- min(round((NCOL(expSampled)/2)),5) ## Parametri fissi da valutare
    pcCov <- elasticnet::spca(covmat, K =k, para = rep(paraSingle,k), type = "Gram", sparse = "varnum")
    eigenvector <- pcCov$loadings[, seq_len(k), drop=F]
    eigenvector
  })
  eigenvector <- averageReplicateMatrices(resamplesEigenMats, na.rm=T)
  scalee <- scale(exp, scale=FALSE)
  scores <- scalee%*%eigenvector
  colnames(scores) <- paste0("PC", seq_len(k))
  colnames(eigenvector) <- paste0("PC", seq_len(k))
  row.names(eigenvector) <- nms
  sd<-apply(scores, 2, sd)
  return(list(x=scores, sdev=sd, loadings=eigenvector))
}

#' Regular PCA
#'
#' @inheritParams topoCompPCsWithResampling
#'
#' @rdname computePCs
#' @export
#'
compPCsWithResampling <- function(resampling, exp, shrink, k) {
  nms <- colnames(exp)

  resamplesEigenMats <- lapply(resampling, function(sample) {
    expSampled <- exp[sample, ,drop=T]
    covmat <- estimateExprCov(expSampled, shrink) ## Consider collapse with the following line!
    covmat <- makePositiveDefinite(covmat)$m1
    scovmat<-covmat
    pcCov <- base::eigen(scovmat)
    eigenvector <- pcCov$vectors[, seq_len(k), drop=F]
    eigenvector
  })

  eigenvector <- averageReplicateMatrices(resamplesEigenMats, na.rm=T)
  scalee <- scale(exp, scale=FALSE)
  scores <- scalee%*%eigenvector
  colnames(scores) <- paste0("PC", seq_len(k))
  colnames(eigenvector) <- paste0("PC", seq_len(k))
  row.names(eigenvector) <- nms
  sd<-apply(scores, 2, sd)
  return(list(x=scores, sdev=sd, loadings=eigenvector))
}
