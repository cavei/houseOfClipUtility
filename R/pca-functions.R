topoCompPCs <- function(exp, shrink, cliques, k) {
  if (is.null(cliques))
    stop("Cliques argument is needed")
  nms <- colnames(exp)
  covmat <- estimateExprCov(exp, shrink)
  covmat <- makePositiveDefinite(covmat)$m1
  cliquesIdx <- lapply(cliques, function(c) match(c, row.names(covmat)))
  if (any(sapply(cliquesIdx, function(x) {any(is.na(x))})))
    stop("Some genes in the cliques are not present as expression.")
  scovmat <- qpgraph::qpIPF(covmat, cliquesIdx)
  pcCov <- base::eigen(scovmat)
  eigenvector <- pcCov$vectors[, seq_len(k), drop=F]
  scalee <- scale(exp, scale=FALSE)
  npc <- min(dim(exp))
  scores <- scalee%*%eigenvector
  colnames(scores) <- paste0("PC", seq_len(k))
  colnames(eigenvector) <- paste0("PC", seq_len(k))
  row.names(eigenvector) <- nms
  sd<-apply(scores, 2, sd)
  return(list(x=scores, sdev=sd, loadings=eigenvector))
}

sparseCompPCs <- function(exp, shrink, k) {
  nms <- colnames(exp)
  covmat <- estimateExprCov(exp, shrink)
  covmat <- makePositiveDefinite(covmat)$m1
  paraSingle <- min(round((NCOL(exp)/2)),5) ## Parametri fissi da valutare
  pcCov <- elasticnet::spca(covmat, K =k, para = rep(paraSingle,k), type = "Gram", sparse = "varnum")
  eigenvector  <- pcCov$loadings[, seq_len(k), drop=F]
  scalee <- scale(exp, scale=FALSE)
  npc <- min(dim(exp))
  scores <- scalee%*%eigenvector
  colnames(scores) <- paste0("PC", seq_len(k))
  colnames(eigenvector) <- paste0("PC", seq_len(k))
  row.names(eigenvector) <- nms
  sd<-apply(scores, 2, sd)
  return(list(x=scores, sdev=sd, loadings=eigenvector))
}

compPCs <- function(exp, shrink, k) {
  nms <- colnames(exp)
  covmat <- estimateExprCov(exp, shrink) ## Consider collapse with the following line!
  covmat <- makePositiveDefinite(covmat)$m1
  scovmat<-covmat
  pcCov <- base::eigen(scovmat)
  eigenvector <- pcCov$vectors[, seq_len(k), drop=F]
  scalee <- scale(exp, scale=FALSE)
  npc <- min(dim(exp))
  scores <- scalee%*%eigenvector
  colnames(scores) <- paste0("PC", seq_len(k))
  colnames(eigenvector) <- paste0("PC", seq_len(k))
  row.names(eigenvector) <- nms
  sd<-apply(scores, 2, sd)
  return(list(x=scores, sdev=sd, loadings=eigenvector))
}

computePCs <- function(exp, shrink, method=c("regular", "topological", "sparse"),
                       cliques=NULL, maxPCs) {
  k<- min(FactoMineR::estim_ncp(exp,scale=FALSE,ncp.min=1)$ncp, maxPCs)
  switch(method,
         regular     = compPCs(exp=exp, shrink=shrink, k=k),
         topological = topoCompPCs(exp=exp, shrink=shrink, cliques=cliques, k=k),
         sparse      = sparseCompPCs(exp=exp, shrink=shrink, k=k))
}
