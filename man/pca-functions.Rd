\name{pca-functions}
\alias{computePCs}
\alias{topoCompPCs}
\alias{sparseCompPCs}
\alias{compPCs}

\title{Multiple way to compute PCA.}
\usage{
computePCs(exp, shrink, method=c("regular", "topological", "sparse"),
cliques=NULL, maxPCs)
topoCompPCs(exp, shrink, cliques, k)
sparseCompPCs(exp, shrink, k)
compPCs(exp, shrink, k)
}
\arguments{
  \item{exp}{expression matrix}
  \item{shrink}{whether to always shrink or not; boolean}
  \item{method}{either regular, topological, sparse}
  \item{cliques}{if NULL, no topology.}
  \item{maxPCs}{mac pcs returned.}
  \item{k}{the number of PC to consider}
}
\description{
Performs different types of PCA. Three methods are allowed so far: regular implement simple PCA; topological implements qpIPF penalization in the covariance matrix; sparse implements a sparce PCA analysis.
}
