% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methylation_dataset_generator.R
\name{makeTheMethylationDataset}
\alias{makeTheMethylationDataset}
\alias{makeUniformMethylationDataset}
\title{Simulate a methylation dataset with survival annotation}
\usage{
makeTheMethylationDataset(pathway, impactThisGenes, impactMeans, ann,
  omicName = "x")

makeUniformMethylationDataset(pathway, ann, omicName = "x")
}
\arguments{
\item{pathway}{a Pathway}

\item{impactThisGenes}{the DEG in the aberrant signal chain}

\item{impactMeans}{list of the same length of impactThisGenes with lower and upper bounds for simulation}

\item{ann}{survival annotations}

\item{omicName}{the name of the omics to consider}
}
\value{
list
  \item{exprs}{expression}
  \item{annotation}{the survival annotation}
  \item{graph}{a graphNEL}
  \item{chain}{the chain of selectP}
}
\description{
To simulate a dataset we need a real pathway.

To simulate a dataset we need a real pathway.
}
