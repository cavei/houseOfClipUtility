\name{superExactTest}
\alias{runSupertest}
\alias{plotSupertest}
\alias{createSignificanceMask}
\alias{computeSupertestFrequencies}
\alias{plotSupertestFrequencies}
\alias{minOrNA}
\alias{summarizeOmicsResByMinPvalue}

\title{
  Super Exact Test utilitues.
}
\description{
  Given a list of sets, perform multiple set intersection and given the background perform enrichment. Many function to handle and plot this kind of data.
}
\usage{
runSupertest(list, universeSize, thr = 0.05)
plotSupertest(mset)
createSignificanceMask(list, thr=0.05)
computeSupertestFrequencies(elementsIntersections)
plotSupertestFrequencies(supertestFreq)
minOrNA(x)
summarizeOmicsResByMinPvalue(col, mat)
}
\arguments{
  \item{list}{a list of set to be intersected.}
  \item{thr}{optional threshold value}
  \item{universeSize}{the universe size.}
  \item{mset}{mset class object}
  \item{elementsIntersections}{object produced by summary(mset)}
  \item{supertestFreq}{object produced by computeSupertestFrequencies}
  \item{x}{a numeric vector}
  \item{col}{col indexes or names, min by row}
  \item{mat}{a numeric matrix}
}
\value{
  a mset class object.
}
\author{
Paolo Martini
}
