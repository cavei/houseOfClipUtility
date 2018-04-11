\name{makePositiveDefinite}
\alias{makePositiveDefinite}
\title{Return a matrix that is positve and definite.}
\description{
Given a covariance matrix, it retruns a positive and defined matrix.
}
\usage{
makePositiveDefinite(m1, m2=NULL, m3=NULL, threshold=0.1)
}

\arguments{
  \item{m1}{a matrix}
  \item{m2}{an optional matrix; default=NULL}
  \item{m3}{an optional matrix; default=NULL}
  \item{threshold}{default=0.1}
}
\details{
This is a function meant to be used only for clip legacy.
}
\value{
A list; list(m1=m1n, m2=m2n, m3=m3n, correction=correction, value=value)
}
\author{
Paolo Martini
}
