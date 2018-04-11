\name{covarianceEstimateFunctions}
\alias{estimateExprCov}
\alias{estimateCov}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions for covariance estimates
}
\description{
Estimate covariance matrices.
}
\usage{
estimateExprCov(expr, shrink)
estimateCov(exp1, exp2, shrink)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expr}{a numeric matrix}
  \item{exp1}{a numeric matrix}
  \item{exp2}{a numeric matrix}
  \item{shrink}{logical: should the covariance be shrinked or not.}
}
\details{
Function for internal usage by clip legacy R packages.
}
\value{Covariances extimates.}
\author{
Paolo Martini
}
