#' Estimate Single Covariance Matrix
#'
#' For internal use only. Estimate Covariance from one matrix
#'
#' @param expr a numeric matrix
#' @inheritParams estimateCov
#'
#' @return a covariance matrix
#' @importFrom stats cov
#' @importFrom corpcor cov.shrink
#' @export
#'
estimateExprCov <- function(expr, shrink) {
  if (shrink)
    unclass(corpcor::cov.shrink(expr, verbose=FALSE))
  else
    cov(expr)
}

#' Estimate Covariance
#'
#' For internal use only. Estimate 2 covariance matrixes
#'
#' @param exp1 first numeric matrix
#' @param exp2 second numeric matrix
#' @param shrink logical wheter to shrink the matrix
#'
#' @return a list
#'   \item{s1}{covariance of exp1}
#'   \item{s2}{covariance of exp2}
#'   \item{s}{combined covariance of exp1 and exp2}
#'
#' @export
#'
estimateCov <- function(exp1, exp2, shrink) {
  ncl1 <- NROW(exp1)
  ncl2 <- NROW(exp2)

  cov1 <- estimateExprCov(exp1, shrink)
  cov2 <- estimateExprCov(exp2, shrink)

  s <- (cov1*(ncl1-1) + cov2*(ncl2-1)) / (ncl1 + ncl2 - 2)
  list(s1=cov1, s2=cov2, s=s)
}
