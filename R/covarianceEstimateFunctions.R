estimateExprCov <- function(expr, shrink) {
  if (shrink)
    unclass(cov.shrink(expr, verbose=FALSE))
  else
    cov(expr)
}

estimateCov <- function(exp1, exp2, shrink) {
  ncl1 <- NROW(exp1)
  ncl2 <- NROW(exp2)

  cov1 <- estimateExprCov(exp1, shrink)
  cov2 <- estimateExprCov(exp2, shrink)

  s <- (cov1*(ncl1-1) + cov2*(ncl2-1)) / (ncl1 + ncl2 - 2)
  list(s1=cov1, s2=cov2, s=s)
}
