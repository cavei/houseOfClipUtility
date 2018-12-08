#' Simulate Survival Data
#'
#' Create survival censored data annotations
#'
#' @param beta the beta values for covariates (see siple.surv.sim)
#' @param x the x values for covariates (see siple.surv.sim)
#' @param patientNum the number of patients to simulate
#' @param fUpTime the length of the followup
#' @param seed random seed setting
#' @param row.names.pref a prefix to row number to simulate real samples names
#'
#' @importFrom survsim simple.surv.sim
#'
#' @export
simulateData <- function(beta, x, patientNum=300, fUpTime=1000, seed=NULL, row.names.pref="P_") {
  ## argument beta: the effet of the covariate.
  dist.ev <- "weibull" # evento Weibull
  anc.ev <- 2.5
  beta0.ev <- 3
  dist.cens <- "weibull" # censored Weibull
  anc.cens <- 2.5
  beta0.cens <- 4
  z <-  list(c("unif", "0.98", "1")) # errore Uniform
  # x <- list(c("bern", 0.6)) # bernulli

  if (!identical(names(beta), names(x)))
    stop("names beta and names x must be equal")

  if (!is.list(beta))
    beta = list(beta)

  if (length(beta) != length(x))
    stop("invalid beta and x specification. Beta and x must be equal in length")

  if (!is.null(seed))
    set.seed(seed)

  survAnnotation <- survsim::simple.surv.sim(n=patientNum, foltime=fUpTime,
                                             dist.ev=dist.ev, anc.ev=anc.ev, beta0.ev=beta0.ev,
                                             dist.cen=dist.cens, anc.cens=anc.cens, beta0.cens=beta0.cens,
                                             z=z, beta=beta, x=x)
  if (!is.null(names(x))){
    cov_col_idx = setdiff(seq_len(ncol(survAnnotation)), 1:5)
    names(survAnnotation)[cov_col_idx] <- names(x)
  }
  row.names(survAnnotation) <- paste0(row.names.pref, row.names(survAnnotation))
  survAnnotation
}


checkGeneInPaths <- function(gene, clipped) {
  if (NROW(clipped)==0)
    return(0)

  f <- sapply(seq_len(NROW(clipped)), function(i) {
    genes <- unlist(strsplit(clipped[i,11], ";"))
    if (gene %in% genes) {
      i
    } else {
      0
    }
  })
  if (all(f == 0))
    return(0)

  return(paste(f[f!=0], collapse=";"))
}

checkChainInBestPath <- function(chain, clipped) {
  if (NROW(clipped)==0)
    return(c("0","NULL","0","NULL"))

  m <- max(as.numeric(clipped[,"maxScore"]))
  bestPathIdx <- which(as.numeric(clipped[,"maxScore"])==m)

  tpExcluded <- Inf
  fp <- Inf

  for (i in bestPathIdx) {
    genes <- unlist(strsplit(clipped[bestPathIdx,11], ";"))
    local.tpEcluded <- length(setdiff(chain, genes))
    local.fp <- length(setdiff(genes, chain))

    if (local.tpEcluded < tpExcluded) {
      tpExcluded <- local.tpEcluded
      vpExcluded <- setdiff(chain, genes)
      vpExcludedC <- paste(setdiff(chain, genes), collapse=";")
      falseP     <- setdiff(genes, chain)
      falsePC    <- paste(setdiff(genes, chain), collapse=";")
    }
    if ((local.tpEcluded == tpExcluded) && (local.fp < fp)) {
      tpExcluded <- local.tpEcluded
      fp <- local.fp

      vpExcluded <- setdiff(chain, genes)
      vpExcludedC <- paste(setdiff(chain, genes), collapse=";")
      falseP     <- setdiff(genes, chain)
      falsePC    <- paste(setdiff(genes, chain), collapse=";")
    }
  }


  return(c(length(vpExcluded), vpExcludedC, length(falseP), falsePC))
}

