#' Ids converter
#'
#' For internal use only. Convert ids from to ids to.
#'
#' @param from ids to be converted
#' @param keytype the type of the starting id
#' @param to the arriving type
#' @param annDB \code{annotationDbi} package
#' @param graphiteStyle logical whether to expect graphite style ids
#'
#' @return vector of translated ids
#' @importFrom AnnotationDbi select
#'
#' @export
#'
convertIds <- function(from, keytype="ENTREZID", to="SYMBOL", annDB="org.Hs.eg.db", graphiteStyle=FALSE) {
  requireNamespace(annDB, character.only = T)
  if(graphiteStyle) {
    suff <- paste0(keytype,":")
    from <- gsub(suff, "", from)
  }
  trans <- AnnotationDbi::select(get(annDB), keys=from, columns = to, keytype=keytype)[[to]]
  if(graphiteStyle) {
    trans <- paste0(to,":", trans)
  }
  trans
}

#' Invert The Dictionary
#'
#' For internal use only. From key to value, to value to key.
#'
#' @param dict a list with keys to values
#'
#' @return a list with value to keys
#'
#' @export
#'
reverseDict <- function(dict) {
  ids <- names(dict)
  df <- lapply(seq_along(ids), function(i){
    data.frame(src=rep(ids[i], length(dict[i])), dest=unname(dict[[i]]), stringsAsFactors = F)
  })
  df <- do.call(rbind, df)
  tapply(df$src, df$dest, function(x) x, simplify=F)
}

#' Map ids using a dict.
#'
#' For internal use only.
#'
#' @param keys a vector of cluster Ids
#' @param dict a vector with keys (names) to values
#'
#' @return a list with value to keys
#'
#' @examples mapFromDict("a", c(a = "5"))
#' @export
#'
mapFromDict <- function(keys, dict) {
  cls <- intersect(names(dict), keys)
  unlist(dict[keys])
}

#' Create fake expression based on variance
#'
#' @param n the number of element (a.k.a. genes)
#' @param nfeatures the number of features (a.k.a. samples)
#'
#' @return dummy expression matrix
#'
#' @examples randomExpression(4)
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats rnorm
#' @export
randomExpression <- function(n, nfeatures=10) {
  set.seed(1234)
  sigma <- matrix(rnorm(n^2), ncol=n)
  sigma[upper.tri(sigma)] <- t(sigma)[upper.tri(sigma)]
  diag(sigma) <- rev(seq_len(n))
  sigma <- makePositiveDefinite(sigma)$m1
  mvtnorm::rmvnorm(n = nfeatures, mean=seq_len(n), sigma = sigma)
}

#' Forse random nodes degress to be even
#'
#' @param degrees the nodes degrees, a named vector
#'
#' @return even nodes degrees
#'
#' @export
#'
makeDegreesEven <- function(degrees) {
  if (sum(degrees)%%2 == 0)
    return(degrees)
  idx <- sample(seq_len(degrees),1)
  degrees[idx] + 1
  degrees
}
