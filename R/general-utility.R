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

#' Create fake expression
#'
#' @param component_x the number of component in x
#' @param component_y the number of component in y
#' @param dimension the final dimension of the expression matrix
#'
#' @return dummy expression matrix
#'
#' @examples randomExpression()
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats rnorm
#' @export
randomExpression <- function(component_x=2, component_y=2, dimension=c(10,20)) {
  set.seed(1234)
  if (component_x > dimension[1] | component_x > 20)
    stop("Too much component_x")
  if (component_y > dimension[2] | component_y > 5)
    stop("Too much component_y")
  if (any(dimension<1))
    stop("both dimension must be positive")

  means <- sample(1:20, component_x)
  features <- seq_len(dimension[1])
  setsLength <- sapply(split(features, features%%component_x), length)
  means <- rep(means, times=setsLength)

  vars <- abs(rnorm(component_y*2, mean = 2, sd = 1))
  n = dimension[1]*dimension[2]
  noise <- matrix(sample(vars,n, replace=T), ncol=dimension[1])
  exprs <- t(mvtnorm::rmvnorm(dimension[2], mean=means)*noise)
}
