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
#' @param dict a list with keys to values
#'
#' @return a list with value to keys
#'
#' @export
#'
mapFromDict <- function(keys, dict) {
  cls <- intersect(names(dict), keys)
  unlist(dict[keys])
}
