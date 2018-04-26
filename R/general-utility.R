convertIds <- function(from, keytype="ENTREZID", to="SYMBOL", annDB="org.Hs.eg.db", graphiteStyle=FALSE) {
  requireNamespace(annDB, character.only = T)
  if(graphiteStyle) {
    suff <- paste0(keytype,":")
    from <- gsub(suff, "", from)
  }
  trans <- select(get(annDB), keys=from, columns = to, keytype=keytype)[[to]]
  if(graphiteStyle) {
    trans <- paste0(to,":", trans)
  }
  trans
}

reverseDict <- function(dict) {
  ids <- names(dict)
  df <- lapply(seq_along(ids), function(i){
    data.frame(src=rep(ids[i], length(dict[i])), dest=unname(dict[[i]]), stringsAsFactors = F)
  })
  df <- do.call(rbind, df)
  tapply(df$src, df$dest, function(x) x, simplify=F)
}

convertClusterIds <- function(clsId, dict) {
  cls <- intersect(names(dict), clsId)
  unlist(dict[clsId])
}
