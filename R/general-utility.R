convertIds <- function(from, keytype="ENTREZID", to="SYMBOL", annDB="org.Hs.eg.db") {
  requireNamespace(annDB, character.only = T)
  # entrez <- gsub("ENTREZID:", "", from)
  select(get(annDB), keys=from, columns = to, keytype=keytype)[[to]]
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
