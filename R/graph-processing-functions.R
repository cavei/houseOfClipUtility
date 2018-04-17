mergeGraphitePathways <- function(pathways, id="path1", title="path1", species="hsapiens",
                                  database="custom", timestamp=Sys.Date()) {
  invisible(lapply(pathways, checkmate::assertClass, classes="Pathway"))

  protEdges <- unique(do.call(rbind, lapply(pathways, function(p) p@protEdges)))
  row.names(protEdges)<-NULL
  metabolEdges<- unique(do.call(rbind, lapply(pathways, function(p) p@metabolEdges)))
  row.names(metabolEdges) <- NULL
  mixedEdges <- unique(do.call(rbind, lapply(pathways, function(p) p@mixedEdges)))
  row.names(mixedEdges) <- NULL
  graphite::buildPathway(id, title, species, database, protEdges, metabolEdges, mixedEdges, timestamp)
}


subPathway <- function(ids, pathway, and=TRUE) {
  invisible(checkmate::assertClass(pathway, classes="Pathway"))
  protEdges <- makeSubset(ids, pathway@protEdges, and)
  metabolEdges <- makeSubset(ids, pathway@metabolEdges, and)
  mixedEdges <- makeSubset(ids, pathway@mixedEdges, and)
  graphite::buildPathway(pathway@id, pathway@title, pathway@species,
                         pathway@database, protEdges, metabolEdges,
                         mixedEdges, Sys.Date())
}

makeSubset <- function(ids, edgeDF, and=TRUE) {
  src <- paste(edgeDF$src_type, edgeDF$src, sep=":")
  dest <- paste(edgeDF$dest_type, edgeDF$dest, sep=":")
  if (and) {
    select <- src %in% ids & dest %in% ids
  } else {
    select <- src %in% ids | dest %in% ids
  }
  edgeDF[select, , drop=F]
}
