#' Merge graphite Pathways
#'
#' For internal use only. Merge a list of Pathways from Graphite
#'
#' @param pathways a list of pathways.
#' @param id string for newly generated pathway id
#' @param title string for newly generated pathway title
#' @param species the species associated to the pathway
#' @param database the name of the databse
#' @param timestamp a date
#'
#' @return a \code{Pathway} object
#'
#' @export
#'
mergeGraphitePathways <- function(pathways, id="path1", title="path1", species="hsapiens",
                                  database="custom", timestamp=Sys.Date()) {
  invisible(lapply(pathways, checkmate::assertClass, classes="Pathway"))

  protEdges <- unique(do.call(rbind, lapply(pathways, function(p) p@protEdges)))
  row.names(protEdges)<-NULL
  protPropEdges <- unique(do.call(rbind, lapply(pathways, function(p) p@protPropEdges)))
  row.names(protPropEdges)<-NULL
  metabolEdges<- unique(do.call(rbind, lapply(pathways, function(p) p@metabolEdges)))
  row.names(metabolEdges) <- NULL
  metabolPropEdges<- unique(do.call(rbind, lapply(pathways, function(p) p@metabolPropEdges)))
  row.names(metabolPropEdges) <- NULL
  mixedEdges <- unique(do.call(rbind, lapply(pathways, function(p) p@mixedEdges)))
  row.names(mixedEdges) <- NULL
  customPath <- graphite::buildPathway(id, title, species, database, protEdges, metabolEdges, mixedEdges, timestamp)
  customPath@metabolPropEdges <- metabolPropEdges
  customPath@protPropEdges <- protPropEdges
  customPath
}

#' Subset of a graphite Pathway
#'
#' For internal use only. Extract sub network from a Graphite Pathway.
#'
#' @param ids nodes ids
#' @param pathway a graphite Pathway object
#' @param and logical keep the edges between ids
#'
#' @return a \code{Pathway} object
#'
#' @export
#'
subPathway <- function(ids, pathway, and=TRUE) {
  invisible(checkmate::assertClass(pathway, classes="Pathway"))
  protEdges <- makeSubset(ids, pathway@protEdges, and)
  metabolEdges <- makeSubset(ids, pathway@metabolEdges, and)
  mixedEdges <- makeSubset(ids, pathway@mixedEdges, and)
  protPropEdges <- makeSubset(ids, pathway@protPropEdges, and)
  metabolPropEdges <- makeSubset(ids, pathway@metabolPropEdges, and)
  customPath <- graphite::buildPathway(pathway@id, pathway@title, pathway@species,
                         pathway@database, protEdges, metabolEdges,
                         mixedEdges, Sys.Date())
  customPath@metabolPropEdges <- metabolPropEdges
  customPath@protPropEdges <- protPropEdges
  customPath
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
