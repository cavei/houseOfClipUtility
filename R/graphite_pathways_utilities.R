#' Filter direct edges
#'
#' Filter all directed edges
#'
#' @param graph a Pathway graph from graphite
#' @param reverse invert selection
#'
#' @return a Pathway with selected edges
#'
#' @importFrom checkmate assertClass
#' @importFrom graphite edges buildPathway
#' @rdname filterDirected
#'
#' @export
filterDirected <- function(graph, reverse=FALSE) {
  checkmate::assertClass(graph, "Pathway")
  eds <- list(
    proteins=graph@protEdges,
    matabolites=graph@metabolEdges,
    mixed=graph@mixedEdges,
    protPropEdges = graph@protPropEdges,
    metabolPropEdges = graph@metabolPropEdges
  )

  newEds <- lapply(eds, function(ed){
    pss <- ed$direction
    sel <- pss == "directed"
    if (reverse)
      sel <- !sel
    ed[sel, ,drop=F]
  })

  customPath <- graphite::buildPathway(graph@id,graph@title,graph@species, graph@database,
                                       proteinEdges=newEds$proteins, metaboliteEdges = newEds$matabolites, mixedEdges = newEds$mixed)
  customPath@metabolPropEdges <- newEds$metabolPropEdges
  customPath@protPropEdges <- newEds$protPropEdges
  customPath
}

#' Create a path from pathway
#'
#' @param start startig node
#' @param end end node
#' @param mode the direction to follow, out mean from source to dest
#'
#' @inheritParams filterDirected
#'
#' @examples
#'   if (require(graphite)) {
#'     p <- pathways("hsapiens", "kegg")[["Cell cycle"]]
#'     createPath(p, "ENTREZID:25", "ENTREZID:7029")
#'   }
#'
#' @importFrom igraph all_shortest_paths graph.data.frame
#' @importFrom checkmate assertClass
#'
#' @rdname filterDirected
#' @export
createPath <- function(graph, start, end, mode="out") {
  checkmate::assertClass(graph, "Pathway")
  eds <- rbind(graph@protEdges, graph@metabolEdges, graph@mixedEdges)
  eds <- data.frame(src=paste(eds$src_type,eds$src, sep=":"),
                    dest=paste(eds$dest_type,eds$dest, sep=":"), stringsAsFactors = F)
  g <- graph.data.frame(eds, directed=T, vertices=graphite::nodes(graph, "mixed"))
  igraph::all_shortest_paths(g, start, end, mode=mode)$res
}

#' Convert a chain in an edgeList
#'
#' @param chain vector of node names
#'
#' @return edgeList
#'
#' @rdname filterDirected
#' @export
fromChain2edgeList <- function(chain) {
  if (length(chain) < 2)
    stop("Chain must be at least of length 2")
  if (length(chain) == 2)
    return(list(chain))
  lapply(seq(2,length(chain)), function(i) {
    c(chain[i-1], chain[i])
  })
}

