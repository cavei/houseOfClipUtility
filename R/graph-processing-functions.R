#' Extract the maximal cliques
#'
#' For internal use only. Extract the cliques.
#'
#'
#' @param dag a Directed Aciclic Graph
#' @param root a node to use as root
#'
#' @return list of nodes cliques
#'
#' @examples
#'   graph <- gRbase::dag(c("me","ve"),c("me","al"),c("me","me"),
#'     c("ve","al"),c("al","an"),
#'     c("al","st"),c("an","st"))
#'   extractCliquesFromDag(graph)
#'
#' @importFrom methods as
#' @importFrom gRbase is.DAG moralize triangulate rip
#' @importFrom checkmate assertClass
#' @rdname graph-processing
#' @export
#'
extractCliquesFromDag <- function(dag, root=NULL) {
  checkmate::assertClass(dag, "graphNEL")
  if (sum(diag(as(dag,"matrix")))!=0){
    dag <- removeSelfLoops(dag) # Done
  }

  if (gRbase::is.DAG(dag)) {
    moral <- gRbase::moralize(dag)
  } else {
    moral <- mmmoralize(dag)
  }

  tg    <- gRbase::triangulate(moral)
  ripped <- gRbase::rip(tg, root=root)
  if (length(ripped)==0)
    return(NULL)
  ripped$cliques
}

#' Remove self loops from a graphNEL
#'
#' Remove the self loops that a present in the graph graphNEL object
#'
#' @param graph a graphNEL object
#'
#' @return a graphNEL object
#'
#'#' @rdname  graph-processing
#' @examples
#'   graph <- gRbase::dag(c("me","ve"),c("me","al"),c("me","me"),
#'     c("ve","al"),c("al","an"),
#'     c("al","st"),c("an","st"))
#'   removeSelfLoops(graph)
#'
#' @importClassesFrom graph graphNEL
#' @importFrom checkmate assertClass
#' @export
#'
removeSelfLoops <- function(graph){
  checkmate::assertClass(graph, "graphNEL")
  edgeL <- graph@edgeL
  for (i in 1:length(edgeL)) {
    pos <- match(i,edgeL[[i]]$edges)
    if (!(is.na(pos)))
      edgeL[[i]]$edges <- edgeL[[i]]$edges[-pos]
  }
  graph@edgeL <- edgeL
  return(graph)
}

#' Moralize
#'
#' For internal use only. Force Moralization
#'
#' @inheritParams removeSelfLoops
#'
#' @importClassesFrom graph graphNEL
#' @importFrom gRbase graphNEL2M moralizeMAT coerceGraph
#' @rdname graph-processing
#'
mmmoralize <- function(graph) {
  m <- gRbase::graphNEL2M(graph)
  m <- gRbase::moralizeMAT(m)
  gRbase::coerceGraph(m, "graphNEL")
}

