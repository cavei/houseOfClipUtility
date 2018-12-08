#' Plot graphNEL with some nodes boerde highlighted
#'
#' @param graph graphNEL
#' @param path sequence of genes in a path
#' @param singles vector of singlets genes
#' @param path_color border color of the genes in the path
#' @param singles_color border color of the singlets
#'
#' @importFrom graph nodes
#' @importFrom checkmate assert
#' @importFrom graphics plot
#' @importFrom stats rbinom
#'
#' @return NULL
#'
#' @export
plotPathInGraphNEL <- function(graph, path=NULL, singles=NULL, path_color="red", singles_color="green") {

  checkmate::assertClass(graph, "graphNEL")

  cols <- rep("grey", length(graph::nodes(graph)))
  names(cols) <- graph::nodes(graph)

  if (!is.null(singles))
    cols[names(cols) %in% singles] <- singles_color
  if (!is.null(path))
    cols[names(cols) %in% path] <- path_color

  nAttrs<-list()
  nAttrs$color <- cols
  plot(graph, nodeAttrs = nAttrs)
}
