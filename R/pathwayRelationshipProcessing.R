getPathFathers <- function(pathway, hierarchyGraph, ord=3, plot=F) {
  assertClass(hierarchyGraph, "igraph")

  if (!(pathway %in% names(V(hierarchyGraph)))){
    warning(paste0("Id ", pathway, " is not in the hierarchy."))
    return(pathway)
  }

  mm <- make_ego_graph(hierarchyGraph, ord, nodes = pathway, mode="in")
  mmlist <- ego(hierarchyGraph, ord, nodes = pathway, mode="in")

  if (plot)
    plot(mm[[1]])

  chain <- as_ids(mmlist[[1]])
  parents <- chain[-1]
  if (length(parents)==0)
    return(chain)

  dis <- distances(mm[[1]])
  idx <- which.max(dis[pathway, ])
  return(colnames(dis)[idx])
}

id2pathwayName <- function(path2fathers, pathwayDict) {
  lapply(path2fathers, function(x) {
    unlist(unname(pathwayDict[x]))
  })
}
