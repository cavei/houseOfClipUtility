#' Simulate an expression dataset from survival annotation
#'
#' To simulate a dataset we need a real pathway.
#'
#' @param pathway a Pathway
#' @param pathBoundaries the start and end point of the disregulated chain
#' @param deInChain the DEG in the aberrant signal chain
#' @param deInPath the DEG in Pathway
#' @param w the index of the path to choose
#' @param ann survival annotations
#' @param omicName the name of the omics to consider
#' @param augmentedMu magnetude of mean increment
#' @param strong should the increment be really marked
#'
#' @return list
#'   \item{exprs}{expression}
#'   \item{annotation}{the survival annotation}
#'   \item{graph}{a graphNEL}
#'   \item{chain}{the chain of selectP}
#'
#' @importFrom checkmate assertClass
#' @importFrom graphite pathwayGraph
#' @importFrom graph nodes randomNodeGraph
#' @importFrom simPATHy simPATHy
#'
#' @rdname simulate_expression_dataset
#' @export
#'
makeTheExpressionDataset <- function(pathway, pathBoundaries, deInChain, deInPath, w, ann, omicName="x",
                           augmentedMu=2, strong=FALSE) {

  assertClass(pathway, "Pathway")

  if (!(omicName %in% names(ann))) {
    stop(paste0(omicName, " does not apper in ann data.frame."))
  }

  # grp   <- filterDirected(pathway)
  grp <- pathway
  graph <- graphite::pathwayGraph(grp)

  path <- createPath(grp, pathBoundaries$start, pathBoundaries$end)

  allPaths <- unique(do.call(rbind,lapply(path,function(x) {names(x)})))
  selectP <- allPaths[w,]
  nd = graph::nodes(graph)

  argList <- computeSimPATHyParameters(selectP, nd, deInChain=deInChain,
                                       deInPath=deInPath, augmentedMu, strong=strong)

  n1 = sum(ann[[omicName]]==1)
  n2 = sum(ann[[omicName]]==0)
  ex <- simPATHy(graph, path=argList$path,
                 min=argList$min, max=argList$max, prob = argList$prob,
                 mu1=argList$mu1, mu2=argList$mu2,
                 n1=n1, n2=n2)

  sel_cl1 = ann[[omicName]]==1
  patients_cl1 <- row.names(ann)[sel_cl1]
  patients_cl2 <- row.names(ann)[!sel_cl1]
  data <- t(ex$dataset)
  colnames(data) <- c(patients_cl1, patients_cl2)
  data <- data[, row.names(ann), drop=F]

  annotations <- data.frame(status=ann$status, days=ann$stop,
             class=ann[[omicName]], row.names=row.names(ann), stringsAsFactors=FALSE)

  return(list(exprs=data, annotation=annotations, graph=graph, chain=selectP))
}

#' Make a dataset with no association between chains and survival
#'
#' @inheritParams makeTheDataset
#'
#' @importFrom graphite pathwayGraph
#' @importFrom checkmate assertClass
#' @importFrom simPATHy simPATHy
#'
#' @rdname simulate_expression_dataset
#' @export
makeUniformExpressionDataset <- function(pathway, ann, omicName="x") {
  assertClass(pathway, "Pathway")

  if (!(omicName %in% names(ann))) {
    stop(paste0(omicName, " does not apper in ann data.frame."))
  }

  grp <- pathway
  graph <- graphite::pathwayGraph(grp)
  nd <- graphite::nodes(grp)

  n1 = sum(ann[[omicName]]==1)
  n2 = sum(ann[[omicName]]==0)
  mu1 <- sample(5:12, length(nd), replace = T)
  names(mu1) <- nd

  # ex <- simPATHy(graph, mu1=mu1, mu2=mu1,n1=n1, n2=n2)
  ex <- tryCatch(simPATHy(graph, mu1=mu1, mu2=mu1,n1=n1, n2=n2),
                 error=function(e) {warning("simPATHy failure");
                   simPATHy(graph, mu1=mu1, mu2=mu1,n1=n1, n2=n2)})

  sel_cl1 = ann[[omicName]]==1
  patients_cl1 <- row.names(ann)[sel_cl1]
  patients_cl2 <- row.names(ann)[!sel_cl1]
  data <- t(ex$dataset)
  colnames(data) <- c(patients_cl1, patients_cl2)
  data <- data[, row.names(ann), drop=F]

  annotations <- data.frame(status=ann$status, days=ann$stop,
                            class=ann[[omicName]], row.names=row.names(ann), stringsAsFactors=FALSE)

  return(list(exprs=data, annotation=annotations, graph=graph, chain=NULL))
}

#' Compute SimPATHy Parameters
#'
#' @param chain a chain o disregulation
#' @param nodes all the other nodes
#' @param deInChain the differentially expressed of the chain
#' @param deInPath the differentially expressed in Pathway
#' @param augmentedMu the magnetude of the mean inclement increment for DEGs
#' @param strong should be the increment really marked
#'
#' @return a list with the parameters
#'   \item{path}{the path}
#'   \item{min}{min}
#'   \item{max}{max}
#'   \item{prob}{prob}
#'   \item{mu1}{mu1}
#'   \item{mu2}{mu2}
#'
#' @importFrom stats runif
#' @export
computeSimPATHyParameters <- function(chain, nodes, deInChain=NULL, deInPath=NULL, augmentedMu=2, strong=FALSE) {

  p <- fromChain2edgeList(chain)
  pLen <- length(p)

  min <- sample(1:2, pLen, replace=T)
  max <- min+runif(pLen,0,1)
  prb <- rep(1,pLen)

  if (strong)
    if (all(deInChain %in% chain)) {
      idx <- which(deInChain %in% chain)
      prb[idx] <- 0
    }

  mu1 <- sample(5:12, length(nodes), replace = T)
  names(mu1) <- nodes
  mu2 <- mu1

  if (!is.null(deInChain))
    if (all(deInChain %in% nodes)) {
      mu2[deInChain] <- mu2[deInChain] + augmentedMu
    } else {
      stop("some deInChain are not nodes")
    }

  if (!is.null(deInPath))
    if (all(deInPath %in% nodes)) {
      mu2[deInPath] <- mu2[deInPath] + augmentedMu*70/100
    } else {
      stop("some deInPath are not nodes")
    }

  return(list(path=p, min=min, max=max, prob=prb,mu1=mu1, mu2=mu2))
}
