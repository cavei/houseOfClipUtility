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
#'
#' @export
filterDirected <- function(graph, reverse=FALSE) {
  checkmate::assertClass(graph, "Pathway")
  eds <- list(
    proteins=graphite::edges(graph, "proteins", stringsAsFactors=FALSE),
    matabolites=graphite::edges(graph, "metabolites", stringsAsFactors=FALSE),
    mixed=graphite::edges(graph, "mixed", stringsAsFactors=FALSE)
  )
  newEds <- lapply(eds, function(ed){
    pss <- ed$direction
    sel <- pss == "directed"
    if (reverse)
      sel <- !sel

    data.frame(src_type=ed$src_type[sel], src=ed$src[sel],
               dest_type=ed$dest_type[sel], dest=ed$dest[sel],
               direction=ed$direction[sel], type=ed$type[sel],
               stringsAsFactors = FALSE)
  })

  graphite::buildPathway(graph@id,graph@title,graph@species, graph@database,
               proteinEdges=newEds$proteins, metaboliteEdges = newEds$matabolites, mixedEdges = newEds$mixed)
}

#' Simulate Survival Data
#'
#' Create survival censored data annotations
#'
#' @param beta the amount of effect of the covariate
#' @param patientNum the number of patients to simulate
#' @param fUpTime the length of the followup
#' @param seed random seed setting
#'
#' @inheritParams simple.surv.sim return
#'
#' @importFrom survsim simple.surv.sim
#'
#' @export
simulateData <- function(beta, patientNum=300, fUpTime=1000, seed=1234) {
  ## argument beta: the effet of the covariate.
  dist.ev <- "weibull" # evento Weibull
  anc.ev <- 2.5
  beta0.ev <- 5
  dist.cens <- "weibull" # censored Weibull
  anc.cens <- 2.5
  beta0.cens <- 5.1
  z <-  list(c("unif", "0.98", "1")) # errore Uniform
  x <- list(c("bern", 0.6)) # bernulli

  if (!is.list(beta))
    beta = list(beta)

  ## set.seed(seed)
  survsim::simple.surv.sim(n=patientNum, foltime=fUpTime,
                  dist.ev=dist.ev, anc.ev=anc.ev, beta0.ev=beta0.ev,
                  dist.cen=dist.cens, anc.cens=anc.cens, beta0.cens=beta0.cens,
                  z=z, beta=beta, x=x)
}

#' Create a path from pathway
#'
#' @param graph Pathway graphite
#' @param start startig node
#' @param end end node
#' @param mode the direction to follow, out mean from source to dest
#'
#' @inheritParams all_shortest_paths return
#'
#' @examples
#'   if (require(graphite)) {
#'     p <- pathways("hsapiens", "kegg")[["Cell cycle"]]
#'     createPath(p, "ENTREZID:25", "ENTREZID:7029")
#'   }
#'
#' @importFrom igraph all_shortest_paths graph.data.frame
#' @importFrom checkmate assertClass
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

#' Simulate a dataset with survival annotation
#'
#' To simulate a dataset we need a real pathway.
#'
#' @param pathway a Pathway
#' @param pathBoundaries the start and end point of the disregulated chain
#' @param deInChain the DEG in the aberrant signal chain
#' @param deInPath the DEG in Pathway
#' @param w the index of the path to choose
#' @param ann survival annotations
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
#' @importFrom graph nodes
#' @importFrom simPATHy simPATHy
#'
#' @export
#'
makeTheDataset <- function(pathway, pathBoundaries, deInChain, deInPath, w, ann, augmentedMu=2, strong=FALSE){
  ## p <- c("EGFR tyrosine kinase inhibitor resistance")
  ## pathBoundaries <- list(start="3082",end="6774")
  ## deInChain de in the interesting chain
  ## deInPath de contained in the pathway
  ## w the index of the shortest path
  ## ann the annotation of the patients

  assertClass(pathway, "Pathway")
  grp   <- filterDirected(pathway)
  graph <- graphite::pathwayGraph(grp)

  path <- createPath(grp, pathBoundaries$start, pathBoundaries$end)

  allPaths <- unique(do.call(rbind,lapply(path,function(x) {names(x)})))
  selectP <- allPaths[w,]
  nd = graph::nodes(graph)

  argList <- computeSimPATHyParameters(selectP, graph::nodes(graph), deInChain=deInChain, deInPath=deInPath, augmentedMu, strong=strong)

  n1 = sum(ann$x==1)
  n2 = sum(ann$x==0)

  ex <- simPATHy(graph, path=argList$path,
                 min=argList$min, max=argList$max, prob = argList$prob,
                 mu1=argList$mu1, mu2=argList$mu2,
                 n1=n1, n2=n2)

  classes <- sapply(ann$x, function(x){
    if (x==1) {
      return("cl1")
    } else {
      return("cl2")
    }
  })

  annForData <- data.frame(pid=paste("P", ann$nid, sep="_"), status=ann$status, days=ann$stop, class=classes, stringsAsFactors=FALSE)

  data <- t(ex$dataset)
  cnames <- c(annForData[annForData$class=="cl1","pid"], annForData[annForData$class=="cl2","pid"])
  colnames(data) <- cnames

  row.names(annForData)<-annForData$pid

  annotations <- annForData[colnames(data),]
  return(list(exprs=data, annotation=annotations, graph=graph, chain=selectP))
}

#' Make a dataset with no association between chains and survival
#'
#' @param pathway a Pathway
#' @param ann survival annotations
#'
#' @inheritParams makeTheDataset return
#'
#' @importFrom graphite pathwayGraph
#' @importFrom checkmate assertClass
#' @importFrom simPATHy simPATHy
#' @export
makeUniformDataset <- function(pathway, ann){
  ## p <- c("EGFR tyrosine kinase inhibitor resistance")

  assertClass(pathway, "Pathway")
  grp   <- filterDirected(pathway)
  graph <- graphite::pathwayGraph(grp)

  n1 = sum(ann$x==1)
  n2 = sum(ann$x==0)
  mu1 <- sample(5:12, length(nodes), replace = T)

  ex <- simPATHy(graph, mu1=mu1, mu2=mu1,
                 n1=n1, n2=n2)

  classes <- sapply(ann$x, function(x){
    if (x==1) {
      return("cl1")
    } else {
      return("cl2")
    }
  })

  annForData <- data.frame(pid=paste("P", ann$nid, sep="_"), status=ann$status, days=ann$stop, class=classes, stringsAsFactors=FALSE)

  data <- t(ex$dataset)
  cnames <- c(annForData[annForData$class=="cl1","pid"], annForData[annForData$class=="cl2","pid"])
  colnames(data) <- cnames

  row.names(annForData)<-annForData$pid

  annotations <- annForData[colnames(data),]
  return(list(exprs=data, annotation=annotations, graph=graph, chain=NULL))
}

checkGeneInPaths <- function(gene, clipped) {
  if (NROW(clipped)==0)
    return(0)

  f <- sapply(seq_len(NROW(clipped)), function(i) {
    genes <- unlist(strsplit(clipped[i,11], ";"))
    if (gene %in% genes) {
      i
    } else {
      0
    }
  })
  if (all(f == 0))
    return(0)

  return(paste(f[f!=0], collapse=";"))
}

checkChainInBestPath <- function(chain, clipped) {
  if (NROW(clipped)==0)
    return(c("0","NULL","0","NULL"))

  m <- max(as.numeric(clipped[,"maxScore"]))
  bestPathIdx <- which(as.numeric(clipped[,"maxScore"])==m)

  tpExcluded <- Inf
  fp <- Inf

  for (i in bestPathIdx) {
    genes <- unlist(strsplit(clipped[bestPathIdx,11], ";"))
    local.tpEcluded <- length(setdiff(chain, genes))
    local.fp <- length(setdiff(genes, chain))

    if (local.tpEcluded < tpExcluded) {
      tpExcluded <- local.tpEcluded
      vpExcluded <- setdiff(chain, genes)
      vpExcludedC <- paste(setdiff(chain, genes), collapse=";")
      falseP     <- setdiff(genes, chain)
      falsePC    <- paste(setdiff(genes, chain), collapse=";")
    }
    if ((local.tpEcluded == tpExcluded) && (local.fp < fp)) {
      tpExcluded <- local.tpEcluded
      fp <- local.fp

      vpExcluded <- setdiff(chain, genes)
      vpExcludedC <- paste(setdiff(chain, genes), collapse=";")
      falseP     <- setdiff(genes, chain)
      falsePC    <- paste(setdiff(genes, chain), collapse=";")
    }
  }


  return(c(length(vpExcluded), vpExcludedC, length(falseP), falsePC))
}
