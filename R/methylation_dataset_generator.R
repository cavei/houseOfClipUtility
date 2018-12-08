#' Simulate a methylation dataset with survival annotation
#'
#' To simulate a dataset we need a real pathway.
#'
#' @param pathway a Pathway
#' @param impactThisGenes the DEG in the aberrant signal chain
#' @param impactMeans list of the same length of impactThisGenes with lower and upper bounds for simulation
#' @param ann survival annotations
#' @param omicName the name of the omics to consider
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
#' @rdname simulate_methylation_dataset
#' @export
#'
makeTheMethylationDataset <- function(pathway, impactThisGenes, impactMeans, ann, omicName="x") {
  checkmate::assertClass(pathway, "Pathway")

  if (!(omicName %in% names(ann))) {
    stop(paste0(omicName, " does not apper in ann data.frame."))
  }

  if (length(impactThisGenes)!=length(impactMeans))
    stop("impactThisGenes length must be equal to impactMeans")
  if (!identical(impactThisGenes, names(impactMeans)))
    stop("impactMeans must be a list with names equal to impactThisGenes and value up_mean, low_mean")

  nd = graphite::nodes(pathway)
  missing_genes <- setdiff(impactThisGenes, nd)
  if (length(missing_genes) != 0) {
    stop(paste0("some of impactThisGenes are missing: ", paste(missing_genes, collapse = ', ')))
  }

  met <- create_random_uniform_beta(nd, nrow(ann), col_names=row.names(ann))

  for (i in impactThisGenes) {
    met[i,] <- create_bimodal_beta(ann[[omicName]], impactMeans[[i]][1], impactMeans[[i]][2])
  }

  if (!identical(colnames(met), row.names(ann)))
    stop("something wrong 340")

  annotations <- data.frame(status=ann$status, days=ann$stop,
                            class=ann[[omicName]], row.names=row.names(ann), stringsAsFactors=FALSE)

  return(list(exprs=met, annotation=annotations, graph=graphite::pathwayGraph(pathway), chain=NULL))
}

#' Simulate a flat methylation dataset with survival annotation
#'
#' To simulate a dataset we need a real pathway.
#'
#' @inheritParams makeTheMethylationDataset
#'
#' @importFrom checkmate assertClass
#' @importFrom graphite pathwayGraph
#' @importFrom graph nodes randomNodeGraph
#' @importFrom simPATHy simPATHy
#'
#' @rdname simulate_methylation_dataset
#' @export
#'
makeUniformMethylationDataset <- function(pathway, ann, omicName="x") {
  checkmate::assertClass(pathway, "Pathway")
  if (!(omicName %in% names(ann))) {
    stop(paste0(omicName, " does not apper in ann data.frame."))
  }

  nd = graphite::nodes(pathway)
  met <- create_random_uniform_beta(nd, nrow(ann), col_names=row.names(ann))

  if (!identical(colnames(met), row.names(ann)))
    stop("something wrong 340")

  annotations <- data.frame(status=ann$status, days=ann$stop,
                            class=ann[[omicName]], row.names=row.names(ann), stringsAsFactors=FALSE)

  return(list(exprs=met, annotation=annotations, graph=graphite::pathwayGraph(pathway), chain=NULL))
}

create_bimodal_beta <- function(classes, up_mean=0.8, low_mean=0.4, offset_up=0.1, offset_low=0.1) {
  up <- clip_value(c(up_mean - offset_up, up_mean + offset_up))
  low <- clip_value(c(low_mean - offset_low, low_mean + offset_low))

  up_values <- runif(sum(classes==1), up[1], up[2])
  low_values <- runif(sum(classes==0), low[1], low[2])
  profile <- rep(0, length(classes))
  profile[classes==1] <- up_values
  profile[classes==0] <- low_values
  profile
}

create_random_uniform_beta <- function(genes, n_patients, col_names=NULL) {
  genes_mean_beta <- runif(length(genes))
  offsets <- runif(length(genes), 0.1, 0.2)

  upper <- clip_value(genes_mean_beta + offsets, 0, 1)
  lower <- clip_value(genes_mean_beta - offsets, 0, 1)

  met_flat_profiles <- lapply(seq_along(genes), function(i) {
    runif(n_patients, lower[i], upper[i])
  })

  met <- do.call(rbind, met_flat_profiles)
  row.names(met) <- genes
  colnames(met) <- col_names
  met
}
