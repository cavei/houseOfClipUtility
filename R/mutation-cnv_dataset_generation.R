#' Simulate a mutation dataset with survival annotation
#'
#' To simulate a dataset we need a real pathway.
#'
#' @param pathway a Pathway
#' @param impacted_genes the genes of a module with the aberrant signal chain
#' @param patients_fractions the fraction of patients with n mutated genes
#' @param ann survival annotations
#' @param omicName the name of the omics to consider
#' @param mutation_rate the basel mutation rate
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
#' @rdname simulate_mutation_dataset
#' @export
#'
makeTheMutationDataset <- function(pathway, impacted_genes, patients_fractions, ann, omicName="x", mutation_rate=0.001) {
  checkmate::assertClass(pathway, "Pathway")

  if (!(omicName %in% names(ann))) {
    stop(paste0(omicName, " does not apper in ann data.frame."))
  }

  nd = graphite::nodes(pathway)
  missing_genes <- setdiff(impacted_genes, nd)
  if (length(missing_genes) != 0) {
    stop(paste0("some of impactThisGenes are missing: ", paste(missing_genes, collapse = ', ')))
  }
  n_patients = nrow(ann)

  mut <- create_random_uniform_mutatios(nd, n_patients, mutation_rate=mutation_rate, col_names=row.names(ann))

  mut_binomial <- create_bimodal_mutations(ann[[omicName]], impacted_genes,
                                           patients_fractions, mutation_rate=mutation_rate)

  mut[row.names(mut_binomial), ] <- mut_binomial

  if (!identical(colnames(mut), row.names(ann)))
    stop("something wrong 340")

  annotations <- data.frame(status=ann$status, days=ann$stop,
                            class=ann[[omicName]], row.names=row.names(ann), stringsAsFactors=FALSE)

  return(list(exprs=mut, annotation=annotations, graph=graphite::pathwayGraph(pathway), chain=NULL))
}

#' Simulate a flat mutation dataset with survival annotation
#'
#' To simulate a dataset we need a real pathway.
#'
#' @inheritParams makeTheMutationDataset
#'
#' @importFrom checkmate assertClass
#' @importFrom graphite pathwayGraph
#' @importFrom graph nodes randomNodeGraph
#' @importFrom simPATHy simPATHy
#'
#' @rdname simulate_mutation_dataset
#' @export
#'
makeUniformMutationsDataset <- function(pathway, ann, omicName="x", mutation_rate=0.001) {
  checkmate::assertClass(pathway, "Pathway")

  if (!(omicName %in% names(ann))) {
    stop(paste0(omicName, " does not apper in ann data.frame."))
  }

  nd = graphite::nodes(pathway)
  n_patients = nrow(ann)
  mut <- create_random_uniform_mutatios(nd, n_patients, mutation_rate=mutation_rate, col_names=row.names(ann))

  if (!identical(colnames(mut), row.names(ann)))
    stop("something wrong 340")

  annotations <- data.frame(status=ann$status, days=ann$stop,
                            class=ann[[omicName]], row.names=row.names(ann), stringsAsFactors=FALSE)

  return(list(exprs=mut, annotation=annotations, graph=graphite::pathwayGraph(pathway), chain=NULL))
}

create_bimodal_mutations <- function(classes, genes=letters[1:5], patients_fractions = c(0.05,0.1,0.1,0.2),
                                     mutation_rate=0.0005) {

  if (!(length(genes)-1 == length(patients_fractions))) {
    stop("something wrong 322")
  }

  cl1 <- sum(classes==1)
  numeric_order <- order(classes, decreasing = T)
  patients_fake_names = paste(seq_along(classes), sep="_", classes)
  ordered_patinets_fake_names <- patients_fake_names[numeric_order]

  patientsPerCategory <- splitPatients(cl1, patients_fractions)
  maxRunningSum <- length(genes)
  maxSum=rep(rev(seq_len(maxRunningSum)), times=patientsPerCategory)

  mutated_data <- sapply(maxSum, function(x) {
    patient_mutation_profile <- rep(0, length(genes))
    patient_mutation_profile[seq_len(x)] <- 1
    patient_mutation_profile <- sample(patient_mutation_profile)
  })

  non_mutated_data <- create_random_uniform_mutatios(genes, n_patients = sum(classes==0))

  data <- cbind(mutated_data, non_mutated_data)
  colnames(data) <- ordered_patinets_fake_names
  data <- data[, patients_fake_names]
  data
}

splitPatients <- function(n, fractions = c(0.05,0.1,0.1,0.2)) {
  blocks <- ceiling(n*fractions)
  last_block <- n - sum(blocks)
  c(blocks, last_block)
}

create_random_uniform_mutatios <- function(genes, n_patients, mutation_rate=0.0005, col_names=NULL) {
  mutation_random_profiles <- lapply(seq_along(genes), function(i) {
    rbinom(n_patients, 2, mutation_rate)
  })
  mut <- do.call(rbind, mutation_random_profiles)
  row.names(mut) <- genes
  colnames(mut) <- col_names
  mut
}

# sample(size = 100,x = seq(-2,2,1),prob = c(0.05,0.15,0.5,0.20,0.1), replace = T)
