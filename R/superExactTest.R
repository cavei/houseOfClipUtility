#' run supertest
#'
#' For internal use only. Perform 'supertest' analysis.
#'
#' @param list a list, each element a list of pvalues
#' @param universeSize the size of all the pvalues computed
#' @param thr threshold to consider significant
#'
#' @return a 'msets' object
#'
#' @importFrom SuperExactTest supertest
#' @importFrom graphics plot
#' @export
#'
runSuperTest <- function(list, universeSize, thr=0.05) {
  echo <- "given a list the function produce the intersection amongst all the sets."
  tmp <- lapply(list, function(x) {
    namesAreNULL <- !is.null(names(x))
    stopifnot(namesAreNULL)
    names(which(x <= thr))
  })
  SuperExactTest::supertest(tmp, n=universeSize)
}

#' Plot supertest
#'
#' For internal use only. Plot 'supertest' analysis.
#'
#' @param mset a 'msets' object
#'
#' @return NULL
#'
#' @importFrom checkmate assertClass
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom SuperExactTest plot.msets
#' @export
#'
plotSuperTest <- function(mset) {
  assertClass(mset, "msets")

  plot(mset, color.on="#409ec3",color.off="white",
       heatmapColor=grDevices::colorRampPalette(RColorBrewer::brewer.pal(9,"OrRd"))(100))
}

#' Create a significance matrix
#'
#' For internal use only. Create a data.frame of TRUE/FALSE
#'
#' @inheritParams runSuperTest list thr
#'
#' @return a 'data.frame'. Each line a feature, columns are covariates: TRUE/FALSE significant
#' @examples
#' l <- list(a=c(0.03, 0.06, 1), b=c(0.6, 0.01, 0.5))
#' createSignificanceMask(l)
#'
#' @rdname runSuperTest
#' @export
#'
createSignificanceMask <- function(list, thr=0.05) {
  len <- length(list[[1]])
  if (!all(sapply(list, length) == len))
    stop("list must have elements equal in length")
  tmp <- lapply(list, function(x) {
    ifelse(x <= thr, 1, 0)
  })
  data.frame(do.call(cbind, tmp), stringsAsFactors = F)
}

#' Compute supertest frequences
#'
#' For internal use only. Compute supertest frequences
#'
#' @param elementsIntersections a summary(msets) object
#' @param sep the field separator of the class
#'
#' @return a 'data.frame' with the frequences
#'
#' @export
#'
computeSupertestFrequencies <- function(elementsIntersections, sep=".") {
  frequencies <- c()
  fathers <- c()
  intersectionClass<-c()
  splitClasses <- c()

  for(i in seq_along(elementsIntersections)){
    class<-elementsIntersections[[i]]
    subClasses <- unlist(strsplit(class,", "))
    spClasses <-unlist(lapply(subClasses, function(subClass) strsplit(subClass,sep,fixed = T)[[1]][1]))
    spClasses <- gsub("^ ", "",spClasses)
    frequencies<-c(frequencies,table(spClasses))
    fathers<-c(fathers,names(table(spClasses)))
    intersectionClass<-c(intersectionClass,rep(names(elementsIntersections[i]),length(table(spClasses))))
  }
  data.frame(category=fathers, frequencies=frequencies, class=intersectionClass, stringsAsFactors = F)
}

#' Compute frequences in list
#'
#' For internal use only. Compute  frequences
#'
#' @inheritParams computeSupertestFrequencies
#'
#' @rdname computeSupertestFrequencies
#'
#' @export
#'
computeFreqs <- function(elementsIntersections) {
  freques <- lapply(names(elementsIntersections), function(x) {
    counts <- table(elementsIntersections[[x]])
    data.frame(category=names(counts), frequencies=as.numeric(counts), class=x, stringsAsFactors = F)
  })
  do.call(rbind, freques)
}

#' Plot supertest frequences
#'
#' For internal use only. Plot supertest frequences
#'
#' @param supertestFreq a data.frame created from 'computeSupertestFrequencies'
#' @param manualColors optional vector of colors
#' @param minSize the minimal fontsize. Maximal frequencies will be added for each class
#' @param maxSize the maximal fontsize dimension, all values above are clipped
#' @param width the numerbe of character to wrap the labels
#'
#' @return NULL
#' @examples
#' df <- data.frame(category=c("talk", "too","mutch", "dear"),
#'   frequencies=c(1,2,1,3),
#'   class=rep("Mut",4), stringsAsFactors = FALSE)
#' plotSupertestFrequencies(df)
#'
#' @importFrom ggplot2 ggplot aes element_text coord_polar geom_point geom_polygon geom_path theme theme_bw labs element_blank
#' @importFrom stats reorder
#' @export
#'
plotSupertestFrequencies <- function(supertestFreq, manualColors=NULL, minSize=4, maxSize=20, width=20){
  if (!all(colnames(supertestFreq) %in% c("category", "frequencies", "class")))
    stop("supertestFreq dataframe must contain columns category, frequencies and class")
  if (!is.null(manualColors)) {
    g <- ggplot2::ggplot(supertestFreq, aes(y = frequencies,
                                            # x = suppressWarnings(stats::reorder(category, class)),
                                            x = factor(category),
                                            group = class, colour=class))
    g <- g + ggplot2::scale_colour_manual(values=manualColors)
  } else {
    g <- ggplot2::ggplot(supertestFreq, aes(y = frequencies,
                                       # x = suppressWarnings(stats::reorder(category, class)),
                                       x = factor(category),
                                       group = class, colour = class))
  }
  size <- tapply(seq_along(supertestFreq$frequencies), factor(supertestFreq$category), function(idx) max(supertestFreq$frequencies[idx]))
  size <- as.numeric(size)+minSize
  size[size > maxSize] <- maxSize
  g + ggplot2::coord_polar() +
    ggplot2::geom_point(stat='identity') +
    ggplot2::geom_polygon(fill=NA)+
    ggplot2::geom_path() +
    ggplot2::labs(x = NULL) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.border = element_blank(), axis.line.x = element_blank(), axis.line.y = element_blank()) +
    ggplot2::theme(axis.text.x=element_text(size=size)) +
    ggplot2::scale_x_discrete(labels=function(x) lapply(strwrap(x, width = width, simplify = FALSE), paste, collapse="\n"))
}

#' Minimum or NA
#'
#' For internal use only. Get back minimum or NA.
#'
#' @param x a numeric
#'
#' @return a numeric. The minimum or NA
#'
#' @examples
#'   minOrNA(c(1,5,0.1,NA))
#'   minOrNA(c(NA,NA,NA))
#'
#' @export
#'
minOrNA <- function(x) {
  if (all(is.na(x)))
    return(NA)
  min(x, na.rm=T)
}

#' Summarize Omics Covaraites By Min Pvalue
#'
#' For internal use only. for each line extrac 'col' and get the minimum.
#'
#' @param col columns to extract from the line
#' @param mat the matrix to be summarized (were to extract lines and 'col')
#'
#' @return a summarized version of the matrix.
#'
#' @examples
#'   summarizeOmicsResByMinPvalue(2:3, mat=matrix(c(1,2,4,1,2,5), nrow=2))
#'
#' @export
#'
summarizeOmicsResByMinPvalue <- function(col, mat) {
  apply(as.matrix(mat[, col, drop=F]), 1, minOrNA)
}
