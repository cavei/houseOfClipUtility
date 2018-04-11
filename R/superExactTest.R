runSupertest <- function(list, universeSize, thr=0.05) {
  echo <- "given a list the function produce the intersection amonst all the sets."
  tmp <- lapply(list, function(x) {
    namesAreNULL <- !is.null(names(x))
    stopifnot(namesAreNULL)
    names(which(x <= thr))
  })
  supertest(tmp, n=universeSize)
}

plotSupertest <- function(mset) {
  assertClass(mset, "msets")
  plot(mset, color.on="#409ec3",color.off="white",
       heatmapColor=colorRampPalette(brewer.pal(9,"OrRd"))(100))
}

createSignificanceMask <- function(list, thr=0.05) {
  len <- length(list[[1]])
  if (!all(sapply(list, length) == len))
    stop("list must have elements equal in length")
  tmp <- lapply(list, function(x) {
    ifelse(x <= thr, 1, 0)
  })
  data.frame(do.call(cbind, tmp), stringsAsFactors = F)
}

computeSupertestFrequencies <- function(elementsIntersections) {
  frequencies <- c()
  fathers <- c()
  intersectionClass<-c()
  splitClasses <- c()

  for(i in seq_along(elementsIntersections)){
    class<-elementsIntersections[[i]]
    subClasses <- unlist(strsplit(class,", "))
    spClasses <-unlist(lapply(subClasses, function(subClass) strsplit(subClass,"-")[[1]][1]))
    spClasses <- gsub("^ ", "",spClasses)
    frequencies<-c(frequencies,table(spClasses))
    fathers<-c(fathers,names(table(spClasses)))
    intersectionClass<-c(intersectionClass,rep(names(elementsIntersections[i]),length(table(spClasses))))
  }
  data.frame(category=fathers, frequencies=frequencies, class=intersectionClass, stringsAsFactors = F)
}

plotSupertestFrequencies <- function(supertestFreq){
  if (!all(colnames(supertestFreq) %in% c("category", "frequencies", "class")))
    stop("supertestFreq dataframe must contain columns category, frequencies and class")
  ggplot(supertestFreq, aes(y = frequencies, x = suppressWarnings(reorder(category, class)), group = class, colour = class)) +
  coord_polar() +
  geom_point(stat='identity') +
  geom_polygon(fill=NA)+
  geom_path() +
  theme(axis.text.x=element_text(size=3, angle=45))+
  labs(x = NULL)+
  theme_bw()
}

minOrNA <- function(x) {
  if (all(is.na(x)))
    return(NA)
  min(x, na.rm=T)
}

summarizeOmicsResByMinPvalue <- function(col, mat) {
  apply(as.matrix(mat[, col, drop=F]), 1, minOrNA)
}
