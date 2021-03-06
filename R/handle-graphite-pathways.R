#' Map Pathways ID from Graphite
#'
#' For internal use only. Retrieve pathway id and names from Pathways object.
#'
#' @param pathways a PathwayList object
#' @param pathwayNames in not NULL, a subset of pathway to extract
#'
#' @return a data frame, id and pathway name
#'
#' @examples
#'   if (require(graphite)){
#'     pathways <- pathways("hsapiens", "kegg")
#'     mapPathwaysIDfromGraphite(pathways, names(pathways)[1:10])
#'   }
#' @importFrom checkmate assertClass
#'
#' @export
#'
mapPathwaysIDfromGraphite <- function(pathways, pathwayNames=NULL) {
  assertClass(pathways, "PathwayList")
  if (is.null(pathwayNames)) {
    pathwayNames <- names(pathways)
  }
  l <- lapply(pathwayNames, function(p) {
    if (!(p %in% names(pathways))) {
      warning(paste0("No id found for ", p ))
      return(NULL)
    }
    data.frame(id=pathways[[p]]@id, pname=p, stringsAsFactors = F)
  })
  do.call(rbind, l)
}
