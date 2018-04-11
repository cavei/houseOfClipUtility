\name{pathwayRelationshipProcessing}
\alias{getPathFathers}
\alias{id2pathwayName}
\title{Extract Pathway Parents.}
\description{
  Given a pathway name extract the father or parents up to a desire order. other Function allows to further process.
}
\usage{
getPathFathers(pathway, hierarchyGraph, ord = 3, plot = F)
id2pathwayName(path2fathers, pathwayDict)
}
\arguments{
  \item{pathway}{a pathway name}
  \item{hierarchyGraph}{an igraph with the heirarchy}
  \item{ord}{the distance of the father}
  \item{plot}{logical: plot or not}
  \item{path2fathers}{list of pathway ids}
  \item{pathwayDict}{a dictionary, id 2 pathway name}
}
\value{
  A character vector woth the names of the fathers.
}
\author{
Paolo Martini
}
