\name{handle-graphite-pathways}
\alias{mapPathwaysIDfromGraphite}
\title{
  Extract pathway ID from Graphite.
}
\description{
  Given a list of pathways from graphite R package this function maps the pathway id.
}
\usage{
mapPathwaysIDfromGraphite(pathways, pathwayNames = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pathways}{a list of pathways from graphite}
  \item{pathwayNames}{an optional vector of pathway names}
}
\details{maps the ids from graphite pathway list.}
\value{
a data frame with pathway id and pathway name.
}
\author{
Paolo Martini
}
