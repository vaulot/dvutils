% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/maps.R
\name{map_world}
\alias{map_world}
\title{Background map using the maps package}
\usage{
map_world(
  color_continents = "grey80",
  color_borders = "white",
  resolution = "coarse"
)
}
\arguments{
\item{color_continents}{color for the continents (grey80 by default)}

\item{color_borders}{color for the borders (white by default)}

\item{resolution}{"coarse" for global maps, "low" for regional maps}
}
\value{
A map as a ggplot object
}
\description{
Background map using the maps package
}
\examples{
map_world() + geom_point(data=cultures_one_year, aes(x=Longitude, y=Latitude), fill="blue", size=2, shape=21)
}
