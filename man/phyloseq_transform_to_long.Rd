% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phyloseq.R
\name{phyloseq_transform_to_long}
\alias{phyloseq_transform_to_long}
\title{Transform a phyloseq object into a long data frame}
\usage{
phyloseq_transform_to_long(ps)
}
\arguments{
\item{ps}{phyloseq object}
}
\value{
Long dataframe with the metadata but without zero, nor NA
}
\description{
Transform a phyloseq object into a long data frame
}
\examples{
df <- phyloseq_transform_to_long(ps)
}
