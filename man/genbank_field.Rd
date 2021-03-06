% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genbank.R
\name{genbank_field}
\alias{genbank_field}
\title{Read a single field from a set of existing GenBank files}
\usage{
genbank_field(accession, directory, field_name)
}
\arguments{
\item{accession}{Character vector of accession numbers}

\item{directory}{Character, directory name which must end by "/"}

\item{field_name}{must be unquoted !!!}
}
\value{
Vector with the field information.
}
\description{
Return a vector containing the value of a GenBank field for a vector of accession numbers. The GenBank files must exist
}
\examples{
gb_organism <- genbank_field(c("JX015376", "JQ768406", "LT621940"), "genbank/", organism)
}
