% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genbank.R
\name{genbank_download_parse}
\alias{genbank_download_parse}
\title{Download and parse sequences from GenBank}
\usage{
genbank_download_parse(
  accession,
  directory,
  sequence_keep = TRUE,
  store_file = FALSE
)
}
\arguments{
\item{accession}{Character vector of accession numbers}

\item{directory}{Character, directory name which must end by "/"}

\item{sequence_keep}{Logical, if FALSE the sequence is not returned in the file data frame}

\item{store_file}{Logical, if FALSE the file is not stored and the sequence is not returned in the file data frame}
}
\value{
Data frame with the metadata information.
}
\description{
If the GenBank file does not exist it fetches and writes it.
If the file exists it just reads it. It extracts the metadata and return into a a data frame.

Note: use only sequence_keep = TRUE for simple gene sequences, else it may crash
}
\examples{
gb_metadata <- genbank_download_parse(accession = c("JX015376", "JQ768406", "LT621940"), directory = "genbank/", sequence_keep=TRUE)
}
