% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rcc.R
\name{rcc_export}
\alias{rcc_export}
\title{Export tables from MySQL database}
\usage{
rcc_export(export_directory = "C:/Databases/web_export/")
}
\arguments{
\item{export_directory}{By defaut "C:/Databases/web_export/"}
}
\value{
TRUE if succesful
}
\description{
This function exports the data from the rcc database to
3 files to be imported into the RCC web site version 1.0
\itemize{
\item web_sequences.csv
\item web_strains.csv
\item web_strains_products.csv
}
}
\examples{
rcc_export()
}
