% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pr2_management.R
\name{pr2_taxo_check}
\alias{pr2_taxo_check}
\title{Check taxonomy}
\usage{
pr2_taxo_check(pr2_taxo, taxo_levels = pr2.env$taxo_levels, dir_taxo = "")
}
\arguments{
\item{pr2_taxo}{dataframe - should contain eight or nine columns kingdom -> species or domain -> species}

\item{taxo_levels}{\itemize{
\item vector with taxo levels
}}

\item{dir_taxo}{character - full path where file are saved (without final //)}
}
\value{
Write several files (see description)
}
\description{
This check the taxonomy for the following conditions:
\itemize{
\item No taxon has 2 different parents
\item Not taxon appears at different ranks
}
}
\note{
The input data frame does not need to be the whole pr2 taxo, can be taxo extracted from an Excel file.
Several files are produced:
\itemize{
\item \emph{taxo_list_(rank_number)_(rank_name).txt} (7 files)
\itemize{
\item name of each taxon
\item taxa_number - number of species in this taxon
\item level of taxon - from 1 to 8 (1 for kingdom or domain)
\item id of the parent taxon (at rank i-1)
\item id of the taxon (from 1 to n for each rank)
}
\item \emph{taxo_list_2_parents_(rank_number)_(rank_name).txt} (7 files).  This is an important file because it list all the taxa that have 2 parents.  When all the files are completely empty, then the taxonomic list is OK.
\itemize{
\item name of the taxon
\item number of parents
}
\item \emph{taxo_level_duplicates.txt}. List the taxa that are found at different ranks.  Fro example Cryptophyta can be in the Class column for one species and in the Division levels for all other species.
\itemize{
\item name
\item number of times found
\item taxa number
\item level at which found
\item id
\item id of the parent
}
\item \emph{taxo_species_duplicates}. List of species that appear at least twice in the species list.
}
}
\examples{
pr2_taxo_check (taxo_list, "C:/Daniel/database")
}
