% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db.R
\name{db_get_query}
\alias{db_get_query}
\title{Read a query into a dataframe.}
\usage{
db_get_query(database_info, query)
}
\arguments{
\item{query}{Character - SQL Query}
}
\value{
The value returned by the query
}
\description{
Read a query into a dataframe.
}
\examples{
taxo <- db_get_query(pr2_db, taxo_query)
}
