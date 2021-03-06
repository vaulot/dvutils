% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pr2_management.R
\name{pr2_extract_one_taxo_level}
\alias{pr2_extract_one_taxo_level}
\title{Extract one level of taxonomy}
\usage{
pr2_extract_one_taxo_level(
  pr2_taxo,
  level_above_name,
  level_name,
  level_number,
  taxo_above
)
}
\arguments{
\item{pr2_taxo}{dataframe - taxonomy with 8 coumns kigdom -> species (but can be more)}

\item{level_above_name}{character - name of the taxo level above (e.g. "order")}

\item{level_name}{character - name the level considered (e.g. "family")}

\item{taxo_above}{dataframe -  the same structure than the output of the function but for the level above (e.g. the order if we are the family level)}
}
\description{
This function returns a data frame that contains:
}
\details{
\itemize{
\item name = name of the taxon
\item taxa_number = number of taxa downstream of this taxon
\item level = number of the level (e.g. genus = level number 7)
\item parend_id = id of the first parent taxon  (sometimes one taxon has more than one taxon)
\item id = id of the taxon (just the row name since they are ordered alphabetically )
}
}
\note{
Here is an example at the genus level
name	            taxa_number	    level	    parent_id	    id
Ascampbelliella	    1	              7	        51	          1
Choreotrichida_XX	  1	              7	        1	            2
Cyttarocylididae_X	1	              7	        52	          3
Cyttarocylis	      4	              7	        52	          4
Petalotricha	      1	              7	        52	          5
Codonaria	          2	              7	        53	          6
Codonella	          2	              7	        53	          7
}
\examples{
genus_df <- pr2_extract_one_taxo_level(pr2_taxo, "family", "genus", 7, family_df)
}
