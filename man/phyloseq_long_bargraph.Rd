% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phyloseq.R
\name{phyloseq_long_bargraph}
\alias{phyloseq_long_bargraph}
\title{Do a bar graph of top taxo_level or asvs based on the long version of a phyloseq file}
\usage{
phyloseq_long_bargraph(
  df,
  n_bars = 30,
  title = "",
  text_scaling = 0.75,
  use_asv = TRUE,
  taxo_level = species,
  taxo_level_fill = division,
  taxo_colors_fill = structure(c("green", "orange", "red", "blue", "brown"), .Names =
    c("Chlorophyta", "Cryptophyta", "Rhodophyta", "Haptophyta", "Ochrophyta"))
)
}
\arguments{
\item{df}{Data frame obtained from a phyloseq file using the function phyloseq_transform_to_long}

\item{n_bars}{numbers of bars to plot}

\item{title}{Title for the treemap}

\item{text_scaling}{Scaling for the text of the graph}

\item{use_asv}{If TRUE use asvs, if FALSE use the taxo_level selected}

\item{taxo_level}{Taxonomic level to use when use_asv = FALSE (should be between species up to class)}

\item{division_colors}{Colors to be used for the different divisions}
}
\value{
Plot the bargraph
Returns a list made of 2 elements
\itemize{
\item gg = the bargraph as ggplot
\item df =  a dataframe with the summarized data
}
}
\description{
Plot the bar graph and returns a list with a ggplot and a df with the summary of the data
}
\examples{
my_list <- phyloseq_long_asv_bargraph(phyloseq_long, title= "Antartica", n_bars=10, use_asv = FALSE, taxo_level= genus, text_scaling=1)
}
