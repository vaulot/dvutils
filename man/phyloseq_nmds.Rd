% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phyloseq.R
\name{phyloseq_nmds}
\alias{phyloseq_nmds}
\title{Do a bar graph of top taxo_level or asvs based on the long version of a phyloseq file}
\usage{
phyloseq_nmds(
  ps,
  title = "",
  env_parameters = c("temperature", "salinity", "Chla", "NO3", "PO4", "Si"),
  classifying_parameters = c("fraction_name", "season"),
  sample_color = fraction_name,
  sample_shape = season,
  sample_label = sample_label,
  taxo_level = division,
  taxo_colors = structure(c("green", "orange", "red", "blue", "brown"), .Names =
    c("Chlorophyta", "Cryptophyta", "Rhodophyta", "Haptophyta", "Ochrophyta"))
)
}
\arguments{
\item{ps}{Phyloseq file}

\item{title}{Title for NMDS plots}

\item{env_parameters}{Environment parameters to use for the envfit}

\item{classifying_parameters}{Parameters used to color and shape the samples points - Unquoted}

\item{sample_color}{Parameter used to color the samples (e.g. fraction_name) - Unquoted}

\item{sample_shape}{Parameter used for shape of the samples (e.g. season) - Unquoted}

\item{sample_label}{Parameter used to label the samples (e.g. sample_label) - Unquoted}

\item{taxo_level}{Taxonomic level to color the OTUs (e.g. division) - Unquoted}

\item{taxo_colors}{Colors to be used to color the OTUs
Plot the bargraph
Returns a list made of 2 elements (ggplots)
\itemize{
\item gg_samples = samples with env parameters as ggplot
\item gg_taxa =  taxa colored
}}
}
\description{
Plot the bar graph and returns a list with a ggplot and a df with the summary of the data
}
\examples{
print(phyloseq_long_asv_bargraph(phyloseq_nmds, title= "Antartica")
}
