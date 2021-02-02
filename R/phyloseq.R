#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import stringr
#' @import Biostrings
#' @import dada2
#' @import phyloseq


# phyloseq_import_mothur : Create a phyloseq file from mothur database file ------------------------------

#' @title Create a phyloseq file from mothur database file
#' @description
#' Write out 3 files:
#' * phyloseq object as a rds file: extension .rds
#' * fasta sequences (unaligned): extension .fas
#' * fasta sequences (aligned): extension .aligned.fas
#' @param file_mothur File created by mothur using the database function
#' @param file_samples A tsv file with the sample metadata.  There must be a column "sample_id" with the ID of the samples which should exactly matching the names of the columns in the mothur database file
#' @param n_taxo_levels Number of taxonomic levels (e.g. 8 if using PR2)
#' @return
#' TRUE if no problem
#' @examples
#' phyloseq_import_mothur("mothur.database", "mothur_samples.tsv", 6)
#' @export
#' @md
phyloseq_import_mothur <- function(file_mothur, file_samples, n_taxo_levels){

# Read the samples files
  samples_metadata <- read_tsv(file_samples)
  otu_database <- read_tsv(file_mothur)

## Create the samples, otu and taxonomy tables

  # 1. samples table : row names are labelled by pcr_codes
  samples_df <- samples_metadata
  row.names(samples_df) <- samples_df$sample_id

  # 2. otu table : select only OTU column and abundance columns (Use a regular expression...)
  otu <- otu_database %>% select(-repSeqName, -repSeq, -OTUConTaxonomy)
  row.names(otu) <- otu$OTUNumber
  otu <- otu %>% select(-OTUNumber)

  # 3. Taxonomy table

  #  Separate the taxonomy and removing the bootstrap if present
  taxo_columns <- str_c("Taxo", 1:n_taxo_levels)
  remove_bootstrap <- function(value) {str_replace_all(value,"\\(\\d+\\)","")}

  tax <- otu_database %>% select(OTUNumber, OTUConTaxonomy)  %>%
                          separate(OTUConTaxonomy, ";", into=taxo_columns) %>%
                          mutate_all(remove_bootstrap) %>%
                          column_to_rownames(var = "OTUNumber")

## Create and save to phyloseq object

# Transform into matrixes
  otu_mat <- as.matrix(otu)
  tax_mat <- as.matrix(tax)

# Transform to phyloseq object and save to Rdata file
  OTU = phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE)
  TAX = phyloseq::tax_table(tax_mat)
  samples = phyloseq::sample_data(samples_df)

  phyloseq_mothur <- phyloseq::phyloseq(OTU, TAX, samples)
  saveRDS(phyloseq_mothur, file = str_c(fs::path_file(file_mothur),".rds") )


## Create fasta file with all sequences removing the alignement

  df <-  otu_database %>%  mutate(sequence_unaligned = str_replace_all(repSeq, "(-|\\.)",""),sequence = repSeq )
  seq_out <- Biostrings::DNAStringSet(df$sequence_unaligned)
  names(seq_out) <- otu_database$OTUNumber
  Biostrings::writeXStringSet(seq_out, str_c(fs::path_file(file_mothur),".fas"), width = 20000)
  seq_out <- Biostrings::DNAStringSet(df$sequence)
  names(seq_out) <- otu_database$OTUNumber
  Biostrings::writeXStringSet(seq_out, str_c(fs::path_file(file_mothur),".aligned.fas"), width = 20000)

  return(TRUE)
}

# phyloseq_filter_autotrophic_taxa : Filter a phyloseq table for autotrophs ------------------------------

#' @title Filter a phyloseq table keeping only autotrophic taxa
#' @description
#' @param ps phyloseq object
#' @return
#' ps object with only autotrophs taxa
#' @examples
#' ps_auto <- phyloseq_filter_autotrophic_taxa(ps)
#' @export
#' @md
phyloseq_filter_autotrophic_taxa <- function(ps) {

  ps <- subset_taxa(ps,  (division %in% c("Chlorophyta", "Cryptophyta", "Rhodophyta",
                                          "Haptophyta", "Ochrophyta")) |
                          ((division ==  "Dinoflagellata") & (class != "Syndiniales")) |
                          (class == "Filosa-Chlorarachnea"))
}

# phyloseq_filter_abundant_taxa : Filter a phyloseq table for abundant taxa ------------------------------

#' @title Filter a phyloseq table keeping only abundant taxa
#' @description
#' @param ps phyloseq object
#' @param fraction_min minimum fraction that an OTU must represent in any given samples (not overall!) - 0.10 corresponds to 10%
#' @return
#' ps object with only abundant taxa
#' @examples
#' # Will return the otus that more that represent more than 10% of total in any given sample. It also normalize the phyloseq to the median of the remaining OTUs
#' ps_abundant <- phyloseq_filter_abundant_taxa(ps)
#' # Idem but with lower threshold (5%)
#' ps_abundant <- phyloseq_filter_abundant_taxa(ps, 0.05)
#' @export
#' @md
phyloseq_filter_abundant_taxa <- function(ps, fraction_min=0.10) {

     total_per_sample <- max(sample_sums(ps))
     ps <- ps %>%
       filter_taxa(function(x) sum(x > total_per_sample*fraction_min) > 0, TRUE) %>%
       phyloseq_normalize_median()
}

# phyloseq_normalize_median : Normalize to the median number of reads ------------------------------

#' @title Normalize to the median number of reads
#' @description
#' @param ps phyloseq object
#' @return
#' ps object normalized
#' @examples
#' ps_norm <- phyloseq_normalize_median(ps)
#' @export
#' @md
phyloseq_normalize_median <- function (ps) {
  ps_median = median(sample_sums(ps))
  normalize_median = function(x, t=ps_median) (if(sum(x) > 0){ t * (x / sum(x))} else {x})
  ps = transform_sample_counts(ps, normalize_median)
  cat(str_c("\n========== \n") )
  print(ps)
  cat(sprintf("\n==========\nThe median number of reads used for normalization is  %.0f", ps_median))
  return(ps)
}

# phyloseq_normalize_percent : Normalize to the percentage for each sample ------------------------------

#' @title Normalize as percent of total reads
#' @description
#' @param ps phyloseq object
#' @return
#' ps object normalized
#' @examples
#' ps_percent <- phyloseq_normalize_percent(ps)
#' @export
#' @md
phyloseq_normalize_percent <- function (ps) {
  normalize_percent = function(x) (if(sum(x) > 0){ x / sum(x)} else {x})
  ps = transform_sample_counts(ps, normalize_percent)
  return(ps)
}


# phyloseq_transform_to_long : Transform a phyloseq object into a long data frame ------------------------------

#' @title Transform a phyloseq object into a long data frame
#' @description
#' @param ps phyloseq object
#' @return
#' Long dataframe with the metadata but without zero, nor NA
#' @examples
#' df <- phyloseq_transform_to_long(ps)
#' @export
#' @md
  phyloseq_transform_to_long <- function(ps) {
    otu_df <- as.data.frame(ps@otu_table@.Data, stringsAsFactors = FALSE) %>%
      rownames_to_column(var = "asv_code") %>%
      pivot_longer(cols = -asv_code,
                   names_to = "file_code",
                   values_to = "n_reads",
                   values_drop_na = TRUE) %>%
      filter(n_reads != 0) %>%
      filter(!is.na(n_reads))

  # See https://github.com/joey711/phyloseq/issues/983
    taxo_df <- as.data.frame(ps@tax_table@.Data, stringsAsFactors = FALSE) %>%
      rownames_to_column(var = "asv_code")

    otu_df <- left_join(taxo_df, otu_df)

    metadata_df <- data.frame(sample_data(ps)) %>%
      rownames_to_column(var = "file_code")

    otu_df <- left_join(otu_df, metadata_df, by = c("file_code"))

    return(otu_df)

  }

# phyloseq_long_treemap : Do a treemap based on the long version of a phyloseq file ------

#' @title Do a treemap based on the long version of a phyloseq file
#' @description Plot the treemaps and returns a list with a ggplot and a df with the summary of the data
#' @param df Data frame obtained from a phyloseq file using the function phyloseq_transform_to_long
#' @param group1 first grouping level (do not quote)
#' @param group2 second grouping level (do not quote)
#' @param title Title for the treemap
#' @param colors If NULL then the default viridis palette is used.  If named vectors then the colors
#' @param label_group1 If true, the higher level is labeled.  If false, only boundaries are marked
#' @return
#' Plot the treemap
#' Returns a list made of 2 elements
#'   * gg = the treemap as ggplot
#'   * df =  a dataframe with the summarized data
#' @examples
#' my_list <- phyloseq_long_treemap(phyloseq_long, division, class, "Singapore strait", colors=
#'                                  structure(c("black", "white"),.Names="Mamiellophyceae", "Dinophyceae"))
#' @export
#' @md
phyloseq_long_treemap <- function(df, group1, group2, title, colors=NULL, label_group1 = TRUE) {

 df <- df %>%
   group_by({{group1}}, {{group2}}) %>%
   summarise(n_reads=sum(n_reads, na.rm = TRUE))

 gg <- ggplot(df, aes(area = n_reads,
                             fill = {{group2}},
                             label = {{group2}},
                             subgroup = {{group1}})) +
    ggtitle(title) +
    treemapify::geom_treemap() +
    treemapify::geom_treemap_subgroup_border() +
    treemapify::geom_treemap_text(colour = "black", place = "topleft", reflow = T,
                                  padding.x =  grid::unit(3, "mm"),
                                  padding.y = grid::unit(3, "mm") )  +
    theme(legend.position="none", plot.title = element_text(size = 16, face = "bold"))

 if (label_group1){
     gg <- gg +
        treemapify::geom_treemap_subgroup_text(place = "centre", grow = T, alpha = 0.5, colour =
                                               "white", fontface = "italic", min.size = 0)
 }

 if (is.null(colors)){
    gg <- gg + scale_fill_viridis_d()
 } else {
    gg <- gg + scale_fill_manual(values = colors)
 }

 print(gg)
 treemap_list <- list(gg = gg, df=df)
 return(treemap_list)
 }

# phyloseq_long_asv_bargraph : Do a bar graph of top asvs based on the long version of a phyloseq file ------

#' @title Do a bar graph of top taxo_level or asvs based on the long version of a phyloseq file
#' @description Plot the bar graph and returns a list with a ggplot and a df with the summary of the data
#' @param df Data frame obtained from a phyloseq file using the function phyloseq_transform_to_long
#' @param n_bars numbers of bars to plot
#' @param title Title for the treemap
#' @param text_scaling Scaling for the text of the graph
#' @param use_asv If TRUE use asvs, if FALSE use the taxo_level selected
#' @param taxo_level Taxonomic level to use when use_asv = FALSE (should be between species up to class)
#' @param division_colors Colors to be used for the different divisions
#' @return
#' Plot the bargraph
#' Returns a list made of 2 elements
#'   * gg = the bargraph as ggplot
#'   * df =  a dataframe with the summarized data
#' @examples
#' my_list <- phyloseq_long_asv_bargraph(phyloseq_long, title= "Antartica", n_bars=10, use_asv = FALSE, taxo_level= genus, text_scaling=1)
#' @export
#' @md

phyloseq_long_bargraph <- function(df, n_bars=30, title="", text_scaling = 0.75,
                                   use_asv = TRUE,
                                   taxo_level= species,
                                   taxo_level_fill = division,
                                   taxo_colors_fill = structure(c("green", "orange", "red", "blue", "brown"),
                                                      .Names=c("Chlorophyta", "Cryptophyta", "Rhodophyta","Haptophyta", "Ochrophyta"))) {
# Define labels for taxa
  if (use_asv) {
    df <- df %>%
      mutate(bar_label = case_when (kingdom == "Eukaryota" ~ str_c(asv_code, species, sep="-"),
                                    TRUE ~  str_c(asv_code, family, sep="-") ))
  } else {
    df <- df %>%
      mutate(bar_label = {{taxo_level}})
  }

  # Clean up the labels
  df <- df %>%
    mutate(bar_label = str_replace_all(bar_label, c("_" = " ",  # Underscocre by space
                                                    "X+" = "",  # One or more X by nothing
                                                    " +" = " ", # One or more space by a single space
                                                    "Radial-centric-basal-" = ""
                                                    )))

  df <- df %>%
    group_by({{taxo_level_fill}}, bar_label) %>%
    summarize(n_reads = sum(n_reads, na.rm = TRUE)) %>%
    arrange(desc(n_reads)) %>%
    filter(n_reads > 0) %>%
    ungroup()

  gg <- ggplot(top_n(df,n=n_bars, wt=n_reads)) +
    geom_col(aes(x=reorder(bar_label, n_reads), y=n_reads, fill={{taxo_level_fill}})) +
    coord_flip() +
    theme_bw() +
    theme(axis.text.x = element_text(size = 16*text_scaling, angle = 0, hjust = 1, vjust = 1)) +
    theme(axis.text.y = element_text(size = 16*text_scaling, angle = 0, hjust = 0, vjust = 0)) +
    theme(legend.title = element_text(size = 24*text_scaling)) +
    theme(legend.text = element_text(size = 16*text_scaling)) +
    xlab("") + ylab("Number of reads") +
    scale_fill_manual(values = taxo_colors_fill, drop = FALSE) +
    ggtitle(title)  +
    theme(plot.title = element_text(size=22*text_scaling, hjust = 0.5)) +
    # theme(axis.text=element_text(size=14), legend.text = element_text(size=16)) +
    theme(legend.position = "top", legend.box = "vertical") +
    guides(fill = guide_legend(title.position="top",
                               ncol = 4, byrow = TRUE))

   print(gg)
   treemap_list <- list(gg = gg, df=df)
   return(treemap_list)

}

# phyloseq_nmds : Do NMDS of a phyloseq file ------

#' @title Do a bar graph of top taxo_level or asvs based on the long version of a phyloseq file
#' @description Plot the bar graph and returns a list with a ggplot and a df with the summary of the data
#' @param ps Phyloseq file
#' @param title Title for NMDS plots
#' @param env_parameters Environment parameters to use for the envfit
#' @param classifying_parameters Parameters used to color and shape the samples points - Unquoted
#' @param sample_color Parameter used to color the samples (e.g. fraction_name) - Unquoted
#' @param sample_shape Parameter used for shape of the samples (e.g. season) - Unquoted
#' @param sample_label Parameter used to label the samples (e.g. sample_label) - Unquoted
#' @param taxo_level Taxonomic level to color the OTUs (e.g. division) - Unquoted
#' @param taxo_colors Colors to be used to color the OTUs
#' Plot the bargraph
#' Returns a list made of 2 elements (ggplots)
#'   * gg_samples = samples with env parameters as ggplot
#'   * gg_taxa =  taxa colored
#' @examples
#' print(phyloseq_long_asv_bargraph(phyloseq_nmds, title= "Antartica")
#' @export
#' @md

phyloseq_nmds <- function(ps,
                          title="",
                          env_parameters = c("temperature", "salinity", "Chla", "NO3", "PO4", "Si"),
                          classifying_parameters = c("fraction_name", "season"),
                          sample_color = fraction_name,
                          sample_shape = season,
                          sample_label = sample_label,
                          taxo_level = division,
                          taxo_colors = structure(c("green", "orange", "red", "blue", "brown"),
                                      .Names=c("Chlorophyta", "Cryptophyta", "Rhodophyta","Haptophyta", "Ochrophyta"))) {


      nmds.ord <- ordinate(ps, "NMDS", "bray")

      # Fit environmental parameters
       env_var <- sample_variables(ps)
       env_matrix <- get_variable(ps,env_parameters )
       env_fit <- vegan::envfit(nmds.ord, env = env_matrix, perm = 999, na.rm=TRUE)
       env_arrows <- data.frame(env_fit$vectors$arrows*sqrt(env_fit$vectors$r)) %>% rownames_to_column(var="parameter")


      nmds_samples <- data.frame(nmds.ord[["points"]],
                                 get_variable(ps, classifying_parameters) ) %>%
                      rownames_to_column(var="sample")

      # Factor to move the labels
      nudge_x <- max(nmds_samples$MDS1)*0.08
      nudge_y <- max(nmds_samples$MDS2)*0.08
      xy_max = max(c(nmds_samples$MDS1, nmds_samples$MDS2))*1.5
      xy_min = min(c(nmds_samples$MDS1, nmds_samples$MDS2))*1.5
      factor <- 3 # for vectors for euks

      print(factor)

    gg1 <- plot_ordination(ps, nmds.ord, type="samples",
                           color=quo_text(enquo(sample_color)),
                           shape=quo_text(enquo(sample_shape)),
                           title=title) +
           geom_point(aes(shape={{sample_shape}}, color={{sample_color}}), size=3.5)  +
           scale_color_viridis_d() +
           # scale_shape_manual(values = c(15,16)) +
           geom_text(aes(label=sample_label, color={{sample_color}}),
                     nudge_x=nudge_x, nudge_y=nudge_y,
                     check_overlap = TRUE, size=2) +
           theme_bw() +
           geom_segment(data=env_arrows, aes(x=0,xend=NMDS1*factor, y=0,yend=NMDS2*factor), inherit.aes = FALSE,
                        arrow = arrow(length = unit(0.5, "cm")), colour = "black")  +
           geom_text(data=env_arrows, aes(x=NMDS1*factor, y=NMDS2*factor, label=parameter), inherit.aes = FALSE ,
                     hjust=-0.2, vjust=-0.2, size=3)
      print(nmds.ord)
      print(gg1)

      # The following lines can be used if you want to avoid using the pjhyloseq functions to plot the data.
      # Notes : - must use inherit.aes = FALSE to add some extra layers
      #         - the saved plots have a different scale for the added layer than the displaued plot can figure out
           # ggplot()+
           # coord_fixed() +
           # xlim(xy_min, xy_max) + ylim(xy_min, xy_max) +
           # geom_point(data=nmds_samples, aes(x=MDS1, y=MDS2, shape=strait, color=monsoon), size=5)  +
           # geom_text(data=nmds_samples, aes(x=MDS1, y=MDS2, label=location_label, color=monsoon),
           #           nudge_x=nudge_x, nudge_y=nudge_y,
           #           check_overlap = FALSE, size=2) +
           # ggtitle(ps_list$title[[i]]) +

    gg2 <- plot_ordination(ps, nmds.ord, type="taxa", color=quo_text(enquo(taxo_level)), title=title) +
       scale_color_manual(values = taxo_colors) +
       geom_point(size=3) +
       theme_bw() +
       geom_segment(data=env_arrows, aes(x=0,xend=NMDS1*factor, y=0,yend=NMDS2*factor), inherit.aes = FALSE,
                    arrow = arrow(length = unit(0.5, "cm")), colour = "black")  +
       geom_text(data=env_arrows, aes(x=NMDS1*factor, y=NMDS2*factor, label=parameter), inherit.aes = FALSE ,
                 hjust=-0.2, vjust=-0.2, size=3)

    print(gg2)


  return(list(gg_samples = gg1, gg_taxa = gg2))
  }

