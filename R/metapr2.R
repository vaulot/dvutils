#' @import dplyr
#' @import stringr
#' @import openxlsx

# metapr2_export_asv -------------------------------------------------------
#' @title Exports the metapr2 database
#' @description
#' Exports a range of file and format.
#'    * fasta file with the full taxonomy or just the genus level
#'    * excel file with the full table of the asv and metadata
#'    * phyloseq file (better when only selecting a single data set)
#'
#' Returns a list with 3 elements
#'    * df = dataframe with all the asv and the stations and read numbers
#'    * ps = phyloseq object
#'    * fasta = df with 2 columns seq_name and sequence
#' @param taxo_level The taxonomic level for selection (do not quote), e.g. class or genus
#' @param taxo_name  The name of the taxonomic level selected, e.g. Chlorophyta
#' @param boot_level The taxonomic level for bootstrap filtering (do not quote), e.g. class_boot or genus_boot
#' @param boot_min  Minimum bootstrap value at the class level, 0 if you want to get all asvs
#' @param directory  Directory where the files are saved (finish with /)
#' @param dataset_id_selected  Integer vector, e.g. 23 or c(21, 23) or 21:23
#' @param export_xls  If TRUE, an xlsx file is produced containing the final data frame
#' @param export_phyloseq  If TRUE, a phyloseq file is produced
#' @param taxonomy_full If TRUE, the fasta file contains the full taxonomy (8 levels), if false only contains the species
#' @return
#' A data frame with asv and metadata for each sample where found.
#' @examples
#' # Export all the asv in a single fasta
#'   metapr2_export_asv()
#'
#' # Export as specific data set as a phyloseq file
#'   metapr2_export_asv(dataset_id_selected = 23, export_phyloseq = TRUE)
#'
#' # Export a specific genus as a fasta file and an excel file
#'   asv_set <- metapr2_export_asv(taxo_level = genus, taxo_name="Pseudohaptolina",
#'                                 taxonomy_full= FALSE, boot_min = 100, export_xls = TRUE)
#' @export
#' @md
#'
metapr2_export_asv <- function(taxo_level = kingdom, taxo_name="Eukaryota",
                               boot_level = class_boot, boot_min = 0,
                               directory = "C:/daniel.vaulot@gmail.com/Databases/_metaPR2/export/",
                               dataset_id_selected = c(1:100),
                               export_xls=FALSE,  export_phyloseq = FALSE, taxonomy_full = TRUE) {

# Define a variable to hold the data set id as character
  dataset_id_char <- case_when ((100 %in% dataset_id_selected) ~ "all",
                                TRUE ~ str_c(dataset_id_selected, collapse = "_"))

# Read the database and already filter for selected records
  metapr2_db <- db_info("metapr2_local")
  metapr2_db_con <- db_connect(metapr2_db)

  asv_set <- tbl(metapr2_db_con, "metapr2_asv") %>%
     filter({{taxo_level}} == taxo_name) %>%
     collect()

  asv_set <- asv_set %>%
     filter(dataset_id %in% dataset_id_selected) %>%
     filter({{boot_level}} >= boot_min)

# Need to use !! before the local varibale
  # Error: Cannot embed a data frame in a SQL query.
  #  If you are seeing this error in code that used to work, the most likely cause is a change dbplyr 1.4.0. Previously `df$x` or
  # `df[[y]]` implied that `df` was a local variable, but now you must make that explict with `!!` or `local()`, e.g., `!!df$x` or
  # `local(df[["y"]))

  metapr2_asv_abundance <- tbl(metapr2_db_con, "metapr2_asv_abundance") %>%
    filter(asv_code %in% !!asv_set$asv_code) %>%
    collect()

  metapr2_samples <- tbl(metapr2_db_con, "metapr2_samples") %>% collect()

  metapr2_metadata <- tbl(metapr2_db_con, "metapr2_metadata") %>% collect()

  metapr2_datasets <- tbl(metapr2_db_con, "metapr2_datasets") %>% collect()

  db_disconnect(metapr2_db_con)

  if ( taxonomy_full) {
    asv_fasta <- asv_set %>%
      rename(seq_name = asv_code)
    fasta_write(asv_fasta, str_c(directory, "metapr2_asv_set_", dataset_id_char ,"_", taxo_name, ".fasta"), compress = FALSE, taxo_include = TRUE)
  } else {
    asv_fasta <- asv_set %>%
       mutate(seq_name = str_c(asv_code, species, sep="|") )
    fasta_write(asv_fasta, str_c(directory, "metapr2_asv_set_", dataset_id_char ,"_", taxo_name, ".fasta"), compress = FALSE, taxo_include = FALSE)
  }

 asv_set <-  asv_set %>%
   left_join(metapr2_asv_abundance) %>%
   left_join(metapr2_samples) %>%
   left_join(metapr2_metadata) %>%
   left_join(select(metapr2_datasets, -gene, -gene_region), by = c("dataset_id" = "dataset_id"))

 if (export_xls) openxlsx::write.xlsx(asv_set, str_c(directory, "metapr2_asv_set_", dataset_id_char ,"_", taxo_name, ".xlsx"))



   ## Create the samples, otu and taxonomy tables

  # 1. samples table : row names are labelled by file_code
  samples_df <- asv_set %>%
    select(file_code:metadata_remark, -n_reads) %>%
    distinct(file_code, .keep_all = TRUE) %>%
    column_to_rownames(var = "file_code")


  # 2. otu table :
  otu <- asv_set %>%
    select(asv_code, file_code, n_reads) %>%
    spread(key=file_code, value = n_reads, fill=0) %>%
    column_to_rownames(var = "asv_code")

  # 3. Taxonomy table

  tax <-  asv_set %>%
    select(asv_code, kingdom:species) %>%
    distinct(asv_code, .keep_all = TRUE) %>%
    column_to_rownames(var = "asv_code")

## Create and save to phyloseq object

# Transform into matrixes
  otu_mat <- as.matrix(otu)
  tax_mat <- as.matrix(tax)

# Transform to phyloseq object and save to Rdata file
  OTU = phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE)
  TAX = phyloseq::tax_table(tax_mat)
  samples = phyloseq::sample_data(samples_df)

  phyloseq_asv <- phyloseq::phyloseq(OTU, TAX, samples)

 if (export_phyloseq) {
  saveRDS(phyloseq_asv, file = str_c(directory, "phyloseq_metapr2_asv_set_", dataset_id_char ,"_", taxo_name, ".rda") )

 }

 return(list(df=asv_set, ps=phyloseq_asv, fasta=asv_fasta))

  }
