#' @import dplyr
#' @import tidyr
#' @import stringr
#' @import openxlsx
#' @import tibble

# metapr2_export_asv -------------------------------------------------------
#' @title Exports the metapr2 database
#' @description
#' Exports a range of file and format.
#'    * fasta file with the full taxonomy or just the genus level
#'    * excel file with the full table of the asv and metadata
#'    * phyloseq file (better when only selecting a single data set)
#'    * list of samples with metadata
#'
#' Returns a list with 4 elements
#'    * df = dataframe with all the asv and the stations and read numbers
#'    * ps = phyloseq object  - if phyloseq = TRUE
#'    * fasta = df with 2 columns seq_name and sequence
#'    * samples= sample_list (list of sammples with metadata)
#'
#' NOTE: if all export_long_xls, export_wide_xls, export_phyloseq are false, abundances are not exported
#' @param taxo_level The taxonomic level for selection (do not quote), e.g. class or genus
#' @param taxo_name  The name of the taxonomic level selected, e.g. "Chlorophyta", can be a vector c("Chlorophyta", "Haptophyta")
#' @param boot_level The taxonomic level for bootstrap filtering (do not quote), e.g. class_boot or genus_boot
#' @param boot_min  Minimum bootstrap value at the class level, 0 if you want to get all asvs
#' @param directory  Directory where the files are saved (finish with /)
#' @param dataset_id_selected  Integer vector, e.g. 23 or c(21, 23) or 21:23
#' @param export_long_xls  If TRUE, an xlsx file is produced containing the final long data frame
#' @param export_wide_xls  If TRUE, an xlsx file is produced containing the final wide data frame
#' @param export_sample_xls  If TRUE, an xlsx file is produced containing the sample list
#' @param export_phyloseq  If TRUE, a phyloseq file is produced and a phyloseq object producted
#' @param export_fasta  If TRUE, a fasta is produced
#' @param taxonomy_full If TRUE, the fasta file contains the full taxonomy (8 levels), if false only contains the species
#' @param use_hash If TRUE, the asvs with identical has will be merged and called by their hash value (sequence_hash)
#' @param sum_reads_min This is the minimum number of reads (summed over the datasets selected) for an asv to be included
#' @return
#' A list with 4 elements
#' @examples
#' # Export all the asv in a single fasta
#'   metapr2_export_asv()
#'
#' # Export as specific data set as a phyloseq file
#'   metapr2_export_asv(dataset_id_selected = 23, export_phyloseq = TRUE)
#'
#' # Export a specific genus as a fasta file and an excel file
#'   asv_set <- metapr2_export_asv(taxo_level = genus, taxo_name=c("Pseudohaptolina","Haptolina"),
#'                                 export_fasta=TRUE, taxonomy_full= FALSE,
#'                                 boot_level = genus_boot, boot_min = 100,
#'                                 export_long_xls = TRUE)
#' @export
#' @md
#'
metapr2_export_asv <- function(taxo_level = kingdom, taxo_name="Eukaryota",
                               boot_level = class_boot, boot_min = 0,
                               directory = "C:/daniel.vaulot@gmail.com/Databases/_metaPR2/export/",
                               dataset_id_selected = c(1:100),
                               export_long_xls=FALSE,
                               export_wide_xls=FALSE,
                               export_sample_xls=FALSE,
                               export_phyloseq = FALSE,
                               export_fasta=FALSE,
                               taxonomy_full = TRUE,
                               use_hash = FALSE,
                               sum_reads_min = 0) {

# Define a variable to hold the data set id as character
  dataset_id_char <- case_when ((100 %in% dataset_id_selected) ~ "all",
                                TRUE ~ str_c(dataset_id_selected, collapse = "_"))


# Read the database and already filter for selected records
  metapr2_db <- db_info("metapr2_google")
  metapr2_db_con <- db_connect(metapr2_db)

# Get the asv filtered by taxonomy
  asv_set <- tbl(metapr2_db_con, "metapr2_asv") %>%
     filter({{taxo_level}} %in% taxo_name) %>%
     collect()

# Create a label for the taxons selected to be used for the files
  if(length(taxo_name) > 1){
    taxo_name_label <- str_c(str_sub(taxo_name, 1, 5), collapse="_")
  } else {
    taxo_name_label <- taxo_name
    }

# Filter the asv based on the datasets selected and the minimum bootstrap
  asv_set <- asv_set %>%
     filter(dataset_id %in% dataset_id_selected) %>%
     filter({{boot_level}} >= boot_min)

# Need to use !! before the local variable
  # Error: Cannot embed a data frame in a SQL query.
  #  If you are seeing this error in code that used to work, the most likely cause is a change dbplyr 1.4.0. Previously `df$x` or
  # `df[[y]]` implied that `df` was a local variable, but now you must make that explict with `!!` or `local()`, e.g., `!!df$x` or
  # `local(df[["y"]))

  metapr2_asv_abundance <- tbl(metapr2_db_con, "metapr2_asv_abundance") %>%
    filter(asv_code %in% !!asv_set$asv_code) %>%
    collect()

  metapr2_samples <- tbl(metapr2_db_con, "metapr2_samples") %>%
    filter(!is.na(file_code)) %>%   # Remove all samples that have not been processed
    collect()

  metapr2_metadata <- tbl(metapr2_db_con, "metapr2_metadata") %>% collect()

  metapr2_datasets <- tbl(metapr2_db_con, "metapr2_datasets") %>% collect()

  db_disconnect(metapr2_db_con)

# Merging together asvs with same hash value

  if (use_hash){
    asv_fasta <- asv_set %>%
        group_by(sequence_hash) %>%
        slice(1) %>%
        mutate(asv_code = sequence_hash) %>%
        ungroup()
  } else {
    asv_fasta <- asv_set
  }


  sample_list <- metapr2_samples %>%
   left_join(metapr2_metadata) %>%
   filter(dataset_id %in% dataset_id_selected) %>%
   select_if(~!all(is.na(.))) # This removes all the column that are empty

  asv_set <-  asv_set %>%
   left_join(metapr2_asv_abundance) %>%
   left_join(metapr2_samples) %>%
   left_join(metapr2_metadata) %>%
   left_join(select(metapr2_datasets, -gene, -gene_region), by = c("dataset_id" = "dataset_id")) %>%
   select_if(~!all(is.na(.))) %>%  # This removes all the column that are empty
   select(-contains(c("dada2", "primer")), -dataset_path, -paper_doi, -paper_reference, -web_link, -asv_id)

# Merging together asvs with same hash value

  if (use_hash) {
    asv_set <- asv_set %>%
        mutate(asv_code = sequence_hash)
  }

# Remove the low abundance ASVs

  asv_set_number <- asv_set %>%
    count(asv_code, wt=n_reads, sort = TRUE) %>%
    filter(n >= sum_reads_min)

  asv_set <- asv_set %>%
    filter(asv_code %in% asv_set_number$asv_code)
  asv_fasta <- asv_fasta %>%
    filter(asv_code %in% asv_set_number$asv_code)

# Export fasta file

    if (taxonomy_full) {
      asv_fasta <- asv_fasta %>%
        rename(seq_name = asv_code)
      if (export_fasta) fasta_write(asv_fasta, str_c(directory, "metapr2_asv_set_", dataset_id_char ,"_", taxo_name_label, ".fasta"),
                                    compress = FALSE, taxo_include = TRUE)
    } else {
      asv_fasta <- asv_fasta %>%
         mutate(seq_name = str_c(asv_code, species, sep="|") )
      if (export_fasta) fasta_write(asv_fasta, str_c(directory, "metapr2_asv_set_", dataset_id_char ,"_", taxo_name_label, ".fasta"),
                                    compress = FALSE, taxo_include = FALSE)
    }

# Export wide Excel file

  if (export_wide_xls){
      if (use_hash) {
        # If use_hash, much remove the bootstrap values that may be different for the different asvs with same hash tag
          asv_set_wide <- asv_set %>%
            select(asv_code, kingdom:species, sequence, file_code, n_reads)
      } else {
          asv_set_wide <- asv_set %>%
            select(asv_code, kingdom:species_boot, sequence, sequence_hash, file_code, n_reads)
      }

    asv_set_wide <- asv_set_wide %>%
      pivot_wider(names_from=file_code, values_from = n_reads, values_fill=list(n_reads=0))

    openxlsx::write.xlsx(asv_set_wide, str_c(directory, "metapr2_wide_asv_set_", dataset_id_char ,"_", taxo_name_label, ".xlsx"))
  }

# Export long excel file

  if (export_long_xls)   openxlsx::write.xlsx(asv_set, str_c(directory, "metapr2_long_asv_set_",
                                                             dataset_id_char ,"_", taxo_name_label, ".xlsx"))

  if (export_sample_xls) openxlsx::write.xlsx(sample_list, str_c(directory, "metapr2_samples_asv_set_", dataset_id_char , ".xlsx"))

# Export Phyloseq file

 if (export_phyloseq) {
       ## Create the samples, otu and taxonomy tables

      # 1. samples table : row names are labelled by file_code
      samples_df <- asv_set %>%
        select(file_code:last_col(), -n_reads) %>%
        distinct(file_code, .keep_all = TRUE) %>%
        column_to_rownames(var = "file_code")


      # 2. otu table :
      otu <- asv_set %>%
        select(asv_code, file_code, n_reads) %>%
        pivot_wider(names_from=file_code, values_from = n_reads, values_fill=list(n_reads=0)) %>%
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

      saveRDS(phyloseq_asv, file = str_c(directory, "phyloseq_metapr2_asv_set_", dataset_id_char ,"_", taxo_name_label, ".rda") )

     }

# Return from function

   if (export_phyloseq) {
     return(list(df=asv_set, ps=phyloseq_asv, fasta=asv_fasta, samples=sample_list))
   } else{
     return(list(df=asv_set, fasta=asv_fasta, samples=sample_list))
  }
  }


# metapr2_export_datasets -------------------------------------------------------
#' @title Exports the metapr2 database
#' @description
#' Exports a data frame with all the available data sets
#' @return
#' A data frame
#' @examples
#' # Export all datasets
#'   datasets <- metapr2_export_datasets()
#' @export
#' @md
#'
metapr2_export_datasets <- function() {

# Read the database

  metapr2_db <- db_info("metapr2_google")
  metapr2_db_con <- db_connect(metapr2_db)
  metapr2_datasets <- tbl(metapr2_db_con, "metapr2_datasets") %>% collect()
  db_disconnect(metapr2_db_con)

  return(metapr2_datasets)
  }

# metapr2_export_reads_total -------------------------------------------------------
#' @title Exports the total number of reads in each dataset
#' @description
#' Exports  the total number of reads for each file_code
#'
#' Data are save in an excel file
#' @param directory  Directory where the files are saved (finish with /)
#' @param dataset_id_selected  Integer vector, e.g. 23 or c(21, 23) or 21:23
#' @return
#' A data frame with 2 columns:
#' * file_code
#' * reads_total
#' @examples
#' # Export data for all datasets
#'   df <- metapr2_export_reads_total()
#' @export
#' @md
#'
metapr2_export_reads_total <- function(directory = "C:/daniel.vaulot@gmail.com/Databases/_metaPR2/export/",
                               dataset_id_selected = c(1:100)) {

# Define a variable to hold the data set id as character
  dataset_id_char <- case_when ((100 %in% dataset_id_selected) ~ "all",
                                TRUE ~ str_c(dataset_id_selected, collapse = "_"))


# Read the database and already filter for selected records
  metapr2_db <- db_info("metapr2_google")
  metapr2_db_con <- db_connect(metapr2_db)

# Get the asv filtered by datasets
  asv_set <- tbl(metapr2_db_con, "metapr2_asv")%>%
     filter(dataset_id %in% dataset_id_selected) %>%
     collect()

# Need to use !! before the local variable
  # Error: Cannot embed a data frame in a SQL query.
  #  If you are seeing this error in code that used to work, the most likely cause is a change dbplyr 1.4.0. Previously `df$x` or
  # `df[[y]]` implied that `df` was a local variable, but now you must make that explict with `!!` or `local()`, e.g., `!!df$x` or
  # `local(df[["y"]))

  metapr2_asv_abundance <- tbl(metapr2_db_con, "metapr2_asv_abundance") %>%
    filter(asv_code %in% !!asv_set$asv_code) %>%
    collect()

  metapr2_samples <- tbl(metapr2_db_con, "metapr2_samples") %>% collect()

  db_disconnect(metapr2_db_con)


  asv_set <-  asv_set %>%
   left_join(metapr2_asv_abundance) %>%
   left_join(metapr2_samples)

  asv_set_total <- asv_set %>%
    filter(!is.na(file_code)) %>%  # Remove empty file codes
    group_by(file_code) %>%
    summarise(reads_total = sum(n_reads))


    openxlsx::write.xlsx(asv_set_total, str_c(directory, "metapr2_reads_total_", dataset_id_char, ".xlsx"))

    return(asv_set_total)

  }

