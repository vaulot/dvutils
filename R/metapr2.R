#' @import dplyr
#' @import tidyr
#' @import stringr
#' @import tibble
#' @importFrom rio export

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
#' @param gene if left empty "" then the whole sequence is exported, if "LSU" or "SSU" look in barrnap table to extract seqence
#' @param taxo_level The taxonomic level for selection (do not quote), e.g. class or genus
#' @param taxo_name  The name of the taxonomic level selected, e.g. "Chlorophyta", can be a vector c("Chlorophyta", "Haptophyta")
#' @param boot_level The taxonomic level for bootstrap filtering (do not quote), e.g. class_boot or genus_boot
#' @param boot_min  Minimum bootstrap value at the class level, 0 if you want to get all asvs
#' @param assigned_with Program used for assignement - "dada2" or "decipher"
#' @param reference_database Reference database used - "pr2_4.14.0" or "pr2_4.12.0"
#' @param directory  Directory where the files are saved (finish with /)
#' @param dataset_id_selected  Integer vector, e.g. 23 or c(21, 23) or 21:23
#' @param filter_samples Character expression for filter, e.g. "DNA_RNA == 'DNA'" (use single quotes inside double quotes)
#' @param filter_metadata  Character expression for filter, e.g. "substrate == 'sediment trap material'" (use single quotes inside double quotes)
#' @param need_accession  If TRUE, only samples with NCBI_SRA_accession not NULL are exported (this can be NA if not available)
#' @param export_long_xls  If TRUE, an xlsx file is produced containing the final long data frame
#' @param export_wide_xls  If TRUE, an xlsx file is produced containing the final wide data frame
#' @param export_sample_xls  If TRUE, an xlsx file is produced containing the sample list
#' @param export_phyloseq  If TRUE, a phyloseq file is produced and a phyloseq object producted
#' @param export_fasta  If TRUE, a fasta is produced
#' @param export_fasta_sum_reads  If TRUE, sum of reads of ASV is added to FASTA file following the VSEARCH format
#' @param taxonomy_full If TRUE, the fasta file contains the full taxonomy (8 levels), if false only contains the species
#' @param use_hash If TRUE, the asvs with identical has will be merged and called by their hash value (sequence_hash)
#' @param sum_reads_min This is the minimum number of reads (summed over the datasets selected) for an asv to be included
#' @param sample_reads_min This is the minimum number of reads for a sample to be included
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
metapr2_export_asv <- function(gene="",
                               taxo_level = domain,
                               taxo_name="Eukaryota",
                               boot_level = class_boot,
                               boot_min = 0,
                               assigned_with = "dada2",
                               reference_database = "pr2_5.1.0",
                               directory = "C:/daniel.vaulot@gmail.com/Databases/_metaPR2/export/",
                               dataset_id_selected = c(1:500),
                               need_accession = TRUE,
                               filter_samples = NULL,
                               filter_metadata = NULL,
                               export_long_xls=FALSE,
                               export_wide_xls=FALSE,
                               export_sample_xls=FALSE,
                               export_phyloseq = FALSE,
                               export_fasta=FALSE,
                               export_fasta_sum_reads=FALSE,
                               taxonomy_full = TRUE,
                               use_hash = FALSE,
                               sum_reads_min = 0,
                               sample_reads_min = 1000) {

  taxo_levels <- c("domain", "supergroup", "division", "subdivision", "class", "order", "family", "genus", "species")
  taxo_levels_boot <- str_c(taxo_levels, "_boot")

  # Use correct values for gene
  if (!(gene %in% c("", "SSU", "LSU", "5S"))) gene = ""

  # For SSU and LSU, do not use the hash values for asv_code
  if (gene !="") use_hash = FALSE

  # Cannot export full taxonomy when adding number of reads
  if (export_fasta_sum_reads) taxonomy_full = FALSE


  # Define a variable to hold the data set id as character
  dataset_id_char <- case_when ((500 %in% dataset_id_selected) ~ "all",
                                length(dataset_id_selected) > 5 ~ "several",
                                TRUE ~ str_c(dataset_id_selected, collapse = "_"))


  # Read the database and already filter for selected records -------------------
  metapr2_db <- db_info("metapr2_google")
  metapr2_db_con <- db_connect(metapr2_db)

  # Get all the asv (filtration is now done latter)

  asv_set <- tbl(metapr2_db_con, "metapr2_asv") %>%
    filter(dataset_id %in% dataset_id_selected) %>%
    filter(is.na(chimera)) %>%
    collect()

  # Use the new assignements ---------------------------------------------------

  # Assign default values if errors
  if(!(assigned_with %in% c("dada2", "decipher"))) {assigned_with = "dada2"}

   table_assignment = str_c("metapr2_asv", assigned_with, reference_database, sep= '_' )
   asv_set_updated <- tbl(metapr2_db_con, table_assignment) %>%
        collect()


    # Merge the assignments

    asv_set <- asv_set %>%
      select(-any_of(c(taxo_levels, taxo_levels_boot))) %>%
      left_join(asv_set_updated, by = c("sequence_hash" = "sequence_hash"))


  # Create a label for the taxons selected to be used for the files -------------------
  if(length(taxo_name) > 1){
    taxo_name_label <- str_c(str_sub(taxo_name, 1, 5), collapse="_")
  } else {
    taxo_name_label <- taxo_name
  }

  cat("asv_set done\n")

  # Filter the asv based on the datasets selected and the minimum bootstrap -------------------
  # asv_set <- asv_set %>%
  #    filter(dataset_id %in% dataset_id_selected)

  # Need to use !! before the local variable
  # Error: Cannot embed a data frame in a SQL query.
  #  If you are seeing this error in code that used to work, the most likely cause is a change dbplyr 1.4.0. Previously `df$x` or
  # `df[[y]]` implied that `df` was a local variable, but now you must make that explict with `!!` or `local()`, e.g., `!!df$x` or
  # `local(df[["y"]))


  # ASV abundance ---

  metapr2_asv_abundance <- tbl(metapr2_db_con, "metapr2_asv_abundance") %>%
    filter(asv_code %in% !!asv_set$asv_code) %>%
    collect()

  cat("asv_abundance done\n")

  # Samples ---

  metapr2_samples <- tbl(metapr2_db_con, "metapr2_samples") %>%
    filter(!is.na(file_code)) %>%   # Remove all samples that have not been processed
    collect()

  # Filter samples using the custom filter
  if (!is.null(filter_samples)) {
    metapr2_samples <- metapr2_samples %>%
      filter(eval(rlang::parse_expr(filter_samples)))
  }

  # Filter samples that do not have Accession numbers
  if (need_accession) {
    metapr2_samples <- metapr2_samples %>%
      filter(!is.na(NCBI_SRA_accession))
  }

  # Filter samples that do not have enough reads
  file_codes_enough_reads <-  metapr2_asv_abundance %>%
    count(file_code, wt = n_reads) %>%
    filter(n >= sample_reads_min) %>%
    pull(file_code)

  metapr2_samples <- metapr2_samples %>%
    filter(file_code %in% file_codes_enough_reads)

  cat("samples done\n")

  # Metadata ---

  metapr2_metadata <- tbl(metapr2_db_con, "metapr2_metadata") %>%
    collect()

  # Filter metadata
  if (!is.null(filter_metadata)) {
    metapr2_metadata <- metapr2_metadata %>%
      filter(eval(rlang::parse_expr(filter_metadata)))
  }

  metapr2_datasets <- tbl(metapr2_db_con, "metapr2_datasets") %>%
    collect()

  metapr2_asv_barrnap <- tbl(metapr2_db_con, "metapr2_asv_barrnap") %>%
    collect()

  db_disconnect(metapr2_db_con)

  # Merging together asvs with same hash value -------------------

  if (use_hash){
    asv_fasta <- asv_set %>%
      group_by(sequence_hash) %>%
      dplyr::slice(1) %>%
      mutate(asv_code = sequence_hash) %>%
      select(-dataset_id) %>%
      ungroup()
  } else {
    asv_fasta <- asv_set
  }


  sample_list <- metapr2_samples %>%
    # Use a inner_join here so that if no metadata then sample is not found in sample list
    # 2024-05-24 - Reverse to left-join to include samples without metadata
    left_join(metapr2_metadata) %>%
    filter(dataset_id %in% dataset_id_selected) %>%
    left_join(metapr2_datasets, by = c("dataset_id" = "dataset_id")) %>%
    select(-contains(c("dada2", "primer", "web")), -dataset_path) %>%
    # This removes all the column that are empty
    select_if(~!all(is.na(.)))

  asv_set <-  asv_set %>%
    left_join(metapr2_asv_abundance) %>%
    # Use a inner_join here so that if no sample then sample is not found in asv_set
    # 2024-05-24 - Reverse to left-join to include samples without metadata
    left_join(metapr2_samples) %>%
    # Use a inner_join here so that if no metadata then sample is not found in asv_set
    # 2024-05-24 - Reverse to left-join to include samples without metadata
    left_join(metapr2_metadata) %>%
    left_join(metapr2_datasets, by = c("dataset_id" = "dataset_id")) %>%
    select(-contains(c("dada2", "primer", "paper", "web")), -dataset_path, -asv_id) %>%
    # Next line removes all the column that are empty
    select_if(~!all(is.na(.)))

  # Merging together asvs with same hash value -------------------

  if (use_hash) {
    asv_set <- asv_set %>%
      mutate(asv_code = sequence_hash)
  }

  # Compute abundance of ASVs --------------------------------------

  asv_set_number <- asv_set %>%
    count(asv_code, wt=n_reads, sort = TRUE, name = "sum_reads_asv") %>%
    filter(sum_reads_asv >= sum_reads_min)

  # Remove the low abundance ASVs ---------------------------------

  asv_set <- asv_set %>%
    filter(asv_code %in% asv_set_number$asv_code)

  asv_fasta <- asv_fasta %>%
    filter(asv_code %in% asv_set_number$asv_code) %>%
    left_join(asv_set_number) %>%
    arrange(desc(sum_reads_asv))

  # Compute the abundance of reads per sample for all and for photosynthetic groups ---

  sample_all_reads <- asv_set %>%
    count(file_code, wt=n_reads, name="reads_corrected_total")

  sample_photo_reads <- asv_set %>%
    filter((division %in% c(
      "Chlorophyta", "Cryptophyta", "Rhodophyta",
      "Haptophyta")
      ) |
      ((subdivision %in% "Gyrista") & str_detect(class, "phyceae")) |
      (class %in% "Chlorarachniophyceae")) %>%
    count(file_code, wt=n_reads, name="reads_corrected_photo")

  # Add corrected counts to asv_set and sample_list -------------------

  asv_set <- asv_set %>%
    left_join(sample_all_reads) %>%
    left_join(sample_photo_reads)

  sample_list <- sample_list %>%
    left_join(sample_all_reads) %>%
    left_join(sample_photo_reads) %>%
    relocate(contains("_corrected_"), .after = reads_total) %>%
    mutate(reads_corrected_photo = replace_na(reads_corrected_photo, 0)) %>%
    filter(!is.na(reads_corrected_total))


  # Filter the ASVs for taxonomic group and minimum bootstrap -------------------

  asv_set <- asv_set %>%
    filter({{taxo_level}} %in% taxo_name) %>%
    filter({{boot_level}} >= boot_min)

  asv_fasta <- asv_fasta %>%
    filter({{taxo_level}} %in% taxo_name) %>%
    filter({{boot_level}} >= boot_min)

  # File names -------------------
  generate_file_name <- function (type, ext = "xlsx") str_c(directory, type, dataset_id_char ,
                                                   "_", taxo_name_label ,
                                                   "_", if_else(gene!="", str_c(gene, "_"),""),
                                                   Sys.Date(),".", ext)


  file_asv_fasta <- generate_file_name("metapr2_asv_set_",  "fasta")
  file_asv_xlsx <-  generate_file_name("metapr2_asv_set_")
  file_wide <- generate_file_name("metapr2_wide_asv_set_")
  file_long <- generate_file_name("metapr2_long_asv_set_")
  file_samples <- generate_file_name("metapr2_samples_asv_set_")
  file_phyloseq <- generate_file_name("metapr2_phyloseq_asv_set_", "rds")


  # Export fasta file ------------

  if (taxonomy_full) {
    asv_fasta <- asv_fasta %>%
      dplyr::mutate(seq_name = asv_code)
  } else {
    asv_fasta <- asv_fasta %>%
      mutate(seq_name = str_c(asv_code, species, sep="|") )
  }

  if (export_fasta_sum_reads) {
    asv_fasta <- asv_fasta %>%
      mutate(seq_name = str_c(seq_name,';size=',sum_reads_asv))
  }

  if (gene != "") {
    metapr2_asv_barrnap <- metapr2_asv_barrnap |>
      filter(gene == gene) |>
      slice_head(by = asv_code) |>
      select(asv_code, start, end)

    asv_fasta <- left_join(asv_fasta, metapr2_asv_barrnap) |>
      filter(!is.na(start)) |>
      mutate(sequence =  str_sub(sequence, start, end))
  }

  # Export asv file
  if (export_fasta){
    fasta_write(asv_fasta,file_asv_fasta, compress = FALSE, taxo_include = taxonomy_full)
    rio::export(asv_fasta, file_asv_xlsx, overwrite = TRUE)
  }


  # Export wide Excel file

  if (export_wide_xls){
    if (use_hash) {
      # If use_hash, much remove the bootstrap values that may be different for the different asvs with same hash tag
      asv_set_wide <- asv_set %>%
        select(asv_code, domain:species, sequence, file_code, n_reads)
    } else {
      asv_set_wide <- asv_set %>%
        select(asv_code, domain:species_boot, sequence, sequence_hash, file_code, n_reads)
    }

    asv_set_wide <- asv_set_wide %>%
      pivot_wider(names_from=file_code,
                  values_from = n_reads,
                  values_fill=list(n_reads=0),
                  values_fn = mean)

    rio::export(asv_set_wide, file_wide, overwrite = TRUE)
  }

  # Export long excel file

  if (export_long_xls)   rio::export(asv_set, file_long, overwrite = TRUE)

  if (export_sample_xls) rio::export(sample_list, file_samples, overwrite = TRUE)

  # Export Phyloseq file

  if (export_phyloseq) {
    ## Create the samples, otu and taxonomy tables

    # 1. samples table : row names are labeled by file_code
    samples_df <- asv_set %>%
      select(dataset_id, file_code:last_col(), -n_reads) %>%
      distinct(file_code, .keep_all = TRUE) %>%
      column_to_rownames(var = "file_code")


    # 2. otu table :
    otu <- asv_set %>%
      select(asv_code, file_code, n_reads) %>%
      pivot_wider(names_from=file_code,
                  values_from = n_reads,
                  values_fill=list(n_reads=0),
                  values_fn = mean) %>%
      column_to_rownames(var = "asv_code")

    # 3. Taxonomy table

    tax <-  asv_set %>%
      select(asv_code, domain:species) %>%
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

    saveRDS(phyloseq_asv, file = file_phyloseq )

  }

  # Return from function

  if (export_phyloseq) {
    return(list(df=asv_set, ps=phyloseq_asv, fasta=asv_fasta, samples=sample_list))
  } else{
    return(list(df=asv_set, fasta=asv_fasta, samples=sample_list))
  }
}


# metapr2_export_datasets -------------------------------------------------------
#' @title Exports the metapr2_datasets table
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

# metapr2_export_metadata -------------------------------------------------------
#' @title Exports the metapr2_metadata table
#' @description
#' Exports a data frame with all the available metadata
#' @return
#' A data frame
#' @examples
#' # Export all datasets
#'   metadata <- metapr2_export_metadata()
#' @export
#' @md
#'
metapr2_export_metadata <- function() {

# Read the database

  metapr2_db <- db_info("metapr2_google")
  metapr2_db_con <- db_connect(metapr2_db)
  metapr2_metadata <- tbl(metapr2_db_con, "metapr2_metadata") %>% collect()
  db_disconnect(metapr2_db_con)

  return(metapr2_metadata)
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
                               dataset_id_selected = c(1:500)) {

# Define a variable to hold the data set id as character
  dataset_id_char <- case_when ((500 %in% dataset_id_selected) ~ "all",
                                length(dataset_id_selected) > 5 ~ "several",
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
   inner_join(metapr2_asv_abundance) %>%
   inner_join(metapr2_samples)

  asv_set_total <- asv_set %>%
    filter(!is.na(file_code)) %>%  # Remove empty file codes
    group_by(file_code) %>%
    summarise(reads_total = sum(n_reads))


    export(asv_set_total, str_c(directory, "metapr2_reads_total_", dataset_id_char,"_", Sys.Date() , ".xlsx"))

    return(asv_set_total)

  }

# metapr2_export_qs -------------------------------------------------------
#' @title Exports the metapr2 database in qs file for shiny app
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
#' @param set_type   "public" (41 sets), "basic" (5 sets only), "all"
#' @param directory  Directory where the files are saved (finish with /)
#' @param do_cluster Should the ASVs be clustered ?
#'
#' @return
#' TRUE if successful
#' @examples
#' # Export public data set
#'
#'   metapr2_export_qs("public")
#'
#' @export
#' @md
#'
metapr2_export_qs <- function(set_type = "public 3.0",
                              do_cluster = FALSE,
                              directory = "data/") {

# Constants  ----------------------------------


  # List to store all global variables
  global <- list()

  # global$taxo_levels <- c("kingdom", "supergroup", "division", "class", "order", "family", "genus", "species", "asv_code") - PR2 version 4.14.0
  global$taxo_levels <- c("domain", "supergroup", "division",  "subdivision", "class", "order", "family", "genus", "species", "asv_code")
  global$traits <- c("mixoplankton", "HAB_species", "ecological_function")

  # All samples are normalized to 100 with 3 decimals, so that it corresponds to a percent
  global$n_reads_tot_normalized = 100

  # Column to removes as well as contains("_boot")

  cols_to_remove = c("asv_id" , "chimera", "sequence_hash", "asv_remark",
                     "sample_id", "file_name", "metadata_id",
                     "NCBI", "sample_name",
                     "sample_concentration","sample_sorting", "sample_sorting_cells",
                     "metadata_code", "replicate",
                     "fraction_name_original", "fraction_min", "fraction_max",
                     "sample_remark", "metadata_code_original",
                     "station_id_num", "year", "time", "O2", "fluorescence",
                     "season", "day_length", "depth_range",
                     "substrate_description", "substrate_description_detailed",
                     "experiment_time" ,
                     "pH" , "ice_type" , "ice_thickness" ,
                     "Secchi_depth",
                     "Chla_0.2_3 um" , "PAR_pct" ,
                     "bact_ml" , "peuk_ml" , "neuk_ml", "crypto_ml",
                     "NO2", "metadata_remark",
                     "dataset_path", "sample_number",
                     "lat_min","lat_max","long_min","long_max",
                     "primers_removed ", "dada2_truncLen","dada2_minLen",
                     "dada2_maxLen","dada2_bigdata","dada2_truncQ",
                     "dada2_maxEE","dada2_max_number_asvs")



  need_accession = FALSE

  # Data sets obsolete - NOT USED  ----------------------------------
  # For class - 5 sets

  if (set_type == "basic"){
    datasets_selected <- metapr2_export_datasets() %>%
      filter(dataset_id %in% c(1, 34, 35, 205, 206))
    do_cluster <- false
    need_accession = TRUE
    # sub_dir = "sets_basic"
  }

  # 41 sets
  if (set_type == "public 1.0"){
    datasets_selected <- metapr2_export_datasets() %>%
      filter(metapr2_version == "1.0")
    cluster_version = "1.0"
    need_accession = TRUE
    # sub_dir = "sets_public"
  }

  # 59 sets
  if (set_type == "public 2.0"){
    datasets_selected <- metapr2_export_datasets() %>%
      filter(metapr2_version %in% c("1.0", "2.0"))
    cluster_version = "2.0"
    need_accession = TRUE
    # sub_dir = "sets_public"
  }


  # 41 sets + green edge
  if (set_type == "green-edge"){
    datasets_selected <- metapr2_export_datasets() %>%
      filter(metapr2_version == "1.0" | dataset_id %in% 21:23)
    cluster_version = "1.0+GE"
    # sub_dir = "sets_public"
  }

  # 59 sets + green edge
  if (set_type == "green-edge2"){
    datasets_selected <- metapr2_export_datasets() %>%
      filter(metapr2_version %in% c("1.0", "2.0") | dataset_id %in% 21:23)
    cluster_version = "2.0+GE"
    # sub_dir = "sets_public"
  }


  # Data sets selected ----------------------------------

  # Nansen legacy cruise
  if (set_type == "nansen"){
    datasets_selected <- metapr2_export_datasets() %>%
      filter(dataset_id %in% c(398))
    do_cluster <- FALSE
    # sub_dir = "sets_basic"
  }

  # PacBio + Illumina TS series + Mahwash dataset
  if (set_type == "pacbio"){
    datasets_selected <- metapr2_export_datasets() %>%
      filter(dataset_id %in% c(393:395))
    do_cluster <- FALSE
    # sub_dir = "sets_basic"
  }

  # Version 3.0 - xx sets
  if (set_type == "public 3.0"){
    datasets_selected <- metapr2_export_datasets() %>%
      filter(metapr2_version %in% c("1.0", "2.0", "3.0"))
    cluster_version = "3.0"
    need_accession = TRUE
    # sub_dir = "sets_public"
  }


  # All - xx sets
  if (set_type == "all"){
    datasets_selected <- metapr2_export_datasets() %>%
      filter(!is.null(metapr2_version),
             processing_pipeline_metapr2 == "dada2",   # To remove Tara Swarm
             gene == "18S rRNA")                       # To remove 16S plastid
    cluster_version = "3.0 all"
    need_accession = FALSE
    # sub_dir = "sets_all"
  }




  cat("Datasets: ", nrow(datasets_selected))

  #
  # # All datasets
  # datasets_selected <- metapr2_export_datasets()

  # Read metaPR2 datasets ----------------------------------


  asv_set <- metapr2_export_asv(
    taxo_level = domain,
    taxo_name = "Eukaryota",
    dataset_id_selected = datasets_selected$dataset_id,
    filter_samples = "fraction_name != 'femto' &
                   !is.na(file_code)  &
                   file_code !='' &
                   !is.na(DNA_RNA)",
    filter_metadata =  "is.na(experiment_name)",
    need_accession = need_accession,
    export_fasta = FALSE,
    export_wide_xls = FALSE,
    export_sample_xls = FALSE,
    export_phyloseq = FALSE,
    directory = directory,
    taxonomy_full = TRUE,
    boot_min = 75,
    boot_level = supergroup_boot,
    use_hash = TRUE,
    sum_reads_min = 100,
    sample_reads_min = 1000
  )

  # Shorten asv_code -------------------

  asv_set$df <- asv_set$df %>%
    mutate(asv_code = str_sub(asv_code, 1,10) )

  asv_set$fasta <- asv_set$fasta %>%
    mutate(asv_code = str_sub(asv_code, 1,10))


  # Cluster -------------------

  if (do_cluster) {
    metapr2_db <- dvutils::db_info("metapr2_google")
    metapr2_db_con <- dvutils::db_connect(metapr2_db)

    clusters <- tbl(metapr2_db_con, "metapr2_asv_clusters") %>%
      collect()

    dvutils::db_disconnect(metapr2_db_con)

    clusters <- clusters %>%
      filter(record_type == "H",
             pct_sim == 100,
             metapr2_version == cluster_version) %>%
      rename(asv_code = hash_value,
             asv_code_centroid = hash_value_centroid) %>%
      mutate(asv_code = str_sub(asv_code,1,10),
             asv_code_centroid = str_sub(asv_code_centroid,1,10)) %>%
      select(-record_type, -pct_sim, -metapr2_version)

    # Check whether any duplicate ASVs
    clusters_duplicated <- clusters %>%
      count(asv_code) %>%
      filter(n> 1)

    cat("Duplicate ASVs in cluster:", nrow(clusters_duplicated), "\n")

    # In the fasta table, remove all ASVs that have a centroid
    # Need also to summarize the number of reads per cluster

    asv_set$fasta_cluster <- asv_set$fasta %>%
      left_join(clusters) %>%
      mutate(asv_code = case_when(!is.na(asv_code_centroid )~ asv_code_centroid,
                                  TRUE ~ asv_code))

    fasta_cluster_sum_reads_asv <- asv_set$fasta_cluster %>%
      group_by(asv_code) %>%
      summarize(sum_reads_asv = sum(sum_reads_asv)) %>%
      ungroup()

    asv_set$fasta_cluster <- asv_set$fasta_cluster %>%
      select(-sum_reads_asv) %>%
      left_join(fasta_cluster_sum_reads_asv) %>%
      filter(is.na(asv_code_centroid)) %>%
      select(-asv_code_centroid)

    cat("Are the 2 sums identical (df_fasta): ",
        sum(asv_set$df_fasta$sum_reads_asv) == sum(asv_set$df_df$sum_reads_asv),
        "\n")

    # In the df table, rename ASVs that have a centroid with the code of the centroid

    asv_set$df_cluster <- asv_set$df %>%
      left_join(clusters) %>%
      mutate(asv_code = case_when(!is.na(asv_code_centroid )~ asv_code_centroid,
                                  TRUE ~ asv_code))  %>%
      select(-asv_code_centroid) %>%
      group_by(dataset_id, file_code, asv_code) %>%
      summarize(n_reads = sum(n_reads)) %>%
      ungroup()

    # Check that no centroid left in cluster


    cat("Are the 2 sums identical (df_cluster): ",
        sum(asv_set$df_cluster$n_reads),
        "Are the 2 sums identical (df): ",
        sum(asv_set$df$n_reads), "\n")

    cat(glue::glue("Number of ASVs before clustering (df): {length(unique(asv_set$df$asv_code))}"), "\n")
    cat(glue::glue("Number of ASVs after clustering (df): {length(unique(asv_set$df_cluster$asv_code))}"), "\n")
    cat(glue::glue("Number of ASVs before clustering (fasta): {length(unique(asv_set$fasta$asv_code))}"), "\n")
    cat(glue::glue("Number of ASVs after clustering (fasta): {length(unique(asv_set$fasta_cluster$asv_code))}"), "\n")

    # Finalize by replacing the df and fasta tables by the one with clusters

    asv_set$df <- asv_set$df_cluster

    asv_set$fasta <- asv_set$fasta_cluster

    asv_set$df_cluster <- NULL
    asv_set$fasta_cluster <- NULL
  }


  # Summarize information for each data set ----------------------------------

  dataset_samples <- asv_set$df %>%
    select(dataset_id, file_code) %>%
    distinct() %>%
    group_by(dataset_id) %>%
    count(name = "sample_number") %>%
    ungroup()

  dataset_asv <- asv_set$df %>%
    select(dataset_id, asv_code) %>%
    distinct() %>%
    group_by(dataset_id) %>%
    count(name = "asv_number") %>%
    ungroup()

  dataset_reads <- asv_set$df %>%
    select(dataset_id, file_code, n_reads) %>%
    group_by(dataset_id, file_code) %>%
    summarize(n_reads = sum(n_reads, na.rm=TRUE)) %>%
    ungroup() %>%
    group_by(dataset_id) %>%
    summarize(n_reads_mean = round(mean(n_reads,  na.rm=TRUE),0)) %>%
    ungroup()


  # Normalize total mumber of sample reads to 100 ------------------------------
  # Only use the first 10 characters for asv_code (non ambiguous)

  asv_set$df <- asv_set$df %>%
    select(file_code, asv_code, n_reads) %>%
    group_by(file_code) %>%
    mutate(n_reads_pct = round(n_reads/sum(n_reads)*global$n_reads_tot_normalized, 3)) %>%
    ungroup()

  # Do not transform the phyloseq data for Alpha and Beta diversity analyses
  # asv_set$ps  = phyloseq::transform_sample_counts(asv_set$ps,
  #                                                 function(x) round( x / sum(x)*global$n_reads_tot_normalized,3 ))

  asv_set$samples <- asv_set$samples %>%
    select(-any_of(cols_to_remove),
           -(processing_pipeline_original:contact_email)
    )

  asv_set$datasets <- datasets_selected %>%
    select(-any_of(cols_to_remove)) %>%
    left_join(dataset_samples) %>%
    left_join(dataset_asv) %>%
    left_join(dataset_reads)

  cat("Datasets: ", nrow(asv_set$datasets))

  pr2_traits <- dvutils::pr2_traits_merge(trait_types = global$traits)  %>%
    select(any_of(c(global$taxo_levels, global$traits)))

  asv_set$fasta <- asv_set$fasta %>%
    select(- seq_name, -any_of(cols_to_remove)) %>% # Remove , -contains("_boot")
    left_join(pr2_traits)






  # The next line is only necessary for exporting as independant phyloseq file
  # asv_set_phyloseq <- asv_set$ps

  # Compute derived data -----------------------------------------------------------

  # Metadata groups

  global$gene_regions <- sort(unique(asv_set$samples$gene_region))
  global$DNA_RNAs <-  sort(unique(asv_set$samples$DNA_RNA))
  global$projects <-  sort(unique(asv_set$samples$project))

  # Reordering the variables ------------------------------------------------


  depth_level_ordered <- c("under ice", "surface", "euphotic",
                           "mesopelagic", "bathypelagic",
                           "composite", "bottom")
  fraction_name_ordered <- c("pico", "pico-nano", "nano",
                             "nano-micro", "micro",
                             "meso" ,"total" )
  substrate_ordered <- c("water", "ice", "sediment trap material", "sediment trap blank",
                         "epibiota", "tissue", "coral tissue", "coral skeleton",
                         "biofilm", "sediment", "sand","soil" )

  ecosystems_ordered <- c( "oceanic", "coastal","estuarine","freshwater lakes",
                           "freshwater rivers","terrestrial")

  asv_set$samples <- asv_set$samples %>%
    mutate(substrate = stringr::str_replace(substrate, "first year ice", "ice"),
           depth_level = forcats::fct_relevel(depth_level, depth_level_ordered),
           fraction_name = forcats::fct_relevel(fraction_name, fraction_name_ordered),
           substrate = forcats::fct_relevel(substrate, substrate_ordered),
           ecosystem = forcats::fct_relevel(ecosystem, ecosystems_ordered)
    )


  # Reorder for the check boxes ---------------------------------------------


  update_order <- function(variable) {
    values <- dplyr::arrange(asv_set$samples, .data[[variable]])
    values <- dplyr::pull(values, .data[[variable]]) %>%
    unique() %>%
    as.character()
  }

  global$depth_levels <- update_order("depth_level")
  global$fraction_names <- update_order("fraction_name")
  global$substrates <- update_order("substrate")
  global$ecosystems <- update_order("ecosystem")


  # Taxonomy structure --------------------------------------------------------

  global$pr2_taxo <- asv_set$fasta %>%
    select(any_of(c(global$taxo_levels, global$traits))) %>%
    distinct() %>%
    arrange(across(any_of(global$taxo_levels)))

  # Colors ------------------------------------------------------------------

  # Supergroups ---

  pr2_colors <- dvutils::pr2_colors_read()

  supergroup_colors <- pr2_colors %>%
    filter(palette_name=="daniel",
           taxo_level == "supergroup")

  global$supergroup_colors <- structure(supergroup_colors$color_hex,
                                        .Names=supergroup_colors$taxo_name)

  # Ecological function ---

  ecological_function_colors <- pr2_colors %>%
    filter(palette_name=="daniel",
           taxo_level == "ecological_function")

  global$ecological_function_colors <- structure(ecological_function_colors$color_name,
                                           .Names=ecological_function_colors$taxo_name)

  trophic_group_colors <- pr2_colors %>%
    filter(palette_name=="daniel",
           taxo_level == "trophic_group")

  global$trophic_group_colors <- structure(trophic_group_colors$color_name,
                                                 .Names=trophic_group_colors$taxo_name)


  # Authentification (move to data_initialize) ------------------------


  # Save data using qs --------------------------------------------------------

  if (do_cluster){
      qs::qsave(asv_set, str_c(directory, "asv_set_cluster.qs"))
  }
  else {
      qs::qsave(asv_set, str_c(directory, "asv_set.qs"))
  }

  qs::qsave(global, str_c(directory,"global.qs"))

  # Object sizes ------------------------------------------------------------

  obj_size <- function(x) {
    cat("Object:",deparse(substitute(x)), "- size: ", round(pryr::object_size(x)/10**6, 2), " Mb \n")
  }

  obj_size(asv_set)
  obj_size(asv_set$df)
  obj_size(asv_set$samples)
  obj_size(asv_set$fasta)

  obj_size(global)
}

