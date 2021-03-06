#' @import dplyr
#' @import dbplyr
#' @import stringr
#' @import readxl
#' @import readr
#' @import tibble
# Notes:
# biocLite("Biostrings")
# library(Biostrings) # To manipulate DNA sequences
# -- The problem of Biostrings is that it writes the file as Linux (LF only),
# -- so one solution is to use sequinr (writes LF/CR)
# -- the other solution is to replace line feed by linefeed + CR [this the chosen solution]
# library(seqinr)     # Used to write fasta files - NOT USED

#  To do:
#  * Use the dbplyr library to build queries...


# source("C:/Users/vaulot/Google Drive/Scripts/R library/dv_function_db.R")
# source("C:/Users/vaulot/Google Drive/Scripts/R library/dv_function_genbank.R")
# source("C:/Users/vaulot/Google Drive/Scripts/R library/dv_function_misc.R")

# pr2_read -------------------------------------------------------
#' @title Reads the whole PR2 database into a data frame
#' @description
#' The database access codes is obtained with the function db_info("pr2_google").
#'
#' **Important**: This gets the whole database, including sequences that have been removed.  These sequences must be filtered out before export using pr2  <- pr2 %>% filter (is.na(removed_version)). However the taxo entries that have been removed are already removed.
#' @param taxo_levels_number integer - number of taxonomic levels, 8 or 9
#' @return
#' A data frame with all the columns from pr2.
#' @examples
#' pr2 <- pr2_read(taxo_levels_number = 8)
#' @export
#' @md
#'
pr2_read <- function(taxo_levels_number = 9) {

  # This the level set not used (i.e. 8 if use 9 and conversely)
  taxo_levels_number_not_used = ifelse(taxo_levels_number == 8, 9, 8)

  # Ending for the taxonomic levels ("_8" or "_9")
  ending = str_c("_",taxo_levels_number) # Ending of taxonomic fields
  ending_removed = str_c("_",taxo_levels_number_not_used) # Ending of taxonomic fields

  taxo_levels<- list()
  taxo_levels[[8]] = c("kingdom", "supergroup", "division", "class", "order", "family", "genus", "species")
  taxo_levels[[9]] = c("domain", "supergroup", "division", "subdivision", "class", "order", "family", "genus", "species")

  print("Using pr2_google")

# Info for the Google dabase
  pr2_db <- db_info("pr2_google")

  pr2_db_con <- db_connect(pr2_db)

  pr2_main <- tbl(pr2_db_con, "pr2_main") %>%
    collect() %>%
    rename_with(~ str_replace(.,ending, ""), contains(ending)) %>%
    select(-contains(ending_removed))%>%
    select(-contains("edited")) %>%
    select(-contains("_id"))


  pr2_seq <- tbl(pr2_db_con, "pr2_sequences") %>%
    collect()

  pr2_taxo <- tbl(pr2_db_con, "pr2_taxonomy") %>%
    filter (is.na(taxo_removed_version)) %>%
    collect()

  pr2_taxo <- pr2_taxo %>%
    rename_with(~ str_replace(.,ending, ""), contains(ending)) %>%
    select(-contains(ending_removed)) %>%
    select(any_of(taxo_levels[[taxo_levels_number]])) %>%
    distinct()   # Very important because some levels are replicated because 2 species_8 may correspond to the same species_9



  pr2_metadata <- tbl(pr2_db_con, "pr2_metadata") %>%
    collect() %>%
    select(-contains("_id"))

  pr2_countries <- tbl(pr2_db_con, "pr2_countries") %>%
    collect()

  pr2_assign_bayes <- tbl(pr2_db_con, "pr2_assign_bayes") %>%
    collect()

  pr2_assign_silva <- tbl(pr2_db_con, "pr2_assign_silva") %>%
    collect() %>%
    select(pr2_accession, silva_taxonomy)

# Join the tables and keep only sequences that are not removed

  pr2 <- pr2_main %>%
    left_join(pr2_taxo, by = c("species"="species")) %>%
    left_join(pr2_seq) %>%
    left_join(pr2_metadata) %>%
    left_join(pr2_countries) %>%
    left_join(pr2_assign_bayes) %>%
    left_join(pr2_assign_silva) %>%
    relocate(!!taxo_levels[[taxo_levels_number]], .after = pr2_accession) %>%
    select_if(~!all(is.na(.)))  # remove empty columns



  db_disconnect(pr2_db_con)


  return(pr2)

}

# pr2_taxo_read -------------------------------------------------------
#' @title Reads the PR2 taxonomic database into a data frame
#' @description
#' The database access codes is obtained with the function db_info("pr2_google").
#'
#' **Important**: This gets the current taxonomic  database and does not include taxa that have been removed.
#' @return
#' A data frame with all the columns from pr2_taxo.
#' @examples
#' pr2_taxo <- pr2_taxo_read()
#' @export
#' @md
#'
pr2_taxo_read <- function() {

  print("Using pr2_google")


# Read the PR2 full tables including removed records
  pr2_db_con <- db_connect(db_info("pr2_google"))

  pr2_taxo <- tbl(pr2_db_con, "pr2_taxonomy") %>%
    filter (is.na(taxo_removed_version))%>%
    collect() %>%
    select(-contains("removed"))

  db_disconnect(pr2_db_con)

  return(pr2_taxo)

}


# pr2_sequence_label ------------------------------
#' @title Create a simple label for a sequence
#' @description
#' Create a label to be used as header a fasta file based on the metadata from the table *pr2_metadata*. Only apply to a single character not to a vector.
#' @param sample_type character - one of "culture", "environmental","isolate", "species"
#' @param strain character - strain name
#' @param clone character - clone name
#' @param specimen_voucher character - voucher
#' @return character
#' @examples
#' fasta_label <- pr2_sequence_label(sample_type="culture",strain="RCC115")
#' @export
#' @md
pr2_sequence_label <- function(sample_type, strain, clone, specimen_voucher) {

  # Note : must use str_replace_na() to make sure that no value is NA, because if a string is NA then str_c will return NA
    sequence_label <- case_when(is.na(sample_type) ~ "",
                      sample_type == "culture" ~ str_c("strain_",str_replace_na(strain, replacement="")),
                      sample_type == "environmental" ~ str_c("clone_",str_replace_na(clone, replacement="")),
                      sample_type == "isolate" ~ "",
                      sample_type == "specimen" ~ str_c("specimen_",str_replace_na(specimen_voucher, replacement="")),
                      TRUE ~ "")
    return(sequence_label)
}

# pr2_export -----------------------------------------
#' @title Export the PR2 database (one file)
#' @description
#' This will save the pr2 database in variety of format. The files are compressed as .gz
#' * "fasta" - the description line is taylored for different applications (usearch, mothur, dada2, blast). for the application dada2 2 kinds of files can be produced one for otu assignement and one for species
#' * "taxo" - the file is saved with the mothur format
#' * "metadata" - only the metadata
#' * "merged" - the whole database in a single file
#'
#' Notes
#' * The data are NOT filtered in any way, this should be done before saving
#' * The following fields are needed to create the sequence labels : pr2_sample_type, gb_strain, gb_clone, gb_specimen_voucher
#' @param pr2_select data frame - the pr2 database or an extract of the pr2 database
#' @param file_name character - full path of file where to save
#' @param file_type character - one of "fasta", "taxo", "metadata", "merged", "merged_excel"
#' @param file_format character - one of "UTAX", "fasta_taxo_short","fasta_taxo_long","dada2", "dada2_species","mothur"
#' @return
#' Write the files in compressed format (.gz)
#' @examples
#' pr2_export(pr2, "C:/Daniel/myfile.fas", "fasta", "mothur", 9)
#' @export
#' @md
pr2_export <- function(pr2_select, file_name, file_type="fasta", file_format="fasta_taxo_long", taxo_levels_number = 9 ) {
  # These lines are for testing and should be commented out
    # pr2 <- pr2_env_not_updated
    # file_name <- "pr2_env_not_updated.fasta"
    # level_select="kingdom"
    # taxo_select="Eukaryota"
    # file_type="fasta"
    # file_format="fasta_taxo_long"

  pr2_seq_out <- Biostrings::DNAStringSet(pr2_select$sequence)  # Store the sequence in a  DNAString - not used anymore

  # Use sequinr (sequences are stored as list)
  # pr2_seq_out <- as.list(pr2_select$sequence)

  # Remove any space from Strain, Clone name, Specimen_voucher
  pr2_select <- pr2_select %>%
     mutate (gb_strain = str_replace_all(gb_strain, " ", "_"),
             gb_clone = str_replace_all(gb_clone, " ", "_"),
             gb_specimen_voucher = str_replace_all(gb_specimen_voucher, " ", "_"),
              # Create a nice label based on the type of sequence
              # if the organelle is empty replace by empty string to avoid NA...
             sequence_label = pr2_sequence_label (pr2_sample_type,
                                                  gb_strain,
                                                  gb_clone,
                                                  gb_specimen_voucher),
             organelle = str_replace_na(organelle, replacement = ""),
             name_taxo_short = str_c(pr2_accession,
                              gene,
                              organelle,
                              sequence_label,
                              species,
                              sep="|"),
             name_dada2_species = str_c(genbank_accession,
                                      str_replace(species,"_"," "),
                                      sep=" ")
     )


    if (taxo_levels_number == 9) {
      pr2_select <- pr2_select %>%
          mutate ( name_UTAX = str_c(pr2_accession,";tax=k:",
                                  # kingdom,",d:",
                                    domain,",d:",
                                    supergroup,",p:",
                                    division,",c:",
                                    class,",o:",
                                    order,",f:",
                                    family,",g:",
                                    genus,",s:",
                                    species,
                                    sep=""),
                   name_taxo_long = str_c(pr2_accession,
                                    gene,
                                    organelle,
                                    sequence_label,
                                    domain,
                                    supergroup,
                                    division,
                                    subdivision,
                                    class,
                                    order,
                                    family,
                                    genus,
                                    species,
                                    sep="|"),
                   name_dada2 = str_c(domain,
                                    supergroup,
                                    division,
                                    subdivision,
                                    class,
                                    order,
                                    family,
                                    genus,
                                    species,
                                    "",
                                    sep=";")
                   )
    } else {
      pr2_select <- pr2_select %>%
          mutate ( name_UTAX = str_c(pr2_accession,";tax=k:",
                                    kingdom,",d:",
                                    supergroup,",p:",
                                    division,",c:",
                                    class,",o:",
                                    order,",f:",
                                    family,",g:",
                                    genus,",s:",
                                    species,
                                    sep=""),
                   name_taxo_long = str_c(pr2_accession,
                                    gene,
                                    organelle,
                                    sequence_label,
                                    kingdom,
                                    supergroup,
                                    division,
                                    class,
                                    order,
                                    family,
                                    genus,
                                    species,
                                    sep="|"),
                   name_dada2 = str_c(kingdom,
                                    supergroup,
                                    division,
                                    class,
                                    order,
                                    family,
                                    genus,
                                    species,
                                    "",
                                    sep=";")
                   )
  }

  switch(file_type,

    fasta = {

    switch(file_format,
         UTAX =  {names(pr2_seq_out) <- pr2_select$name_UTAX},
         fasta_taxo_short= {names(pr2_seq_out) <- pr2_select$name_taxo_short},
         fasta_taxo_long= {names(pr2_seq_out) <- pr2_select$name_taxo_long},
         dada2= {names(pr2_seq_out) <- pr2_select$name_dada2},
         dada2_species= {names(pr2_seq_out) <- pr2_select$name_dada2_species},
         mothur= {names(pr2_seq_out) <- pr2_select$pr2_accession
         }  )
  # When using DNA strings only LF added (Unix convention)
  # Set width to max possible so that no carriage return (causes problems with dada2)
    Biostrings::writeXStringSet(pr2_seq_out, file_name, compress=TRUE, width = 20000)
  # Convert to dos
    # unix2dos(file_name)
    # next line is when using seqinr
    # write.fasta(sequences=pr2_seq_out, names=names(pr2_seq_out), file.out=file_name, as.string = TRUE)
    },
    taxo = {
      # mothur taxonomy file terminates with semi column at end of line
      pr2_select <- pr2_select %>% mutate(taxo_string = str_c(domain, supergroup, division, subdivision,
                                                                       class, order, family, genus, species,"", sep=";")  )
      pr2_select_taxo <- pr2_select %>% select(pr2_accession, taxo_string)
      write.table(pr2_select_taxo, gzfile(file_name),
                  col.names = FALSE, row.names=FALSE,
                  sep="\t", quote= FALSE)

    },
    metadata = {
      pr2_select_metadata <- pr2_select %>% select(pr2_accession,
                                                   genbank_accession,
                                                   starts_with("gb_"),
                                                   starts_with("pr2_"),
                                                   starts_with("eukref_"),
                                                   starts_with("silva_"),
                                                   pubmed_id,
                                                   metadata_remark)
      # write.table(pr2_select_metadata, gzfile(file_name),
      #             col.names = TRUE, row.names=FALSE,
      #             sep="\t", quote= FALSE, na="")
      write_tsv(pr2_select_metadata, gzfile(file_name))
    },
    merged = {
      write.table(pr2_select, gzfile(file_name),
                  col.names = TRUE, row.names=FALSE,
                  sep="\t", quote= FALSE, na="")
    },
    merged_excel = {

      if (taxo_levels_number == 9) {
      pr2_select <- pr2_select %>%
        arrange(domain, supergroup, division, subdivision, class, order, family, genus, species)
      } else {
      pr2_select <- pr2_select %>%
        arrange(kingdom, supergroup, division, class, order, family, genus, species)
      }

      pr2_select %>% select(-starts_with("name")) %>%
        openxlsx::write.xlsx(file_name, zoom = 90, firstRow = TRUE, firstCol = TRUE)

    }
  )
}

# pr2_export SQLite -----------------------------------------
#' @title Export the PR2 database to a SQLite file
#' @description
#' This will save the pr2 database in SQLite databas as 4 tables (pr2_main, pr2_sequences, pr2_metadata, pr2_taxo)
#' @param file_name character - full path of file where to save
#' @return
#' TRUE is succesful
#' @examples
#' pr2_export_sqlite("pr2_google.sqlite")
#' @export
#' @md

pr2_export_sqlite <- function(file_name) {

# Read the PR2 full tables including removed records
  pr2_db_con <- db_connect(db_info("pr2_google"))
  pr2_main <- tbl(pr2_db_con, "pr2_main") %>%
    filter (is.na(removed_version))%>%
    collect()

  pr2_seq <- tbl(pr2_db_con, "pr2_sequences") %>%
    collect()

  pr2_taxo <- tbl(pr2_db_con, "pr2_taxonomy") %>%
    filter (is.na(taxo_removed_version))%>%
    collect()

  pr2_metadata <- tbl(pr2_db_con, "pr2_metadata")%>%
    collect()

  pr2_countries <- tbl(pr2_db_con, "pr2_countries")%>%
    collect()

  pr2_assign_silva <- tbl(pr2_db_con, "pr2_assign_silva")%>%
    collect()

  pr2_assign_bayes <- tbl(pr2_db_con, "pr2_assign_bayes")%>%
    collect()

  db_disconnect(pr2_db_con)

# Write to local database
  pr2_db_con <- db_connect_sqlite(file_name)

  copy_to(pr2_db_con, pr2_main, name = "pr2_main", indexes= list("pr2_accession", "genbank_accession", "species_9"),
          temporary=FALSE, overwrite = TRUE)
  copy_to(pr2_db_con, pr2_seq, name = "pr2_sequences", indexes= list("pr2_accession"),
          temporary=FALSE, overwrite = TRUE)
  copy_to(pr2_db_con, pr2_metadata, name = "pr2_metadata", indexes= list( "genbank_accession"),
          temporary=FALSE, overwrite = TRUE)
  copy_to(pr2_db_con, pr2_taxo, name = "pr2_taxonomy", indexes= list( "species_9"),
          temporary=FALSE, overwrite = TRUE)
  copy_to(pr2_db_con, pr2_countries, name = "pr2_countries", indexes= list( "pr2_country"),
          temporary=FALSE, overwrite = TRUE)
  copy_to(pr2_db_con, pr2_assign_silva, name = "pr2_assign_silva", indexes= list( "pr2_accession"),
          temporary=FALSE, overwrite = TRUE)
  copy_to(pr2_db_con, pr2_assign_bayes, name = "pr2_assign_bayes", indexes= list( "pr2_accession"),
          temporary=FALSE, overwrite = TRUE)

  db_disconnect(pr2_db_con)

  R.utils::gzip(file_name, overwrite = TRUE, remove=TRUE)

  return(TRUE)

}

# pr2_export_all -----------------------------------------
#' @title Export the PR2 database (all files)
#' @description
#' This will save the pr2 database in all the format necessary for the repository on GitHub. It calls \code{link{pr2_export}}
#' @param pr2.env list containing the name of files
#' @return
#' Write the all the files in compressed format (.gz)
#' @examples
#' pr2_export_all(pr2.env)
#' @export
pr2_export_all <- function(pr2.env) {

  taxo_levels_number = pr2.env$taxo_levels_number

  pr2 <- pr2_read(taxo_levels_number)
  pr2_taxonomy <- pr2_taxo_read()

  pr2_chimera <- pr2 %>%
    filter(!is.na(chimera)) %>%
    select_if(~!all(is.na(.)))  # remove empty columns

  pr2_unassigned <- pr2 %>%
    filter(is.na(removed_version)) %>%
    filter(is.na(species)) %>%
    select_if(~!all(is.na(.)))  # remove empty columns



# Filter out sequences that have been removed and sequences without species names
  pr2 <- pr2 %>%
    filter(is.na(removed_version)) %>%
    filter(!is.na(species)) %>%
    select(-contains("chimera")) %>%
    select(-contains("bayes")) %>%
    select(-contains("boot")) %>%
    select_if(~!all(is.na(.)))  # remove empty columns


# All sequences (merged file only)

  print("Exporting merged file")

  pr2_export(pr2, file_name=pr2.env$file_merged_excel, file_type="merged_excel")
  # pr2_export(pr2, file_name=pr2.env$file_metadata, file_type="metadata")

# 18S nuclear only

  print("Exporting 18S")

  pr2_18S <- filter(pr2, gene == "18S_rRNA")
  pr2_export(pr2_18S, file_name=pr2.env$file_18S_fasta_UTAX, file_type="fasta", file_format="UTAX", taxo_levels_number = taxo_levels_number )

  pr2_export(pr2_18S, file_name=pr2.env$file_18S_fasta_taxo_long, file_type="fasta", file_format="fasta_taxo_long", taxo_levels_number = taxo_levels_number)
  # pr2_export(pr2, file_name=pr2.env$file_fasta_taxo_short, file_type="fasta", file_format="fasta_taxo_short", taxo_levels_number = taxo_levels_number)

  pr2_export(pr2_18S, file_name=pr2.env$file_18S_fasta_mothur, file_type="fasta", file_format="mothur", taxo_levels_number = taxo_levels_number)
  pr2_export(pr2_18S, file_name=pr2.env$file_18S_taxo_mothur, file_type="taxo", taxo_levels_number = taxo_levels_number)

  pr2_export(pr2_18S, file_name=pr2.env$file_18S_fasta_dada2, file_type="fasta", file_format="dada2", taxo_levels_number = taxo_levels_number)

# 18S V4 - dada2 only

  # print("Exporting 18S V4")
  #
  # # 1. pcr in silico
  # # 2. remove amplicon that are NA
  # # 3. replace sequence by amplicon
  # pr2_18S_V4 <- data.frame(amplicon = pcr_sequences(pr2.env$primer_V4_fwd,
  #                                                   pr2.env$primer_V4_rev,
  #                                                   pr2_18S$sequence,
  #                                                   pr2.env$primer_V4_mismatches))  %>%
  #   bind_cols(pr2_18S) %>%
  #   filter(!is.na(amplicon)) %>%
  #   mutate(sequence=amplicon)
  #
  # pr2_export(pr2_18S_V4, file_name=pr2.env$file_18S_V4_fasta_dada2, file_type="fasta", file_format="dada2", taxo_levels_number = taxo_levels_number)

# 16S plastid only

  print("Exporting 16S plastid")

  pr2_16S <- filter(pr2, gene == "16S_rRNA")
  pr2_export(pr2_16S, file_name=pr2.env$file_16S_fasta_UTAX, file_type="fasta", file_format="UTAX", taxo_levels_number = taxo_levels_number)

  pr2_export(pr2_16S, file_name=pr2.env$file_16S_fasta_taxo_long, file_type="fasta", file_format="fasta_taxo_long", taxo_levels_number = taxo_levels_number)

  pr2_export(pr2_16S, file_name=pr2.env$file_16S_fasta_mothur, file_type="fasta", file_format="mothur", taxo_levels_number = taxo_levels_number)
  pr2_export(pr2_16S, file_name=pr2.env$file_16S_taxo_mothur, file_type="taxo", taxo_levels_number = taxo_levels_number)

  pr2_export(pr2_16S, file_name=pr2.env$file_16S_fasta_dada2, file_type="fasta", file_format="dada2", taxo_levels_number = taxo_levels_number)

# R package

  print("Exporting R package")


  pr2.env$pr2database_directory

  save(pr2,  file=str_c(pr2.env$pr2database_directory, "data/pr2.rda"))
  save(pr2_taxonomy,  file=str_c(pr2.env$pr2database_directory, "data/pr2_taxonomy.rda"))
  save(pr2_chimera,  file=str_c(pr2.env$pr2database_directory, "data/pr2_chimera.rda"))
  save(pr2_unassigned,  file=str_c(pr2.env$pr2database_directory, "data/pr2_unassigned.rda"))

# SQLite

  print("Exporting SQLite")
  pr2_export_sqlite(pr2.env$file_sqlite)

}

# pr2_taxo_check   ---------------
#' @title Check taxonomy
#' @description
#' This check the taxonomy for the following conditions:
#' * No taxon has 2 different parents
#' * Not taxon appears at different ranks
#'
#' @note
#' The input data frame does not need to be the whole pr2 taxo, can be taxo extracted from an Excel file.
#' Several files are produced:
#'
#' * *taxo_list_(rank_number)_(rank_name).txt* (7 files)
#'      * name of each taxon
#'      * taxa_number - number of species in this taxon
#'      * level of taxon - from 1 to 8 (1 for kingdom or domain)
#'      * id of the parent taxon (at rank i-1)
#'      * id of the taxon (from 1 to n for each rank)
#' * *taxo_list_2_parents_(rank_number)_(rank_name).txt* (7 files).  This is an important file because it list all the taxa that have 2 parents.  When all the files are completely empty, then the taxonomic list is OK.
#'      * name of the taxon
#'      * number of parents
#' * *taxo_level_duplicates.txt*. List the taxa that are found at different ranks.  Fro example Cryptophyta can be in the Class column for one species and in the Division levels for all other species.
#'      * name
#'      * number of times found
#'      * taxa number
#'      * level at which found
#'      * id
#'      * id of the parent
#' * *taxo_species_duplicates*. List of species that appear at least twice in the species list.
#'
#' @param pr2_taxo dataframe - should contain eight or nine columns kingdom -> species or domain -> species
#' @param taxo_levels - vector with taxo levels
#' @param dir_taxo character - full path where file are saved (without final //)
#' @return
#' Write several files (see description)
#' @examples
#' pr2_taxo_check (taxo_list, "C:/Daniel/database")
#' @export
#' @md
pr2_taxo_check <- function(pr2_taxo, taxo_levels = pr2.env$taxo_levels, dir_taxo="") {

  # Initialize variables

    taxo_levels_number <- length(taxo_levels)
    taxo_list <- list()
    taxo_list_2_parents <- list()

  # Extract all taxo levels

  # Level 1 - Kingdom or Domain
    i <- 1
    # print(taxo_levels[[1]])
    taxo <- pr2_taxo %>%
            group_by (!!as.name(taxo_levels[[1]])) %>%
            summarise(taxa_number=n()) %>%
            rename(name=!!as.name(taxo_levels[[1]]))
    taxo$level<-i
    taxo$id<-1:nrow(taxo)
    taxo$parent_id<-0
    taxo_list[[1]] <- taxo


  # Level 2 to 8 or 9

    for (i in c(2:taxo_levels_number)) {
      taxo_list[[i]] <-pr2_extract_one_taxo_level(pr2_taxo,
                                              taxo_levels[[i-1]],
                                              taxo_levels[[i]],
                                              i,
                                              taxo_list[[i-1]])
      write_tsv(taxo_list[[i]], path=str_c(dir_taxo,"/","taxo_list_", i, "_", taxo_levels[[i]], ".txt"))
    }

  #  Find taxo levels that have more than one parent
  # Note : it is necessary to examine each file one by one going in descending order and rerun the script each time

    for (i in c(2:taxo_levels_number)) {
      taxo_list_2_parents[[i]] <- taxo_list[[i]] %>%
                                  group_by (name) %>%
                                  summarize (n_parents=n()) %>%
                                  filter (n_parents>1)
      write_tsv(taxo_list_2_parents[[i]],
                  path=str_c(dir_taxo,"/","taxo_list_2_parents_", i, "_", taxo_levels[[i]], ".txt"))
    }

  #  Find taxo names that are found at different levels
     taxo_list_all <- taxo_list[[1]]
     for (i in c(2:taxo_levels_number)) {
      taxo_list_all <- rbind (taxo_list_all, taxo_list[[i]])
    }
     taxo_list_duplicate<- taxo_list_all %>%
                           group_by (name) %>%
                           summarize (n_found=n()) %>%
                           filter ((n_found>1))
     taxo_list_duplicate<-left_join(taxo_list_duplicate, taxo_list_all)
     write_tsv(taxo_list_duplicate, path=str_c(dir_taxo,"/","taxo_level_duplicates", ".txt"))

  #  Find species that are duplicated in the species list
     taxo_species_duplicate<- pr2_taxo %>%
                             group_by (species) %>%
                             summarize (n_found=n()) %>%
                             filter (n_found>1)
     write_tsv(taxo_species_duplicate, path=str_c(dir_taxo,"/","taxo_species_duplicates", ".txt"))

}


# pr2_extract_one_taxo_level -------------------------------------------------------------------------
#' Extract one level of taxonomy
#'
#' This function returns a data frame that contains:
#'
#' * name = name of the taxon
#' * taxa_number = number of taxa downstream of this taxon
#' * level = number of the level (e.g. genus = level number 7)
#' * parend_id = id of the first parent taxon  (sometimes one taxon has more than one taxon)
#' * id = id of the taxon (just the row name since they are ordered alphabetically )
#' @note
#' Here is an example at the genus level
#'   name	            taxa_number	    level	    parent_id	    id
#' Ascampbelliella	    1	              7	        51	          1
#' Choreotrichida_XX	  1	              7	        1	            2
#' Cyttarocylididae_X	1	              7	        52	          3
#' Cyttarocylis	      4	              7	        52	          4
#' Petalotricha	      1	              7	        52	          5
#' Codonaria	          2	              7	        53	          6
#' Codonella	          2	              7	        53	          7
#' @param pr2_taxo dataframe - taxonomy with 8 coumns kigdom -> species (but can be more)
#' @param level_above_name character - name of the taxo level above (e.g. "order")
#' @param level_name character - name the level considered (e.g. "family")
#' @param taxo_above dataframe -  the same structure than the output of the function but for the level above (e.g. the order if we are the family level)
#' @examples
#' genus_df <- pr2_extract_one_taxo_level(pr2_taxo, "family", "genus", 7, family_df)
#' @export
#' @md
  pr2_extract_one_taxo_level <- function(pr2_taxo, level_above_name, level_name, level_number, taxo_above) {
    taxo <-  pr2_taxo %>%
                group_by(!!as.symbol(level_above_name), !!as.symbol(level_name)) %>%
                summarise(n=n()) %>%
                transmute(name= !!as.symbol(level_name), number_of_taxa=n)
    taxo$level<-level_number
    taxo <-merge(taxo, taxo_above, all.x, by.x=level_above_name, by.y="name")
    taxo <- taxo %>%
      transmute(name, taxa_number=number_of_taxa, level=level.x, parent_id = id) %>%
      arrange(name) %>%
      rowid_to_column( var = "id")

    return(taxo)
  }

# pr2_taxo_list -----------------------------------------

#' @title Build a list of taxa
#' @description
#' Construct a list of taxa from the pr2 database
#' @param pr2 - data frame - whole or filtered pr2 data base
#' @param level_select - taxonomic rank
#' @return
#' A character vector
#' @examples
#' pr2_taxo_list(pr2, "class")
#' @export
#' @md
  pr2_taxo_list <- function(pr2, level_select) {

  taxo_list <- pr2 %>%
                distinct_(level_select) %>%
                arrange_(level_select)

  return(taxo_list)
}

# pr2_taxo_X -----------------------------------------

#' @title Add _X for taxon that are repeated on same line
#' @description
#' If more than 2 ranks on a given line have the same name, rename with Mytaxon_X, Mytaxon_XX etc.. up to the species name which becomes Mytaxon_XX_sp.
#' @param pr2_taxo Data frame containing columns from kingdom to species (or else see the next parameter).  It is not necessary that the columns be ordered.
#' @param taxo_levels A vector of ordered taxonomic ranks such as c("kingdom", "supergroup", "division", "class", "order", "family", "genus", "species")
#' @return
#' A corrected data frame
#' @examples
#' pr2_taxo_clean <- pr2_taxo_X(pr2_taxo, pr2.env$taxo_levels)
#' @export
#' @md

  pr2_taxo_X <- function(pr2_taxo, taxo_levels) {

  taxo_levels_number <- length(taxo_levels)

  for (i in 1:nrow(pr2_taxo)) {
    stem = ""
    for (j in 1:(taxo_levels_number-1) ){
      # Check first if the next level is equal to previous level which is stored in stem
      if (pr2_taxo[i,taxo_levels[j]]==stem) {
              # If true then add _XX where the number of X is proportional to the difference in number of columns since the start
              pr2_taxo[i,taxo_levels[j]] <-  str_c(stem, "_", str_dup("X", j-stem_level))
              # If we are at the end (species level) then use _XX_sp. there
              if (j == taxo_levels_number-1) {
                pr2_taxo[i,taxo_levels[j+1]] <-  str_c(stem,"_", str_dup("X", j-stem_level), "_sp.")
                }
            }
      else {
        # If the level is not repeated then use it to check the next column
        stem  <-  pr2_taxo[i,taxo_levels[j]]
        stem_level <- j
      }
    }
  }
  return(pr2_taxo)
}


# pr2_build_taxons   ---------------
#' @title Build taxon table (long format)
#' @description
#' This constructs a taxon table (long format) from a taxonomy table (wide format)
#' If add_old_fields = TRUE then the following columns are appended
#' old_species_id, taxon_edited_version, taxon_edited_by,  taxon_removed_version, taxon_remark,  taxon_reference
#'
#' @note
#' The input data frame does not need to be the whole pr2 taxo, can be taxo extracted from an Excel file.
#'
#' @param pr2_taxonomy dataframe - should contain eight columns kingdom -> species (this can be modified by using the pr2.env list)
#' @param pr2.env list containing the taxon levels
#' @param add_old_fields logical - If add_old_fields = TRUE then the following columns are appended: old_species_id, taxon_edited_version, taxon_edited_by,  taxon_removed_version, taxon_remark,  taxon_reference
#' @return
#' Data frame pr2_taxons with columns "taxon_id" , "taxon_name", "taxon_level", "parent_id", "old_species_id",
#' "taxon_edited_version", "taxon_edited_by", "taxon_removed_version" ,"taxon_remark", "taxon_reference"
#' @examples
#' pr2_taxons <- pr2_buid_taxons(pr2_taxonomy)
#' @export
#' @md
pr2_buid_taxons <- function(pr2_taxonomy, pr2.env, add_old_fields=FALSE) {

  # Initialize variables

    taxo_levels_number <- length(pr2.env$taxo_levels)
    taxo_list <- list()
    taxo_list_2_parents <- list()

  # Extract all taxo levels

  # Level 1 - Kingdom
    i <- 1
    taxo <- pr2_taxonomy %>%
            group_by (any_of(pr2.env$taxo_levels[1])) %>%
            summarise(number_of_taxa=n())
    taxo<-taxo %>%
      transmute(name=kingdom, taxa_number=number_of_taxa)
    taxo$level<-i
    taxo$id<-1:nrow(taxo)
    taxo$parent_id<-0
    taxo_list[[1]] <- taxo


  # Level 2 to 8

    for (i in c(2:taxo_levels_number)) {
      taxo_list[[i]] <-pr2_extract_one_taxo_level(pr2_taxonomy,
                                              pr2.env$taxo_levels[[i-1]],
                                              pr2.env$taxo_levels[[i]],
                                              i,
                                              taxo_list[[i-1]])
     }

  # Reduce the 8 lists to a single dataframe
    pr2_taxons <- taxo_list %>% reduce (bind_rows)

  # Get a serial number for the taxons
    pr2_taxons <- pr2_taxons %>% rowid_to_column(var="taxon_id")

  # Now replace the parent_id by the new id
    pr2_taxons_1 <- pr2_taxons %>% mutate(parent_level=level-1)

    pr2_taxons_2 <- pr2_taxons %>% select(parent_id_new=taxon_id, level, id)

    pr2_taxons_3 <- left_join(pr2_taxons_1, pr2_taxons_2, by = c("parent_level" = "level", "parent_id"="id"))

  # Change the names
    pr2_taxons_3 <- pr2_taxons_3 %>% select(taxon_id, taxon_name = name, taxon_level = level, parent_id = parent_id_new)

  # Left join to recover the old information from the pr2_taxonomy table

    if(add_old_fields==TRUE) {
    pr2_taxo_1 <- pr2_taxo %>% select(old_species_id = taxo_id,
                                      taxon_name=species,
                                      taxon_edited_version = taxo_edited_version,
                                      taxon_edited_by = taxo_edited_by,
                                      taxon_removed_version = taxo_removed_version,
                                      taxon_remark = taxo_remark,
                                      taxon_reference = reference)
    pr2_taxons_3 <- pr2_taxons_3 %>% left_join(pr2_taxo_1)
    }


    return(pr2_taxons_3)

}

# pr2_build_taxonomy   ---------------
#' @title Buikd taxonomy table (wide format) from taxon table (long format)
#' @description
#' This constructs a taxonomy table from a taxon table
#'
#' @note
#' The input data frame does not need to be the whole pr2 taxons, can be taxo extracted from an Excel file.
#'
#' Programming note:
#' This function uses the following trick for mutating the colums with the taxo_level string
#'     rank_name = **sym**(pr2.env$taxo_levels[i])
#'     mutate (**!!** rank_name **:=** taxon_name)
#' @param pr2_taxons dataframe - should contain at least three columns taxon_id, taxon_name, taxon_level, parent_id
#' @param pr2.env list containing taxo_levels
#' @return
#' Data frame pr2_taxonomy with 9 columns
#' @examples
#' pr2_taxo <- pr2_buil_taxonomy(pr2_taxons, pr2.env)
#' @export
#' @md
pr2_build_taxonomy <- function(pr2_taxons, pr2.env) {


  taxo_levels_number <- length(pr2.env$taxo_levels)
   taxo_list <- list()
   for (i in c(taxo_levels_number:1)){
      # i=7
    # Needs to transform the taxo_level string into an expression
       rank_name = sym(pr2.env$taxo_levels[i])

    # Select one level of taxonomy and keep its ID and its Parent ID
       taxo_list[[i]] <- pr2_taxons %>%
         filter(taxon_level == i) %>%
         dplyr::select(taxon_id, taxon_name, parent_id) %>%
         mutate (!! rank_name := taxon_name) %>%
         select(-taxon_name)

    # For the last level (species) just take all the species. For the other levels do a join with next level (e.g. if genus, join to species)
       if (i == taxo_levels_number) {
         pr2_taxonomy <-  taxo_list[[i]]
       } else {
         pr2_taxonomy <- pr2_taxonomy %>%
           left_join(taxo_list[[i]], by = c("parent_id"="taxon_id") ) %>%
           select(-parent_id) %>%
           rename(parent_id = parent_id.y)
       }
   }

    pr2_taxonomy <- pr2_taxonomy %>%
      select_( .dots = c("taxon_id", pr2.env$taxo_levels))

    return(pr2_taxonomy)

}

