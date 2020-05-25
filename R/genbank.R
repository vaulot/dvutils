#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import lubridate
#' @import stringr
#' @import fs
#' @import genbankr
#' @import rentrez
#' @import reutils



# Notes :
#   Genbank taxonomy can be downloaded from : ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/

# genbank_search ---------------------------------------------
#'@title Search GenBank
#'@description
#'Return sequences information from an ENTREZ search query
#'@param query Character string for ENTREZ search
#'@param db Character string for database to be used
#'@param seq_max Maximum number of sequences to obtain
#'@return
#'Data frame with the following columns
#' * genbank_accession
#' * gb_definition
#' * gb_date
#' * gb_organism
#'@examples
#'genbank_search(query = "28S[TITL] AND rRNA[TITL] AND Chlorophyta[ORGN]", seq_max = 500)

#'@export
#'@md

genbank_search <-function(query, db="nuccore", seq_max = 2000) {

  # The next line is to prevent curl errors : https://github.com/ropensci/rentrez/issues/127
  httr::set_config(httr::config(http_version = 0))

  # Use API key to go much faster
  set_entrez_key("97ce6407215b5d1b6f5ee3ce8a6703793608")

  GB_entrez <- rentrez::entrez_search(db=db, term=query, use_history=TRUE)
  print(GB_entrez)

  seq_max = min(GB_entrez$count, seq_max)
  seq_step = 50
  seq_info <- list()
  for( seq_start in seq(1,seq_max, seq_step)){
     recs <- entrez_summary(db="nuccore", web_history=GB_entrez$web_history, retmax=seq_step, retstart=seq_start)
     genbank_accession <- extract_from_esummary(recs, "caption")
     gb_definition <- extract_from_esummary(recs, "title")
     gb_date <- extract_from_esummary(recs, "createdate")
     gb_organism <- extract_from_esummary(recs, "organism")
     gb_info_fields <- extract_from_esummary(recs, "subtype")
     gb_info_data <- extract_from_esummary(recs, "subname")
     seq_info[[as.character(seq_start)]] <- data.frame(genbank_accession, gb_definition,
                                                       gb_date, gb_organism,
                                                       gb_info_fields, gb_info_data)
     cat(seq_start+seq_step-1, "sequences obtained\r")
  }

  seq_df <- purrr::reduce(seq_info, dplyr::bind_rows)
  return(seq_df)
}


# genbank_download ---------------------------------------------
#'@title Download sequences from GenBank
#'@description
#'Write GenBank files from an vector of accession numbers
#'@param accession Character vector of accession numbers
#'@param directory Character, directory name which must end by "/"
#'@return
#'TRUE if successful
#'@examples
#'genbank_download(c("JX015376", "JQ768406", "LT621940"), "genbank/")

#'@export
#'@md

genbank_download <-function(accession,directory) {

  # The next line is to prevent curl errors : https://github.com/ropensci/rentrez/issues/127
  httr::set_config(httr::config(http_version = 0))

  # Use API key to go much faster
  set_entrez_key("97ce6407215b5d1b6f5ee3ce8a6703793608")

  # Create the directory if does not exist (directory must end by "/")
  if (!dir_exists(directory))
    { dir_create(directory, recursive = TRUE)}

  for (i in 1:length(accession)) {
      print(accession[i])
      GB_entry <- rentrez::entrez_search(db="nuccore", term=str_c(accession[i],"[ACCN]"))

	    # if the Accession has been found write the file
      if (length(GB_entry$ids) > 0) {
        GB_seq <- rentrez::entrez_fetch(db="nuccore", id=GB_entry$ids, rettype="gb")
		    GB_file <- str_c(directory,accession[i],".gb")
        write(GB_seq, file = GB_file )
		    print(str_c("file written to: ", GB_file))
      }
  }
  return(TRUE)
}

# genbank_download_parse  ---------------------------------------------
#'@title Download and parse sequences from GenBank
#'@description
#'If the GenBank file does not exist it fetches and writes it.
#'If the file exists it just reads it. It extracts the metadata and return into a a data frame.
#'
#'Note: use only sequence_keep = TRUE for simple gene sequences, else it may crash
#'@param accession Character vector of accession numbers
#'@param directory Character, directory name which must end by "/"
#'@param sequence_keep Logical, if FALSE the sequence is not returned in the file data frame
#'@param store_file Logical, if FALSE the file is not stored and the sequence is not returned in the file data frame
#'@return
#'Data frame with the metadata information.
#'@examples
#'gb_metadata <- genbank_download_parse(accession = c("JX015376", "JQ768406", "LT621940"), directory = "genbank/", sequence_keep=TRUE)

#'@export
#'@md

genbank_download_parse <-function(accession,directory, sequence_keep=TRUE, store_file=FALSE) {

  # The next line is necessary to read the date correctly
  Sys.setlocale("LC_TIME", "C")

  # The next line is to prevent curl errors : https://github.com/ropensci/rentrez/issues/127
  httr::set_config(httr::config(http_version = 0))

  # Use API key to go much faster
  rentrez::set_entrez_key("97ce6407215b5d1b6f5ee3ce8a6703793608")

  # Create the directory if does not exist (directory must end by "/")
  if (!fs::dir_exists(directory))
    { fs::dir_create(directory, recursive = TRUE)}

  # Create empty dataframe
  metadata <- data.frame()
  metadata_list <- list()

  for (i in 1:length(accession)) {

      print(str_c("i = ",i, " - ",accession[i]))
      GB_file <- str_c(directory,accession[i],".gb")

      GB_entrez <- rentrez::entrez_search(db="nuccore", term=str_c(accession[i],"[ACCN]"))
      print(GB_entrez)

      if(length(GB_entrez$ids) == 0) {
        metadata_one_row <- data.frame (
  		      genbank_accession = accession[i],
  		      gb_definition = "Genbank entry does not exist anymore")
      }

      # The following line checks whether the GenBank file exists - if it is there no need to rewrite
      # and go to GenBank to get it

      if (!(fs::file_exists(GB_file))) {
            # if the Accession has been found write the file and read the metadata
            if (length(GB_entrez$ids) > 0) {
                GB_seq <- rentrez::entrez_fetch(db="nuccore", id=GB_entrez$ids, rettype="gb")
                write(GB_seq, file = GB_file )
		            print(str_c("file written to: ", GB_file))
            }
      }
		  # There are some cases where the accession number either
      #  - does not exist anymore in GenBank (old entry deleted)
      #  - not yet (new submission not published yet)
      if ((fs::file_exists(GB_file))) {

		    GB_entry <- tryCatch({genbankr::readGenBank(GB_file, partial=TRUE)},
		                         error=function(e) NA)  # Catch error in read the GenBank file...
		    if (typeof(GB_entry)=="S4") {
  		    GB_meta <- sources(GB_entry) # Could use also GB_entry@sources$isolate for example
  		    # print(GB_meta)
  		    # print(GB_entry@definition)
  		    # print(GB_entry@locus)

  		    if (sequence_keep) {
            one_sequence = as.character(getSeq(GB_entry))
          } else {
            one_sequence = NA
          }

          metadata_one_row <- data.frame (
  		      genbank_accession = accession[i],
  		      gb_definition = GB_entry@definition,
  		      gb_organism = ifelse(is.null(GB_meta$organism), NA, GB_meta$organism),
  		      gb_organelle = ifelse(is.null(GB_meta$organelle), NA, GB_meta$organelle),
  		      gb_strain = ifelse(is.null(GB_meta$strain), NA, GB_meta$strain),
  		      gb_isolate = ifelse(is.null(GB_meta$isolate), NA, GB_meta$isolate),
  		      gb_clone = ifelse(is.null(GB_meta$clone), NA, GB_meta$clone),
  		      gb_specimen_voucher = ifelse(is.null(GB_meta$specimen_voucher), NA, GB_meta$specimen_voucher),
  		      gb_collected_by = ifelse(is.null(GB_meta$collected_by), NA, GB_meta$collected_by),
  		      gb_lat_lon = ifelse(is.null(GB_meta$lat_lon), NA, GB_meta$lat_lon),
  		      gb_note = ifelse(is.null(GB_meta$note), NA, GB_meta$note),
  		      gb_culture_collection = ifelse(is.null(GB_meta$culture_collection), NA, GB_meta$culture_collection),
  		      gb_isolation_source = ifelse(is.null(GB_meta$isolation_source), NA, GB_meta$isolation_source),
  		      gb_host = ifelse(is.null(GB_meta$host), NA, GB_meta$host),
  		      gb_environmental_sample = ifelse(is.null(GB_meta$environmental_sample), NA, GB_meta$environmental_sample),
  		      gb_collection_date = ifelse(is.null(GB_meta$collection_date), NA, GB_meta$collection_date),
  		      gb_country = ifelse(is.null(GB_meta$country), NA, GB_meta$country),
  		      gb_date = lubridate::as_date(str_sub(GB_entry@locus, start =-11),format="%d-%B-%Y", tz = "UTC"),
  		      gb_locus = str_sub(GB_entry@locus, start =-15, end=-13),
  		      sequence = one_sequence
          )

		    } else {
		       # If genbankr cannot read the file, then use rentrez which has some limitations (e.g. no locus)
           # cat("GB_entrez$ids : ", GB_entrez$ids) # DEBUG

		       gb_summary <- entrez_summary(db="nuccore", id=GB_entrez$ids)

		       # cat("gb_summary$subname : ", gb_summary$subname) # DEBUG
           # cat("gb_summary$subtype : ", gb_summary$subtype) # DEBUG

		       GB_meta <- as.list(unlist(str_split(gb_summary$subname, pattern="[|]")))
           names(GB_meta) <- as.list(unlist(str_split(gb_summary$subtype, pattern="[|]")))
           print(str_c("date : ", gb_summary$createdate))

		       metadata_one_row <- data.frame (
  		      genbank_accession = accession[i],
  		      gb_definition = gb_summary$title,
  		      gb_organism = gb_summary$organism,
  		      gb_strain = gb_summary$strain,
  		      gb_organelle = ifelse(is.null(GB_meta$organelle), NA, GB_meta$organelle),
  		      gb_isolate = ifelse(is.null(GB_meta$isolate), NA, GB_meta$isolate),
  		      gb_clone = ifelse(is.null(GB_meta$clone), NA, GB_meta$clone),
  		      gb_specimen_voucher = ifelse(is.null(GB_meta$specimen_voucher), NA, GB_meta$specimen_voucher),
  		      gb_collected_by = ifelse(is.null(GB_meta$collected_by), NA, GB_meta$collected_by),
  		      gb_lat_lon = ifelse(is.null(GB_meta$lat_lon), NA, GB_meta$lat_lon),
  		      gb_note = ifelse(is.null(GB_meta$note), NA, GB_meta$note),
  		      gb_culture_collection = ifelse(is.null(GB_meta$culture_collection), NA, GB_meta$culture_collection),
  		      gb_isolation_source = ifelse(is.null(GB_meta$isolation_source), NA, GB_meta$isolation_source),
  		      gb_host = ifelse(is.null(GB_meta$host), NA, GB_meta$host),
  		      gb_environmental_sample = ifelse(is.null(GB_meta$environmental_sample), NA, GB_meta$environmental_sample),
  		      gb_collection_date = ifelse(is.null(GB_meta$collection_date), NA, GB_meta$collection_date),
  		      gb_country = ifelse(is.null(GB_meta$country), NA, GB_meta$country),
  		      gb_date = lubridate::as_date(gb_summary$createdate,format="%Y/%m/%d", tz = "UTC"),
  		      gb_locus = "",
  		      sequence = NA
  		      )
		      }


          if (!is.na(metadata_one_row$gb_lat_lon)) {
            # Separate the lat_lon field
            metadata_one_row <- metadata_one_row %>% separate(col=gb_lat_lon,
                                                              into=c("latitude","lat_NS", "longitude", "long_EW"),
                                                              sep=" ", remove = FALSE)

            # If the data are correcly formatted compute pr2_lat and long
            if (str_detect(metadata_one_row$lat_NS,"[NS]")){

            # compute the two fields
            metadata_one_row <- metadata_one_row %>%
              mutate(pr2_longitude=as.numeric(longitude)*ifelse((long_EW == "E"), 1, -1)) %>%
              mutate(pr2_latitude=as.numeric(latitude)*ifelse((lat_NS == "N"), 1, -1))
            }

            # Remove fields that are not used
            metadata_one_row <- metadata_one_row%>%
              select(-latitude, - longitude, -lat_NS, -long_EW)
          }

          metadata_one_row$pr2_sample_type <- NULL

          if (!is.na(metadata_one_row$gb_isolate))
            {metadata_one_row$pr2_sample_type <-"isolate" }
          else
            {if (!is.na(metadata_one_row$gb_strain)|!is.na(metadata_one_row$gb_culture_collection))
              { metadata_one_row$pr2_sample_type <-"culture" }
          }
          if (metadata_one_row$gb_locus == "ENV")
          {metadata_one_row$pr2_sample_type <- "environmental"}

          if (sequence_keep ==TRUE) metadata_one_row$sequence_length <- str_length(metadata_one_row$sequence)

  		    print(str_c("i = ",i, " - ",accession[i], " - gb file read"))
  		    # print(metadata)
      }
   metadata_list[[i]] <- metadata_one_row
   if (!store_file) file.remove(GB_file)

  }

  # Only combine rows if they are not empty
  if(length(metadata_list)>0){
    metadata <- purrr::reduce(metadata_list, bind_rows)
  }
  return(metadata)
}

# genbank_download_parse_rentrez  ---------------------------------------------
#'@title Download and parse sequences from GenBank using rentrez only
#'@description
#'Does not write the GenBank file.
#'Note: use only sequence_keep = TRUE for simple gene sequences, else it may crash
#'@param accession Character vector of accession numbers
#'@param sequence_keep Logical, if FALSE the sequence is not returned in the file data frame
#'@return
#'Data frame with the metadata information.
#'@examples
#'gb_metadata <- genbank_download_parse(accession = c("JX015376", "JQ768406", "LT621940"), sequence_keep=TRUE)

#'@export
#'@md
genbank_download_parse_rentrez <-function(accession, sequence_keep=TRUE) {

  # accession <- c("AH006017","AH006986") # for testing

  # The next line is necessary to read the date correctly
  Sys.setlocale("LC_TIME", "C")

  # The next line is to prevent curl errors : https://github.com/ropensci/rentrez/issues/127
  httr::set_config(httr::config(http_version = 0))

  # Use API key to go much faster
  rentrez::set_entrez_key("97ce6407215b5d1b6f5ee3ce8a6703793608")

  acc_max = length(accession)
  acc_step = 50

  metadata_list <- list()

  for (i in seq(1, acc_max , acc_step)) {

    # i = 1 # for testing
    i_max = min(i+acc_step-1, acc_max)

    query = str_c(str_c(accession[i:i_max], "[ACCN]"), collapse = "|")
    GB_entrez <- rentrez::entrez_search(db="nuccore", term=query, use_history=TRUE)
    print(GB_entrez)

    if (sequence_keep){
      sequences <- entrez_fetch(db="nuccore", web_history=GB_entrez$web_history, rettype="fasta")
      tmp = tempfile()
      writeLines(sequences, tmp)
      sequences <- readDNAStringSet(tmp)
      sequences <- data.frame(gb_sequence=sequences) %>%
                  rownames_to_column(var = "seq_name") %>%
                  tidyr::separate(col = seq_name, into = "genbank_accession", sep = "[. ]", extra = "drop" )

    }

    recs <- entrez_summary(db="nuccore", web_history=GB_entrez$web_history)
    genbank_accession <- extract_from_esummary(recs, "caption")
    gb_definition <- extract_from_esummary(recs, "title")
    gb_date <- extract_from_esummary(recs, "createdate")
    gb_organism <- extract_from_esummary(recs, "organism")
    gb_info_fields <- extract_from_esummary(recs, "subtype")
    gb_info_data <- extract_from_esummary(recs, "subname")

    for (j in 1:length(recs)) {
        # j = 1 # for testing

        GB_meta <- as.list(unlist(str_split(gb_info_data[j], pattern="[|]")))
        names(GB_meta) <- as.list(unlist(str_split(gb_info_fields[j], pattern="[|]")))

        metadata_one_row <- data.frame (
        		      genbank_accession = genbank_accession[j],
        		      gb_definition = gb_definition[j],
        		      gb_organism = gb_organism[j],
        		      gb_strain =  ifelse(is.null(GB_meta$strain), NA, GB_meta$strain),
        		      gb_organelle = ifelse(is.null(GB_meta$organelle), NA, GB_meta$organelle),
        		      gb_isolate = ifelse(is.null(GB_meta$isolate), NA, GB_meta$isolate),
        		      gb_clone = ifelse(is.null(GB_meta$clone), NA, GB_meta$clone),
        		      gb_specimen_voucher = ifelse(is.null(GB_meta$specimen_voucher), NA, GB_meta$specimen_voucher),
        		      gb_collected_by = ifelse(is.null(GB_meta$collected_by), NA, GB_meta$collected_by),
        		      gb_lat_lon = ifelse(is.null(GB_meta$lat_lon), NA, GB_meta$lat_lon),
        		      gb_note = ifelse(is.null(GB_meta$note), NA, GB_meta$note),
        		      gb_culture_collection = ifelse(is.null(GB_meta$culture_collection), NA, GB_meta$culture_collection),
        		      gb_isolation_source = ifelse(is.null(GB_meta$isolation_source), NA, GB_meta$isolation_source),
        		      gb_host = ifelse(is.null(GB_meta$host), NA, GB_meta$host),
        		      gb_environmental_sample = ifelse(is.null(GB_meta$environmental_sample), NA, GB_meta$environmental_sample),
        		      gb_collection_date = ifelse(is.null(GB_meta$collection_date), NA, GB_meta$collection_date),
        		      gb_country = ifelse(is.null(GB_meta$country), NA, GB_meta$country),
        		      gb_date = lubridate::as_date(gb_date[j],format="%Y/%m/%d", tz = "UTC"),
        		      gb_locus = ""
        		      ) %>%
        tidyr::separate(col=gb_lat_lon,
                        into=c("latitude","lat_NS", "longitude", "long_EW"),
                        sep=" ", remove = FALSE, extra = "drop") %>%
        mutate(pr2_longitude = case_when(str_detect(long_EW,"[EW]") ~ as.numeric(longitude)*ifelse((long_EW == "E"), 1, -1)) ) %>%
        mutate(pr2_latitude= case_when (str_detect(lat_NS,"[NS]") ~ as.numeric(latitude)*ifelse((lat_NS == "N"), 1, -1)  ) ) %>%
        select(-latitude, - longitude, -lat_NS, -long_EW) %>%
        mutate(pr2_sample_type = case_when(!is.na(gb_isolate) ~ "isolate",
                                           !is.na(gb_strain)|!is.na(gb_culture_collection) ~ "culture",
                                           gb_locus == "ENV" ~ "environmental"))


        if (sequence_keep) {
           metadata_one_row <- left_join( metadata_one_row, sequences) %>%
          mutate(sequence_length = str_length(gb_sequence))
        }


       metadata_list[[i+j-1]] <-  metadata_one_row
    }

   cat(length(recs), " sequences obtained\n")

  }

  # Only combine rows if they are not empty
  if(length(metadata_list)>0){
    metadata <- purrr::reduce(metadata_list, bind_rows)
  }

  return(metadata)
}

# genbank_features  ---------------------------------------------
#'@title Read features of sequences from GenBank
#'@description
#'Return a dataframe containing all features for a vector of accession numbers. The GenBank files must exist.
#'@param accession Character vector of accession numbers
#'@param directory Character, directory name which must end by "/"
#'@return
#'Data frame with the features information.
#'@examples
#'gb_features <- genbank_features(c("JX015376", "JQ768406", "LT621940"), "genbank/")

#'@export
#'@md

genbank_features <-function(accession,directory) {

  # Create empty dataframe
  features <- data.frame()

  for (i in 1:length(accession)) {
      print(str_c("i = ",i, " - ",accession[i]))
      GB_file <- str_c(directory,accession[i],".gb")
      if ((file_exists(GB_file))) {
		    GB_entry <- genbankr::readGenBank(GB_file, partial=TRUE)
		    if (length(otherFeatures(GB_entry)) > 0 ) {
        features_one_accession <- as.data.frame (otherFeatures(GB_entry))
        features_one_accession$genbank_accession <- accession[i]
		    features <- bind_rows(features, features_one_accession)
        }
      }
  }
  return(features)
}


# genbank_field  ---------------------------------------------------------

#'@title Read a single field from a set of existing GenBank files
#'@description
#'Return a vector containing the value of a GenBank field for a vector of accession numbers. The GenBank files must exist
#'@param accession Character vector of accession numbers
#'@param directory Character, directory name which must end by "/"
#'@param field_name must be unquoted !!!
#'@return
#'Vector with the field information.
#'@examples
#'gb_organism <- genbank_field(c("JX015376", "JQ768406", "LT621940"), "genbank/", organism)

#'@export
#'@md

genbank_field <-function(accession, directory, field_name) {

  field_value <- character()
  for (i in 1:length(accession)) {
    print(str_c("i = ",i, " - ",accession[i]))
	  GB_file <- str_c(directory,accession[i],".gb")
	  field_value[i] <- NA
    if (file_exists(GB_file)) {
        GB_entry <- genbankr::readGenBank(GB_file)                                    # Get the GenBank file
		    GB_meta <- sources(GB_entry)                                        # Store the metadata associated with the
		    field_value[i] <- ifelse(is.null(eval(substitute(GB_meta$field_name))), NA, eval(substitute(GB_meta$field_name)))
	     }
  }
      return(field_value)
}

# genbank_locus ---------------------------------------------------------

#'@title Read locus from a set of existing GenBank files
#'@description
#'Return a vector containing the value of a GenBank locus for a vector of accession numbers. The GenBank files must exist
#'@param accession Character vector of accession numbers
#'@param directory Character, directory name which must end by "/"
#'@return
#'Vector with the locus information.
#'@examples
#'gb_locus <- genbank_locus(c("JX015376", "JQ768406", "LT621940"), "genbank/")

#'@export
#'@md

genbank_locus <-function(accession, directory) {

  field_value <- character()

  for (i in 1:length(accession)) {
    print(str_c("i = ",i, " - ",accession[i]))
	  GB_file <- str_c(directory,accession[i],".gb")
	  field_value[i] <-NA
    if (file_exists(GB_file)) {
        GB_entry <- genbankr::readGenBank(GB_file)                                    # Get the GenBank file
		    GB_meta <- sources(GB_entry)                                        # Store the metadata associated with the locus
		    field_value[i] <- str_sub(GB_entry@locus, start =-15, end=-12)

	     }
  }
      return(field_value)
}

# genbank_taxonomy ---------------------------------------------------------
#'@title Read NCBI taxo
#'@description
#'Reads directly from NCBI information taxo_lineage from set of taxo_id.
#'
#'Note : the full NCBI taxonomy has been uploaded in the PR2 database, so this function is not really necessary.
#'@param taxo_id  A vector of genbank taxonomy ids
#'@return
#'A data frame with 2 columns.
#'@examples
#'gb_tax <- genbank_taxonomy(c("1230134", "1230316", "1905175"))
#'@export
#'@md

genbank_taxonomy <-function(taxo_id) {

  taxo_lineage <- character()

  for (i in 1:length(taxo_id)) {
      taxo_lineage[i] <- NA
    # The commented line is not necessary since we have the id already
    # If necessary we will the use GB_entry$ids in the following fetch
	    # GB_entry <- entrez_search(db="taxonomy", str_c(taxo_id[i],[taxID]"))

    # Directly getting the entry since we have the taxonomy
      GB_tax <- rentrez::entrez_fetch(db="taxonomy", id=taxo_id[i], rettype="xml")

    # entrez_fetch returns an XML string from which we extract the lineage
      taxo_lineage[i] <- str_replace(str_extract(GB_tax,"<Lineage>[a-zA-Z;\\s]+"),"<Lineage>","")
  }
  # Binds the 2 columns in a dataframe
    taxo <- cbind(taxo_id, taxo_lineage)
    colnames(taxo) <- c("taxo_id", "taxo_lineage")
    return(taxo)

}
