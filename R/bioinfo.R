#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import stringr
#' @import Biostrings
#' @import dada2

# fasta_write : Write fasta file with taxo ----------------------------------------------

#' @title Write a fasta file with the taxonomy
#'
#' @description
#' Write a fasta file from a set of sequences
#' Option : add to the definition line the the taxonomy separated by separator character (e.g. |)
#'
#' >Otu0001|Alveolata|Dinophyta|Syndiniales|Dino-Group-I|Dino-Group-I-Clade-1|Dino-Group-I-Clade-1_X|Dino-Group-I-Clade-1_X_sp.
#'
#' AGCTCCAATAGCGTATATTAAAGTTGTTGCGGTTAAAAAGCTCGTAGTTGGA...
#' @param df The data frame with the otu names, the taxonomy and the sequences. It should have the following columns (with exactly these names)
#'
#'       * seq_name : the sequence name
#'       * supergroup: species
#'       * sequence
#' @param file_name Character, where to save the fasta file
#' @param compress If TRUE produces a gz file
#' @param taxo_include If TRUE then add taxo information which must be provided
#' @param taxo_separator Character used to separate the different taxonomic levels
#' TRUE if it terminates OK
#'
#' @examples
#' fasta_write(df,"otu_taxo.fasta", compress=FALSE, include_taxo=TRUE, taxo_separator=";")
#' @md
#' @export

fasta_write <- function(df,file_name, compress=FALSE, taxo_include=TRUE, taxo_separator="|") {

# First remove the gaps (can be - or .)
  df <-  df %>%  mutate(sequence = str_replace_all(sequence, "(-|\\.)",""))

  seq_out <- DNAStringSet(df$sequence)

  if (taxo_include==TRUE) {
      names(seq_out) <- str_c(df$seq_name,
                                df$supergroup,
                                df$division,
                                df$class,
                                df$order,
                                df$family,
                                df$genus,
                                df$species,
                                sep=taxo_separator)
  }
  else { names(seq_out) <- df$seq_name
  }

  writeXStringSet(seq_out, file_name, compress=compress, width = 20000)

  return(TRUE)
}

# fasta_read  : Read a fasta file into a data frame --------------------------------

#' @title Read a fasta file into a data frame
#'
#' @description
#' Read a fasta file from a set of sequences into a data frame
#' @param file_name Character, where to save the fasta file
#' @param compress If TRUE reads a gz file
#' @examples
#' df <- fasta_read(df,"otu_taxo.fasta", compress=FALSE)
#' @md
#' @export
fasta_read <- function(file_name, compress=FALSE) {

  # file_name <- "C:/daniel.vaulot@gmail.com/Scripts/R/dvutils/tests/testthat/mothur.database.fas"

  df <- Biostrings::readDNAStringSet(file_name)

  df <- data.frame(sequence=df) %>%
                rownames_to_column(var = "seq_name")
}


# fasta_filter : Filter a fasta file based on length an ambiguities ----------------------

#' @title Filter a fasta file
#'
#' @description
#' Read a fasta file and write back a fasta file keeping only sequences with less than max_ambig and between min and max length.
#' The file name is changed from xxx.fas to xxx.filtered.fas.
#' Prints also the number of sequences in and out.
#'
#' @param file_name Name the input fasta file
#' @param min_length Minimum length of sequences to keep
#' @param max_length Maximum length of sequences to keep
#' @param type "DNA" or "AA"
#' @param max_ambig Maximum number of ambiguites
#'
#' @return
#' TRUE if it terminates OK
#' @examples
#' fasta_filter("protein.faa", min_length=100, max_length=10000, type="AA", max_ambig=5)
#' @section To do:
#' Nothing =)
#' @md
#' @export

fasta_filter <- function(file_name, min_length=100, max_length=10000, type="AA", max_ambig=5){

# Read the file
  if (type=="AA") {seq_in <- readAAStringSet(file_name)}

# Transform to data frame
  df <- data.frame (names=names(seq_in),sequence=as.character(seq_in), length=nchar(seq_in) )
  print(str_c("Number of sequences initially: ", nrow(df)))

# Remove ambiguities
  if (type=="AA") {df <- df %>% filter(str_count(sequence,"X")<= max_ambig)}

# Filter based on lengths
  df <- df %>% filter((length >= min_length)& (length <= max_length))

# Construct the Biostring set
  seq_out <- AAStringSet(df$sequence)
  names(seq_out) <- df$names
  print(str_c("Number of sequences after filtration: ", nrow(df)))

# Write fasta file
  file_name_out <- str_c(fs::path_ext_remove(file_name),".filtered.",fs::path_ext(file_name))
  writeXStringSet(seq_out, file_name_out, compress=FALSE, width = max_length)
  return(TRUE)
}

# fastq_subsample : Function to sample fastq files ----------------------
#' @title Subsample fastq files
#'
#' @description
#' Subsample fastq files for the first n_seq sequences (! not a random sample)
#' @param fastq_path File directories
#' @param n_seq Number of sequences to subsample
#' @param random If TRUE, then sample a range around n_seq (gaussian with sd=0.1*nseq)
#'
#' @return
#' TRUE if it terminates OK
#' @examples
#' fastq_subsample("C:Mydata", n_seq=1000, random=FALSE)
#' @section To do:
#' Implement compressed files
#' @md
#' @export
#'
fastq_subsample <- function(fastq_path, n_seq=1000, random=FALSE) {

  # fastq_path <- "C:/Data Biomol/RNA/Tags/CEE 2014 Andres/fastq"
  # fastq_path <- "C:/Data Biomol/RNA/Tags/CARBOM/fastq/tutorial"

   fns <- sort(list.files(fastq_path, full.names = TRUE))
   fnsR1 <- fns[str_detect( basename(fns),"R1.fastq")]
   fnsR2 <- fns[str_detect( basename(fns),"R2.fastq")]

 for(i in 1:length(fnsR1)) {

   n_sampled = n_seq
   if (random) n_sampled = round(rnorm(1,n_seq,0.05*n_seq))

   print(fnsR1[i])

  # Process R1
    # Create connection to be able to close it after
     file_in <- FastqFile(fnsR1[i])

    # Sample the number of sequences needed
     sampler <- FastqStreamer(con=fnsR1[i], n=n_sampled)

    # Draw an instance of the sequences
     sample<-yield(sampler)

    # Write the Fastq file
     file_out <- str_replace(fnsR1[i],"R1", "R1.subsample")
     print(file_out)
     writeFastq(sample, file=file_out, compress=FALSE)
     close(file_in)

  # Process R2
     file_in <- FastqFile(fnsR2[i])
     sampler <- FastqStreamer(con=fnsR2[i], n=n_sampled)
     set.seed(123) # The seed need to be reset before each yield to have the matching pairs
     sample<-yield(sampler)
     file_out <- str_replace(fnsR2[i],"R2", "R2.subsample")
     print(file_out)
     writeFastq(sample, file=file_out, compress=FALSE)
     close(file_in)

  }

}

# XStringSet_to_df : Transforms a DNA string set into a data frame --------------------------------

#' @title Transforms a DNA or AA String set into a data frame
#' @description
#' Very simple helper function. The data frame has two columns seq_name and sequence
#' @param seq_set DNAStringSet or AA
#' @return
#' A data frame
#' @examples
#' df.fasta <- DNAStrinSet_to_df(DNAset_fasta)
#' @export
#'
XStringSet_to_df <- function(seq_set){

  df <- data.frame(seq_set) %>%
                rownames_to_column(var = "seq_name") %>%
                rename(sequence=seq_set)
}

# kmer : Get all kmers in a sequence --------------------------------

#' @title Get all kmers in a sequence
#' @description
#' Helper function. The input is a set of sequences in a df and it returne a data frame with all kmer in column "sequence".
#'
#' Note: Remove duplicates.
#' @param seq Character, the sequence to aanalyze
#' @param kmer_width kmer size
#' @return
#' A data frame
#' @examples
#' kmer("CATATGCTTGTCTCAAAGTTAAGCCA", kmer_width = 8)
#' @export
#'

kmer  <- function(seq, kmer_width = 15) {

  n_pos <- str_length(seq)- kmer_width - 1
  i_pos <- c(1:n_pos)

  df_kmer <- map_chr(i_pos, ~ str_sub(seq, start=., end = . + kmer_width-1))
  df_kmer <- data.frame(df_kmer)
  colnames(df_kmer) <- "sequence"
  df_kmer <- df_kmer %>% distinct(kmer_seq)

  return(df_kmer)
}

# kmer_set : Get a set of kmer common to set of sequences --------------------------------

#' @title Get all kmers in a set of sequence
#' @description
#' From a set of sequences with name provided, it returns a data frame with the kmer common to
#' all sequences from the list.  It contains 4 columns: sequence, ref_sequence, start, end.
#'
#' The start and end correspond to the position of the kmers on the reference sequences.
#' @param seq_name Vector containing the sequence names (eg accession)
#' @param sequence Vector containing the actual sequences
#' @param kmer_width kmer size
#' @return
#' A data frame with 4 columns
#' @section To do:
#' This may break down if the same kmer is found twice...
#' @examples
#'
#' @export
#'

  kmer_set <- function(seq_name, sequence, kmer_width = 8) {

# Build data frame
  seq_set <- data.frame(sequence=sequence, seq_name=seq_name)
# Remove any sequence that contains anything else than ATGC
  seq_set <- seq_set  %>% mutate(sequence=str_to_upper(sequence)) %>%
                          filter(!str_detect(sequence, "[^ATGC]"))

# Take as reference sequence the shortes one
  seq_ref <- seq_set %>%  top_n(n=1,wt=dplyr::desc(str_length(sequence)))

# Construct the set of kmers from the reference sequence
  kmer_set <- kmer(seq_ref$sequence,kmer_width = kmer_width)

# Map the kmers over the other sequences.  Produce a matrix of 0 and 1 with the reference sequence as rows and the kmers as column
  kmer_set_map.list <- map(as.character(kmer_set$sequence), ~ vcountPattern(pattern=.x,subject=DNAStringSet(seq_set$sequence)))

# transform to a data frame with the column as the kmer sequences and the line the input sequences
  kmer_set_map <- data.frame(kmer_set_map.list)
  colnames(kmer_set_map) <- kmer_set$sequence
  kmer_set_map$seq_name <-seq_set$seq_name

# Looking for all kmers that are absent at least in 1 sequence.
# Transform the matrix in the long form and screen for any kmer which has at least one zero,
# which means that it is absent in one of the sequence from the set
  kmer_rejected <- kmer_set_map %>% tidyr::gather(key="sequence", value="kmer_found", -seq_name) %>%
                                  filter(kmer_found==0) %>%
                                  distinct(sequence)
# By difference the kmer selected are the one which are not rejected
  kmer_selected <- dplyr::setdiff(kmer_set,kmer_rejected)

# If no kmer return NA
  if (nrow(kmer_selected) > 0) {

# Next lines is to construct a DNA String, not used at this point
  kmer_selected.DNAString <- DNAStringSet(kmer_selected$sequence)
  names(kmer_selected.DNAString) <- kmer_selected$sequence

# Compute the position of the oligomers on the reference sequence (the shortest one selected)
  # This is very tricky because map returns a list
  # It is necessary to use the purr functon map to extract the elements of the list (Nice trick)
  kmer_pos <- map(as.character(kmer_selected$sequence), ~ str_locate(seq_ref$sequence, .x) )

# Build the final df
  kmer_out <- data.frame(sequence = kmer_selected$sequence,
                         ref_seq_name=seq_ref$seq_name,
                         start=map_int(kmer_pos, 1),
                         end=map_int(kmer_pos, 2))

  return(kmer_out)

  } else {
  return(NA)
  }
  }




# dada2_assign : Assign sequences from a fasta file to a taxonomy ------------------------------

#' @title Assign sequences using dada2 wang assigner
#' @description
#' Write out 2 files with the taxonomy and the bootstrap values. Extnesions are .dada2.taxo and .dada2.boot
#'
#' The reference file should be a fasta fiel with taxonomy separated by ";"
#'
#' >Level1;Level2;Level3;Level4;Level5;Level6;
#'
#' ACCTAGAAAGTCGTAGATCGAAGTTGAAGCATCGCCCGATGATCGTCTGAAGC
#' @param seq_file_name Name of the file to assign
#' @param ref_file_name Name of the reference sequence file. By default the latest PR2 database
#' @param tax_levels Vector with the taxa levels
#' @return
#' A data frame with the sequence names, the assignements and the bootstrap values
#' @examples
#' @export
#' @md
#'
dada2_assign <- function(seq_file_name,
                         ref_file_name="C:/daniel.vaulot@gmail.com/Databases/_PR2/versions/4.10.0/pr2_version_4.10.0_dada2.fasta.gz",
                         tax_levels=c("kingdom", "supergroup", "division", "class", "order",
                                      "family", "genus", "species")){

# It is necessary to read the sequences to get the names because dada2 takes the sequecne themselves as names.
  seq_in <- readDNAStringSet(seq_file_name)
  seq_names <- names(seq_in)

  taxa <- assignTaxonomy(seqs=seq_file_name,
                         refFasta=ref_file_name,
                         taxLevels = tax_levels,
                         minBoot = 0, outputBootstraps = TRUE,
                         verbose = TRUE)


  boot <- data.frame(taxa$boot) %>% rename_all(funs(str_c(.,"_boot")))
  dada2_result <- bind_cols(data.frame(seq_name=seq_names) , data.frame(taxa$tax)) %>% bind_cols(boot)

  write_tsv(dada2_result, filename_change_ext(seq_file_name,"dada2.taxo"), na="")

  return(dada2_result)

    }

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

  ps <- subset_taxa(ps, (division %in% c("Chlorophyta", "Dinophyta", "Cryptophyta",
                                                 "Haptophyta", "Ochrophyta", "Cercozoa")) &
                                              !(class %in% c("Syndiniales", "Sarcomonadea")))
}

# phyloseq_filter_abundant_taxa : Filter a phyloseq table for abundant taxa ------------------------------

#' @title Filter a phyloseq table keeping only abundant taxa
#' @description
#' @param ps phyloseq object
#' @param fraction_min minimum fraction that an OTU must represent in any given samples (not overall!)
#' @return
#' ps object with only abundant taxa
#' @examples
#' # Will return the otus that more that represent more than 10% of total in any given sample
#' ps_abundant <- phyloseq_filter_abundant_taxa(ps)
#' # Idem but with lower threshold (5%)
#' ps_abundant <- phyloseq_filter_abundant_taxa(ps, 0.05)
#' @export
#' @md
phyloseq_filter_abundant_taxa <- function(ps, fraction_min=0.10) {

  ps <- filter_taxa(ps, function(x) sum(x > total*fraction_min) > 0, TRUE)
}

# get_primer_position : get primer position ----------------------------------------------

#' @title Get primer position on sequence set
#'
#' @description
#' Returns start end end of the primer.  If there are more than one matches, return the first match and the number of matches
#' @param primer_seq primer sequence
#' @param target_seq target sequences (vector)
#' @param orientation "fwd" or "rev"
#' @param mismatches Number of mismatches allowed
#' @return
#' data frame with four columns:
#' - $id - row of sequences
#' - $start - start of match (NA if no match)
#' - $end - end of match (NA if no match)
#' - $n_matches - number of matches (NA if no match)
#' @examples
#' get_primer_position("ATT",c("ATTTTCGGG", "AGTTTCGGG"), orientation="fwd", mismatches=0)
#' @md
#' @export
get_primer_position<- function(primer_seq, target_seq, orientation="fwd", mismatches=0){

  # Strings treated as strings
  options(stringsAsFactors = FALSE)

  # For debugging
  # primer_seq = "TTT"
  # target_seq = "ATTTTCGGG"
  # target_seq = c("ATTTTCGGG", "AGGGAAAA")
  # orientation="fwd"
  # mismatches=0

  # Transform primer in DNA String
  primer_seq <-  DNAString(primer_seq)
  if (orientation == "rev") primer_seq <-  reverseComplement(primer_seq)

  # Create a data frame to hold the sequence position
  df.seq <- data.frame(id=as.character(1:length(target_seq)))

  # Transform target sequences into DNA String set
  target_seq <- DNAStringSet(target_seq)
  names(target_seq) <- df.seq$id

  pos <- vmatchPattern(primer_seq, target_seq, max.mismatch=mismatches, min.mismatch=0,
                      with.indels=FALSE, fixed=FALSE, algorithm="auto")

  # Need to unlist the MIndex position
  list.pos <- unlist(pos)

  # Create data frames from the list
  df.pos <- data.frame(id = list.pos@NAMES,
                        start=list.pos@start,
                        end=list.pos@start + list.pos@width - 1)

  # Compute number of matches
  df.matches <- df.pos %>%
    group_by(id) %>%
    summarize(n_matches=n())

  # Retain only the first match
  df.pos <- df.pos %>%
    group_by(id) %>%
    summarise(start = min(start),
              end = min(end))

  # Merge with the sequence list to identify sequences with n_matches
  df.pos <- df.seq %>%
    left_join(df.pos) %>%
    left_join(df.matches)

 return(df.pos)
 }

# pcr_sequences : in silico amplification of sequences ----------------------------------------------

#' @title In silico amplification
#'
#' @description
#' Returns start end end of the primer.  If there are more than one matches, return the first match and the number of matches
#' @param fwd_seq forward primer sequence
#' @param rev_seq reverse primer sequence
#' @param target_seq target sequences (vector)
#' @param mismatches Number of mismatches allowed
#' @return
#' vector with the amplicon or NA if no matches
#' @examples
#' pcr_sequences("ATT","CCC", c("ATTTTCGGG", "AGTTTCGGG"), mismatches=0)
#' @md
#' @export
pcr_sequences <- function(fwd_seq, rev_seq, target_seq, mismatches=0){

  # Strings treated as strings
  options(stringsAsFactors = FALSE)

  # For debugging
  # fwd_seq = "TTT"
  # rev_seq = "CCC"
  # target_seq = "ATTTTCGGG"
  # target_seq = c("ATTTTCGGG", "AGGGAAAA")
  # mismatches=0

  fwd_pos <- get_primer_position(fwd_seq, target_seq, orientation="fwd", mismatches = mismatches) %>%
    select(fwd_start = start)
  rev_pos <- get_primer_position(rev_seq, target_seq, orientation="rev", mismatches = mismatches) %>%
    select(rev_end = end)

  df <- data.frame(sequence = target_seq) %>%
    bind_cols(fwd_pos) %>%
    bind_cols(rev_pos) %>%
    mutate(amplicon = case_when ((is.na(fwd_pos)|is.na(rev_pos)) ~ NA_character_,
                                 TRUE ~ stringr::str_sub(sequence, fwd_start, rev_end)))


 return(df$amplicon)
 }
