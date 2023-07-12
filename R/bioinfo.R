#' @import dplyr
#' @import tidyr
#' @import stringr
#' @import tibble
#' @import Biostrings

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

  seq_out <- Biostrings::DNAStringSet(df$sequence)

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

  Biostrings::writeXStringSet(seq_out, file_name, compress=compress, width = 20000)

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

  df <- data.frame(sequence=df, seq_name=names(df))
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
  if (type=="AA") {seq_in <- Biostrings::readAAStringSet(file_name)}

# Transform to data frame
  df <- data.frame (names=names(seq_in),sequence=as.character(seq_in), length=nchar(seq_in) )
  print(str_c("Number of sequences initially: ", nrow(df)))

# Remove ambiguities
  if (type=="AA") {df <- df %>% filter(str_count(sequence,"X")<= max_ambig)}

# Filter based on lengths
  df <- df %>% filter((length >= min_length)& (length <= max_length))

# Construct the Biostring set
  seq_out <- Biostrings::AAStringSet(df$sequence)
  names(seq_out) <- df$names
  print(str_c("Number of sequences after filtration: ", nrow(df)))

# Write fasta file
  file_name_out <- str_c(fs::path_ext_remove(file_name),".filtered.",fs::path_ext(file_name))
  Biostrings::writeXStringSet(seq_out, file_name_out, compress=FALSE, width = max_length)
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
     file_in <- ShortRead::FastqFile(fnsR1[i])

    # Sample the number of sequences needed
     sampler <- ShortRead::FastqStreamer(con=fnsR1[i], n=n_sampled)

    # Draw an instance of the sequences
     sample <- ShortRead::yield(sampler)

    # Write the Fastq file
     file_out <- str_replace(fnsR1[i],"R1", "R1.subsample")
     print(file_out)
     ShortRead::writeFastq(sample, file=file_out, compress=FALSE)
     close(file_in)

  # Process R2
     file_in <- ShortRead::FastqFile(fnsR2[i])
     sampler <- ShortRead::FastqStreamer(con=fnsR2[i], n=n_sampled)
     set.seed(123) # The seed need to be reset before each yield to have the matching pairs
     sample<-ShortRead::yield(sampler)
     file_out <- str_replace(fnsR2[i],"R2", "R2.subsample")
     print(file_out)
     ShortRead::writeFastq(sample, file=file_out, compress=FALSE)
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

# dada2_export -----------------------------------------
#' @title Export a data frame to dada2 AssignTaxonomy compatible format
#' @description
#' This will save  data frame to dada2 AssignTaxonomy compatible format
#' @param df data frame - should contain 8 or 9 columns with taxo levels and 1 column with sequences
#' @param file_name character - full path of file where to save
#' @return
#' Write the file in compressed format (.gz)
#' @examples
#' pr2_export_dada2(my_sequences, "C:/Daniel/myfile.fas.gz", 8)
#' @export
#' @md
dada2_export <- function(df, file_name, taxo_levels_number = 9 ) {

  seq_out <- Biostrings::DNAStringSet(df$sequence)  # Store the sequence in a  DNAString - not used anymore

  # Remove any space from Strain, Clone name, Specimen_voucher

    if (taxo_levels_number == 9) {
      df <- df %>%
          mutate ( name_dada2 = str_c(domain,
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
      df <- df %>%
          mutate ( name_dada2 = str_c(kingdom,
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

  names(seq_out) <- df$name_dada2

  # When using DNA strings only LF added (Unix convention)
  # Set width to max possible so that no carriage return (causes problems with dada2)

    Biostrings::writeXStringSet(seq_out, file_name, compress=TRUE, width = 20000)
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
                         ref_file_name="C:/daniel.vaulot@gmail.com/Databases/_PR2/versions/4.12.0/pr2_version_4.12.0_18S_dada2.fasta.gz",
                         tax_levels=c("kingdom", "supergroup", "division", "class", "order",
                                      "family", "genus", "species")){

# It is necessary to read the sequences to get the names because dada2 takes the sequecne themselves as names.
  seq_in <- Biostrings::readDNAStringSet(seq_file_name)
  seq_names <- names(seq_in)

  taxa <- dada2::assignTaxonomy(seqs=seq_file_name,
                         refFasta=ref_file_name,
                         taxLevels = tax_levels,
                         minBoot = 0, outputBootstraps = TRUE,
                         verbose = TRUE)


  boot <- data.frame(taxa$boot) %>%
    rename_all(funs(str_c(.,"_boot")))
  dada2_result <- data.frame(seq_name=seq_names) %>%
    bind_cols(data.frame(taxa$tax)) %>%
    bind_cols(boot)

  readr::write_tsv(dada2_result, filename_change_ext(seq_file_name,"dada2.taxo"), na="")

  return(dada2_result)

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
  primer_seq <-  Biostrings::DNAString(primer_seq)
  if (orientation == "rev") primer_seq <-  reverseComplement(primer_seq)

  # Create a data frame to hold the sequence position
  df.seq <- data.frame(id=as.character(1:length(target_seq)))

  # Transform target sequences into DNA String set
  target_seq <- Biostrings::DNAStringSet(target_seq)
  names(target_seq) <- df.seq$id

  pos <- Biostrings::vmatchPattern(primer_seq, target_seq, max.mismatch=mismatches, min.mismatch=0,
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
                                 TRUE ~ str_sub(sequence, fwd_start, rev_end)))


 return(df$amplicon)
 }


# seq_reverse_complement : Reverse complement ----------------------------------------------

#' @title Reverse complement string
#'
#' @description
#' Reverse complement a set of sequences
#' @param seq vector of sequences
#' @return
#' Reverse complement vector of strings
#' @examples
#' seq_reverse_complement (c("ATTTTCGGG", "ATTTTCGGG"))
#' @md
#' @export
seq_reverse_complement <- function(seqs){

  return(as.character(reverseComplement(DNAStringSet(seqs))))

}

# seq_hash : Sequence hash (vectorial) -----------------------------------------

#' @title Sequence hash (vectorial)
#'
#' @description
#' Get the hash value for a sequences
#' @param seq vector of sequences
#' @return
#' Hash value as a vector
#' @examples
#' seq_hash (c("ATTTTCGGG", "ATTTTTGGG"))
#' @md
#' @export
seq_hash <- function(seqs){

  purrr::map_chr(seqs, digest::sha1)

}
