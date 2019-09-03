#' @import readr
#' @import dplyr
#' @import stringr
#' @import pr2database

# options(stringsAsFactors=FALSE)

#' @title Process a tabular output from blastn (BLAST+)
#' @description
#' The BLAST file must originate from blastn with the follwing output format option:
#'
#' -outfmt"6 qseqid sseqid sacc stitle sscinames staxids sskingdoms sblastnames pident
#'           slen length mismatch gapopen qstart qend sstart send evalue bitscore"
#'
#' It is very important that the columns are in this precise order and no column is missing.\cr
#'
#' For the formating options see:
#' * [NCBI](https://www.ncbi.nlm.nih.gov/books/NBK279684/)
#' * [Metagenomics wiki](http://www.metagenomics.wiki/tools/blast/blastn-output-format-6).
#'
#' What does the function do :
#' 1. Group all GenBank accession
#' 2. Obtain taxonomy from GenBank (note the GenBank taxonomy is now in the PR2 database after downloading from \url{ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/})
#' 3. Check if the sequence is in PR2 and get the PR2 taxo (this is done with the pr2 database package)
#' 4. Detect whether the Subject sequence is uncultured or not
#' 5. Merge back into the BLAST file
#' 6. Compute a summary with best hit, best hit to PR2, best hit to cultured, taxo consensus (identity>96%), contradiction at division level (identity>90%)
#'
#'The **modified BLAST output file** includes additional columns
#' * kingdom -> species : PR2 taxonomy for those accession numbers that are present in PR2
#' * hit_rank : the rank of the hit based on decreasing % identity and decreasing bit scores
#' * uncultured : TRUE if the hit corresponds to an uncultured item
#' * hit_lineage : GenBank taxonomy of the hit
#'
#' The **summary file** contains several set of columns
#' 1. The top hit (column with prefix hit_top_)
#' 2. The top hit for which a PR2 sequence is available (columns starting with hit_pr2_)
#' 3. The top hit corresponding to a culture or an isolate (columns starting with hit_cul)
#' 4. A "consensus" taxonomy based on all the hits with more than 98\% identity, keeping the most frequent hit at the genus level and using the sum of the bit scores to decide if there are some ties.
#' 5. Contradiction between hits >90\% identity at the division level
#' @param file_name The name of the BLAST file with full path
#' @return
#' TRUE if the function has been successful.
#'
#' The modified table is saved by changing the name of the file by replacing the extension by **_pr2.tsv**.
#'
#' The summary table is saved by changing the name of the file by replacing the extension by **_summary.tsv**.
#' @section To do:
#'
#' @section Column names:
#' The columns for the Blast are named as follows. For the summary a prefix is added
#'
#'       hit_top_ / hit_pr2_ / hit_cult_
#'
#'       query_id, hit_id, hit_acc, hit_title, hit_sci_names
#'
#'       hit_tax_ids, hit_super_kingdoms, hit_blast_names,
#'
#'       pct_identity, hit_length, alignment_length, mismatches,
#'
#'       gap_opens, query_start, query_end, hit_start, hit_end,
#'
#'       evalue, bit_score
#' @section Programming notes:
#' The following functions must be used with libary qualifier **dplyr::** because they are also in the plyr library
#' * ungroup
#' * desc
#' * rename
#'
#' Uses the pr2database package for faster access (much faster !!)
#' @examples
#' blast_18S_reformat("C:/BLAST_output.tsv")
#' @export
#' @md

blast_18S_reformat <- function(file_name){

#  For testing
#  file_name = "C:/Data Biomol/RNA/Tags/BIOSOPE/blast/dada2_seq_18S.tsv"

# Define the columns of the blast file

  blast_columns<-data.frame(name=c("query_id","hit_id","hit_acc","hit_title","hit_sci_names",
                                   "hit_tax_ids","hit_super_kingdoms","hit_blast_names",
                                   "pct_identity","hit_length","alignment_length","mismatches",
                                   "gap_opens","query_start","query_end","hit_start","hit_end",
                                   "evalue","bit_score"),
                           outfmt=c("qseqid","sseqid","sacc","stitle","sscinames",
                                    "staxids","sskingdoms","sblastnames",
                                    "pident","slen","length","mismatch",
                                    "gapopen","qstart","qend","sstart","send",
                                    "evalue","bitscore"),
                           stringsAsFactors=FALSE
                  )

# Read the file
  # some sequences have # so cannot use comment = "#"
  # some fields have N/A when missing information
  blast<- read_tsv(file_name, col_names = blast_columns$name, na=c("N/A",""), guess_max=10000)

# order by otu, decreasing percent, increasing evalue, and add a rank variable
  blast<- blast %>% arrange(query_id, dplyr::desc(pct_identity),evalue, dplyr::desc(bit_score)) %>%
                    group_by(query_id) %>%
                    mutate(hit_rank=row_number()) %>%
                    ungroup()

# Keep a single taxo_id per line (some are of the form "100272;588879")
  blast <- blast %>%
           mutate(hit_tax_ids=str_replace_all(hit_tax_ids,";(\\d)+",""))

# List distinct taxonomy id numbers to search PR2
# Note : Pull is necessary to extract a vector from a tibble
  blast_taxo_id <- blast %>%
               select(hit_tax_ids) %>%
               distinct()

# Get the lineage from GenBank
  # The taxonomy table is now in the pr2 database.
  # Using dplyr to create SQL query

# 1. Create the database connection
# Use the local PR2 database
  pr2_db_con <- db_connect(db_info("pr2_local"))

# 2. Connect to table gb_taxonomy
  gb_taxo <-tbl(pr2_db_con, "gb_taxonomy")

# 3. Join the blast taxo_id and the genbank taxo_id.
  # The order of the join is VERY important because the right side is sent to the server.
  # If we had a left_join then the whole genbank database would have been sent to the server !
  # Do not forget to collect to get the database
  blast_taxo_id <- right_join(gb_taxo, blast_taxo_id, by=c("tax_id"="hit_tax_ids"), copy=TRUE) %>%
                   collect()
  db_disconnect(pr2_db_con)

# 4. Some accession numbers have \r at the end... I am not sure where comes from
  blast_taxo_id <- blast_taxo_id %>% mutate(tax_id=str_replace(tax_id,"\\r", ""))

# Get the PR2 sequences and taxonomy

# List distinct accession numbers to search PR2
  blast_acc <- blast %>%
               select(hit_acc, hit_title) %>%
               distinct()

# Merge blast with pr2
  data(pr2)
  # Remove entries in PR2 corresponding to the same GenBank (they create duplicated lines in the summary output)
  pr2 <- pr2 %>% select(genbank_accession, kingdom:genus, species) %>% distinct()
  blast_acc <- left_join(blast_acc, pr2,  by=c("hit_acc"="genbank_accession"))

# Rename the taxonomy columns to pr2_xxx (very very pedestrian...)
  # blast_acc <-  blast_acc %>%
  #               dplyr::rename(kingdom_pr2=kingdom , supergroup_pr2=supergroup , division_pr2=division ,
  #                             class_pr2=class , order_pr2=order , family_pr2=family ,
  #                             genus_pr2=genus , species_pr2=species)

# Detect uncultured
  blast_acc <- blast_acc %>% mutate(uncultured=str_detect(hit_title, "Uncultured|uncultured|eukaryote clone"))

# Do final merging
  blast_new <- left_join(blast, blast_acc)
  blast_new <- left_join(blast_new, blast_taxo_id, by=c("hit_tax_ids"="tax_id"))

# Do some summaries ---

  # Count the number of hit lines per query and build a table
  blast_count <-  blast_new %>%
                    group_by(query_id) %>%
                    summarize(hit_number=n())

  # Extract the top hit line
  blast_first <- blast_new %>%
                    filter(hit_rank==1) %>%
                    rename_all(funs(str_c("hit_top_",.))) %>%
                    dplyr::rename(query_id=hit_top_query_id)

  # Extract the top hit for cultured (must use desc() because if not then). Also MUST ungroup before renaming !!!
  blast_first_uncultured <- blast_new %>%
                    filter(uncultured==FALSE) %>%
                    group_by(query_id) %>%
                    top_n(1, dplyr::desc(hit_rank)) %>%
                    ungroup()%>%
                    rename_all(funs(str_c("hit_cult_",.))) %>%
                    dplyr::rename(query_id=hit_cult_query_id)

  # Extract the top hit corresponding to a pr2 sequence (must use desc() because if not then).
  # Also MUST ungroup before renaming !!!
  blast_first_pr2 <- blast_new %>%
                    filter(!is.na(kingdom)) %>%
                    group_by(query_id) %>%
                    top_n(1, dplyr::desc(hit_rank)) %>%
                    ungroup()%>%
                    rename_all(funs(str_c("hit_pr2_",.))) %>%
                    dplyr::rename(query_id=hit_pr2_query_id)

  # Compute consensus taxo at the genus from the number of hits that have more than 98% similiarity
  # The sum of the bit scores is also computed in case of ties, only the top sum is kept
  # Finally if there is still a tie, the first row is kept...

  # One could use contains and a regular expression to select the columns. I left the explicit ranks for clarity
  taxo_levels = c("kingdom", "supergroup", "division", "class", "order", "family", "genus", "species")
  taxo_regex <- str_c(taxo_levels[1:7], collapse="|")
  taxo_regex <- str_c("(",taxo_regex, ")")

  blast_taxo_consensus <- blast_new %>%
                    filter(!is.na(kingdom)&(pct_identity >=96)) %>%
                    select(query_id, kingdom, supergroup, division, class, order, family, genus, bit_score)%>%
                    group_by(query_id, kingdom, supergroup, division, class, order, family, genus) %>%
                    summarize(n_hits=n(),sum_bit_score = sum(bit_score) )%>%
                    group_by(query_id) %>%
                    top_n(1, wt=n_hits) %>%
                    top_n(1, wt=sum_bit_score) %>%
                    distinct(query_id, .keep_all=TRUE) %>%
                    ungroup()%>%
                    rename_all(funs(str_c("consensus_",.)))%>%
                    dplyr::rename(query_id=consensus_query_id)

  # Check otus with hits belonging to two different divisions (only if > 90% similiarity)

  blast_taxo_contradict <- blast_new %>%
                    filter(!is.na(kingdom)&(pct_identity >=90)) %>%
                    select(query_id, division, pct_identity)%>%
                    group_by(query_id, division) %>%
                    summarize(n_hits=n(),mean_pct_identity = mean(pct_identity) )%>%
                    group_by(query_id) %>%
                    filter(n()>1) %>%
                    arrange(dplyr::desc(n_hits)) %>%
                    summarize(taxo_contradict=str_c("Division: ",division," (N hits>90% = ",n_hits,
                                                    " - mean id = ",sprintf("%.2f",mean_pct_identity), ")",
                                                    collapse = " - "))
  # write_tsv(blast_taxo_contradict, "junk.txt")



  # Assemble now everything together
  blast_summary <-  left_join(blast_count, blast_first) %>%
                    left_join(blast_first_pr2) %>%
                    left_join(blast_first_uncultured) %>%
                    left_join(blast_taxo_consensus) %>%
                    left_join(blast_taxo_contradict)

# Save the files
  file_name_out <- str_c(fs::path_ext_remove(file_name),
                                  "_pr2.",
                                  fs::path_ext(file_name))
  write_tsv(blast_new, path=file_name_out, na="")

  file_name_out <- str_c(fs::path_ext_remove(file_name),
                                  "_summary.",
                                  fs::path_ext(file_name))
  write_tsv(blast_summary, path=file_name_out, na="")

  return(TRUE)
}

#' @title Write a summary for a tabular output from blastn (BLAST+)
#' @description
#' The BLAST file must originate from blastn with the follwing output format option:
#'
#' -outfmt"6 qseqid sseqid sacc stitle sscinames staxids sskingdoms sblastnames pident
#'           slen length mismatch gapopen qstart qend sstart send evalue bitscore"
#'
#' It is very important that the columns are in this precise order and no column is missing.\cr
#'
#' For the formating options see:
#' * [NCBI](https://www.ncbi.nlm.nih.gov/books/NBK279684/)
#' * [Metagenomics wiki](http://www.metagenomics.wiki/tools/blast/blastn-output-format-6).
#'
#' What does the function do :
#' 1. Group all GenBank accession
#' 2. Obtain taxonomy from GenBank (note the GenBank taxonomy is now in the PR2 database after downloading from \url{ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/})
#' 5. Merge back into the BLAST file
#' 6. Compute a summary with best hit,
#'
#' The **summary file** contains several set of columns
#' 1. The top hit (column with prefix hit_top_)
#' @param file_name The name of the BLAST file with full path
#' @return
#' TRUE if the function has been successful.
#'
#' The summary table is saved by changing the name of the file by replacing the extension by **_summary.tsv**.
#' @section To do:
#' @section Column names:
#' The columns for the Blast are named as follows. For the summary a prefix is added
#'
#'       query_id, hit_id, hit_acc, hit_title, hit_sci_names
#'
#'       hit_tax_ids, hit_super_kingdoms, hit_blast_names,
#'
#'       pct_identity, hit_length, alignment_length, mismatches,
#'
#'       gap_opens, query_start, query_end, hit_start, hit_end,
#'
#'       evalue, bit_score
#' @section Programming notes:
#' The following functions must be used with libary qualifier **dplyr::** because they are also in the plyr library
#' * ungroup
#' * desc
#' * rename
#'
#' Uses the local version of the PR2 database for faster access (much faster !!)
#' @examples
#' blast_reformat("C:/BLAST_output.txt")
#' @export
#' @md

blast_summary <- function(file_name){

#  For testing
#  file_name = "C:/Data Biomol/RNA/Tags/BIOSOPE/blast/dada2_seq_18S.tsv"

# Define the columns of the blast file

  blast_columns<-data.frame(name=c("query_id","hit_id","hit_acc","hit_title","hit_sci_names",
                                   "hit_tax_ids","hit_super_kingdoms","hit_blast_names",
                                   "pct_identity","hit_length","alignment_length","mismatches",
                                   "gap_opens","query_start","query_end","hit_start","hit_end",
                                   "evalue","bit_score"),
                           outfmt=c("qseqid","sseqid","sacc","stitle","sscinames",
                                    "staxids","sskingdoms","sblastnames",
                                    "pident","slen","length","mismatch",
                                    "gapopen","qstart","qend","sstart","send",
                                    "evalue","bitscore"),
                           stringsAsFactors=FALSE
                  )

# Read the file
  # some sequences have # so cannot use comment = "#"
  # some fields have N/A when missing information
  blast<- read_tsv(file_name, col_names = blast_columns$name, na=c("N/A",""), guess_max=10000)

# order by otu, decreasing percent, increasing evalue, and add a rank variable
  blast<- blast %>% arrange(query_id, dplyr::desc(pct_identity),evalue, dplyr::desc(bit_score)) %>%
                    group_by(query_id) %>%
                    mutate(hit_rank=row_number()) %>%
                    ungroup()

# Keep a single taxo_id per line (some are of the form "100272;588879")
  blast <- blast %>%
           mutate(hit_tax_ids=str_replace_all(hit_tax_ids,";(\\d)+",""))

# List distinct taxonomy id numbers to search PR2
# Note : Pull is necessary to extract a vector from a tibble
  blast_taxo_id <- blast %>%
               select(hit_tax_ids) %>%
               distinct()

# Get the lineage from GenBank
  # The taxonomy table is now in the pr2 database.
  # Using dplyr to create SQL query

# 1. Create the database connection
# Use the local PR2 database
  pr2_db_con <- db_connect(db_info("pr2_local"))

# 2. Connect to table gb_taxonomy
  gb_taxo <-tbl(pr2_db_con, "gb_taxonomy")

# 3. Join the blast taxo_id and the genbank taxo_id.
  # The order of the join is VERY important because the right side is sent to the server.
  # If we had a left_join then the whole genbank database would have been sent to the server !
  # Do not forget to collect to get the database
  blast_taxo_id <- right_join(gb_taxo, blast_taxo_id, by=c("tax_id"="hit_tax_ids"), copy=TRUE) %>%
                   collect()
  db_disconnect(pr2_db_con)

# 4. Some accession numbers have \r at the end... I am not sure where comes from
  blast_taxo_id <- blast_taxo_id %>% mutate(tax_id=str_replace(tax_id,"\\r", ""))


# Do final merging
  blast_new <- left_join(blast, blast_taxo_id, by=c("hit_tax_ids"="tax_id"))

# Do some summaries ---

  # Count the number of hit lines per query and build a table
  blast_count <-  blast_new %>%
                    group_by(query_id) %>%
                    summarize(hit_number=n())

  # Extract the top hit line
  blast_first <- blast_new %>%
                    filter(hit_rank==1)

# Save the files
  file_name_out <- str_c(fs::path_ext_remove(file_name),
                                  ".summary.",
                                  fs::path_ext(file_name))
  write_tsv(blast_first, path=file_name_out, na="")

  return(TRUE)
}

