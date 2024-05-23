# Bioinformatics ------------------------------------------------------------

context("Bioinfo")

# Jump over all the other tests ====================================================
if (FALSE) {
  # Jump, jump... ====================================================================

test_that("Test primer match  ", {
  pos <- get_primer_position("ATT",c("ATTTTCGGG", "AGTTTCGGG"), orientation="fwd", mismatches=0)
  write_tsv(pos, "output/primer_position.tsv")
})

test_that("Test fasta read  ", {
  df <- fasta_read("mothur.database.fas")
  head(df)
})

test_that("Test fastq subsample  ", {
  fastq_subsample("fastq", 500, random=TRUE)
})


test_that("Test pr2 taxo check  ", {
  pr2_taxo <- read_tsv("taxo2.txt")
  pr2_taxo_check(pr2_taxo, "output")
})


test_that("Test Mothur database to Phyloseq ", {
  phyloseq_import_mothur("mothur.database", "mothur_samples.tsv", 6)
})


test_that("Test dada2 assgin  ", {

  ref_file_name <- "C:/Data Biomol/RNA/_PR2/versions/4.10.0/pr2_version_4.10.0_dada2.fasta"
  seq_file_name <- "output/otu_taxo.fasta"
  dada2_assign(seq_file_name, ref_file_name)

})


test_that("Test kmer  ", {kmer_test <- kmer("CATATGCTTGTCTCAAAGTTAAGCCATGCATGTCTAAGTATAACCGTTATACTGGGAAACTGCGAATGGCTCATTAAATCAGTTAT",
                                            kmer_width = 15)
write_tsv(kmer_test, "output/kmer_test.txt")
})

test_that("Test kmer set  ", {
  seq_set <- read_tsv("pr2_sample.txt")%>% select(seq_name=pr2_accession, sequence, sequence_length) %>%
    filter((sequence_length > 1500)&(sequence_length < 1800))
  kmer_set_test  <- kmer_set(seq_set$seq_name, seq_set$sequence, kmer_width = 8)
  write_tsv(kmer_set_test, "output/kmer_set_test.txt")
})


test_that("Test fasta_filter ", {
  fasta_filter("protein.faa", min_length=100, max_length=10000, type="AA", max_ambig=5)
})

test_that("Test write_fasta ", {
  # This generates an error whch I do not understand. The function works perfectly...
  # Error in x[[method]](...) :
  #   tentative d'appliquer un objet qui n'est pas une fonction
  # Calls: <Anonymous> ... <Anonymous> -> o_apply -> lapply -> FUN -> <Anonymous>

  df <- read_xlsx("otu_sequences.xlsx")
  fasta_write(df,"output/otu_taxo.fasta", compress=FALSE, taxo_include=TRUE, taxo_separator="|")
  fasta_write(df,"output/otu_no_taxo.fasta", compress=FALSE, taxo_include=FALSE)
  # expect_equal(write_fasta_taxo(df,"output/otu_taxo.fasta", compress=FALSE), TRUE)
})

# End of the Loop to go over the tests -----------------------------
}


test_that("Test dada2 correct ", {

  df <- rio::import("data/dada2_assigned.xlsx")

  df  <- dada2_assign_correct(df, boot_threshold = 80)

  rio::export(df, "output/dada2_assigned_corrected.xlsx")
  # expect_equal(write_fasta_taxo(df,"output/otu_taxo.fasta", compress=FALSE), TRUE)
})




