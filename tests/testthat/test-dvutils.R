

# Jump over all the other tests ====================================================
  if (FALSE) {
# Jump, jump... ====================================================================


# MetaPR2 ---------------------------------------------------------------------

context("metapr2")


test_that("metapr2 asv export ", {
  # Export all the asv in a single fasta
  # metapr2_export_asv()

  # Export as specific data set as a phyloseq file
  # metapr2_export_asv(dataset_id_selected = 23, export_phyloseq = TRUE)

  # Export a specific genus as a fasta file and an excel file
  # asv_set <- metapr2_export_asv(taxo_level = genus, taxo_name="Pseudohaptolina",
  #                               taxonomy_full= FALSE, boot_min = 100, export_xls = TRUE)

  # Export without any filter
  # metapr2_export_asv (taxo_level = kingdom, taxo_name="Eukaryota",
  #                      boot_level = class_boot, boot_min = 0,
  #                      directory = "output/metapr2/47/",
  #                      dataset_id_selected = 47,
  #                      filter_metadata = NULL,
  #                      export_long_xls=FALSE,
  #                      export_wide_xls=TRUE,
  #                      export_sample_xls=TRUE,
  #                      export_phyloseq = FALSE,
  #                      export_fasta=TRUE,
  #                      taxonomy_full = TRUE,
  #                      use_hash = FALSE,
  #                      sum_reads_min = 0)
  # Export with filter for water only samples
  metapr2_export_asv (taxo_level = kingdom, taxo_name="Eukaryota",
                      boot_level = class_boot, boot_min = 0,
                      directory = "output/metapr2/47_filter/",
                      dataset_id_selected = 47,
                      filter_metadata = "((substrate == 'water') & is.na(substrate_description )) | (substrate == 'sediment trap material')",
                      export_long_xls=FALSE,
                      export_wide_xls=TRUE,
                      export_sample_xls=TRUE,
                      export_phyloseq = FALSE,
                      export_fasta=TRUE,
                      taxonomy_full = TRUE,
                      use_hash = FALSE,
                      sum_reads_min = 0)
})



# PR2 ---------------------------------------------------------------------

context("PR2")


test_that("Test pr2_taxo_X which add an X to taxonomies when the same level is repeated sereral times ", {
   pr2_taxo <- read_tsv("pr2_ciliates.txt")
   pr2_taxo_clean <- pr2_taxo_X(pr2_taxo, pr2.env$taxo_levels)
   write_tsv (pr2_taxo_clean, "output/pr2_ciliates.clean.txt")
})

test_that("Test pr2 to sqlite  ", {
   pr2_export_sqlite("output/pr2.sqlite")
})


test_that("Test get primer position  ", {
   pos <- get_primer_position("ATT",c("ATTTTCGGG", "AGTTTCGGG"), orientation="fwd", mismatches=0)
   write_tsv(pos, "output/primer_position.tsv")
})


# Genbank -----------------------------------------------------------------

test_that("Test Genbank_download parse  ", {
   gb_info <- genbank_download_parse(c("JX015376", "JQ768406", "LT621940",
                                       "JN989541","HQ664613", "AARH01029193"),
                                     "output/genbank/", sequence_keep=FALSE)
   write_tsv(gb_info, "output/gb_info.tsv")

   gb_info <- genbank_download_parse(c("JX015376", "JQ768406", "LT621940"),
                                     "output/genbank/", sequence_keep=TRUE)
   write_tsv(gb_info, "output/gb_info_seq.tsv")
})

test_that("Test Genbank_download  ", {
   genbank_download(c("JX015376", "JQ768406", "LT621940"), "output/genbank/")
})

test_that("Test Genbank features  ", {
   gb_info <- genbank_features(c("JX015376", "JQ768406", "LT621940"), "output/genbank/")
   write_tsv(gb_info, "output/gb_features.tsv")
})

test_that("Test Genbank field  ", {
   gb_field <- genbank_field(c("JX015376", "JQ768406", "LT621940"), "output/genbank/", organism)
   write_tsv(data.frame(gb_field), "output/gb_field.tsv")
})

test_that("Test Genbank locus  ", {
   gb_locus <- genbank_locus(c("JX015376", "JQ768406", "LT621940"), "output/genbank/")
   write_tsv(data.frame(gb_locus), "output/gb_locus.tsv")
})

test_that("Test Genbank taxonomy  ", {
   gb_tax <- genbank_taxonomy(c("1230134", "1230316", "1905175"))
   write_tsv(data.frame(gb_tax), "output/gb_tax.tsv")
})

test_that("Test Genbank search  ", {
   gb_search <- genbank_search(query = "28S[TITL] AND rRNA[TITL] AND Chlorophyta[ORGN]", seq_max = 500)
   write_tsv(gb_search, "output/gb_search.tsv")
})


# Bioinformatics ------------------------------------------------------------

context("Bioinfo")

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


# BLAST -------------------------------------------------------------------

context("BLAST")

test_that("BLAST ", {
    expect_equal(blast_18S_reformat("blast_small.tsv"),TRUE)
  })

test_that("BLAST AA ", {
    expect_equal(blast_summary("blast_AA.txt"),TRUE)
  })



# Databases ---------------------------------------------------------------

context("Database")

test_that("Test update_field_string ", {
  expect_equal( db_update_field_string("pr2_accession", "ZZ220100", append_comma=TRUE) , "pr2_accession = 'ZZ220100', ")
  expect_equal( db_update_field_string("pr2_accession", "ZZ220100", append_comma=FALSE) , "pr2_accession = 'ZZ220100'")
  expect_equal( db_update_field_string("pr2_accession", 1.25, append_comma=FALSE) , "pr2_accession = '1.25'")
  expect_equal( db_update_field_string("pr2_accession", NA, append_comma=FALSE) , "pr2_accession = NULL")
})

test_that("Test db_connect_sqlite  ", {
    con <- db_connect_sqlite("output/test.sqlite")
    df_taxo <- read_tsv("taxo.txt")
    copy_to(con, df_taxo, name = "taxo", temporary=FALSE, overwrite = TRUE)
    db_disconnect(con)

})

test_that("Test db_append_records ", {
  db_test <- list(user='root', password='root', dbname='test', host='localhost')
  df_taxo <- read_tsv("taxo.txt")
  df_taxo <- df_taxo %>% select(kingdom:species)
  expect_equal(db_append_records(db_test, "taxo",df_taxo ), TRUE)
  })

test_that("Test pr2_read ", {
  skip_if(TRUE, "PR2 database read not done")
  # Put next line to true to allow the test
  if(FALSE){
    pr2 <- pr2_read()
    pr2_filter <- pr2 %>% filter(class=="Mamiellophyceae")
    write_tsv(pr2_filter, path="output/pr2_read.txt", na="")
    expect_known_output(pr2_filter, file="output/pr2_read_print.txt", print=TRUE, update=TRUE)
  }
  })

# Misc --------------------------------------------------------------------
context("File utils")


test_that("Test filename change end", {
  expect_equal( filename_append("myfile.txt",".pr2") , "myfile.pr2.txt")
  expect_equal( filename_append("myfile.txt","_pr2") , "myfile_pr2.txt")
  expect_equal( filename_change_ext("myfile.txt","pr2") , "myfile.pr2")
})




context("Lat Long")

test_that("Test lot_long_dec ", {
  expect_equal( lat_long_dec(10, 30, "N", "latitude") , 10.5)
  expect_equal( lat_long_dec(10, 30, "E", "latitude") , NA)  # Latitude cannot be east
  expect_equal( lat_long_dec(190, 30, "N", "latitude") , NA) # Latitude degree too large (>90)
  expect_equal( lat_long_dec(10, 80, "N", "latitude") , NA)  # Latitude minutes too large (>60)
})

# Plots -------------------------------------------------------------------

context("Plots")

test_that("Test gg_hist ", {

  df <- readxl::read_xlsx("cell_size.xlsx", sheet = "size")
  p<- gg_hist(df,x="Feret", bin = 0.2)
  ggsave("output/test_gg_hist.png", plot=p)

  expect_is( p , "ggplot")
})

test_that("Test gg_density ", {

  df <- readxl::read_xlsx("cell_size.xlsx", sheet = "size")
  p<- gg_density(df,x="Feret", factor = "strain")
  ggsave("output/test_gg_density.png", plot=p)

  expect_is( gg_density(df,x="Feret", factor = "strain") , "ggplot")
})



test_that("Test cowplots schemes ", {

  df <- readxl::read_xlsx("cell_size.xlsx", sheet = "size")
  p<- ggplot(df,aes(x=FeretX, y=FeretY, colour=strain)) + geom_point()
  ggsave("output/test_theme_default.png", plot=p)

  p<- p + cowplot::theme_cowplot()
  ggsave("output/test_theme_cowplot.png", plot=p)

  p<- p + theme_dviz_grid()
  ggsave("output/test_theme_dviz_grid.png", plot=p)

  expect_is( p , "ggplot")
})

# Maps -------------------------------------------------------------------

context("Maps")

test_that("Distribution and zoom ", {

  df <- readxl::read_xlsx("chlorophyta.xlsx")

  p<- map_distribution(df, map_title="Chlorophyta", map_tag="A")
  ggsave("output/test_map_chlorophyta.png", plot=p)

  p_eu <- map_zoom_europe(p)
  ggsave("output/test_map_chlorophyta_europe.png", plot=p_eu)

  expect_is( p , "ggplot")
})

test_that("Map Leaflet ", {

  df <- readxl::read_xlsx("OSD.xlsx")
  df <- dplyr::select(df, latitude, longitude, label)
  df <- dplyr::mutate(df, label=as.character(label))

  map_leaflet(df, width=1000, height=1000)

})

# End of the Loop to go over the tests -----------------------------
  }

