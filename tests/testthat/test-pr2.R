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

