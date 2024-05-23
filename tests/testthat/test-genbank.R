# Bioinformatics ------------------------------------------------------------

context("GenBank")

test_that("Test Genbank Search  ", {
  query = "(ITS1[All Fields] OR ITS2[All Fields] OR 28S[All Fields])  AND (protists[filter] AND is_nuccore[filter])"
  pr2_query <- genbank_search(query = query , seq_max = 1000)
  write_tsv(pr2_query, "output/genbank_search.tsv")
})

