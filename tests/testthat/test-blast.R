# BLAST -------------------------------------------------------------------

context("BLAST")

test_that("BLAST ", {
  expect_equal(blast_18S_reformat("blast_small.tsv"),TRUE)
})

test_that("BLAST AA ", {
  expect_equal(blast_summary("blast_AA.txt"),TRUE)
})

