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
