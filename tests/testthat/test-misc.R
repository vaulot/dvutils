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
