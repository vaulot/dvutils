
# Plots -------------------------------------------------------------------

context("Plots")

test_that("Test gg_hist ", {

  df <- readxl::read_xlsx("data/cell_size.xlsx", sheet = "size")
  p<- gg_hist(df,x="Feret", bin = 0.2)
  ggsave("output/test_gg_hist.png", plot=p)

  expect_is( p , "ggplot")
})

test_that("Test gg_density ", {

  df <- readxl::read_xlsx("data/cell_size.xlsx", sheet = "size")
  p<- gg_density(df,x="Feret", factor = "strain")
  ggsave("output/test_gg_density.png", plot=p)

  expect_is( gg_density(df,x="Feret", factor = "strain") , "ggplot")
})



test_that("Test cowplots schemes ", {

  df <- readxl::read_xlsx("data/cell_size.xlsx", sheet = "size")
  p<- ggplot(df,aes(x=FeretX, y=FeretY, colour=strain)) + geom_point()
  ggsave("output/test_theme_default.png", plot=p)

  p<- p + cowplot::theme_cowplot()
  ggsave("output/test_theme_cowplot.png", plot=p)

  p<- p + theme_dviz_grid()
  ggsave("output/test_theme_dviz_grid.png", plot=p)

  expect_is( p , "ggplot")
})
