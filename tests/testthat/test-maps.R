context("Maps")

test_that("Distribution and zoom ", {

  df <- readxl::read_xlsx("data/chlorophyta.xlsx")

  p<- map_distribution(df, map_title="Chlorophyta", map_tag="A")
  ggsave("output/test_map_chlorophyta.png", plot=p)

  p_eu <- map_zoom_europe(p)
  ggsave("output/test_map_chlorophyta_europe.png", plot=p_eu)

  expect_is( p , "ggplot")
})

test_that("Map Leaflet ", {

  df <- readxl::read_xlsx("data/OSD.xlsx")
  df <- dplyr::select(df, latitude, longitude, label)
  df <- dplyr::mutate(df, label=as.character(label))

  map_leaflet(df, width=1000, height=1000)

})
