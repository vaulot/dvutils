#' @import ggplot2
#' @import rworldmap
#' @import viridis
#' @import ggmap
#' @import leaflet


# Get the world map as data frame -----------------------------------------

#' @title Get the world map for the rworldmap package
#' @description
#' @param resolution "coarse" for global maps, "low" for regional maps
#' @return
#' Data frame with world map
#' @examples
#' @export
#' @md
#'
map_get_world <- function(resolution="coarse"){
  worldMap <- getMap(resolution = resolution) # Change to "coarse" for global maps / "low" for regional maps
  world.points <- fortify(worldMap)
  world.points$region <- world.points$id
  world.df <- world.points[,c("long","lat","group", "region")]
  }

# Function to draw a world background map ---------------------------------
#' @title Background map using the maps package
#' @description
#' @param color_continents color for the continents (grey80 by default)
#' @param color_borders color for the borders (white by default)
#' @param resolution "coarse" for global maps, "low" for regional maps
#' @return
#' A map as a ggplot object
#' @examples
#' map_world() + geom_point(data=cultures_one_year, aes(x=Longitude, y=Latitude), fill="blue", size=2, shape=21)
#' @export
#' @md

map_world <- function(color_continents = "grey80", color_borders = "white", resolution = "coarse") {

  # Background map using the maps package
  # world.df <- map_data("world")

  world.df <- map_get_world(resolution)

  map <- ggplot() +
    geom_polygon(data = world.df, aes(x=long, y = lat, group = group), fill=color_continents, color=color_borders) +
    # scale_fill_manual(values= color_continents , guide = FALSE) +
    scale_x_continuous(breaks = (-4:4) * 45) +
    scale_y_continuous(breaks = (-2:2) * 30) +
    xlab("Longitude") + ylab("Latitude") +
    coord_fixed(1.3) +
    theme_bw()
    # species_map <- species_map + coord_map ()  # Mercator projection
    # species_map <- species_map + coord_map("gilbert") # Nice for the poles
  return(map)
  }

# Function to draw a US background map ---------------------------------

  # NOTE - This function can be made more general

#' @title Background US map using the maps package to extract a specific country
#' @description
#' @param color_continents color for the continents (grey80 by default)
#' @param color_borders color for the borders (white by default)
#' @return
#' A map as a ggplot object
#' @examples
#' map_US() + geom_point(data=cultures_one_year, aes(x=Longitude, y=Latitude), fill="blue", size=2, shape=21)
#' @export
#' @md

map_US <- function(color_continents = "grey80", color_borders = "white") {

  world.df <- map_get_world()
  north_america.df <- world.df %>% filter(region %in% c("United States of America","Canada"))

  map <- ggplot() +
    geom_polygon(data = north_america.df, aes(x=long, y = lat, group = group), fill=color_continents, color=color_borders) +
    # scale_fill_manual(values= color_continents, guide = FALSE) +
    scale_x_continuous(breaks = (-4:4) * 45) +
    scale_y_continuous(breaks = (-3:3) * 30) +
    xlab("Longitude") + ylab("Latitude") +
    coord_fixed(1.3, xlim = c(-180, -45),ylim = c(20, 90)) +
    theme_bw()
    # species_map <- species_map + coord_map ()  # Mercator projection
    # species_map <- species_map + coord_map("gilbert") # Nice for the poles
  return(map)
  }

# Function to draw a world background map with Google in color---------------------------------

#' @title Background map using the ggmap package
#' @description
#' @param location an address, longitude/latitude pair (in that order), or left/bottom/right/top bounding box
#' @param zoom map zoom, an integer from 3 (continent) to 21 (building), default value 10 (city). "auto" automatically determines the zoom for bounding box specifications, and is defaulted to 10 with center/zoom specifications. maps of the whole world currently not supported.
#' @param maptype character string providing map theme. options available are "terrain", "terrain-background", "satellite", "roadmap", and "hybrid" (google maps)
#' @return
#' A map as a ggplot object
#' @examples
#' map_world_google() + geom_point(data=cultures_one_year, aes(x=Longitude, y=Latitude), fill="blue", size=2, shape=21)
#' map_world_google(location="Singapore", zoom=11 ,maptype="hybrid")
#' map_world_google(location="Paris", zoom="auto", maptype="terrain")
#' @export
#' @md

map_world_google <- function(location=c(0, 0), zoom = 3, maptype = "satellite") {

  background <- ggmap::get_map(location = location, zoom = zoom,  scale = 2,
                               color = "color", maptype = maptype)
  map_gg <- ggmap::ggmap(background) + xlab("Longitude") + ylab("Longitude")
  return(map_gg)

}

# Function to draw a  map with Leafket in color---------------------------------

#' @title Map with leaflet package (OpenStreetMap)
#' @description
#' @param df Dataframe having at least 3 columns longitude, latitude,label for label to be plotted
#' @param lng_center longitude of the center of the map
#' @param lat_center latitude of the center of the map
#' @param zoom map zoom, an integer from 3 (continent) to 21 (building), default value 10 (city). "auto" automatically determines the zoom for bounding box specifications, and is defaulted to 10 with center/zoom specifications. maps of the whole world currently not supported.
#' @param width width in pixels
#' @param height height in pixels
#' @return
#' A map as a leaflet object
#' @examples
#' map_leaflet(df, lng_center=0, lat_center=45, zoom=4)
#' @export
#' @md

map_leaflet <- function(df, lng_center=0, lat_center=0, zoom = 3, width=500,  height=500) {

  map <- leaflet(width = width, height = height) %>%
    addTiles() %>%
    setView(lng=lng_center, lat=lat_center, zoom=zoom) %>%
    addCircleMarkers(data = df, lat = ~ latitude, lng = ~ longitude,
                     radius = 5,
               label = ~ label,
               labelOptions = labelOptions(textsize = "10px", noHide = T),
               clusterOptions = markerClusterOptions())

}

# Function to draw a map of a dataframe ---------------------------------

#' @title Draw distribution map from data frame
#' @description
#' @param df Dataframe having at least 3 columns z for the quantity to be mapped, longitude, latitude
#' @param z_limits Limits for the legend
#' @param z_breaks Breaks for the legend
#' @param z_min Lower limit below which the location is represented by a cross
#' @param color_not_present Color used to display points where the z is below z_min
#' @param base_size Font size for the theme (can increase to make axes legends and titles bigger)
#' @param map_title Title of the map
#' @param map_tag Tag ("A", "B") to be displayed at top right
#' @param legend_position Position of the legend
#' @param legend_title Title of the legend
#' @return
#' A map as a ggplot object
#' @examples
#' map_distribution(df)
#' @export
#' @md

map_distribution <- function(df,z_limits = c(0,100),
                             z_breaks=c(10,25, 50, 75, 100),
                             z_min=1,
                             color_not_present= "blue",
                             size_not_present=1 ,
                             base_size = 14,
                             map_title="",
                             map_tag="",
                             legend_position =c(0.15,0.25),
                             legend_title="%" ) {

theme_map_distribution <- theme_light(base_size = base_size) +
                          theme ( legend.position = legend_position,
                                  legend.background = element_rect(fill = "transparent"),
                                  legend.text=element_text(size=12),
                                  legend.title =element_text(size=12),
                                  plot.tag.position="topright",
                                  plot.tag= element_text(size=20, face="bold"),
                                  panel.grid.minor = element_blank())

 map <-  map_world()+
    # coord_map(xlim = c(-180, 180),ylim = c(-90, 90)) +
    # scale_x_continuous(breaks = (-4:4) * 45) +
    # scale_y_continuous(breaks = (-2:2) * 30) +
    geom_point(data= filter(df,z < z_min),
               aes(x=longitude, y=latitude),
               color=color_not_present,
               size=size_not_present, shape=3) +
    geom_point(data=filter(df, z >= z_min),
               aes(x=longitude, y=latitude, size=z, color=z, alpha=z) ) +
    scale_size(name = legend_title, range = c(0, 8),limits = z_limits, breaks=z_breaks)+
    scale_alpha_continuous(name=legend_title, range=c(0.5, .9), limits = z_limits, breaks=z_breaks) +
    viridis::scale_color_viridis(option="magma", name=legend_title, limits = z_limits, breaks=z_breaks ) +
    guides( colour = guide_legend()) +
    labs(title = map_title, tag = map_tag) +
    theme_map_distribution

  return(map)
  }


# Function to zoom on Europe ---------------------------------

#' @title Zoom a map on Europe
#' @description
#' @param map Map to be zoomed
#' @return
#' A map as a ggplot object
#' @examples
#' map_zoom_europe(map)
#' @export
#' @md

map_zoom_europe <- function(map){
  map <- map +
    coord_fixed(1.3, xlim = c(-40, 40),ylim = c(30, 70)) +
    scale_x_continuous(breaks = (-4:4) * 10) +
    scale_y_continuous(breaks = (3:7) * 10)
}
