# Latitude, Longitude conversion ------------------------------------------

#' Convert lat and long to decimal
#'
#' @param degree integer
#' @param minute integer
#' @param orientation "E" or "W" for long, "N" or "S" for lat, for the equitor can be anything
#' @param coordinate_type "longitude" or "latitude"
#'
#' @return
#' Decimal coordinates or NA if coordinate_type is wrong
#' @export
#'
#' @examples
#' lat_dec <- lat_long_dec(10, 30, "N", "latitude")

lat_long_dec <- function(degree, minute, orientation, coordinate_type) {
      if (is.na(orientation)) {return(NA)}
      else {
        if (coordinate_type == "longitude") {
          if (orientation == "E") {sign = 1}
          else if (orientation == "W") {sign = -1}
          else {sign=0}
          if ((degree > 180)| (minute >= 60)) {sign=0}
        }
        if (coordinate_type == "latitude") {
          if (orientation == "N") {sign = 1}
          else if (orientation == "S") {sign = -1}
          else {sign = 0}
          if ((degree > 90) | (minute >= 60)) {sign=0}
        }
        if (sign != 0) {return((degree+minute/60)*sign)}
        else {return(NA)}
      }
}

# Latitude, Longitude conversion ------------------------------------------

#' Convert lat and long to decimal
#'
#' @param x character of the form 34° 52' 49"S	or 173° 17' 53"E
#' Must use purrr::map_dbl for mutate
#'
#' @return
#' Decimal coordinates or NA if no coordinate
#' @export
#'
#' @examples
#' lat_dec <- lat_long_string_dec("34° 52' 49\"S")

lat_long_string_dec <- function(x) {

  if(is.na(x)) return(NA)

  # Remove any space
  x <- stringr::str_replace_all(x," ", "")

  # Convert using sp function
  as.numeric(sp::char2dms(from = x,  chd="°", chm="'", chs = "\""))
}
