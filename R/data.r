#' A dataset taken from the global terrorism database (GTD) (http://www.start.umd.edu/gtd/)
#' containing information of terrorism activity 2010--2016
#' @name terrorism
#' @format A data frame with 55215 rows and 8 variables:
#' \describe{
#' \item{nkill}{number of fatalities per attack}
#' \item{latitude}{latitude location}
#' \item{longitude}{longitude location}
#' \item{iyear}{year index 2010--2016}
#' \item{x.coord}{x coordinate from location projected onto a spheree}
#' \item{y.coord}{y coordinate from location projected onto a sphere}
#' \item{z.coord}{z coordinate from location projected onto a sphere}
#' \item{fatal}{binary, 1 = fatalities occured, 0 = no fatalities}
#' }
#' @docType data
#' @usage data(terrorism)
NULL

#' Spatial polygon of countries worldwide
#' @name world
#' @format A Spatial Polygons Data Frame of contries worldwide
#' @docType data
#' @keywords datasets
#' @usage data(world)
NULL
#' Spatial polygon of New Zealand (including the Chatham islands)
#' @name NZ
#' @format A Spatial Polygons Data Frame of New Zealand
#' @docType data
#' @keywords datasets
#' @usage data(NZ)
NULL
#' A dataset taken from the GeoNet Quake (http://quakesearch.geonet.org.nz/)
#' containing earthquake information in New Zealand (and the Chatham islands) from 14-Februrary-2018 to 14-March-2018
#' @name earthquakes
#' @format A data frame containing 382 rows and 6 variables:
#' \describe{
#' \item{latitude}{latitude location}
#' \item{longitude}{longitude location}
#' \item{origintime} {The UTC time of the event's occurrence (in ISO 8601 format)}
#' \item{depth}{The focal depth of the event (km)}
#' \item{magnitude}{The magnitude of the earthquake}
#' \item{publicid}{The Public Id, a unique earthquake reference code}
#' }
#' @docType data
#' @usage data(earthquakes}
NULL
