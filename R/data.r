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
