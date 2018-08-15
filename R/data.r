#' A dataset taken from the global terrorism database (GTD) (http://www.start.umd.edu/gtd/)
#' containing information of terrorism activity 2010--2017
#' @name terrorism
#' @format A data frame with 72366 rows and 16 variables:
#' \describe{
#' \item{iyear}{numeric year 2010--2017}
#' \item{imonth}{numeric month index 1--12}
#' \item{iday}{numeric day 1--31 (zeros are a non-entry)}
#' \item{country}{country}
#' \item{latitude}{latitude location}
#' \item{longitude}{longitude location}
#' \item{sucess}{Logical fatal or not. TRUE = fatal}
#' \item{nkill}{number of fatalities per attack}
#' \item{specificity}{factor which represents the apatial accuracy of the evets: 1 = most accurate, 5 = worst}
#' \item{gname}{character name of attack perpetrators}
#' \item{x.coord}{x coordinate from location projected onto a sphere}
#' \item{y.coord}{y coordinate from location projected onto a sphere}
#' \item{z.coord}{z coordinate from location projected onto a sphere}
#' \item{popdensity}{scaled: number of people per kilometer squared}
#' \item{luminosity}{scaled: luminosity}
#' \item{tt}{scaled: time to nearest city in minutes}
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
#' containing earthquake information in New Zealand (and the Chatham islands) from 26-March-2018 to 26-April-2018
#' @name earthquakes
#' @format A data frame containing 368 rows and 5 variables:
#' \describe{
#' \item{latitude}{latitude location}
#' \item{longitude}{longitude location}
#' \item{origintime}{The UTC time of the event's occurrence}
#' \item{depth}{The focal depth of the event (km)}
#' \item{magnitude}{The magnitude of the earthquake}
#' }
#' @docType data
#' @usage data(earthquakes)
NULL
#' A simulated dataset indicating breeding pair presence of Grus Grus in the wetlands of England.
#' @details Data simulated by Andrea Soriano Redondo (A.Soriano-Redondo@exeter.ac.uk) for the paper
#' "Estimating species distribution in highly dynamic populations using point process models" submitted to Ecography.
#' @name cranes
#' @format A data frame containing 5052 rows and 10 variables:
#' \describe{
#' \item{Wetland_Identity}{numeric wetland ID}
#' \item{Lon}{longitude epicentre location of wetland}
#' \item{Lat}{latitude epicentre location of wetland}
#' \item{Area}{Area of wetland in m^2}
#' \item{Perimiter}{Wetland perimiter in m}
#' \item{Wet_density_buf_NoSea}{Surrounding wetland density}
#' \item{Urb_density_buf_NoSea}{Surrouning urban density}
#' \item{Year}{Year of observation, i.e., 2014 or 2015}
#' \item{mark}{A simulated binary mark indicating presence of a Grus Grus breeding pair at the wetland}
#' \item{PA_ratio}{Wetland perimiter to area ratio}
#' }
#' @docType data
#' @usage data(cranes)
NULL
#' A simulated realisation of a log-Gaussian Cox process 2D
#' @name lgcp2D
#' @format A data frame containing 281 rowa and 2 variables:
#' \describe{
#' \item{x}{x-coordinate of the point pattern}
#' \item{y}{y-coordinate of the point pattern}
#' }
#' @docType data
#' @usage data(lgcp2D)
NULL
#' A cleaned aggregated dataset orginally taken from the global terrorism database (GTD)
#' (http://www.start.umd.edu/gtd/)
#' containing information of terrorism activity worldwide 2010--2017
#' @name terrorism_aggregate
#' @format A data frame containing 10845 rows and 20 variables:
#' \describe{
#' \item{latitude}{latitude location}
#' \item{longitude}{longitude location}
#' \item{iyear}{year index 2010--2016}
#' \item{total}{total number of events}
#' \item{lethal}{number of lethal events}
#' \item{imonth}{numeric month index 1--12}
#' \item{iday}{numeric day 1--31 (zeros are a non-entry)}
#' \item{country}{country}
#' \item{sucess}{Logical fatal or not. TRUE = fatal}
#' \item{nkill}{total number of fatalities}
#' \item{specificity}{factor which represents the apatial accuracy of the evets: 1 = most accurate, 5 = worst}
#' \item{gname}{character name of attack perpetrators}
#' \item{x.coord}{x coordinate from location projected onto a sphere}
#' \item{y.coord}{y coordinate from location projected onto a sphere}
#' \item{z.coord}{z coordinate from location projected onto a sphere}
#' \item{popdensity}{scaled: number of people per kilometer squared}
#' \item{luminosity}{scaled: luminosity}
#' \item{tt}{scaled: time to nearest city in minutes}
#' \item{country.idx}{numeric country index}
#' \item{originalcount}{total number of fatalities}
#' }
#' @docType data
#' @usage data(terrorism_aggregate)
NULL
