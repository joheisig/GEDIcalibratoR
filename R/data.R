#' Small subset of GEDI observations for demo
#'
#' GEDI shots observed over Germany, subset from orbit
#' GEDI02_A_2020181043919_O08756_02_T05629_02_003_01_V002,
#' mainly located in forested areas.
#'
#' @usage data("gedi")
#' @format A sf object with 51 rows and 104 variables:
#' \describe{
#'   \item{rh1-rh99}{relative heights in millimeters}
#'   \item{shot_number}{shot number}
#'   \item{solar_azimuth}{solar azimuth angle}
#'   \item{solar_elevation}{solar elevation angle}
#'   \item{geom}{sf geometry cloumn}
#' }
"gedi"
