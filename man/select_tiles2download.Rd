% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tiles_download.R
\name{select_tiles2download}
\alias{select_tiles2download}
\title{Get ALS tiles by manual selection}
\usage{
select_tiles2download(
  region = "all",
  gedi_obs = NULL,
  quiet = FALSE,
  mode = "click"
)
}
\arguments{
\item{region}{Region to select tiles from. Options are 'NRW', 'TH', 'SN', or 'all' (default). Selecting a region speeds up map display.}

\item{gedi_obs}{Point locations of GEDI observations. Should be a \code{sf} or \code{sfc} object. Not required, but helpful for visual comparison.}

\item{quiet}{If TRUE, hide message with ALS tile names. Default is FALSE.}

\item{mode}{Select tiles with either a "click" (default) or by "draw"-ing a polygon.}
}
\value{
sf object with download links for tiles of ALS point clouds,
digital elevation models, and digital surface models. Used as input to download_tiles().
}
\description{
Interactive version of \code{intersect_tiles2download()}. Opens a map from which one can select ALS tiles for download.
Currently supports open LiDAR archives from three German states: Northrhine-Westfalia (NRW), Thuringia (TH), and Saxony (SN).
}
\examples{
\dontrun{
library(sf)
gedi = system.file("demodata/gedi.gpkg", package = "GEDIcalibratoR") |> st_read()
tiles = select_tiles2download(gedi_obs = gedi)
head(tiles)
}
}
\seealso{
\code{\link[=intersect_tiles2download]{intersect_tiles2download()}}
}
