% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/location_calibration.R
\name{get_best_shift}
\alias{get_best_shift}
\alias{get_no_shift}
\title{Match original GEDI data with no vs. the best-performing shift}
\usage{
get_best_shift(ext_data, shifts)

get_no_shift(ext_data, shifts)
}
\arguments{
\item{ext_data}{Extracted reference values. Output of \code{\link[=cal_extract]{cal_extract()}}.}

\item{shifts}{Output of \code{\link[=shift_stats]{shift_stats()}}.}
}
\value{
sf object with original GEDI measurements and locations shifted according to the best-performing parameters.
}
\description{
Match original GEDI data with no vs. the best-performing shift
}
\examples{
grid_ext = system.file("demodata/grid_ext.rds", package = "GEDIcalibratoR") |> readRDS()
shifts = shift_stats(grid_ext)
get_best_shift(grid_ext, shifts)
get_no_shift(grid_ext, shifts)
}
