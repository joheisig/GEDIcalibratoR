#' Get ALS tiles by intersection
#'
#' Non-interactive version of `select_tiles2download()`. Takes an input of point locations and finds intersecting ALS tiles.
#' Currently supports open LiDAR archives from three German states: Northrhine-Westfalia (NRW), Thuringia (TH), and Saxony (SN).
#'
#' @param region Region to select tiles from. Options are 'NRW', 'TH', 'SN', or 'all' (default).
#' @param gedi_obs Point locations of GEDI observations. Should be a `sf` or `sfc` object.
#' @param quiet If TRUE, hide message with ALS tile names. Default is FALSE.
#'
#' @return sf object with download links for tiles of ALS point clouds,
#' digital elevation models, and digital surface models. Used as input to download_tiles().
#' @export
#' @seealso [GEDIcalibratoR::select_tiles2download()]
#' @examples
#' library(sf)
#' gedi = system.file("demodata/gedi.gpkg", package = "GEDIcalibratoR") |> st_read()
#' tiles = intersect_tiles2download(gedi_obs = gedi)
#' head(tiles)
intersect_tiles2download = function(region = "all", gedi_obs = NULL, quiet = F){
  if (! any(inherits(gedi_obs, c("sf", "sfc"))))
    stop("Please provide GEDI observations as georeferenced point dataset (sf/sfc object).", call. = F)
  if (missing(region)) region = "all"
  if (! region %in%  c("TH","NRW","SN", "all"))
    stop(call. = F,paste(region,"is not a valid region. Please select one of TH, NRW, SN or all."))
  if (! region == "all") als = dplyr::filter(als, state == region)
  if (! sf::st_crs(gedi_obs) == sf::st_crs(als)) als = sf::st_transform(als, sf::st_crs(gedi_obs))

  tiles = als[gedi_obs,] |> unique()

  count_obs = sf::st_intersection(als, sf::st_geometry(gedi_obs)) |>
    sf::st_drop_geometry() |>
    dplyr::group_by(tile_laz) |>
    dplyr::summarise(n_obs = dplyr::n())

  tiles = dplyr::inner_join(tiles, count_obs, by = "tile_laz")

  if (!quiet){
    message(paste(nrow(gedi_obs), "GEDI observations intersect the following", nrow(tiles),"tiles:"))
    message(crayon::yellow(paste(tiles$tile_laz, collapse = "\n")))
  }
  return(tiles)
}

#' Get ALS tiles by manual selection
#'
#' Interactive version of `intersect_tiles2download()`. Opens a map from which one can select ALS tiles for download.
#' Currently supports open LiDAR archives from three German states: Northrhine-Westfalia (NRW), Thuringia (TH), and Saxony (SN).
#'
#' @param region Region to select tiles from. Options are 'NRW', 'TH', 'SN', or 'all' (default). Selecting a region speeds up map display.
#' @param gedi_obs Point locations of GEDI observations. Should be a `sf` or `sfc` object. Not required, but helpful for visual comparison.
#' @param quiet If TRUE, hide message with ALS tile names. Default is FALSE.
#' @param mode Select tiles with either a "click" (default) or by "draw"-ing a polygon.
#'
#' @return sf object with download links for tiles of ALS point clouds,
#' digital elevation models, and digital surface models. Used as input to download_tiles().
#' @export
#' @seealso [GEDIcalibratoR::intersect_tiles2download()]
#' @examples
#' \dontrun{
#' library(sf)
#' gedi = system.file("demodata/gedi.gpkg", package = "GEDIcalibratoR") |> st_read()
#' tiles = select_tiles2download(gedi_obs = gedi)
#' head(tiles)
#' }
select_tiles2download = function(region = "all", gedi_obs = NULL, quiet = FALSE, mode = "click"){
  if (!missing(gedi_obs) & !any(inherits(gedi_obs, c("sf", "sfc"))))
    stop("Please provide GEDI observations as georeferenced point dataset (sf/sfc object).
         GEDI observations are helpful but not required for select_tiles2download() to work.", call. = F)
  if (missing(region)) region = "all"
  if (! region %in%  c("TH","NRW","SN", "all")){
    stop(call. = F,
         message(paste(region,"is not a valid region. Please select one of TH, NRW, SN or all."))
    )}

  if (! region == "all"){
    als = dplyr::filter(als, state == region)
    overview_layers = lapply(overview_layers, function(x) dplyr::filter(x, state == region))
  }

  basemap = mapview::mapview(overview_layers[[1]], zcol="year_laz", layer.name = "ALS Year",hide = T,
                             popup = leafpop::popupTable(overview_layers[[1]], zcol = c(1,2), feature.id = F)) +
    mapview::mapview(overview_layers[[2]], zcol="year_dem", hide = T, legend = T,layer.name = "DEM Year",
                     popup = leafpop::popupTable(overview_layers[[2]], zcol = c(1,2), feature.id = F)) +
    mapview::mapview(overview_layers[[3]], zcol="year_dsm", hide = T, legend = T,layer.name = "DSM Year",
                     popup = leafpop::popupTable(overview_layers[[3]], zcol = c(1,2), feature.id = F))

  if (! missing(gedi_obs)){
    gedi_obs = dplyr::group_by(gedi_obs, file, beam) |>
      dplyr::summarise() |> sf::st_cast("LINESTRING")
    basemap = basemap + mapview::mapview(gedi_obs, color = "orange",
                                         layer.name = "GEDI", popup=F)
  }
  files = mapedit::selectFeatures(als, alpha = 0.2, color="white",
                                  title="Select ALS tiles for download",
                                  label = als$tile_laz, map = basemap, mode)
  if (!quiet){
    message(paste(nrow(files),"tiles selected:"))
    message(crayon::yellow(paste(files$tile_laz, collapse = "\n")))
  }
  return(files)
}


#' Download selected LAZ, DEM, and DSM tiles
#'
#' @param tiles Table containing tiles to download. Output of [GEDIcalibratoR::intersect_tiles2download()]
#' or [GEDIcalibratoR::select_tiles2download()]. Tiles are skipped if they already exist in `dir`.
#' @param dir Directory to download tiles to.
#' @param what Either "LAZ", "DEM", or "DSM" (default).
#' @param setTimeOut If e.g. laz files are large, download performance may suffer from the default
#' timeout setting (60 sec). Adjust this in case you run into difficulties.
#'
#' @return Paths to downloaded files.
#' @export
#'
#' @examples
#' \dontrun{
#' library(sf)
#' library(dplyr)
#'
#' gedi = system.file("demodata/gedi.gpkg", package = "GEDIcalibratoR") |> st_read()
#' tiles = intersect_tiles2download(gedi_obs=gedi)
#'
#' # subset to keep one tile from each region
#' tiles = group_by(tiles, state) |> slice_head(n=1)
#'
#' # download to temporary directory
#' td = tempdir()
#' dsm_tiles = download_tiles(tiles, dir = td, what = "DSM")
#' dem_tiles = download_tiles(tiles, dir = td, what = "DEM")
#' laz_tiles = download_tiles(tiles, dir = td, what = "LAZ", setTimeOut = 300)
#'
#' file.exists(c(dsm_tiles, dem_tiles, laz_tiles))
#'
#' # delete files
#' unlink(c(dsm_tiles, dem_tiles, laz_tiles))
#'}
download_tiles = function(tiles, dir, what = c("DEM","DSM"), setTimeOut = NULL){
  if (missing(tiles)) stop("First select ALS tiles using select_tiles2download() or intersect_tiles2download().")
  if (! any(what %in%  c("LAZ","DEM","DSM"))){
    stop(call. = F, "Invalid file type for Download. Choose from \n(1) LAZ - raw ALS point cloud,
         \n(2) DEM - digital elevation model, or \n(3) DSM - digital surface model.")
  }
  if (! dir.exists(dir)) stop("Directory does not exist!")
  if (! missing(setTimeOut)) options(timeout = setTimeOut)

  dlfiles = as.vector(NULL)
  for (w in what) {
    tiles_w = dplyr::select(tiles, dplyr::ends_with(w))
    names(tiles_w) = gsub(paste0(w,"|_"), "", names(tiles_w), ignore.case = T)
    tiles_w = tiles_w[order(tiles_w$ftype),]
    if (any(grepl("zip|gz", tiles_w$ftype))) tmp = tempfile()

    pb = txtProgressBar(min = 0, max = nrow(tiles_w), initial = 0, style = 3)
    for (t in 1:nrow(tiles_w)) {
      if (tiles_w[t,]$ftype %in% c("laz","tif")){
        download_als_laz(tiles_w[t,], dir)
      } else if (tiles_w[t,]$ftype == "gz") {
        download_als_gz(tiles_w[t,], dir, tmp)
      } else {
        download_als_zip(tiles_w[t,], dir, tmp)
      }
      dlfiles = append(dlfiles, tiles_w[t,]$tile)
      setTxtProgressBar(pb, t)
      close(pb)
    }
  }
  message(crayon::yellow(paste("\n", paste(dlfiles, sep=" -- "))))
  message(paste(length(dlfiles), "tiles were successfully downloaded to", dir, "."))
  return(file.path(dir, dlfiles))
}

#' Download LAZ, DEM, and DSM tiles
#'
#' Internal function used by [GEDIcalibratoR::download_tiles()]
#'
#' @param tile single row from tiles table, an output of [GEDIcalibratoR::intersect_tiles2download()]
#' or [GEDIcalibratoR::select_tiles2download()].
#' @param dir directory to download tiles to.
#' @param tmp temporary file
download_als_zip = function(tile, dir, tmp){
  laz = file.path(dir, tile$tile)
  if (! file.exists(laz)){
    download.file(tile$url, tmp, quiet = F)
    unzip(tmp, exdir = dir, files = basename(laz))
  }
}

#' @rdname download_als_zip
download_als_gz = function(tile, dir, tmp){
  laz = file.path(dir, tile$tile)
  if (! file.exists(laz)){
    download.file(tile$url, tmp, quiet = F)
    R.utils::gunzip(tmp, destname = file.path(dir,basename(laz)))
  }
}

#' @rdname download_als_zip
download_als_laz = function(tile, dir){
  laz = file.path(dir, tile$tile)
  if (! file.exists(laz)){
    download.file(tile$url, laz, quiet = T)
  }
}
