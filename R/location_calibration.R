#' @aliases cal_grid_rect
#'
#' @title Create grid layouts for GEDI location calibration
#' @description [GEDIcalibratoR::cal_grid_rect()] and [GEDIcalibratoR::cal_grid_circle()] both create grids to be used in the location
#' calibration process. The original GEDI observations serves as the center of the grid. GEDI Level 2 data comes with a
#' mean location error of 10 meters. Its relationship with ALS data may be enhanced through calibration.
#'
#' @param gedi GEDI point locations
#' @param steps Sequence (vector) of sample point distances from the original GEDI location to span the grid. Should start with 0.
#' @param n_circle_elements Number of points per circle. Only applies to cal_grid_circle(). Defaults to 8.
#' @param copy_GEDI_vars Which GEDI variables to keep for calibration, e.g. "rh98".
#'
#' @return sf object containing calibration grids for each GEDI observation.
#' @export
#'
#' @examples
#' library(sf)
#' library(dplyr)
#' gedi = system.file("demodata/gedi.gpkg", package = "GEDIcalibratoR") |> st_read()
#'
#' gedi_rect = cal_grid_rect(gedi, steps = c(0, 3, 6, 9, 12),
#'                           copy_GEDI_vars = c("file","id","rh98"))
#' gedi_circ = cal_grid_circle(gedi, steps = c(0, 3, 6, 9, 12),
#'                             n_circle_elements = 12,
#'                             copy_GEDI_vars = c("file","id","rh98"))
#' par(mfrow=c(1,2))
#' filter(gedi_rect) |> slice_head(n = 1) |> st_geometry() |> plot()
#' filter(gedi_circ) |> slice_head(n = 1) |> st_geometry() |> plot()
cal_grid_circle = function(gedi, steps = c(0,5,10), n_circle_elements = 8, copy_GEDI_vars = NULL){
  n_circle_elements = as.integer(round(n_circle_elements))
  if (sf::st_is_longlat(gedi)) stop(.call = F, "Please sf::st_transfrom() your GEDI input to a projected CRS
         with units 'metre', like for example ETRS89 / UTM.")
  if (!is.null(copy_GEDI_vars) & all(copy_GEDI_vars %in% names(gedi))){
    coo = cbind(sf::st_coordinates(gedi), sf::st_drop_geometry(gedi)[copy_GEDI_vars]) |>
      as.data.frame() |> setNames(c("GEDI_x", "GEDI_y", copy_GEDI_vars))
  } else {
    coo = sf::st_coordinates(gedi) |> as.data.frame() |> setNames(c("GEDI_x", "GEDI_y"))
  }
  coo$cid = 1:nrow(coo)

  slice = 360 / n_circle_elements
  angle = seq(0,359, slice)
  grid = expand.grid(steps, angle) |>
    unique() |>
    setNames(c("step","angle")) |>
    as.data.frame() |>
    dplyr::filter(step != 0)
  grid = rbind(c(0,0), grid)
  grid$shift_x = grid$step*sin(grid$angle * (pi/180)) |> round(2)
  grid$shift_y = grid$step*cos(grid$angle * (pi/180)) |> round(2)

  coo = do.call("rbind", replicate(nrow(grid), coo, simplify = FALSE)) |>
    dplyr::arrange(cid)
  grid = do.call("rbind", replicate(nrow(gedi), grid, simplify = FALSE))
  grid = cbind(coo, grid)

  grid = dplyr::mutate(grid, grid_x = GEDI_x + shift_x, grid_y = GEDI_y + shift_y) |>
    sf::st_as_sf(coords=c("grid_x","grid_y"), crs=sf::st_crs(gedi)) |>
    dplyr::select(dplyr::starts_with("GEDI"), dplyr::starts_with("shift"), dplyr::everything(),-cid)
  attr(grid, "grid_type") = "circular"
  return(grid)
}

#' @export
#'
#' @rdname cal_grid_circle
cal_grid_rect = function(gedi, steps = c(0,5,10), copy_GEDI_vars = NULL){
  if (sf::st_is_longlat(gedi))
    stop(.call = F, "Please sf::st_transfrom your GEDI input to a projected CRS with units 'metre', like for example ETRS89 / UTM.")
  if (!is.null(copy_GEDI_vars) & all(copy_GEDI_vars %in% names(gedi))){
    coo = cbind(sf::st_coordinates(gedi), sf::st_drop_geometry(gedi)[copy_GEDI_vars]) |>
      as.data.frame() |> setNames(c("GEDI_x", "GEDI_y", copy_GEDI_vars))
  } else {
    coo = sf::st_coordinates(gedi) |> as.data.frame() |> setNames(c("GEDI_x", "GEDI_y"))
  }
  coo$cid = 1:nrow(coo)
  grid = expand.grid(c(-steps, steps), c(-steps, steps)) |>
    unique() |>
    setNames(c("shift_x","shift_y"))
  coo = do.call("rbind", replicate(nrow(grid), coo, simplify = FALSE)) |>
    dplyr::arrange(cid)
  grid = do.call("rbind", replicate(nrow(gedi), grid, simplify = FALSE))
  grid = cbind(coo, grid)

  grid = dplyr::mutate(grid, grid_x = GEDI_x + shift_x, grid_y = GEDI_y + shift_y) |>
    sf::st_as_sf(coords=c("grid_x","grid_y"), crs=sf::st_crs(gedi)) |>
    dplyr::select(dplyr::starts_with("GEDI"), dplyr::starts_with("shift"), dplyr::everything(),-cid)
  attr(grid, "grid_type") = "rectangular"
  return(grid)
}


#' Merge downloaded raster tiles
#'
#' @param paths paths downloaded raster tiles
#' @param CRS force assign CRS, e.g. "EPSG:25832"
#'
#' @return SpatRaster object from merging one or multiple tiles
merge_tiles = function(paths, CRS){
  r = NULL
  if (length(paths) > 1){
    r = lapply(paths, terra::rast)
    r = do.call(terra::merge, r)
    names(r) = "value"
  } else if (length(paths) == 1){
    r = terra::rast(paths)
    names(r) = "value"
  }
  if (! is.null(r) & ! missing(CRS)) terra::crs(r) = CRS
  return(r)
}

#' Create Relative Height Model from DEM and DSM
#'
#' A Relative Height Model (RHM) is a normalized Digital Surface Model (DSM). It is created by simply deducting the
#' Digital Elevation Model (DEM) from the DSM.
#'
#' @param tile_paths Paths to downloaded DSM and DEM tile. Output of [GEDIcalibratoR::download_tiles()].
#' @param aggregate_res Aggregation factor to decrease spatial detail of RHM.
#' @param overwrite_if_exists Recompute and overwrite existing RHM tiles?
#' @param quiet If TRUE, hide message with tile names. Default is FALSE.
#' @param ... Further arguments passed to terra::aggregate(), e.g. fun = "max"
#'
#' @return paths to RHM tiles
#' @export
#'
#' @examples
#' \dontrun{
#' library(sf)
#' gedi = system.file("demodata/gedi.gpkg", package = "GEDIcalibratoR") |> st_read()
#' tiles = intersect_tiles2download(gedi_obs = gedi)
#' down_tiles = download_tiles(tiles, dir = tempdir())
#' rhm_paths = make_rhm(tile_paths, overwrite_if_exists = T)
#' }
make_rhm = function(tile_paths, aggregate_res = NULL, overwrite_if_exists = F, quiet = F, ...){
  nonex = sum(!file.exists(tile_paths))
  if (nonex == length(tile_paths)) stop(.call=F, "\nNone of the provided tiles seem to exists. Please check and download tiles if necessary.")
  if (nonex > 0) warning(.call=F, paste(nonex, "out of", length(tile_paths), "provided tiles do not seem to exists. Continuing with existing tiles."))
  # check length + message
  dir = dirname(tile_paths[1])
  dsms = grep("dom", basename(tile_paths), value = T)
  dems = grep("dgm", basename(tile_paths), value = T)

  order_dsms = gsub(".?dom\\d\\d?|_","", dsms) |> order()
  order_dgms = gsub("dgm\\d?|_","", dems) |> order()

  dsms = dsms[order_dsms]
  dems = dems[order_dgms]

  rhms = sub("dgm1","rhm", dems)
  rhms = file.path(dir, sub(".xyz",".tif", rhms))

  rhmfiles = as.vector(NULL)
  for (i in 1:length(dems)){
    if(any(!file.exists(rhms[i]), overwrite_if_exists)){
      (em = terra::rast(file.path(dir, dems[i])))
      (sm = terra::rast(file.path(dir, dsms[i])))

      if (terra::crs(em) != terra::crs(sm)) terra::crs(em) = terra::crs(sm) = ""
      rem = terra::res(em)[1]; rsm = terra::res(sm)[1]
      if (rem != rsm) {
        larger = rem > rsm
        ifelse(larger, sm <- terra::resample(sm, em), em <- terra::resample(em, sm))
      }
      if (grepl("ndom", dsms[i])){
        rhm = sm
      } else {
        rhm = sm - em
      }
      if (!is.null(aggregate_res)) rhm = terra::aggregate(rhm, aggregate_res, ...)
      terra::writeRaster(rhm, rhms[i], overwrite=T)
      rhmfiles = append(rhmfiles, rhms[i])
      if (!quiet) message(basename(rhms[i]))
    }
  }
  return(rhmfiles)
}

#' Extract reference raster values for a calibration grid
#'
#' @param tile_paths Paths to DSM, DEM or pre-processed RHM tiles.
#' Output of [GEDIcalibratoR::download_tiles()] or [GEDIcalibratoR::make_rhm()].
#' @param gedi_grid Calibration grid. Output of [GEDIcalibratoR::cal_grid_rect()] or [GEDIcalibratoR::cal_grid_circle()].
#' @param buffer spatial circular buffer around each grid point. Default is 12.5 m to match GEDI's 25 m footprint.
#' @param fun see [terra::extract()]
#' @param quiet If FALSE (default), print progress messages.
#'
#' @return The same input calibration grid (sf object) with an additional column named 'ext_value'.
#' @export
#'
#' @examples
#' \dontrun{
#' library(sf)
#' gedi = system.file("demodata/gedi.gpkg", package = "GEDIcalibratoR") |> st_read()
#'
#' tiles = intersect_tiles2download(gedi_obs = gedi)
#' down_tiles = download_tiles(tiles, dir = tempdir())
#' rhm_paths = make_rhm(tile_paths, overwrite_if_exists = T)
#'
#' gedi_circ = cal_grid_circle(gedi)
#'
#' grid_ext = cal_extract(rhm_paths, gedi_circ, buffer = 12.5)
#' summary(grid_ext)
#' hist(grid_ext$ext_value)
#' }
cal_extract = function(tile_paths, gedi_grid, buffer = 12.5, fun = 'max', quiet = F){
  nonex = sum(!file.exists(tile_paths))
  if (nonex == length(tile_paths)) stop(.call=F, "\nNone of the provided tiles seem to exists. Please check and download tiles if necessary.")
  if (nonex > 0) warning(.call=F,paste(nonex, "out of", length(tile_paths), "provided tiles do not seem to exists. Continuing with existing tiles."))
  if (!is.numeric(buffer) | buffer < 0) stop(.call=F, "Please select a buffer size > 0.")
  gtype = attr(gedi_grid, "grid_type")

  # make individual stacks
  if (!quiet) message("Stacking tiles...")
  nrw = merge_tiles(grep("_nw", tile_paths, value = T), "EPSG:25832")
  th = merge_tiles(grep("_th_", tile_paths, value = T), "EPSG:25832")
  sn = merge_tiles(grep("_sn", tile_paths, value = T), "EPSG:25833")
  bb = merge_tiles(grep("_bb", tile_paths, value = T), "EPSG:25833")

  # extract
  if (!quiet) message("Extracting...")
  gedi32 = gedi33 = enrw = eth = esn = ebb = NULL
  if (any(!is.null(nrw), !is.null(th))){
    gedi32 = sf::st_transform(gedi_grid, sf::st_crs(25832)) |> sf::st_buffer(buffer)
    if (!is.null(nrw)){
      enrw = terra::extract(nrw, terra::vect(gedi32), fun = fun, ID = F)[,1] |> as.vector()
      enrw[is.nan(enrw)] = NA
    }
    if (!is.null(th)){
      eth = terra::extract(th, terra::vect(gedi32), fun = fun, ID = F)[,1] |> as.vector()
      eth[is.nan(eth)] = NA
    }
  }
  if (any(!is.null(sn), !is.null(bb))){
    gedi33 = sf::st_transform(gedi_grid, sf::st_crs(25833)) |> sf::st_buffer(buffer)
    if (!is.null(sn)){
      esn = terra::extract(sn, terra::vect(gedi33), fun = fun, ID = F)[,1] |> as.vector()
      esn[is.nan(esn)] = NA
    }
    if (!is.null(bb)){
      ebb = terra::extract(bb, terra::vect(gedi33), fun = fun, ID = F)[,1] |> as.vector()
      ebb[is.nan(ebb)] = NA
    }
  }
  # merge
  gedi_grid$ext_value = rowSums(cbind(enrw, eth, esn, ebb), na.rm = T)
  gedi_grid = dplyr::select(gedi_grid, dplyr::everything(), dplyr::starts_with("geom"))
  attr(gedi_grid, "grid_type") = gtype
  if (!quiet) message("Complete!")
  return(gedi_grid)
}

#' R-squared and RMSE helper functions
#'
#' @param a actual values
#' @param p predicted values
#'
#' @return numeric
rmse = function(a, p){
  sqrt(mean((a - p)^2, na.rm=T))
}

#' @rdname rmse
r2  = function(a, p){
  cor(a, p, use="complete.obs")^2
}

#' Calculate statistics for each calibration shift
#'
#' @param ext_data Extracted reference values. Output of [GEDIcalibratoR::cal_extract()].
#' @param cal_var GEDI variable in ext_data to use for comparison.
#'
#' @return Table with an entry for each file and shift summarizing calibration measures RMSE and R-squared.
#' @export
#'
#' @examples
#' grid_ext = system.file("demodata/grid_ext.rds", package = "GEDIcalibratoR") |> readRDS()
#' shifts = shift_stats(grid_ext)
#' plot_shift_stats(shifts)
shift_stats = function(ext_data, cal_var = "rh98"){
  gtype = attr(ext_data, "grid_type")
  if (gtype == "rectangular") ext_data = dplyr::select(ext_data, dplyr::everything(),
                                                       var = dplyr::all_of(cal_var)) |>
      dplyr::mutate(x = shift_x, y = shift_y)
  if (gtype == "circular") ext_data = dplyr::select(ext_data, dplyr::everything(),
                                                    var = dplyr::all_of(cal_var)) |>
      dplyr::mutate(x = step, y = angle)

  ext_data = dplyr::mutate(ext_data, diff_ext = var - ext_value) |>
    dplyr::group_by(file, x, y) |>
    dplyr::summarise(n = dplyr::n(),
            mean_diff = mean(abs(diff_ext), na.rm=T),
            sd_diff = sd(diff_ext, na.rm=T),
            rmse = rmse(var, ext_value),
            r2 = r2(var, ext_value)) |>
    dplyr::ungroup()
  if (gtype == "rectangular") ext_data = dplyr::select(ext_data, dplyr::everything(),
                                                       shift_x = x, shift_y = y)
  if (gtype == "circular") ext_data = dplyr::select(ext_data, dplyr::everything(),
                                                    step = x, angle = y)
  attr(ext_data, "cal_var") = cal_var
  attr(ext_data, "grid_type") = gtype
  return(ext_data)
}


#' Comparative scatter plot of calibration shifts
#'
#' @param shifts Output of [GEDIcalibratoR::shift_stats()]
#'
#' @export
#' @examples
#' grid_ext = system.file("demodata/grid_ext.rds", package = "GEDIcalibratoR") |> readRDS()
#' shifts = shift_stats(grid_ext)
#' plot_shift_stats(shifts)
plot_shift_stats = function(shifts){
  gtype = attr(shifts, "grid_type")
  if (gtype == "rectangular") plot_shift_stats_rect(shifts)
  if (gtype == "circular") plot_shift_stats_polar(shifts)
}

#' Comparative scatter plot of calibration shifts
#'
#' @param shifts Output of [GEDIcalibratoR::shift_stats()]
plot_shift_stats_rect = function(shifts){
  no_shift = dplyr::group_by(shifts, file) |>
    dplyr::filter(shift_x == 0 , shift_y == 0)

  min_rmse = dplyr::group_by(shifts, file) |>
    dplyr::filter(rmse == min(rmse))
  min_rmse$diff_rmse = no_shift$rmse - min_rmse$rmse

  ggplot2::ggplot(per_shift, ggplot2::aes(x=shift_x, y=shift_y)) +
    ggplot2::geom_point(ggplot2::aes(size=r2, color=rmse)) +
    ggplot2::scale_size(name="R-squared") +
    ggplot2::scale_color_viridis_c(direction = 1, option="H", name="RMSE [m]") +
    ggplot2::facet_wrap(~file) + ggplot2::coord_equal() +
    ggplot2::geom_point(cex=6, stroke=2, data = min_rmse, color="red",
               fill="transparent", pch=21) +
    ggplot2::geom_label(ggplot2::aes(label=paste0(round(min_rmse$rmse,2), " (",round(diff_rmse,2),")")),
                        data = min_rmse,
                        label.size = 0.2, nudge_y = 2, nudge_x = 1,
                        label.padding = ggplot2::unit(0.15, "lines")) +
    ggplot2::ggtitle("Calibration Results per Shift",
                     subtitle = "Red circle indicates minimum RMSE [m] and difference to original observation.")
}

#' @rdname plot_shift_stats_rect
plot_shift_stats_polar = function(shifts){
  no_shift = dplyr::group_by(shifts, file) |>
    dplyr::filter(step == 0 , angle == 0)

  min_rmse = dplyr::group_by(shifts, file) |>
    dplyr::filter(rmse == min(rmse))
  diff_rmse = min_rmse$rmse - no_shift$rmse

  angles = unique(shifts$angle)
  ggplot2::ggplot(shifts, ggplot2::aes(x=angle, y=step)) +
    ggplot2::geom_point(ggplot2::aes(size=r2, color=rmse)) +
    ggplot2::scale_size(name="R-squared", limits = c(0.1,0.9)) +
    ggplot2::scale_color_viridis_c(direction = 1, option="H", name="RMSE [m]") +
    ggplot2::scale_x_continuous(limits=c(0, 360),breaks = angles) +
    ggplot2::facet_wrap(~file) + ggplot2::coord_polar() +
    ggplot2::geom_point(cex=6, stroke=2, data = min_rmse, color="red",
               fill="transparent", pch=21) +
    ggplot2::geom_label(ggplot2::aes(label=paste0(round(min_rmse$rmse,2), " (",round(diff_rmse,2),")")),
                        data = min_rmse,
                        label.size = 0.2, nudge_y = 2, nudge_x = 1,
                        label.padding = ggplot2::unit(0.15, "lines")) +
    ggplot2::ggtitle("Calibration Results per Shift",
                     subtitle = "Red circle indicates minimum RMSE [m] and difference to original observation.") +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
}


#' Match original GEDI data with no vs. the best-performing shift
#'
#' @param ext_data Extracted reference values. Output of [GEDIcalibratoR::cal_extract()].
#' @param shifts Output of [GEDIcalibratoR::shift_stats()].
#'
#' @return sf object with original GEDI measurements and locations shifted according to the best-performing parameters.
#' @export
#'
#' @examples
#' grid_ext = system.file("demodata/grid_ext.rds", package = "GEDIcalibratoR") |> readRDS()
#' shifts = shift_stats(grid_ext)
#' get_best_shift(grid_ext, shifts)
#' get_no_shift(grid_ext, shifts)
get_best_shift = function(ext_data, shifts){

  if (any(c("geometry", "geom") %in% names(shifts))) shifts = sf::st_drop_geometry(shifts)
  shifts$shift_index = 1:nrow(shifts)
  index = shifts |>
    dplyr::group_by(file) |>
    dplyr::summarize(ind = shift_index[which.min(rmse)])

  best_shift = shifts[index$ind,] |> dplyr::select(-shift_index)

  gtype = attr(ext_data, "grid_type")
  if (gtype == "rectangular") best = dplyr::inner_join(ext_data, best_shift, by = c("file", "shift_x", "shift_y"))
  if (gtype == "circular") best = dplyr::inner_join(ext_data, best_shift, by = c("file", "step", "angle"))
  attr(best, "cal_var") = attr(shifts, "cal_var")
  return(best)
}

#' @rdname get_best_shift
#' @export
get_no_shift = function(ext_data, shifts){
  if (any(c("geometry", "geom") %in% names(shifts))) shifts = sf::st_drop_geometry(shifts)
  gtype = attr(ext_data, "grid_type")
  if (gtype == "rectangular"){
    shifts = dplyr::filter(shifts, shift_x == 0, shift_y == 0)
    no = dplyr::inner_join(ext_data, shifts, by = c("file", "shift_x", "shift_y"))
  }
  if (gtype == "circular"){
    shifts = dplyr::filter(shifts, step == 0, angle == 0)
    no = dplyr::inner_join(ext_data, shifts, by = c("file", "step", "angle"))
  }
  attr(no, "cal_var") = attr(shifts, "cal_var")
  return(no)
}

#' Plot calibration results
#'
#' @param cal Calibrated GEDI locations. Output of [GEDIcalibratoR::get_best_shift()].
#' @param ref Optional. Reference GEDI locations, e.g. output of [GEDIcalibratoR::get_no_shift()].
#' @param label Length 2 character vector. Name of the table to be displayed. e.g. "best shift" and "reference".
#' @param color Length 2 character vector. Colors identifying ref and cal.
#'
#' @export
#'
#' @examples
#' grid_ext = system.file("demodata/grid_ext.rds", package = "GEDIcalibratoR") |> readRDS()
#' shifts = shift_stats(grid_ext)
#' bs = get_best_shift(grid_ext, shifts)
#' ns = get_no_shift(grid_ext, shifts)
#' plot_calibration(bs, ns,  label = c("Best","Ref"), color =  c("forestgreen", "firebrick"))
plot_calibration = function(cal, ref = NULL, label = "", color=NULL){

  cal_var = attr(cal, "cal_var")
  cal$cal_val = cal |> dplyr::pull(cal_var)
  cal$label = label[1]
  r1 = mean(cal$rmse)
  r2 = NULL
  if (!is.null(ref)){
    if (length(label) != 2) stop('Please provide labels for inputs cal and ref in the form c("label1","label2")')
    ref$cal_val = ref |> dplyr::pull(cal_var)
    ref$label = label[2]
    rmse2 = paste("\n", label[2], "RMSE =", round(mean(ref$rmse),2),  "m")
    cal = rbind(cal, ref)
  }
  if (is.null(color)) color = rainbow(length(label))
  p = ggplot2::ggplot(cal) +
    ggplot2::geom_point(ggplot2::aes(x = cal_val, y=ext_value, color = label), alpha=0.4) +
    ggplot2::geom_abline(slope = 1, intercept = 0) +
    ggplot2::labs(x = cal_var, y = "ALS") +
    ggplot2::scale_x_continuous(limits = c(0,50)) +
    ggplot2::scale_y_continuous(limits = c(0,50)) +
    ggplot2::scale_color_manual(values = color) +
    ggplot2::ggtitle("GEDI Location Calibration")

  if (is.null(ref)) p = p +
    ggplot2::geom_label(ggplot2::aes(x=10,y=45, label = paste("RMSE =", round(r1, 2), "m")))
  if (!is.null(ref)) p = p +
    ggplot2::geom_label(ggplot2::aes(x=15,y=45, label = paste(label[1], "RMSE =",round(r1, 2), "m", rmse2)))
  print(p)
}

#' Apply spatial shift to GEDI data
#'
#' @param gedi_obs Original GEDI observations
#' @param shift shift containing one unique shift per file. Preferably output of [GEDIcalibratoR::get_best_shift()]
#'
#' @return Calibrated GEDI observations
#' @export
#'
#' @examples
#' library(sf)
#' gedi = system.file("demodata/gedi.gpkg", package = "GEDIcalibratoR") |> st_read()
#' grid_ext = system.file("demodata/grid_ext.rds", package = "GEDIcalibratoR") |> readRDS()
#' shifts = shift_stats(grid_ext)
#' best_shift = get_best_shift(grid_ext, shifts)
#' gedi_cal = apply_shift(gedi, best_shift)
apply_shift = function(gedi_obs, shift){

  if (length(unique(shift$shift_x)) > length(unique(shift$file))) stop("Too many possible shifts per file. Provide only one per file.")
  gedi_cal = cbind(sf::st_drop_geometry(gedi_obs), sf::st_coordinates(gedi_obs))
  gedi_cal = sf::st_drop_geometry(shift) |> dplyr::group_by(file) |>
    dplyr::summarise(shift_x = min(shift_x), shift_y = min(shift_y)) |>
    dplyr::inner_join(gedi_cal, by="file") |>
    dplyr::mutate(X=X+shift_x, Y=Y+shift_y) |>
    sf::st_as_sf(coords=c("X","Y"), crs = sf::st_crs(gedi_obs)) |>
    dplyr::select(-shift_x, -shift_y)
}
