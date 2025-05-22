#' Modified version of ee_extract function
#'
#' See Issue: https://github.com/r-spatial/rgee/issues/367
#'
#' See [rgee::ee_extract()] for documentation
#'
#' @inheritParams rgee::ee_extract
#'
#' @export
ee_extract <- function(x, y, fun = rgee::ee$Reducer$mean(), scale = NULL, sf = FALSE,
                       via = "getInfo", container = "rgee_backup", lazy = FALSE,
                       quiet = FALSE, ...) {
  rgee:::ee_check_packages("ee_extract", c("geojsonio", "sf"))
  if (!quiet & is.null(scale)) {
    scale <- 1000
    message(sprintf("The image scale is set to %s.", scale))
  }
  if (!any(class(x) %in% rgee:::ee_get_spatial_objects("i+ic"))) {
    stop("x is neither an rgee::ee$Image nor rgee::ee$ImageCollection")
  }
  if (any(class(x) %in% "ee.imagecollection.ImageCollection")) {
    x <- rgee::ee$ImageCollection$toBands(x)
  }
  oauth_func_path <- system.file("python/ee_extract.py", package = "rgee")
  extract_py <- rgee:::ee_source_python(oauth_func_path)
  sp_objects <- rgee:::ee_get_spatial_objects("Table")
  if (!any(class(y) %in% c("sf", "sfc", sp_objects))) {
    stop("y is not a sf, sfc, rgee::ee$Geometry, rgee::ee$Feature or rgee::ee$FeatureCollection object.")
  }
  if (any("sf" %in% class(y))) {
    sf_y <- y
    ee_y <- rgee::sf_as_ee(y[[attr(y, "sf_column")]], quiet = TRUE)
  } else if (any("sfc" %in% class(y))) {
    sf_y <- sf::st_sf(id = seq_along(y), geometry = y)
    ee_y <- rgee::sf_as_ee(y, quiet = TRUE)
  } else if (any(rgee::ee_get_spatial_objects("Table") %in% class(y))) {
    ee_y <- rgee::ee$FeatureCollection(y)
    sf_y <- tryCatch(
      expr = rgee::ee_as_sf(y, quiet = FALSE, maxFeatures = 10000),
      error = function(e) {
        stop("The rgee::ee$FeatureCollection (y) must be not higher than 10 000.")
      }
    )
  }
  ee_add_rows <- function(f) {
    f_prop <- rgee::ee$Feature$get(f, "system:index")
    rgee::ee$Feature(rgee::ee$Feature$set(f, "ee_ID", f_prop))
  }
  ee_y <- rgee::ee$FeatureCollection(ee_y) %>% rgee::ee$FeatureCollection$map(ee_add_rows)
  fun_name <- gsub("Reducer.", "", (rgee::ee$Reducer$getInfo(fun))[["type"]])
  x_ic <- rgee:::bands_to_image_collection(x)
  create_tripplets <- function(img) {
    img_reduce_regions <- img$reduceRegions(
      collection = ee_y, reducer = fun, scale = scale,
      ...
    )
    rgee::ee$FeatureCollection$map(img_reduce_regions, function(f) {
      rgee::ee$Feature$set(f, "imageId", rgee::ee$Image$get(img, "system:index"))
    })
  }
  triplets <- x_ic %>%
    rgee::ee$ImageCollection$map(create_tripplets) %>%
    rgee::ee$ImageCollection$flatten()
  table <- extract_py$table_format(
    triplets, "ee_ID", "imageId",
    fun_name
  )$map(function(feature) {
    rgee::ee$Feature$setGeometry(feature, NULL)
  })
  if (via == "drive") {
    table_id <- basename(tempfile("rgee_file_"))
    ee_user <- ee_exist_credentials()
    dsn <- sprintf("%s/%s.csv", tempdir(), table_id)
    table_task <- ee_init_task_drive_fc(
      x_fc = table, dsn = dsn,
      container = container, table_id = table_id, ee_user = ee_user,
      selectors = NULL, timePrefix = TRUE, quiet = quiet
    )
    if (lazy) {
      prev_plan <- future::plan(future::sequential, .skip = TRUE)
      on.exit(future::plan(prev_plan, .skip = TRUE), add = TRUE)
      future::future(
        {
          ee_extract_to_lazy_exp_drive(
            table_task, dsn,
            quiet, sf, sf_y
          )
        },
        lazy = TRUE
      )
    } else {
      ee_extract_to_lazy_exp_drive(
        table_task, dsn, quiet,
        sf, sf_y
      )
    }
  } else if (via == "gcs") {
    table_id <- basename(tempfile("rgee_file_"))
    ee_user <- ee_exist_credentials()
    dsn <- sprintf("%s/%s.csv", tempdir(), table_id)
    table_task <- ee_init_task_gcs_fc(
      x_fc = table, dsn = dsn,
      container = container, table_id = table_id, ee_user = ee_user,
      selectors = NULL, timePrefix = TRUE, quiet = quiet
    )
    if (lazy) {
      prev_plan <- future::plan(future::sequential, .skip = TRUE)
      on.exit(future::plan(prev_plan, .skip = TRUE), add = TRUE)
      future::future(
        {
          ee_extract_to_lazy_exp_gcs(
            table_task, dsn,
            quiet, sf, sf_y
          )
        },
        lazy = TRUE
      )
    } else {
      ee_extract_to_lazy_exp_gcs(
        table_task, dsn, quiet,
        sf, sf_y
      )
    }
  } else {
    table_geojson <- table %>%
      rgee::ee$FeatureCollection$getInfo() %>%
      rgee::ee_utils_py_to_r()
    class(table_geojson) <- "geo_list"
    table_sf <- geojsonio::geojson_sf(table_geojson)
    sf::st_geometry(table_sf) <- NULL
    table_sf <- table_sf[, order(names(table_sf))]
    table_sf["id"] <- NULL
    table_sf["ee_ID"] <- NULL
    if (isTRUE(sf)) {
      table_geometry <- sf::st_geometry(sf_y)
      table_sf <- sf_y %>%
        sf::st_drop_geometry() %>%
        cbind(table_sf) %>%
        sf::st_sf(geometry = table_geometry)
    } else {
      table_sf <- sf_y %>%
        sf::st_drop_geometry() %>%
        cbind(table_sf)
    }
    table_sf
  }
}
