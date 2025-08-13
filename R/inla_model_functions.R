# functions for inla model workflow


#' Get distance from point to ebird range
#'
#' Distance will be zero inside the species' range.
#'
#' @param sp_dat combined species count and survey data
#' @param sp_code species four letter code
#' @param species_ranges A list containing ebird range data named by species code
#'
#' @returns data frame with distance_from_range column added
#' @export
#'
#' @examples
#' library(dplyr)
#' library(sf)
#' test_dat <- readRDS(system.file("extdata", "analysis_data_test.rds", package = "birdDistribution"))
#' surv_pt <- test_dat$all_surveys %>%
#'   slice(1) %>%
#'   select(Survey_Type, Survey_Duration_Minutes)
#' get_dist_to_range(surv_pt, "WTSP", test_dat$species_ranges)
get_dist_to_range <- function(sp_dat, sp_code, species_ranges) {
  # Extract ebird range for this species (if it exists)

  if (sp_code %in% names(species_ranges)) {
    range <- species_ranges[[sp_code]] %>% sf::st_transform(sf::st_crs(sp_dat))

    if (any(sf::st_geometry_type(sp_dat) == "POLYGON" |
            sf::st_geometry_type(sp_dat) == "MULTIPOLYGON")) {
      start_pt <- sf::st_centroid(sp_dat)
    } else {
      start_pt <- sp_dat
    }

    # Identify distance of each survey to the edge of species range (in km)
    sp_dat$distance_from_range <- ((sf::st_distance(start_pt, range)) %>% as.numeric()) / 1000
  } else {
    sp_dat$distance_from_range <- 0
  }
  sp_dat
}

#' Get QPAD Offsets
#'
#' Use a pre-generated table of EDR and cue rate from QPAD to calculate offsets
#'
#' @param sp_dat combined species count and survey data
#' @param sp_code species four letter code
#' @param offset_table data frame of offset data
#'
#' @returns data.frame with QPAD offset column added. Column is either
#'   `log_QPAD_offset` if `Survey_Duration_Minutes` is included in `sp_dat`,
#'   or `log_offset_5min` if not
#' @export
#'
#' @examples
#' test_dat <- readRDS(system.file("extdata", "analysis_data_test.rds", package = "birdDistribution"))
#'
#' # for a 5 min survey
#' data.frame(survey_id = 1:5) %>%
#'   get_QPAD_offsets("WTSP", test_dat$species_to_model)
#'
#' # based on survey duration
#' data.frame(survey_id = 1:5, Survey_Duration_Minutes = 1:5) %>%
#'   get_QPAD_offsets("WTSP", test_dat$species_to_model)

get_QPAD_offsets <- function(sp_dat, sp_code, offset_table) {
  # ---
  # Generate QPAD offsets for each survey (assumes unlimited distance point counts)
  # ---

  species_offsets <- subset(offset_table, Species_Code_BSC == sp_code)

  if (!species_offsets$offset_exists) {
    if (hasName(sp_dat, "Survey_Duration_Minutes")) {
      sp_dat$log_QPAD_offset <- 0
    } else {
      sp_dat$log_offset_5min <- 0
    }
  }

  if (species_offsets$offset_exists) {
    if (hasName(sp_dat, "Survey_Duration_Minutes")) {
      # Calculate offset for duration of survey from species overall offset value
      # Not using HSS because it is in the INLA model
      A_metres <- pi * species_offsets$EDR^2
      p <- 1 - exp(-sp_dat$Survey_Duration_Minutes * species_offsets$cue_rate)
      sp_dat$log_QPAD_offset <- log(A_metres * p)
    } else {
      # QPAD offsets associated with a 5-minute unlimited distance survey
      sp_dat$log_offset_5min <- species_offsets$log_offset_5min
    }
  }
  sp_dat
}

#' Combine survey and count data for a species, optionally filter the data used
#'
#' Also gets distance to range and QPAD offsets
#'
#' @param analysis_data list containing analysis data, including `full_count_matrix` and `all_surveys`
#' @param sp_code species four letter code
#' @param proj_use projection/coordinate reference system to use
#' @param train_dat_filter a string that will be used to filter the input data, the default will not filter anything
#' @param survey_types character vector of values in the `Survey_Type` column to keep
#'
#' @returns data frame of combined data
#'
#' @export
#'
#' @examples
#' # get test data set
#' test_dat <- readRDS(system.file("extdata", "analysis_data_test.rds", package = "birdDistribution"))
#' prep_sp_dat(test_dat, "WTSP", sf::st_crs(test_dat$all_survey))
#'
#' # With a filter
#' prep_sp_dat(test_dat, "WTSP", sf::st_crs(test_dat$all_survey),
#'             train_dat_filter = "Date_Time < lubridate::ymd('2003-05-01')")
#'
# prep_sp_dat <- function(analysis_data, sp_code, proj_use, train_dat_filter = "TRUE",
#                         survey_types = c("Point_Count", "ARU")) {
#
#   sp_dat <- analysis_data$all_surveys %>%
#     mutate(count = analysis_data$full_count_matrix[[sp_code]]) %>%
#
#     # select types of data, could have multiple in one model
#     subset(Survey_Type %in% survey_types) %>%
#     sf::st_transform(proj_use) %>%
#     filter(!!rlang::parse_expr(train_dat_filter))
#
#   sp_dat <- get_dist_to_range(sp_dat, sp_code, analysis_data$species_ranges)
#
#   # sp_dat <- get_QPAD_offsets(sp_dat, sp_code, analysis_data$species_to_model)
#
#   sp_dat
# }

#' Create a spatial mesh, which is used to fit the residual spatial field
#'
#' @param poly sf polygon of area to make mesh in. Eg the study area
#' @param proj_use projection/coordinate reference system to use. Note this
#'  should have units of kms to avoid an extremely dense mesh
#' @param max.edge The largest allowed triangle edge length. See [fmesher::fm_mesh_2d_inla()] for details
#' @param cutoff The minimum allowed distance between points. See [fmesher::fm_mesh_2d_inla()] for details
#'
#' @returns an inla.spde2 mesh object
#' @export
#'
#' @examples
#'
#' pg <- sf::st_polygon(list(rbind(c(0,0), c(100000,0), c(100000,100000), c(0,100000), c(0,0))))
#' mesh <- pg %>%
#'  make_mesh(proj_use = 9311)
#' plot(mesh$mesh)
#' plot(pg, add = TRUE, border = "red")
make_mesh <- function(poly, proj_use, max.edge = c(70, 100), cutoff = 30) {

  # make a two extension hulls and mesh for spatial model
  hull <- fmesher::fm_extensions(
    poly,
    convex = c(50, 200),
    concave = c(350, 500)
  )
  mesh_spatial <- fmesher::fm_mesh_2d_inla(
    boundary = hull,
    max.edge = max.edge, # km inside and outside
    cutoff = cutoff,
    crs = fmesher::fm_crs(proj_use)
  ) # cutoff is min edge

  mesh_locs <- mesh_spatial$loc[, c(1, 2)] %>% as.data.frame()
  message("Mesh created with ", dim(mesh_locs)[1], " vertices.")

  prior_range <- c(500, 0.5) # 50% chance range is smaller than 500 km
  prior_sigma <- c(0.1, 0.1) # 10% chance sd is larger than 0.1
  INLA::inla.spde2.pcmatern(mesh_spatial,
                            prior.range = prior_range,
                            prior.sigma = prior_sigma
  )
}


#' Fit INLA model to two atlas cyles (labeled as OBBA3 and OBBA2)
#'
#' @param sp_code species four letter code
#' @param analysis_data list containing analysis data, including `full_count_matrix`,
#'  `all_surveys`, `species_ranges`, and `species_to_model`
#' @param proj_use projection/coordinate reference system to use. Note this
#'  should have units of kms to avoid an extremely dense mesh
#' @param study_boundary sf polygon of study area
#' @param covariates a data frame with columns covariate, model, mean, prec, beta
#'  which is used to define the model formula.
#' @param mod_dir directory where the model object will be saved or loaded from
#'  if it already exists
#' @param train_dat_filter a string that will be used to filter the input data,
#'  the default will not filter anything
#' @param save_mod logical. Should the model object be saved?
#' @param file_name_bit suffix attached to the file name eg to identify
#'  cross-validation fold
#' @param bru_verbose level of verbosity from bru. Lower number leads to less
#'   output. See [inlabru::bru_options()]
#'
#' @note Helpful information about setting priors can be found here: https://tutorials.inbo.be/tutorials/r_inla/spatial.pdf
#'
#' @returns An INLA object with the fit model
#' @export
#'
#' @examples
#' # get test data set
#' test_dat <- readRDS(system.file("extdata", "analysis_data_test.rds", package = "birdDistribution"))
#'
#' # Using km for projection units because it helps with INLA model calculations to have smaller numbers
#' AEA_proj <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-106 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs "
#'
#' # takes awhile to run
#' if(FALSE){
#' mod <- fit_inla(
#'   sp_code = test_dat$species_to_model$Species_Code_BSC[1],
#'   analysis_data = test_dat,
#'   proj_use = AEA_proj,
#'   study_boundary = test_area,
#'   covariates = cov_df,
#'   mod_dir = tempdir(),
#'   save_mod = FALSE
#' )
#' }

fit_inla2 <- function(sp_dat,
                      study_boundary,
                      covariates,
                      error_type = "poisson",
                      prior_range_abund = c(500,0.1),  # 50% chance range is smaller than 500 km
                      prior_sigma_abund = c(0.5,0.1),  # 10% chance SD is larger than 0.5
                      prior_range_change = c(500,0.1), # 50% chance range is smaller than 500 km
                      prior_sigma_change = c(0.5,0.1), # 10% chance SD is larger than 0.5
                      bru_verbose = 4) {

  # Create spatial mesh
  hull <- fm_extensions(
    study_boundary,
    convex = c(50, 200),
    concave = c(350, 500)
  )

  mesh_spatial <- fm_mesh_2d_inla(
    boundary = hull,
    max.edge = c(50, 200), # km inside and outside
    cutoff = 10,
    crs = st_crs(sp_dat)
  )

  # Controls residual spatial field for abundance
  matern_abund <- inla.spde2.pcmatern(mesh_spatial,
                                      prior.range = prior_range_abund,
                                      prior.sigma = prior_sigma_abund,
                                      constr = TRUE
  )

  # Controls residual spatial field for change over time
  matern_change <- inla.spde2.pcmatern(mesh_spatial,
                                       prior.range = prior_range_change,
                                       prior.sigma = prior_sigma_change,
                                       constr = TRUE
  )

  # Create mesh to model effect of time since sunrise (HSS)
  sp_dat$Hours_Since_Sunrise <- as.numeric(sp_dat$Hours_Since_Sunrise)
  HSS_range <- range(sp_dat$Hours_Since_Sunrise)
  HSS_meshpoints <- seq(HSS_range[1] - 1, HSS_range[2] + 1, length.out = 21)
  HSS_mesh1D <- INLA::inla.mesh.1d(HSS_meshpoints, boundary = "free")
  HSS_spde <- INLA::inla.spde2.pcmatern(HSS_mesh1D,
                                        prior.range = c(5, 0.1),
                                        prior.sigma = c(2, 0.1)
  )

  # iid random effect for atlas squares
  pc_prec <- list(prior = "pcprec", param = c(0.1, 0.1))

  # Model formulas
  covariates <- covariates %>%
    mutate(
      components = paste0(
        "Beta", beta, "_", covariate, '(1,model=\"', model,
        '\"', ", mean.linear = ", mean, ", prec.linear = ",
        prec, ")"
      ),
      formula = paste0("Beta", beta, "_", covariate, "*", covariate, "^", beta)
    )

  model_components <- as.formula(paste0(
    '~
            Intercept_OBBA2(1)+
            Intercept_OBBA3(1)+
            HSS(main = Hours_Since_Sunrise,model = HSS_spde) +
            range_effect(1,model="linear", mean.linear = -0.046, prec.linear = 10000)+
            kappa(square_atlas, model = "iid", constr = TRUE, hyper = list(prec = pc_prec)) +
            spde_abund(main = geometry, model = matern_abund) +
            spde_change(main = geometry, model = matern_change) +',
    paste0(covariates$components, collapse = " + ")
  ))

  model_formula_OBBA2 <- as.formula(paste0("count ~
                  Intercept_OBBA2 +
                  HSS +
                  kappa +
                  range_effect * distance_from_range +
                  spde_abund +",
                                           paste0(covariates$formula, collapse = " + ")
  ))

  model_formula_OBBA3 <- as.formula(paste0("count ~
                  Intercept_OBBA3 +
                  HSS +
                  kappa +
                  range_effect * distance_from_range +
                  spde_abund + spde_change +",
                                           paste0(covariates$formula, collapse = " + ")
  ))

  # Fit model to both atlas periods
  start <- Sys.time()
  fit_INLA <- NULL
  while (is.null(fit_INLA)) {
    fit_INLA <- inlabru::bru(
      components = model_components,

      inlabru::like(
        family = error_type,
        formula = model_formula_OBBA2,
        data = subset(sp_dat, Atlas == "OBBA2")
      ),

      inlabru::like(
        family = error_type,
        formula = model_formula_OBBA3,
        data = subset(sp_dat, Atlas == "OBBA3")
      ),

      options = list(inla.mode = "experimental",
                     control.compute = list(waic = FALSE, cpo = FALSE),
                     bru_verbose = bru_verbose
      )
    )
    if ("try-error" %in% class(fit_INLA)) fit_INLA <- NULL
  }

  end <- Sys.time()
  runtime_INLA <- difftime(end, start, units = "mins") %>% round(2)
  message(paste0(sp_code, " - ", runtime_INLA, " min to fit model"))

  return(fit_INLA)
}

summarize_posterior <- function(mat, CI_probs = c(0.05, 0.95), prefix = "var") {
  stopifnot(is.matrix(mat))

  mean_vals   <- matrixStats::rowMeans2(mat, na.rm = TRUE)
  median_vals <- matrixStats::rowMedians(mat, na.rm = TRUE)
  sd_vals     <- matrixStats::rowSds(mat, na.rm = TRUE)
  cv_vals     <- sd_vals / median_vals

  lower_vals <- matrixStats::rowQuantiles(mat, probs = CI_probs[1], na.rm = TRUE)
  upper_vals <- matrixStats::rowQuantiles(mat, probs = CI_probs[2], na.rm = TRUE)

  out <- data.frame(
    setNames(list(mean_vals),   paste0(prefix, "_mean")),
    setNames(list(median_vals), paste0(prefix, "_q50")),
    setNames(list(sd_vals),     paste0(prefix, "_sd")),
    setNames(list(cv_vals),     paste0(prefix, "_cv_median")),
    setNames(list(lower_vals),  paste0(prefix, "_lower")),
    setNames(list(upper_vals),  paste0(prefix, "_upper"))
  )

  return(out)
}

#' Use a fit INLA model to generate predictions
#'
#' Generate predictions from a predictor data set and a fitted INLA model.
#' Predictions are summarised using the median, credible interval, and
#' coefficient of variation. The Continuous Ranked Probability Score and the Log
#' score are optionally calculated.
#'
#' @param dat data to make predictions from
#' @param analysis_data list containing analysis data, including `full_count_matrix`,
#'  `all_surveys`, `species_ranges`, and `species_to_model`
#' @param mod fitted INLA model
#' @param sp_code species four letter code
#' @param covariates a data frame with columns covariate, model, mean, prec, beta
#'  which is used to define the model formula.
#' @param do_crps logical. Should the Continuous Ranked Probability Score and
#'  the Log score be calculated from the raw predictions?
#'
#' @returns `dat` with columns added for predictions.
#' @export
#'
#' @examples
#'
#' # get test data set
#' test_dat <- readRDS(system.file("extdata", "analysis_data_test.rds", package = "birdDistribution"))
#'
#' # Using km for projection units because it helps with INLA model calculations to have smaller numbers
#' AEA_proj <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-106 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs "
#'
#' # takes awhile to run
#' if(FALSE){
#'   mod <- fit_inla(
#'     sp_code = test_dat$species_to_model$Species_Code_BSC[1],
#'     analysis_data = test_dat,
#'     proj_use = AEA_proj,
#'     study_boundary = test_area,
#'     covariates = cov_df,
#'     mod_dir = tempdir(),
#'     save_mod = FALSE
#'   )
#'
#'   pred <- predict_inla(
#'     dat = test_dat$ONGrid,
#'     analysis_data = test_dat,
#'     mod = mod,
#'     sp_code = test_dat$species_to_model$Species_Code_BSC[1],
#'     covariates = cov_df,
#'     do_crps = FALSE
#'   )
#' }

predict_inla <- function(mod, grid, pred_formula) {

  start <- Sys.time()

  pred <- inlabru::generate(mod,
                            grid,
                            formula = pred_formula,
                            n.samples = 1000,
                            seed = 123
  )

  # Reformat pred such to a named list, with prediction matrices for each object (n_grid_cells x n_samples)
  pred_vars <- names(pred[[1]])
  preds <- lapply(pred_vars, function(v) sapply(pred, function(x) x[[v]]))
  names(preds) <- pred_vars

  # Convert to count scale
  preds$OBBA3 <- exp(preds$OBBA3)
  preds$OBBA2 <- exp(preds$OBBA2)
  preds$log_change <- log(preds$OBBA3) - log(preds$OBBA2)
  preds$pct_change <- (exp(preds$log_change) - 1) * 100
  preds$abs_change <- preds$OBBA3 - preds$OBBA2

  end <- Sys.time()

  runtime_pred <- difftime(end, start, units = "mins") %>% round(2)
  message(paste0(sp_code, " - ", runtime_pred, " min to generate predictions"))
  return(posterior_list)
}

#' Make maps of INLA model predictions
#'
#' Makes several maps from different prediction outputs
#'
#' @param sp_code species four letter code
#' @param analysis_data list containing analysis data, including `full_count_matrix`,
#'  `all_surveys`, `species_ranges`, and `species_to_model`
#' @param preds data frame of predictions
#' @param proj_use projection/coordinate reference system to use.
#' @param atlas_squares grid of squares to show predictions and observations in.
#' @param bcr_poly polygon of BCR boundaries to use in map.
#' @param study_boundary sf polygon of study area.
#' @param target_raster raster with desired structure eg resolution, crs etc
#' @param map_dir directory where the map images should be saved
#' @param train_dat_filter a string that will be used to filter the input data,
#'  the default will not filter anything
#' @param file_name_bit suffix attached to the file name eg to identify
#'  cross-validation fold
#'
#' @returns Saves maps to output `map_dir`
#' @export
#'
#' @examples
#' test_dat <- readRDS(system.file("extdata", "analysis_data_test.rds", package = "birdDistribution"))
#'
#' # Using km for projection units because it helps with INLA model calculations to have smaller numbers
#' AEA_proj <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-106 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs "
#'
#' # takes awhile to run
#' if(FALSE){
#'   mod <- fit_inla(
#'     sp_code = test_dat$species_to_model$Species_Code_BSC[1],
#'     analysis_data = test_dat,
#'     proj_use = AEA_proj,
#'     study_boundary = test_area,
#'     covariates = cov_df,
#'     mod_dir = tempdir(),
#'     save_mod = FALSE
#'   )
#'
#'   pred <- predict_inla(
#'     dat = test_dat$ONGrid,
#'     analysis_data = test_dat,
#'     mod = mod,
#'     sp_code = test_dat$species_to_model$Species_Code_BSC[1],
#'     covariates = cov_df,
#'     do_crps = FALSE
#'   )
#'
#'   map_raster_out <- terra::rast(
#'     resolution = 10, crs = AEA_proj,
#'     extent = terra::ext(test_area %>% terra::vect()),
#'     vals = 1
#'   )
#'
#'   dir_use <- tempdir()
#'   map_inla_preds(
#'     sp_code = test_dat$species_to_model$Species_Code_BSC[1],
#'     analysis_data = test_dat,
#'     pred,
#'     proj_use = AEA_proj,
#'     study_boundary = test_area,
#'     atlas_squares = atl_sq %>% sf::st_transform(AEA_proj),
#'     bcr_poly = bcr_poly,
#'     target_raster = map_raster_out,
#'     map_dir = dir_use
#'   )
#'
#'   out_maps <- list.files(dir_use, pattern = "png$", full.names = TRUE)
#' }

map_relabund <- function(species_name,
                         obs_dat,
                         grid,
                         preds_summarized,
                         atlas_squares,
                         study_boundary,
                         map_dir = "figures/species_maps/",
                         train_dat_filter = "TRUE",
                         prefix = "OBBA3",
                         plot_obs_data = TRUE,
                         title = "Relative Abundance",
                         subtitle = "Per 5-minute point count",
                         upper_bound = 1,
                         lower_bound = 0.01,
                         res = 1.1) {

  proj_use <- st_crs(obs_dat)

  # Helper to split species name if it's too long
  wrap_species_label <- function(label, max_length = 15) {
    if (nchar(label) <= max_length) return(label)

    words <- strsplit(label, " ")[[1]]
    if (length(words) == 1) return(label)  # Single word, don't split

    # Put everything except the last word on the first line
    paste0(paste(words[-length(words)], collapse = " "), "<br>", words[length(words)])
  }

  # Summarize atlas_squares where species was detected
  sp_detected <- obs_dat %>%
    sf::st_intersection(atlas_squares %>% st_transform(st_crs(obs_dat))) %>%
    as.data.frame() %>%
    group_by(square_id_) %>%
    summarize(
      sp_detected = as.numeric(sum(count) > 0),
      sp_mean_count = mean(count) %>% round(2)
    )

  atlas_squares_species <- atlas_squares %>%
    relocate(geometry, .after = last_col()) %>%
    left_join(sp_detected, by = join_by(square_id_)) # %>% left_join(CL_detected)

  atlas_squares_centroids <- sf::st_centroid(atlas_squares_species)

  # ---- Bind summarized predictions to grid
  q50_col <- names(preds_summarized)[endsWith(names(preds_summarized), "_q50")]
  grid$pred_q50 <- preds_summarized[[q50_col]]

  CV_col <- names(preds_summarized)[endsWith(names(preds_summarized), "_cv_median")]
  grid$pred_CV <- preds_summarized[[CV_col]]

  # ---- Plot median predictions

  # Bounds for plotting and associated labels
  breaks <- 10^seq(log10(lower_bound),log10(upper_bound),length.out = 5) %>% signif(2)

  # Cap values at upper/lower bounds
  grid$pred_capped <- as.numeric(pmax(pmin(grid$pred_q50, upper_bound), lower_bound))

  # Convert sf to SpatVector
  v <- terra::vect(grid)

  # Create a raster template with desired resolution
  r_template <- terra::rast(v, res = res)

  # Rasterize pred_capped values using mean within each cell
  pred_rast <- terra::rasterize(v, r_template, field = "pred_capped", fun = mean)

  # Convert raster to stars for plotting
  pred_rast_stars <- stars::st_as_stars(pred_rast)

  # Set legend labels
  break_labels <- as.character(breaks)
  break_labels[1] <- paste0("<", break_labels[1])
  break_labels[length(break_labels)] <- paste0(">", break_labels[length(break_labels)])

  colscale_q50 <- c(
    "#FEFEFE", "#FBF7E2", "#FCF8D0", "#EEF7C2", "#CEF2B0",
    "#94E5A0", "#51C987", "#18A065", "#008C59", "#007F53", "#006344"
  )
  colpal_q50 <- colorRampPalette(colscale_q50)

  q50_plot <- ggplot() +
    stars::geom_stars(data = pred_rast_stars) +
    scale_fill_gradientn(
      name = paste0(
        "<span style='font-size:20pt; font-weight:bold'>", wrap_species_label(species_name), "</span><br><br>",
        "<span style='font-size:14pt'>", title, "</span><br>",
        "<span style='font-size:7pt'>", subtitle, "</span><br>",
        "<span style='font-size:7pt'>Posterior Median</span>"
      ),
      colors = colpal_q50(10),
      trans = "log10",
      na.value = "transparent",
      breaks = breaks,
      labels = break_labels,
      limits = c(min(breaks) / 1.1, max(breaks) * 1.1)
    ) +
    geom_sf(data = atlas_squares_centroids %>% subset(!is.na(sp_mean_count)), colour = "gray50", size = 0.5, shape = 1, stroke = 0.1) +
    geom_sf(data = atlas_squares_centroids %>% subset(sp_mean_count>0), colour = "black", size = 0.5, shape = 1, stroke = 0.2) +
    geom_sf(data = study_boundary, colour = "black", fill = NA, lwd = 0.3, show.legend = FALSE) +
    coord_sf(clip = "off") +
    theme_void() +
    theme(
      plot.margin = unit(c(0, 0, 0, 0), "cm"),
      legend.title = ggtext::element_markdown(lineheight = .9),
      legend.position = c(1,0.9),
      legend.justification = c(1,1),
      legend.background = element_rect(fill = "transparent", color = "transparent")
    )

  png(paste0(map_dir,"/",species_name,"_",prefix,"_relabund_q50.png"), width = 10, height = 8, units = "in", res = 1000, type = "cairo")
  print(q50_plot)
  dev.off()

  # ---- Plot uncertainty in predictions (width of 90% CRI)

  colscale_uncertainty <- c("#FEFEFE", "#FFF4B3", "#F5D271", "#F2B647", "#EC8E00", "#CA302A")
  colpal_uncertainty <- colorRampPalette(colscale_uncertainty)

  # Cap values at upper/lower bounds
  grid$pred_capped <- as.numeric(pmax(pmin(grid$pred_CV, 1), 0))

  # Convert sf to SpatVector
  v <- terra::vect(grid)

  # Create a raster template with desired resolution
  r_template <- terra::rast(v, res = res)

  # Rasterize pred_capped values using mean within each cell
  pred_rast <- terra::rasterize(v, r_template, field = "pred_capped", fun = mean)

  # Convert raster to stars for plotting
  pred_rast_stars <- stars::st_as_stars(pred_rast)

  # Set legend labels/breaks
  breaks <- seq(0,1,length.out = 5) %>% signif(2)
  break_labels <- as.character(breaks)
  break_labels[length(break_labels)] <- paste0(">", break_labels[length(break_labels)])

  CV_plot <- ggplot() +
    stars::geom_stars(data = pred_rast_stars) +
    scale_fill_gradientn(
      name = paste0(
        "<span style='font-size:20pt; font-weight:bold'>", wrap_species_label(species_name), "</span><br><br>",
        "<span style='font-size:14pt'>", title, "</span><br>",
        "<span style='font-size:7pt'>", subtitle, "</span><br>",
        "<span style='font-size:7pt'>Width of 90% CRI</span>"
      ),
      colors = colpal_uncertainty(11),
      na.value = "transparent",
      breaks = breaks,
      labels = break_labels,
      limits = c(min(breaks) / 1.1, max(breaks) * 1.1)
    ) +
    geom_sf(data = atlas_squares_centroids %>% subset(!is.na(sp_mean_count)), colour = "gray50", size = 0.5, shape = 1, stroke = 0.1) +
    geom_sf(data = atlas_squares_centroids %>% subset(sp_mean_count>0), colour = "black", size = 0.5, shape = 1, stroke = 0.2) +
    geom_sf(data = study_boundary, colour = "black", fill = NA, lwd = 0.3, show.legend = FALSE) +
    coord_sf(clip = "off") +
    theme_void() +
    theme(
      plot.margin = unit(c(0, 0, 0, 0), "cm"),
      legend.title = ggtext::element_markdown(lineheight = .9),
      legend.position = c(1,0.9),
      legend.justification = c(1,1),
      legend.background = element_rect(fill = "transparent", color = "transparent")
    )

  png(paste0(map_dir,"/",species_name,"_",prefix,"_relabund_CV.png"), width = 10, height = 8, units = "in", res = 1000, type = "cairo")
  print(CV_plot)
  dev.off()

}

map_change <- function(species_name,
                       grid,
                       preds_summarized,
                       study_boundary,
                       map_dir = "figures/species_maps/",
                       upper_bound = -1,
                       lower_bound = 1,
                       res = 1.1,
                       change_type = "Percent") {

  proj_use <- st_crs(obs_dat)

  # Helper to split species name if it's too long
  wrap_species_label <- function(label, max_length = 15) {
    if (nchar(label) <= max_length) return(label)

    words <- strsplit(label, " ")[[1]]
    if (length(words) == 1) return(label)  # Single word, don't split

    # Put everything except the last word on the first line
    paste0(paste(words[-length(words)], collapse = " "), "<br>", words[length(words)])
  }

  # ---- Bind summarized predictions to grid
  q50_col <- names(preds_summarized)[endsWith(names(preds_summarized), "_q50")]
  grid$pred_q50 <- preds_summarized[[q50_col]]

  CV_col <- names(preds_summarized)[endsWith(names(preds_summarized), "_cv_median")]
  grid$pred_CV <- preds_summarized[[CV_col]]

  # ---- Plot median predictions

  grid$pred_capped <- as.numeric(pmax(pmin(grid$pred_q50, upper_bound), lower_bound))

  # Convert sf to SpatVector
  v <- terra::vect(grid)

  # Create a raster template with desired resolution
  r_template <- terra::rast(v, res = res)

  # Rasterize pred_capped values using mean within each cell
  pred_rast <- terra::rasterize(v, r_template, field = "pred_capped", fun = mean)

  # Convert raster to stars for plotting
  pred_rast_stars <- stars::st_as_stars(pred_rast)

  # Bounds for plotting and associated labels
  breaks <- seq(lower_bound,upper_bound,length.out = 7)

  # Set legend labels
  if (change_type == "Percent"){
    break_labels <- (100 * (exp(breaks) - 1)) %>% signif(2)
    break_labels <- paste0(break_labels,"%")
    prefix = "pct_change"
    title = "Percent change"
  } else{
    break_labels = breaks %>% signif(2)
    prefix = "abs_change"
    title = "Absolute change"
  }
  break_labels[breaks>0] <- paste0("+",break_labels[breaks>0])
  break_labels[1] <- paste0("< ", break_labels[1])
  break_labels[length(break_labels)] <- paste0("> ", break_labels[length(break_labels)])

  colscale_q50 <- RColorBrewer::brewer.pal(11,"RdBu")
  colpal_q50 <- colorRampPalette(colscale_q50)

  q50_plot <- ggplot() +
    stars::geom_stars(data = pred_rast_stars) +
    scale_fill_gradientn(
      name = paste0(
        "<span style='font-size:20pt; font-weight:bold'>", wrap_species_label(species_name), "</span><br><br>",
        "<span style='font-size:14pt'>", title, "</span><br>",
        "<span style='font-size:7pt'>OBBA2 to OBBA3</span><br>",
        "<span style='font-size:7pt'>Posterior Median</span>"
      ),
      colors = colpal_q50(11),
      na.value = "transparent",
      breaks = breaks,
      labels = break_labels,
      limits = c(min(breaks) * 1.1, max(breaks) * 1.1)
    ) +
    geom_sf(data = study_boundary, colour = "black", fill = NA, lwd = 0.3, show.legend = FALSE) +
    coord_sf(clip = "off") +
    theme_void() +
    theme(
      plot.margin = unit(c(0, 0, 0, 0), "cm"),
      legend.title = ggtext::element_markdown(lineheight = .9),
      legend.position = c(1,0.9),
      legend.justification = c(1,1),
      legend.background = element_rect(fill = "transparent", color = "transparent")
    )

  png(paste0(map_dir,"/",species_name,"_",prefix,"_q50.png"), width = 10, height = 8, units = "in", res = 1000, type = "cairo")
  print(q50_plot)
  dev.off()

  # # ---- Plot uncertainty in predictions (width of 90% CRI)
  #
  # colscale_uncertainty <- c("#FEFEFE", "#FFF4B3", "#F5D271", "#F2B647", "#EC8E00", "#CA302A")
  # colpal_uncertainty <- colorRampPalette(colscale_uncertainty)
  #
  # # Cap values at upper/lower bounds
  # grid$pred_capped <- as.numeric(pmax(pmin(grid$pred_CV, 1), 0))
  #
  # # Convert sf to SpatVector
  # v <- terra::vect(grid)
  #
  # # Create a raster template with desired resolution
  # r_template <- terra::rast(v, res = res)
  #
  # # Rasterize pred_capped values using mean within each cell
  # pred_rast <- terra::rasterize(v, r_template, field = "pred_capped", fun = mean)
  #
  # # Convert raster to stars for plotting
  # pred_rast_stars <- stars::st_as_stars(pred_rast)
  #
  # # Set legend labels/breaks
  # breaks <- seq(0,1,length.out = 5) %>% signif(2)
  # break_labels <- as.character(breaks)
  # break_labels[length(break_labels)] <- paste0(">", break_labels[length(break_labels)])
  #
  # CV_plot <- ggplot() +
  #   stars::geom_stars(data = pred_rast_stars) +
  #   scale_fill_gradientn(
  #     name = paste0(
  #       "<span style='font-size:20pt; font-weight:bold'>", wrap_species_label(species_name), "</span><br><br>",
  #       "<span style='font-size:14pt'>", title, "</span><br>",
  #       "<span style='font-size:7pt'>", subtitle, "</span><br>",
  #       "<span style='font-size:7pt'>Width of 90% CRI</span>"
  #     ),
  #     colors = colpal_uncertainty(10),
  #     na.value = "transparent",
  #     breaks = breaks,
  #     labels = break_labels,
  #     limits = c(min(breaks) / 1.1, max(breaks) * 1.1)
  #   ) +
  #   geom_sf(data = atlas_squares_centroids %>% subset(!is.na(sp_mean_count)), colour = "gray50", size = 0.5, shape = 1, stroke = 0.1) +
  #   geom_sf(data = atlas_squares_centroids %>% subset(sp_mean_count>0), colour = "black", size = 0.5, shape = 1, stroke = 0.2) +
  #   geom_sf(data = study_boundary, colour = "black", fill = NA, lwd = 0.3, show.legend = FALSE) +
  #   coord_sf(clip = "off") +
  #   theme_void() +
  #   theme(
  #     plot.margin = unit(c(0, 0, 0, 0), "cm"),
  #     legend.title = ggtext::element_markdown(lineheight = .9),
  #     legend.position = c(1,0.9),
  #     legend.justification = c(1,1),
  #     legend.background = element_rect(fill = "transparent", color = "transparent")
  #   )
  #
  # png(paste0(map_dir,"/",species_name,"_",prefix,"_relabund_CV.png"), width = 10, height = 8, units = "in", res = 1000, type = "cairo")
  # print(CV_plot)
  # dev.off()

}

#' Build INLA prediction map
#'
#' @param pred_rast Raster of predictions
#' @param title map title describing the metric shown
#' @param subtitle subtitle under title with more details
#' @param subsubtitle subtitle under subtitle eg with units or further description
#' @param samp_grid Points in sampling grid where species was or wasn't detected
#' @param bcr_poly sf polygon for BCR to be added to plot
#' @param col_pal_fn function to create colour palette
#' @param species_label species name which will be added to the map
#' @param levs_nm name of levels of colour ramp
#' @param file_nm file name where map should be saved
#'
#' @returns saves the map to `file_nm`
#'

do_res_plot <- function(preds,
                        species_name,
                        title,
                        subtitle,
                        subsubtitle = "",
                        study_boundary,
                        col_pal_fn,
                        breaks,
                        res = 1,
                        lower_bound = 0.01,
                        upper_bound = 1) {

  # Helper to split species name if it's too long
  wrap_species_label <- function(label, max_length = 15) {
    if (nchar(label) <= max_length) return(label)

    words <- strsplit(label, " ")[[1]]
    if (length(words) == 1) return(label)  # Single word, don't split

    # Put everything except the last word on the first line
    paste0(paste(words[-length(words)], collapse = " "), "<br>", words[length(words)])
  }

  # Cap values at upper/lower bounds
  preds$pred_capped <- pmax(pmin(preds$pred_q50, upper_bound), lower_bound)

  # Set legend labels
  break_labels <- as.character(breaks)
  break_labels[1] <- paste0("<", break_labels[1])
  break_labels[length(break_labels)] <- paste0(">", break_labels[length(break_labels)])

  # Convert sf to SpatVector
  v <- terra::vect(preds)

  # Create a raster template with desired resolution
  r_template <- terra::rast(v, res = res)

  # Rasterize pred_capped values using mean within each cell
  pred_rast <- terra::rasterize(v, r_template, field = "pred_capped", fun = mean)

  # Convert raster to stars for plotting
  pred_rast_stars <- stars::st_as_stars(pred_rast)

  res_plot <- ggplot() +
    stars::geom_stars(data = pred_rast_stars) +
    scale_fill_gradientn(
      name = paste0(
        "<span style='font-size:20pt; font-weight:bold'>", wrap_species_label(species_name), "</span><br><br>",
        "<span style='font-size:14pt'>", title, "</span><br>",
        "<span style='font-size:7pt'>", subtitle, "</span><br>",
        "<span style='font-size:7pt'>", subsubtitle, "</span>"
      ),
      colors = col_pal_fn(10),
      trans = "log10",
      na.value = "transparent",
      breaks = breaks,
      labels = break_labels,
      limits = c(min(breaks) / 1.1, max(breaks) * 1.1)
    ) +
    geom_sf(data = study_boundary, colour = "black", fill = NA, lwd = 0.3, show.legend = FALSE) +
    coord_sf(clip = "off") +
    theme_void() +
    theme(
      plot.margin = unit(c(0, 0, 0, 0), "cm"),
      legend.title = ggtext::element_markdown(lineheight = .9),
      legend.position = c(1,0.7),
      legend.justification = c(1,1),
      legend.background = element_rect(fill = "transparent", color = "transparent")
    )

  # Save and return
  png(file_nm, width = 10, height = 8, units = "in", res = 1000, type = "cairo")
  print(res_plot)
  dev.off()

  return(res_plot)
}

# do_res_plot <- function(preds, title, subtitle, subsubtitle = "", samp_grid, study_boundary,
#                         col_pal_fn, species_label, breaks, file_nm) {
#
#   break_labels <- as.character(breaks)
#   break_labels[1] <- paste0("<",break_labels[1])
#   break_labels[length(break_labels)] <-  paste0(">",break_labels[length(break_labels)])
#
#   preds$pred_q50[preds$pred_q50 > upper_bound] <- upper_bound
#   preds$pred_q50[preds$pred_q50 < lower_bound] <- 0
#
#   res_plot <- ggplot2::ggplot() +
#     ggplot2::geom_sf(data = preds, aes(col = pred_q50), size = 0.1) +
#     ggplot2::scale_color_gradientn(
#       name = paste0("<span style='font-size:13pt'>", title,
#                     "</span><br><span style='font-size:7pt'>", subtitle,
#                     "</span><br><span style='font-size:7pt'>", subsubtitle,
#                     "</span>"),
#       colors = colpal_relabund(10),
#       trans = "log10",
#       na.value = "black",
#       breaks = breaks,
#       labels = break_labels,
#       limits = c(min(breaks)/1.1,max(breaks)*1.1))+
#
#     ggplot2::geom_sf(data = study_boundary,colour="black",fill=NA,lwd=0.3,show.legend = F) +
#
#     ggplot2::coord_sf(clip = "off",xlim = range(as.data.frame(st_coordinates(ONBoundary))$X)) +
#     ggplot2::theme(panel.background = element_blank(),
#                    panel.grid.major = element_blank(),
#                    panel.grid.minor = element_blank(),
#                    axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),
#                    plot.margin = unit(c(0, 0, 0, 0), "cm"))+
#
#     ggplot2::annotate(geom="text",x=400,y=1800, label= paste0(species_name),lineheight = .85,hjust = 0,size=6,fontface =2) +
#     ggplot2::annotate(geom="text",x=410,y=1700, label= "OBBA3",lineheight = .85,hjust = 0,size=5,fontface =2)
#
#   png(file_nm, width = 10, height = 6.5, units = "in", res = 1000, type = "cairo")
#   print(res_plot)
#   dev.off()
#
#   return(res_plot)
# }


#' Rasterize a series of spatial predictions (needed for plotting)
#'
#' @param df dataframe of predictions
#' @param target_raster raster with desired structure eg resolution, crs etc
#' @param column_name column in `df` to rasterize
#' @param lower_bound,upper_bound upper and lower bounds for labels

cut_fn <- function(df = NA,
                   target_raster = NA,
                   column_name = NA,
                   lower_bound = NA,
                   upper_bound = NA) {
  max_val <- upper_bound
  max_val <- ifelse(is.na(max_val), 0, max_val)
  max_lev <- ifelse(max_val > 1.6, 4, ifelse(max_val > 0.8, 4, 3))

  cut_levs <- signif(max_val / (2^((max_lev - 1):0)), 2)
  cut_levs <- unique(cut_levs)
  cut_levs <- ifelse(is.na(cut_levs), 0, cut_levs)

  if (lower_bound %in% cut_levs) cut_levs <- cut_levs[-which(cut_levs == lower_bound)]
  if (lower_bound > min(cut_levs)) cut_levs <- cut_levs[-which(cut_levs < lower_bound)]

  max_lev <- length(cut_levs)

  cut_levs_labs <- c(
    paste0("0-", lower_bound),
    paste(lower_bound, cut_levs[1], sep = "-"),
    paste(cut_levs[-max_lev], cut_levs[-1], sep = "-"),
    paste(cut_levs[max_lev], "+")
  )

  cut_levs <- c(-1, lower_bound, cut_levs, 1000) %>% unique()

  df <- mutate(df, levs = cut(.data[[column_name]], cut_levs, labels = cut_levs_labs))
  # TODO: change this to use terra. It is faster now and easier I think
  tgt <- stars::st_as_stars(target_raster)
  tmp <- stars::st_rasterize(df %>% dplyr::select(levs, geometry),
                             nx = dim(tgt)[1], ny = dim(tgt)[2]
  )

  return(list(raster = tmp, cut_levs = cut_levs))
}


#' Evaluate model performance based on predicted vs observed count
#'
#' Compare predicted median abundance to observed count. Also summarise CRPS and log score if available
#'
#' @param pred prediction data frame
#' @param mod fitted model
#' @param sp_code species code
#' @param analysis_data list containing analysis data, including `full_count_matrix`
#'
#' @returns data frame of performance metrics
#' @export
#'
#' @examples
#'
#' test_dat <- readRDS(system.file("extdata", "analysis_data_test.rds", package = "birdDistribution"))
#'
#' test_dat$all_surveys$Obs_Index <- 1:nrow(test_dat$all_surveys)
#'
#' # Using km for projection units because it helps with INLA model calculations to have smaller numbers
#' AEA_proj <- "+proj=aea +lat_1=50 +lat_2=70 +lat_0=40 +lon_0=-106 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs "
#'
#' # takes awhile to run
#' if(FALSE){
#' # train on old data
#'   mod <- fit_inla(
#'     sp_code = test_dat$species_to_model$Species_Code_BSC[1],
#'     analysis_data = test_dat,
#'     proj_use = AEA_proj,
#'     study_boundary = test_area,
#'     covariates = cov_df,
#'     mod_dir = tempdir(),
#'     save_mod = FALSE,
#'     train_dat_filter = "Date_Time < lubridate::ymd('2010-05-01')"
#'   )
#' # Predict on new data
#'   pred <- predict_inla(
#'     dat = test_dat$all_surveys %>% filter(Date_Time >= lubridate::ymd('2010-05-01')),
#'     analysis_data = test_dat,
#'     mod = mod,
#'     sp_code = test_dat$species_to_model$Species_Code_BSC[1],
#'     covariates = cov_df,
#'     do_crps = TRUE
#'   )
#'
#' evaluate_preds(pred %>% mutate(Crossval_Fold = 1), mod,
#'                sp_code = test_dat$species_to_model$Species_Code_BSC[1],
#'                test_dat)
#' }
#'
evaluate_preds <- function(pred, mod, sp_code, analysis_data) {
  obs_count <- analysis_data$full_count_matrix[pred$Obs_Index, sp_code]

  rmse_pred <- sqrt(mean((obs_count - pred$pred_q50)^2))

  size <- mod$summary.hyperpar$"0.5quant"[1]

  auc_pred <- pROC::auc(response = obs_count > 0, predictor = pred$pObs_5min)

  # log pointwise predictive density, I think...
  lppd_pred <- sum(dnbinom(obs_count, mu = pred$pred_q50, size = size, log = TRUE))

  med_crps <- median(pred$crps)

  med_crps_pres <- median(pred$crps[obs_count > 0])
  med_crps_abs <- median(pred$crps[obs_count == 0])

  med_logs <- median(pred$logs)

  med_logs_pres <- median(pred$logs[obs_count > 0])
  med_logs_abs <- median(pred$logs[obs_count == 0])

  tibble::tibble(
    species = sp_code, fold = unique(pred$Crossval_Fold),
    rmse = rmse_pred, auc = as.numeric(auc_pred), lppd = lppd_pred,
    med_crps = med_crps, med_logs = med_logs,
    med_crps_pres = med_crps_pres, med_logs_pres = med_logs_pres,
    med_crps_abs = med_crps_abs, med_logs_abs = med_logs_abs
  )
}
