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
#' surv_pt <- test_dat$all_surveys %>% slice(1) %>% select(Survey_Type, Survey_Duration_Minutes)
#' get_dist_to_range(surv_pt, "WTSP", test_dat$species_ranges)
get_dist_to_range <- function(sp_dat, sp_code, species_ranges){
  # Extract ebird range for this species (if it exists)

  if (sp_code %in% names(species_ranges)){

    range <- species_ranges[[sp_code]] %>% sf::st_transform(sf::st_crs(sp_dat))

    if(any(sf::st_geometry_type(sp_dat) == "POLYGON" |
       sf::st_geometry_type(sp_dat) == "MULTIPOLYGON")){
      start_pt <- sf::st_centroid(sp_dat)
    } else {
      start_pt <- sp_dat
    }

    # Identify distance of each survey to the edge of species range (in km)
    sp_dat$distance_from_range <- ((sf::st_distance(start_pt, range)) %>% as.numeric())/1000
  } else{
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
#' @returns
#' @export
#'
#' @examples
get_QPAD_offsets <- function(sp_dat, sp_code, offset_table){
  # ---
  # Generate QPAD offsets for each survey (assumes unlimited distance point counts)
  # ---

  species_offsets <- subset(offset_table, Species_Code_BSC == sp_code)

  if (!species_offsets$offset_exists){
    if(hasName(sp_dat, "Survey_Duration_Minutes")){
      sp_dat$log_QPAD_offset <- 0
    } else {
      sp_dat$log_offset_5min <- 0
    }
  }

  if (species_offsets$offset_exists){
    if(hasName(sp_dat, "Survey_Duration_Minutes")){
      # Calculate offset for duration of survey from species overall offset value
      # Not using TSS because it is in the INLA model
      A_metres <- pi*species_offsets$EDR^2
      p <- 1-exp(-sp_dat$Survey_Duration_Minutes*species_offsets$cue_rate)
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
#' @param analysis_data list containing analysis data, including `full_count_matrix` and `all_surveys`
#' @param sp_code species four letter code
#' @param proj_use projection/coordinate reference system to use
#' @param train_dat_filter a string that will be used to filter the input data, the default will not filter anything
#' @param survey_types character vector of values in the `Survey_Type` column to keep
#'
#' @returns
#'
#' @examples
prep_sp_dat <- function(analysis_data, sp_code, proj_use, train_dat_filter = "TRUE",
                        survey_types = c("Point_Count","ARU_SPT","ARU_SPM")){
  sp_dat <- analysis_data$all_surveys %>%
    mutate(count = analysis_data$full_count_matrix[,sp_code]) %>%
    # select types of data, could have multiple in one model
    subset(Survey_Type %in% survey_types) %>%
    sf::st_transform(proj_use) %>%
    filter(!!rlang::parse_expr(train_dat_filter))

  sp_dat <- get_dist_to_range(sp_dat, sp_code, analysis_data$species_ranges)

  sp_dat <- get_QPAD_offsets(sp_dat, sp_code, analysis_data$species_to_model)

  sp_dat
}

#' Create a spatial mesh, which is used to fit the residual spatial field
#'
#' @param poly sf polygon of area to make mesh in. Eg the study area
#' @param proj_use projection/coordinate reference system to use. Note this
#'  should have units of kms to avoid an extremely dense mesh
#' @param max.edge The largest allowed triangle edge length. See [fmesher::fm_mesh_2d_inla()] for details
#' @param cutoff The minimum allowed distance between points. See [fmesher::fm_mesh_2d_inla()] for details
#'
#' @returns
#' @export
#'
#' @examples
make_mesh <- function(poly, proj_use, max.edge = c(70000, 100000), cutoff = 30000){
  # make a two extension hulls and mesh for spatial model
  hull <- fmesher::fm_extensions(
    poly,
    convex = c(50000, 200000),
    concave = c(350000, 500000)
  )
  mesh_spatial <- fmesher::fm_mesh_2d_inla(
    boundary = hull,
    max.edge = max.edge, # km inside and outside
    cutoff = cutoff,
    crs = fmesher::fm_crs(proj_use)
  ) # cutoff is min edge
  mesh_locs <- mesh_spatial$loc[,c(1,2)] %>% as.data.frame()
  message("Mesh created with ", dim(mesh_locs)[1], " vertices.")

  prior_range <- c(300000,0.1) # 10% chance range is smaller than 300000
  prior_sigma <- c(0.5,0.1) # 10% chance sd is larger than 0.5
  INLA::inla.spde2.pcmatern(mesh_spatial,
                            prior.range = prior_range,
                            prior.sigma = prior_sigma)
}


#' Fit INLA model
#'
#' Fit a model of bird abundance using INLA
#'
#' @param sp_code species four letter code
#' @param analysis_data list containing analysis data, including `full_count_matrix`,
#'  `all_surveys`, `species_ranges`, and `species_to_model`
#' @param proj_use projection/coordinate reference system to use. Note this
#'  should have units of kms to avoid an extremely dense mesh
#' @param study_poly sf polygon of study area polygon
#' @param covariates a data frame with columns covariate, model, mean, prec, beta
#'  which is used to define the model formula.
#' @param mod_dir directory where the model object will be saved or loaded from
#'  if it already exists
#' @param train_dat_filter a string that will be used to filter the input data,
#'  the default will not filter anything
#' @param save_mod logical. Should the model object be saved?
#' @param file_name_bit suffix attached to the file name eg to identify
#'  cross-validation fold
#'
#' @returns An INLA object with the fit model
#' @export
#'
#' @examples
fit_inla <- function(sp_code, analysis_data, proj_use, study_poly, covariates,
                     mod_dir = "data/derived-data/INLA_results/models/",
                     train_dat_filter = "TRUE", save_mod = TRUE, file_name_bit = "all"){
  message("starting model for: ", sp_code)

  model_file <- paste0(mod_dir, sp_code, "_",
                       file_name_bit, "_mod.rds")

  if(!dir.exists(mod_dir)){
    dir.create(mod_dir)
  }

  if(file.exists(model_file)){
    message("Using saved model")
    return(readRDS(model_file))
  }

  # Prepare data for this species
  sp_dat <- prep_sp_dat(analysis_data, sp_code, proj_use, train_dat_filter)

  matern_coarse <- make_mesh(study_poly, proj_use)

  # FIT MODEL WITH INLA

  # ---
  # Create mesh to model effect of time since sunrise (TSS)
  # ---
  sp_dat$Hours_Since_Sunrise <- as.numeric(sp_dat$Hours_Since_Sunrise)
  TSS_range <- range(sp_dat$Hours_Since_Sunrise)
  TSS_meshpoints <- seq(TSS_range[1]-0.1, TSS_range[2]+0.1, length.out = 11)
  TSS_mesh1D = INLA::inla.mesh.1d(TSS_meshpoints, boundary="free")
  TSS_spde = INLA::inla.spde2.pcmatern(TSS_mesh1D,
                                       prior.range = c(6,0.1),
                                       prior.sigma = c(1,0.1)) # 10% chance sd is larger than 1

  # ---
  # Model formulas
  # ---
  covariates <- covariates %>%
    mutate(components = paste0("Beta", beta, "_", covariate, '(1,model=\"', model,
                               '\"', ", mean.linear = ", mean, ", prec.linear = ",
                               prec, ")"),
           formula = paste0("Beta", beta, "_", covariate, "*", covariate, "^", beta))

  model_components <- as.formula(paste0('~
            Intercept_PC(1)+
            range_effect(1,model="linear", mean.linear = -0.046, prec.linear = 10000)+
            TSS(main = Hours_Since_Sunrise,model = TSS_spde) +
            spde_coarse(main = coordinates, model = matern_coarse) +',
         paste0(covariates$components, collapse = " + ")))

  model_formula_PC = as.formula(paste0('count ~
                  Intercept_PC +
                  log_QPAD_offset +
                  TSS +
                  range_effect * distance_from_range +
                  spde_coarse +',
               paste0(covariates$formula, collapse = " + ")))


  # ---
  # Fit with nbinomial error
  # ---

  PC_sp <- sp_dat %>% as('Spatial')
  start <- Sys.time()
  fit_INLA <- NULL
  while(is.null(fit_INLA)){
    fit_INLA <- inlabru::bru(components = model_components,

                             inlabru::like(family = "nbinomial",
                                           formula = model_formula_PC,
                                           data = PC_sp),

                    options = list(
                      control.compute = list(waic = FALSE, cpo = FALSE),
                      bru_verbose = 4))
    if ("try-error" %in% class(fit_INLA)) fit_INLA <- NULL
  }

  end <- Sys.time()
  runtime_INLA <- difftime( end,start, units="mins") %>% round(2)
  message(paste0(sp_code," - ",runtime_INLA," min to fit model"))

  if(save_mod){
    saveRDS(fit_INLA, model_file)
  }

  return(fit_INLA)
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
predict_inla <- function(dat, analysis_data, mod, sp_code, covariates, do_crps = TRUE){
  dat <- get_dist_to_range(dat, sp_code, analysis_data$species_ranges)

  dat <- get_QPAD_offsets(dat, sp_code, analysis_data$species_to_model)

  # either based on actual survey or assumes 5-minute unlimited distance survey
  offset_var <- str_subset(names(dat), "offset")

  covariates <- covariates %>%
    mutate(formula = paste0("Beta", beta, "_", covariate, "*", covariate, "^", beta))

  # get formula from model object
  mod_form <- mod$bru_info$lhoods[[1]]$formula

  if(offset_var != "log_QPAD_offset"){
    dat <- dat %>% rename(log_QPAD_offset = offset_var)
  }

  # Predictions are on log scale, and do not include variance components
  start2 <- Sys.time()
  pred <- NULL
  pred <- inlabru::generate(mod,
                   as(dat,'Spatial'),
                   formula =  mod_form,
                   n.samples = 1000)

  pred <- exp(pred)

  # Median and upper/lower credible intervals (90% CRI)
  prediction_quantiles = apply(pred,1,function(x) quantile(x,c(0.05,0.5,0.95),na.rm = TRUE))
  dat$pred_q05 <- prediction_quantiles[1,]
  dat$pred_q50 <- prediction_quantiles[2,]
  dat$pred_q95 <- prediction_quantiles[3,]
  dat$pred_CI_width_90 <- prediction_quantiles[3,] - prediction_quantiles[1,]
  dat$CV <- apply(pred,1,function(x) sd(x,na.rm = TRUE)/mean(x,na.rm = TRUE))
  if(do_crps){
    if(is.null(dat$Obs_Index)){
      warning("dat does not contain a column Obs_Index so dat cannot be ",
              "connected to full_count_matrix and crps cannot be calculated")
    } else {
      # CRPS requires comparison to observed value, needs full sample for prediction
      dat$obs_count <- analysis_data$full_count_matrix[dat$Obs_Index, sp_code]
      dat$crps <- scoringRules::crps_sample(dat$obs_count, pred)
      dat$logs <- scoringRules::logs_sample(dat$obs_count, pred)
    }
  }


  # Probability of observing species in 5-minute point count
  size <- mod$summary.hyperpar$'0.5quant'[1] # parameter of negative binomial

  # Probability of detecting species in a 5-minute point count
  #TODO: get family from INLA object so works for multiple
  prob_zero_PC <- dnbinom(0,mu=prediction_quantiles[2,],size=size)
  dat$pObs_5min <- 1-prob_zero_PC

  end2 <- Sys.time()
  runtime_pred <- difftime( end2,start2, units="mins") %>% round(2)
  message(paste0(sp_code," - ",runtime_pred," min to generate predictions"))
  return(dat %>% mutate(species = sp_code))
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
#' @param map_dir directory where the map images should be saved
#' @param train_dat_filter a string that will be used to filter the input data,
#'  the default will not filter anything
#' @param file_name_bit suffix attached to the file name eg to identify
#'  cross-validation fold
#'
#' @returns
#' @export
#'
#' @examples
map_inla_preds <- function(sp_code, analysis_data, preds, proj_use, atlas_squares, bcr_poly, study_poly,
                           map_dir = "data/derived-data/INLA_results/maps/",
                           train_dat_filter = "TRUE", file_name_bit = "all"){
  map_file <- file.path(map_dir, paste0(sp_code,"_",
                     file_name_bit, "_q50.png"))

  if(!dir.exists(map_dir)){
    dir.create(map_dir)
  }


  # Prepare data for this species
  sp_dat <- prep_sp_dat(analysis_data, sp_code, proj_use, train_dat_filter)

  # Summarize atlas_squares where species was detected
  PC_detected <- sp_dat %>%
    sf::st_intersection(atlas_squares)%>%
    as.data.frame() %>%
    group_by(sq_id) %>%
    summarize(PC_detected = as.numeric(sum(count)>0),
              PC_mean_count = mean(count) %>% round(2))

  # CL_detected <-sp_dat %>%
  #   subset(Survey_Type %in% c("Breeding Bird Atlas","Linear transect")) %>%
  #   as.data.frame() %>%
  #   group_by(sq_id) %>%
  #   summarize(CL_detected = as.numeric(sum(count)>0),
  #             CL_mean_count = mean(count))

  atlas_squares_species <- atlas_squares %>%
    relocate(geometry,.after = last_col()) %>%
    left_join(PC_detected, by = join_by(sq_id)) #%>% left_join(CL_detected)

  atlas_squares_centroids <- sf::st_centroid(atlas_squares_species)

  # Label for figure and ebird range limit

  species_name = analysis_data$ON_spcd$CommonName[which(analysis_data$ON_spcd$spcd == sp_code)]
  species_label = analysis_data$ON_spcd$Label[which(analysis_data$ON_spcd$spcd == sp_code)]

  # sf object for ebird range limit (optional - not needed for plotting)

  range <- NA
  if (sp_code %in% names(analysis_data$species_ranges)){
    range <- analysis_data$species_ranges[[sp_code]]  %>%
      sf::st_transform(sf::st_crs(study_poly)) %>%
      sf::st_intersection(study_poly)
  }

  # Plot median prediction

  colscale_relabund <-c("#FEFEFE", "#FBF7E2", "#FCF8D0", "#EEF7C2", "#CEF2B0",
                        "#94E5A0", "#51C987", "#18A065", "#008C59", "#007F53", "#006344")
  colpal_relabund <- colorRampPalette(colscale_relabund)

  lower_bound <- 0.01
  upper_bound <- quantile(preds$pred_q50,0.99,na.rm = TRUE) %>% signif(2)
  if (lower_bound >= (upper_bound/5)) lower_bound <- (upper_bound/5) %>% signif(2)

  target_raster <- terra::rast(resolution = 10, crs = proj_use,
                               extent = terra::ext(study_poly %>% terra::vect()),
                               vals = 1)

  sp_cut <- cut.fn(df = preds,
                   target_raster = target_raster,
                   column_name = "pred_q50",
                   lower_bound = lower_bound,
                   upper_bound = upper_bound)

  raster_q50 <- sp_cut$raster

  # Median of posterior
  plot_q50 <- do_res_plot(raster_q50, "Relative Abundance",
                          "Per 5-minute point count", "(Posterior Median)",
                          atlas_squares_centroids, bcr_poly,
                          colpal_relabund, species_label,
                          "levs",
                          map_file)

  # Plot uncertainty in prediction (width of 90% CRI)

  colscale_uncertainty <- c("#FEFEFE", "#FFF4B3", "#F5D271", "#F2B647", "#EC8E00", "#CA302A")
  colpal_uncertainty <- colorRampPalette(colscale_uncertainty)

  lower_bound <- 0.01
  upper_bound <- quantile(preds$pred_CI_width_90,0.99,na.rm = TRUE) %>% signif(2)
  if (lower_bound >= (upper_bound/5)) lower_bound <- (upper_bound/5) %>% signif(2)

  raster_CI_width_90 <- cut.fn(df = preds,
                               target_raster = target_raster,
                               column_name = "pred_CI_width_90",
                               lower_bound = lower_bound,
                               upper_bound = upper_bound)$raster

  plot_CI_width_90 <- do_res_plot(raster_CI_width_90, "Relative Uncertainty",
                                  "Per 5-minute point count", "Width of 90% CI",
                                  atlas_squares_centroids,  bcr_poly,
                                  colpal_uncertainty, species_label,
                                  "levs",
                                  map_file %>% str_replace("_q50", "_CI_width_90"))

  # Plot uncertainty in prediction (coefficient of variation)

  colscale_uncertainty <- c("#FEFEFE", "#FFF4B3", "#F5D271", "#F2B647", "#EC8E00", "#CA302A")
  colpal_uncertainty <- colorRampPalette(colscale_uncertainty)

  cut_levs <- c(-0.1,0.25,0.5,1,2,5,2000)
  cut_levs_labs <- c("0 to 0.25",
                     "0.25 to 0.5",
                     "0.5 to 1",
                     "1 to 2",
                     "2 to 5",
                     "> 5")

  preds$CV_levs <- cut(as.data.frame(preds)[,"CV"],
                                cut_levs,labels=cut_levs_labs)
  raster_CV <- stars::st_rasterize(preds %>% dplyr::select(CV_levs, geometry),
                                   nx = dim(raster_q50)[1],ny = dim(raster_q50)[2])

  plot_CV <- do_res_plot(raster_CV, "Coef. of Variation",
                         "Per 5-minute point count", "",
                         atlas_squares_centroids, bcr_poly,
                         colpal_uncertainty, species_label,
                         "CV_levs",
                         map_file %>% str_replace("_q50", "_CV"))

  # Plot probability of observing species in a 5-minute point count

  colscale_pObs <- c("#FEFEFE",RColorBrewer::brewer.pal(5,"BuGn")[2:5])
  colpal_pObs <- colorRampPalette(colscale_pObs)

  cut_levs <- c(-0.1,0.01,0.05,0.125,0.5,1)
  cut_levs_labs <- c("0 to 0.01",
                     "0.01 to 0.05",
                     "0.05 to 0.125",
                     "0.125 to 0.50",
                     "0.50 to 1")

  preds$pObs_levs <- cut(as.data.frame(preds)[,"pObs_5min"],
                                  cut_levs,labels=cut_levs_labs)

  raster_pObs = stars::st_rasterize(preds %>% dplyr::select(pObs_levs, geometry),
                                    nx = dim(raster_q50)[1], ny = dim(raster_q50)[2])

  plot_pObs <- do_res_plot(raster_pObs, "Prob. of Observation",
                           "Per 5-minute point count",
                           "(Posterior Median)",
                           atlas_squares_centroids, bcr_poly,
                           colpal_pObs, species_label,
                           "pObs_levs",
                           map_file %>% str_replace("_q50", "_PObs"))

  # Density estimate (per m2) - subtract detectability offset
  species_offsets <- subset(analysis_data$species_to_model, Species_Code_BSC == sp_code)

  log_offset_5min <- 0
  if (species_offsets$offset_exists == TRUE) log_offset_5min <- species_offsets$log_offset_5min

  if (log_offset_5min != 0){

    preds$density_per_ha_q50 <- preds$pred_q50 / exp(log_offset_5min) * 10000

    colscale_relabund <-c("#FEFEFE", "#FBF7E2", "#FCF8D0", "#EEF7C2", "#CEF2B0",
                          "#94E5A0", "#51C987", "#18A065", "#008C59", "#007F53", "#006344")
    colpal_relabund <- colorRampPalette(colscale_relabund)

    lower_bound <- 0.01
    upper_bound <- quantile(preds$density_per_ha_q50,0.99,na.rm = TRUE) %>% signif(2)
    if (lower_bound >= (upper_bound/5)) lower_bound <- (upper_bound/5) %>% signif(2)

    sp_cut <- cut.fn(df = preds,
                     target_raster = target_raster,
                     column_name = "density_per_ha_q50",
                     lower_bound = lower_bound,
                     upper_bound = upper_bound)

    raster_dens <- sp_cut$raster

    # Median of posterior
    plot_dens <- do_res_plot(raster_dens, "Density", "Males per hectare",
                             "(Posterior Median)",
                             atlas_squares_centroids, bcr_poly,
                             colpal_relabund, species_label, "levs",
                             map_file %>% str_replace("_q50", "_density"))
  }
}

#' Build INLA prediction map
#'
#' @param pred_rast Raster of predictions
#' @param title map title describing the metric shown
#' @param subtitle subtitle under title with more details
#' @param subsubtitle subtitle under subtitle eg with units or further description
#' @param samp_grid Points in sampling grid where species was or wasn't detected
#' @param col_pal_fn function to create colour palette
#' @param species_label species name which will be added to the map
#' @param levs_nm name of levels of colour ramp
#' @param file_nm file name where map should be saved
#'
#' @returns saves the map to `file_nm`
#'
#' @examples
do_res_plot <- function(pred_rast, title, subtitle, subsubtitle = "", samp_grid, bcr_poly,
                        col_pal_fn, species_label, levs_nm, file_nm){
  res_plot <- ggplot() +

    stars::geom_stars(data = pred_rast, na.rm = TRUE) +
    scale_fill_manual(name = paste0("<span style='font-size:13pt'>", title,
                                    "</span><br><span style='font-size:7pt'>", subtitle,
                                    "</span><br><span style='font-size:7pt'>", subsubtitle,
                                    "</span>"),
                      values = col_pal_fn(length(levels(pred_rast[[levs_nm]]))),
                      drop = FALSE, na.translate = FALSE)+

    # BCR boundaries
    geom_sf(data = bcr_poly, fill = "transparent", col = "gray20", linewidth = 0.5)+

    # Point count detections and surveyed squares
    geom_sf(data = subset(samp_grid, !is.na(PC_detected)),
            aes(col = as.factor(PC_detected)), size = 0.5, stroke = 0, shape = 16)+
    scale_colour_discrete(type = c("gray70", "black"), guide = NULL)+

    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())+
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    theme(legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(5,10,5,-20),
          legend.title.align=0.5,
          legend.title = ggtext::element_markdown(lineheight=.9,hjust = 1),
          legend.justification = c(1,1),
          legend.position = "inside",
          legend.position.inside = c(0.99,0.95))+
    theme(legend.key = element_rect(fill = "transparent", colour = "transparent"))+

    annotate(geom="text",x=1700000,y=1930000, label= paste0(species_label),lineheight = .85,hjust = 0,size=6,fontface =2) +
    annotate(geom="text",x=2155000,y=530000, label= paste0("Prepared on ",Sys.Date()),size=2,lineheight = .75,hjust = 0,color="gray60")+

    guides(fill = guide_legend(order = 1),
           size = guide_legend(order = 2))

  png(file_nm, width=10, height=6.5, units="in", res=1000, type="cairo")
  print(res_plot)
  dev.off()

  return(res_plot)
}


#' Rasterize a series of spatial predictions (needed for plotting)

cut.fn <- function(df = NA,
                   target_raster = NA,
                   column_name = NA,
                   lower_bound = NA,
                   upper_bound = NA){

  max_val <- upper_bound
  max_val <- ifelse(is.na(max_val), 0, max_val)
  max_lev <- ifelse(max_val > 1.6, 4,ifelse(max_val > 0.8, 4, 3))

  cut_levs <- signif(max_val/(2^((max_lev-1):0)), 2)
  cut_levs <- unique(cut_levs)
  cut_levs <- ifelse(is.na(cut_levs), 0, cut_levs)

  if (lower_bound %in% cut_levs) cut_levs <- cut_levs[-which(cut_levs == lower_bound)]
  if (lower_bound > min(cut_levs)) cut_levs = cut_levs[-which(cut_levs < lower_bound)]

  max_lev <- length(cut_levs)

  cut_levs_labs <- c(paste0("0-",lower_bound),
                     paste(lower_bound, cut_levs[1], sep="-"),
                     paste(cut_levs[-max_lev], cut_levs[-1], sep="-"),
                     paste(cut_levs[max_lev], "+"))

  cut_levs <- c(-1, lower_bound, cut_levs, 1000) %>% unique()

  df <- mutate(df, levs = cut(.data[[column_name]], cut_levs, labels = cut_levs_labs))
  # TODO: change this to use terra. It is faster now and easier I think
  tgt <- stars::st_as_stars(target_raster)
  tmp = stars::st_rasterize(df %>% dplyr::select(levs, geometry),
                            nx = dim(tgt)[1],ny = dim(tgt)[2])

  return(list(raster = tmp,cut_levs = cut_levs))
}


# evaluate model performance based on predicted vs observed count
evaluate_preds <- function(pred, mod, sp_code, analysis_data){
  obs_count <- analysis_data$full_count_matrix[pred$Obs_Index, sp_code]

  rmse_pred <- sqrt(mean((obs_count - pred$pred_q50)^2))

  size <- mod$summary.hyperpar$'0.5quant'[1]

  auc_pred <- pROC::auc(response = obs_count > 0, predictor = pred$pObs_5min)

  # log pointwise predictive density, I think...
  lppd_pred <- sum(dnbinom(obs_count, mu = pred$pred_q50, size = size, log = TRUE))

  med_crps <- median(pred$crps)

  med_crps_pres <- median(pred$crps[obs_count > 0])
  med_crps_abs <- median(pred$crps[obs_count == 0])

  med_logs <- median(pred$logs)

  med_logs_pres <- median(pred$logs[obs_count > 0])
  med_logs_abs <- median(pred$logs[obs_count == 0])

  tibble(species = sp_code, fold = unique(pred$Crossval_Fold),
         rmse = rmse_pred, auc = as.numeric(auc_pred), lppd = lppd_pred,
         med_crps = med_crps, med_logs = med_logs,
         med_crps_pres = med_crps_pres, med_logs_pres = med_logs_pres,
         med_crps_abs = med_crps_abs, med_logs_abs = med_logs_abs)
}
